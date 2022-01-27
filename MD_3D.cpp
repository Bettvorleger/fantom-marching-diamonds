#include <fantom/algorithm.hpp>
#include <fantom/datastructures/interfaces/Field.hpp>
#include <fantom/graphics.hpp>
#include <fantom/register.hpp>
#include <math.h>

#include <fantom-plugins/utils/Graphics/HelperFunctions.hpp>
#include <fantom-plugins/utils/Graphics/ObjectRenderer.hpp>

#include <stdexcept>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>

using namespace fantom;

namespace
{

    class MarchingDiamondsAlgorithm : public VisAlgorithm
    {

    protected:
        std::shared_ptr<const Field<3, Scalar>> field;
        std::shared_ptr<const Grid<3>> grid;
        std::vector<Point3> controlDiamond;

    public:
        struct Options : public VisAlgorithm::Options
        {
            Options(fantom::Options::Control &control)
                : VisAlgorithm::Options(control)
            {
                add<Field<3, Scalar>>("Field", "A 2D scalar field", definedOn<Grid<3>>(Grid<3>::Points));
                add<double>("Isovalue", "The desired isovalue to show.", 1.05);
                add<Color>("Color", "The color of the graphics.", Color(1.0, 0.0, 0.0));
                add<int>("Max Splitting", "Max amount of splitting steps.", 0, &acceptNumber);
            }

            static int acceptNumber(const int &i)
            {
                return std::max(i, 0);
            }
        };

        struct VisOutputs : public VisAlgorithm::VisOutputs
        {
            VisOutputs(fantom::VisOutputs::Control &control)
                : VisAlgorithm::VisOutputs(control)
            {
                addGraphics("Iso-surface");
                addGraphics("Diamond");
            }
        };

        MarchingDiamondsAlgorithm(InitData &data)
            : VisAlgorithm(data)
        {
        }

        virtual void execute(const Algorithm::Options &options, const volatile bool &abortFlag) override
        {
            field = options.get<Field<3, Scalar>>("Field");
            std::shared_ptr<const Function<Scalar>> function = options.get<Function<Scalar>>("Field");
            double isovalue = options.get<double>("Isovalue");
            Color color = options.get<Color>("Color");
            int maxSplit = options.get<int>("Max Splitting");

            //check if input field is set
            if (!field)
            {
                debugLog() << "Input field not set." << std::endl;
                return;
            }

            grid = std::dynamic_pointer_cast<const Grid<3>>(function->domain());

            //check if grid can be used or is correct
            if (!grid)
            {
                throw std::logic_error("Wrong type of grid!");
            }

            //map of edges point indices as key and tetrahedron index as aggregated value
            std::map<std::pair<size_t, size_t>, std::vector<size_t>> tetrahedronEdges;
            size_t tetrCellCount = grid->numCells();

            //vector of point and indices vector for diamonds
            std::map<std::vector<size_t>, std::vector<Point3>> diamonds;

            //aggregates intersection for each triangle by index
            std::unordered_map<size_t, std::vector<Point3>> tetrIntersectPoints;

            //points for final iso segment drawing
            std::vector<VectorF<3>> isoSurfaces;

            //load all cells into vector
            for (size_t i = 0; i < tetrCellCount; i++)
            {
                Cell cell = grid->cell(i);
                std::vector<size_t> pI;

                for (size_t j = 0; j < 4; j++)
                    pI.push_back(cell.index(j));

                std::sort(pI.begin(), pI.end());

                tetrahedronEdges[std::make_pair(pI[0], pI[1])].push_back(i);
                tetrahedronEdges[std::make_pair(pI[0], pI[2])].push_back(i);
                tetrahedronEdges[std::make_pair(pI[0], pI[3])].push_back(i);
                tetrahedronEdges[std::make_pair(pI[1], pI[2])].push_back(i);
                tetrahedronEdges[std::make_pair(pI[1], pI[3])].push_back(i);
                tetrahedronEdges[std::make_pair(pI[2], pI[3])].push_back(i);
            }

            tetrIntersectPoints = marchingDiamonds(tetrahedronEdges, tetrIntersectPoints, diamonds, tetrCellCount, isovalue, maxSplit + 1);

            //add triangle intersection points to final iso segment vector
            for (auto &tr : tetrIntersectPoints)
            {
                if (tr.second.size() >= 3)
                {
                    isoSurfaces.insert(isoSurfaces.end(), {VectorF<3>(tr.second[0]), VectorF<3>(tr.second[1]), VectorF<3>(tr.second[3])});
                }
            }

            /************LOGGING***************/
            std::vector<VectorF<3>> vertGrid;

            size_t k = controlDiamond.size() - 1;

            for (size_t i = 0; i < k; i++)
            {
                vertGrid.insert(vertGrid.end(), {VectorF<3>(controlDiamond[i]), VectorF<3>(controlDiamond[i + 1])});
                vertGrid.insert(vertGrid.end(), {VectorF<3>(controlDiamond[i]), VectorF<3>(controlDiamond[k])});
                vertGrid.insert(vertGrid.end(), {VectorF<3>(controlDiamond[i]), VectorF<3>(controlDiamond[k + 1])});
            }

            if (abortFlag)
            {
                return;
            }

            auto const &system = graphics::GraphicsSystem::instance();
            std::string resourcePath = PluginRegistrationService::getInstance().getResourcePath("utils/Graphics");

            auto bs = graphics::computeBoundingSphere(isoSurfaces);
            std::vector<unsigned int> indices(isoSurfaces.size());
            std::iota(indices.begin(), indices.end(), 0);               //fill indices incrementally 0,1,2,3.....
            auto norm = graphics::computeNormals(isoSurfaces, indices); //calc normals between every surface and the viewer

            std::shared_ptr<graphics::Drawable> surface = system.makePrimitive(
                graphics::PrimitiveConfig{graphics::RenderPrimitives::TRIANGLES}
                    .vertexBuffer("position", system.makeBuffer(isoSurfaces))
                    .vertexBuffer("normal", system.makeBuffer(norm))
                    .indexBuffer(system.makeIndexBuffer(indices))
                    .uniform("color", color)
                    .renderOption(graphics::RenderOption::Blend, true)
                    .boundingSphere(bs),
                system.makeProgramFromFiles(resourcePath + "shader/surface/phong/singleColor/vertex.glsl",
                                            resourcePath + "shader/surface/phong/singleColor/fragment.glsl"));
            setGraphics("Iso-surface", surface);

            std::shared_ptr<graphics::Drawable> diamond = system.makePrimitive(
                graphics::PrimitiveConfig{graphics::RenderPrimitives::LINES}
                    .vertexBuffer("in_vertex", system.makeBuffer(vertGrid))
                    .uniform("u_lineWidth", 3.0f)
                    .uniform("u_color", color)
                    .boundingSphere(bs),
                system.makeProgramFromFiles(resourcePath + "shader/line/noShading/singleColor/vertex.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/fragment.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/geometry.glsl"));
            setGraphics("Diamond", diamond);
        }

        /**
         * @brief Search for diamonds and find intersections inside them. If splitting is neccessary the function is called recursively with new tetrahedron.
         * 
         * @param trianglePoints - all triangle cells as point vectors
         * @param trIntersectPoints - map for aggregating intersection points for segment connection
         * @param diamondPoints - point vector for all found diamonds, first call empty, then reused for splitting
         * @param newtetrahedron - new tetrahedron found through splitting
         * @param trCellCount - triangle cell count changed through splitting tetrahedron
         * @param isovalue - isovalue to be displayed
         * @param splitCount - recursion depth controling splitting steps
         * @return std::unordered_map<size_t, std::vector<Point3>> - map for aggregating intersection points for segment connection
         */
        template <typename T, typename U, typename V>
        std::unordered_map<size_t, std::vector<Point3>> marchingDiamonds(T tetrahedronEdges, U tetrIntersectPoints, V diamonds, size_t tetrCellCount, double isovalue, int splitCount)
        {
            if (splitCount == 0)
            {
                return tetrIntersectPoints;
            }
            tetrIntersectPoints.clear();

            auto evaluator = field->makeEvaluator();
            //get control points of the grid
            const ValueArray<Point3> &points = grid->points();

            int test = 1;
            //FIND DIAMONDS
            for (auto &tEdge : tetrahedronEdges)
            {
                if (tEdge.second.size() >= 4)
                {
                    std::vector<Point3> d;
                    std::vector<std::vector<size_t>> tetrIndices;
                    std::unordered_set<size_t> diamondIndices;

                    for (auto &tEdgeTetr : tEdge.second)
                    {
                        Cell cell = grid->cell(tEdgeTetr);
                        std::vector<size_t> tempTetrI;
                        for (size_t i = 0; i < 4; i++)
                        {
                            if (cell.index(i) != (tEdge.first).first && cell.index(i) != (tEdge.first).second)
                                tempTetrI.push_back(cell.index(i));
                        }
                        tetrIndices.push_back(tempTetrI);
                    }

                    size_t k = 0;
                    std::vector<size_t> visitedInd;
                    for (size_t i = 0; i < tEdge.second.size() - 1; i++)
                    {
                        for (size_t j = 0; j < tEdge.second.size(); j++)
                        {
                            if (k == j || std::find(visitedInd.begin(), visitedInd.end(), j) != visitedInd.end())
                                continue;

                            if (tetrIndices[k][0] == tetrIndices[j][0])
                            {
                                diamondIndices.insert(tetrIndices[k][1]);
                                diamondIndices.insert(tetrIndices[k][0]);
                                diamondIndices.insert(tetrIndices[j][1]);
                                k = j;
                                visitedInd.push_back(k);
                                break;
                            }
                            else if (tetrIndices[k][0] == tetrIndices[j][1])
                            {
                                diamondIndices.insert(tetrIndices[k][1]);
                                diamondIndices.insert(tetrIndices[k][0]);
                                diamondIndices.insert(tetrIndices[j][0]);
                                k = j;
                                visitedInd.push_back(k);
                                break;
                            }
                            else if (tetrIndices[k][1] == tetrIndices[j][0])
                            {
                                diamondIndices.insert(tetrIndices[k][0]);
                                diamondIndices.insert(tetrIndices[k][1]);
                                diamondIndices.insert(tetrIndices[j][1]);
                                k = j;
                                visitedInd.push_back(k);
                                break;
                            }
                            else if (tetrIndices[k][1] == tetrIndices[j][1])
                            {
                                diamondIndices.insert(tetrIndices[k][0]);
                                diamondIndices.insert(tetrIndices[k][1]);
                                diamondIndices.insert(tetrIndices[j][0]);
                                k = j;
                                visitedInd.push_back(k);
                                break;
                            }
                        }
                    }

                    //diamondIndices.erase((tEdge.first).first);
                    //diamondIndices.erase((tEdge.first).second);

                    for (auto &dI : diamondIndices)
                    {
                        d.push_back(points[dI]);
                    }
                    d.push_back(points[(tEdge.first).first]);
                    d.push_back(points[(tEdge.first).second]);

                    if (test > 0)
                    {
                        infoLog() << tEdge.second.size() << std::endl;
                        for (auto &dI : d)
                        {
                            infoLog() << dI << " ";
                        }
                        infoLog() << std::endl;
                        test--;
                        controlDiamond = d;
                    }

                    diamonds[tEdge.second] = d;
                }
            }

            //check if splitting was neccesary
            bool splitting = false;

            //find all intersections on grid, using constant iterator for manipulating the original structure
            for (auto dP = diamonds.cbegin(); dP != diamonds.cend();)
            {
                std::vector<Point3> tetrIntPoints = findIntersections(dP->second, isovalue, evaluator);

                //splitting for more than one intersection found
                if (tetrIntPoints.size() == 2)
                {
                    splitting = true;

                    continue;
                }
                else if (tetrIntPoints.size() == 1)
                {
                    //add the one intersection two all k Tetraeders
                    for (size_t i = 0; i < (dP->first).size(); i++)
                        tetrIntersectPoints[dP->first[i]].push_back(tetrIntPoints[0]);
                }
                ++dP;
            }

            infoLog() << tetrIntersectPoints.size() << std::endl;

            if (splitting)
            {
                return marchingDiamonds(tetrahedronEdges, tetrIntersectPoints, diamonds, tetrCellCount, isovalue, splitCount - 1);
            }
            return tetrIntersectPoints;
        }

    private:
        /**
         * @brief Handle intersection search for given diamond points, isovalue and evaluator
         * 
         * @tparam T 
         * @param dP 
         * @param isovalue 
         * @param evaluator 
         * @return std::vector<Point3> - vector of intersection points, two at most
         */
        template <typename T>
        std::vector<Point3> findIntersections(std::vector<Point3> dP, double isovalue, T &evaluator)
        {
            std::vector<Point3> tetrIntPoints;
            std::vector<double> s;

            size_t k = dP.size() - 2;

            for (size_t i = 0; i < k + 1; i++)
            {
                evaluator->reset(dP[i]);
                s.push_back(norm(evaluator->value()));
            }

            double C = -2 * sin(2 * M_PI / k) + sin(4 * M_PI / k);
            double D = (k - 2) * abs(3 * sin(2 * M_PI / k) - 4 * sin(4 * M_PI / k) - sin(6 * M_PI / k));

            std::vector<double> intersects = calcIntersection(s, C, D, isovalue);

            for (auto &z : intersects)
            {
                if (z >= 0 && z <= 2)
                {
                    tetrIntPoints.push_back(barycentricTransform(dP, C, D, z));
                }
            }

            if (tetrIntPoints.size() > 2)
                infoLog() << "Mehr als 2 Intersections, nicht möglich" << std::endl;

            return tetrIntPoints;
        }

        /**
         * @brief Cubic solving according to Schwarze[1] and Shelbey[2] 
         * 
         * [1] = Jochen Schwarze. Cubic and quartic roots. In Graphics Gems, pages 404–407. Academic Press Professional, Inc., San Diego, CA, USA, 1990.
         * [2] = Shelbey, Samuel, ed. (1975). CRC Standard Mathematical Tables. CRC Press. 
         * 
         * @param s 
         * @param C 
         * @param D 
         * @param E 
         * @param isovalue 
         * @return std::vector<double> 
         */
        std::vector<double> calcIntersection(std::vector<double> s, double C, double D, double isovalue)
        {
            std::vector<double> intersections;

            double sRingSum = 0;
            size_t k = s.size() - 2;

            for (size_t i = 0; i < k; i++)
                sRingSum += s[i];

            double t = k * isovalue - sRingSum;

            double denominator = 4 * abs(C) * (-t) + D * (s[k + 1] - s[k]);

            if (denominator == 0)
                return intersections;

            double a = (16 * abs(C) * t + 6 * D * (s[k] - isovalue)) / denominator;
            double b = (16 * abs(C) * (-t) + 12 * D * (isovalue - s[k])) / denominator;
            double c = (8 * D * (s[k] - isovalue)) / denominator;

            double p = pow(a, 2) / 3 + b;
            double q = c - (2 * pow(a, 3) / 27) - a * b / 3;

            if (p != 0 && (4 * pow(p, 3) + 27 * pow(q, 3)) < 0)
            {
                for (size_t i = 0; i < 3; i++)
                {
                    double y = 2 * sqrt(-p / 3) * cos(1 / 3 * acos((3 * q / 2 * p) * sqrt(-3 / p)) - (2 * M_PI * k) / 3);
                    intersections.push_back(y);
                }
            }

            return intersections;
        }

        std::vector<Point3> findDiamond()
        {
            return std::vector<Point3>();
        }

        Point3 barycentricTransform(std::vector<Point3> diamondPoints, double C, double D, double z)
        {
            Point3 p;
            size_t k = diamondPoints.size() - 2;

            double E = 4 * k * abs(C) * z * pow((2 - z), 2) + D * (pow((2 - z), 3) + pow(z, 3));

            double a1 = 4 * abs(C) * z * pow((2 - z), 2) / E;
            double a2 = D * pow((2 - z), 3) / E;
            double a3 = D * pow(z, 3) / E;

            for (size_t i = 0; i < k; i++)
                p += diamondPoints[i] * a1;
            p += a2 * diamondPoints[k] + a3 * diamondPoints[k + 1];

            return p;
        }

        template <typename T>
        Point3 linearInterpolation(Point3 p0, Point3 p1, double isovalue, T &evaluator)
        {
            Point3 intersectPoint;

            evaluator->reset(p0);
            double s0 = norm(evaluator->value());
            evaluator->reset(p1);
            double s1 = norm(evaluator->value());

            double alpha = (isovalue - s1) / (s0 - s1);

            if (0 <= alpha && alpha <= 1)
            {
                intersectPoint = p1 + alpha * (p0 - p1);
            }
            return intersectPoint;
        }
    };

    AlgorithmRegister<MarchingDiamondsAlgorithm> dummy("Iso-Surface/MarchingDiamonds3D",
                                                       "Show iso-surface in 3D grid for given isovalue.");
} // namespace