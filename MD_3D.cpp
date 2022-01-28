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
        std::map<std::pair<size_t, size_t>, std::vector<size_t>> tetrahedronEdges;
        std::map<size_t, Point3> gridPoints;

        std::vector<Point3> controlDiamond;
        std::vector<VectorF<3>> vertGrid;
        std::vector<VectorF<3>> edge;

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
                addGraphics("Diamond-Edge");
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
            size_t gridPointCount = grid->numPoints();

            //vector of point and indices vector for diamonds
            std::map<std::vector<size_t>, std::pair<std::vector<Point3>, std::pair<size_t, size_t>>> diamonds;

            //aggregates intersection for each triangle by index
            std::unordered_map<size_t, std::vector<Point3>> tetrIntersectPoints;

            //points for final iso segment drawing
            std::vector<VectorF<3>> isoSurfaces;

            const ValueArray<Point3> &points = grid->points();

            //load all cells into vector
            for (size_t i = 0; i < tetrCellCount; i++)
            {
                Cell cell = grid->cell(i);
                std::vector<size_t> pI;

                for (size_t j = 0; j < 4; j++)
                {
                    pI.push_back(cell.index(j));
                    gridPoints[cell.index(j)] = points[cell.index(j)];
                }

                std::sort(pI.begin(), pI.end());

                tetrahedronEdges[std::make_pair(pI[0], pI[1])].push_back(i);
                tetrahedronEdges[std::make_pair(pI[0], pI[2])].push_back(i);
                tetrahedronEdges[std::make_pair(pI[0], pI[3])].push_back(i);
                tetrahedronEdges[std::make_pair(pI[1], pI[2])].push_back(i);
                tetrahedronEdges[std::make_pair(pI[1], pI[3])].push_back(i);
                tetrahedronEdges[std::make_pair(pI[2], pI[3])].push_back(i);
            }

            tetrIntersectPoints = marchingDiamonds(tetrahedronEdges, tetrIntersectPoints, diamonds, tetrCellCount, gridPointCount, isovalue, maxSplit + 1);

            //add triangle intersection points to final iso segment vector
            for (auto &tr : tetrIntersectPoints)
            {
                if (tr.second.size() >= 3)
                {
                    isoSurfaces.insert(isoSurfaces.end(), {VectorF<3>(tr.second[0]), VectorF<3>(tr.second[1]), VectorF<3>(tr.second[2])});
                }
            }

            /************LOGGING****************/
            size_t k = controlDiamond.size() - 2;

            for (size_t i = 0; i < k; i++)
            {
                size_t index = (i + 1) % k;
                vertGrid.insert(vertGrid.end(), {VectorF<3>(controlDiamond[i]), VectorF<3>(controlDiamond[index])});
                vertGrid.insert(vertGrid.end(), {VectorF<3>(controlDiamond[i]), VectorF<3>(controlDiamond[k])});
                vertGrid.insert(vertGrid.end(), {VectorF<3>(controlDiamond[i]), VectorF<3>(controlDiamond[k + 1])});
            }
            edge.insert(edge.end(), {VectorF<3>(controlDiamond[k]), VectorF<3>(controlDiamond[k + 1])});

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

            std::shared_ptr<graphics::Drawable> diamondE = system.makePrimitive(
                graphics::PrimitiveConfig{graphics::RenderPrimitives::LINES}
                    .vertexBuffer("in_vertex", system.makeBuffer(edge))
                    .uniform("u_lineWidth", 3.0f)
                    .uniform("u_color", Color(0.0, 1.0, 0.0))
                    .boundingSphere(bs),
                system.makeProgramFromFiles(resourcePath + "shader/line/noShading/singleColor/vertex.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/fragment.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/geometry.glsl"));
            setGraphics("Diamond-Edge", diamondE);
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
        std::unordered_map<size_t, std::vector<Point3>> marchingDiamonds(T tetrahedronEdges, U tetrIntersectPoints, V diamonds, size_t tetrCellCount, size_t gridPointCount, double isovalue, int splitCount)
        {
            if (splitCount == 0)
            {
                return tetrIntersectPoints;
            }
            tetrIntersectPoints.clear();

            auto evaluator = field->makeEvaluator();

            //find diamonds sharing an edge
            for (auto &tEdge : tetrahedronEdges)
            {
                std::vector<size_t> tetrIndices = tEdge.second;
                std::pair<std::vector<Point3>, std::pair<size_t, size_t>> diamondEl = findDiamond(tEdge);

                if (!diamondEl.first.empty())
                    diamonds[tetrIndices] = diamondEl;
            }

            controlDiamond = (diamonds.begin()->second).first;

            //check if splitting was neccesary
            bool splitting = false;

            std::map<std::vector<size_t>, std::pair<std::vector<Point3>, std::pair<size_t, size_t>>> newDiamondPoints;

            //find all intersections on grid, using constant iterator for manipulating the original structure
            for (auto dP = diamonds.cbegin(); dP != diamonds.cend();)
            {
                std::vector<size_t> tetrIndices = dP->first;
                std::vector<Point3> diamPoints = (dP->second).first;
                std::pair<size_t, size_t> edgeE = (dP->second).second;
                size_t kTetrCount = tetrIndices.size();
                std::vector<Point3> tetrIntPoints = findIntersections(diamPoints, isovalue, evaluator);

                //splitting for more than one intersection found
                if (tetrIntPoints.size() == 2)
                {
                    splitting = true;

                    //new point halfway between the two intersections
                    Point3 newVertex = (tetrIntPoints[0] + tetrIntPoints[1]) / 2;

                    //new triangles and diamonds created through splitting with new vertex, added by expanding cell indices
                    gridPoints[gridPointCount] = newVertex;

                    std::vector<size_t> oldTetrInd = tetrahedronEdges[edgeE];
                    std::vector<size_t> newTetrInd{tetrCellCount, tetrCellCount + 1, tetrCellCount + 2, tetrCellCount + 3};

                    tetrahedronEdges[std::make_pair(edgeE.first, gridPointCount)] = oldTetrInd;
                    tetrahedronEdges[std::make_pair(edgeE.second, gridPointCount)] = newTetrInd;

                    gridPointCount++;
                    tetrCellCount += 4;

                    //erase the two original triangles and found diamond
                    tetrahedronEdges.erase(edgeE);
                    dP = diamonds.erase(dP);

                    continue;
                }
                else if (tetrIntPoints.size() == 1)
                {
                    //add the one intersection two all k Tetraeders
                    for (auto tetrInd : dP->first)
                        tetrIntersectPoints[tetrInd].push_back(tetrIntPoints[0]);
                }
                ++dP;
            }

            infoLog() << tetrIntersectPoints.size() << std::endl;

            if (splitting)
            {
                //return marchingDiamonds(tetrahedronEdges, tetrIntersectPoints, diamonds, tetrCellCount, gridPointCount, isovalue, splitCount - 1);
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
         * @brief Cubic solving according to Schwarze[1], which is based on Cardano’s formula for cubic equations
         * 
         * [1] = Jochen Schwarze. Cubic and quartic roots. In Graphics Gems, pages 404–407. Academic Press Professional, Inc., San Diego, CA, USA, 1990.
         * 
         * @param s - vector of values on grid points
         * @param C,D,E - helper numbers for calculation
         * @param isovalue 
         * @return std::vector<double> 
         */
        std::vector<double> calcIntersection(std::vector<double> s, double C, double D, double isovalue)
        {
            std::vector<double> intersections{999};

            double sRingSum = 0;
            size_t k = s.size() - 2;

            for (size_t i = 0; i < k; i++)
                sRingSum += s[i];

            double t = k * isovalue - sRingSum;

            double denominator = 4 * abs(C) * (-t) + D * (s[k + 1] - s[k] - isovalue);

            if (abs(denominator) < 0.00000000001)
                return intersections;
            intersections.clear();

            double a = (16 * abs(C) * t + D * (6 * s[k] - isovalue)) / denominator;
            double b = (16 * abs(C) * (-t) + 4 * D * (isovalue - 3 * s[k])) / denominator;
            double c = (4 * D * (2 * s[k] - isovalue)) / denominator;

            double p = (1 / 3) * b - (pow(a, 2) / 9);
            double q = (1 / 2) * c + (pow(a, 3) / 27) - ((a * b) / 6);

            double determinant = pow(q, 2) + pow(p, 3);
            double u = cbrt(-q + sqrt(D));
            double v = cbrt(-q - sqrt(D));

            if (determinant > 0)
                intersections.push_back(u + v);
            else if (determinant == 0)
            {
                intersections.push_back(u + v);
                intersections.push_back(-(u + v) / 2);
            }
            else
            {
                double y1 = 2 * sqrt(-p) * cos((1 / 3) * (-q / sqrt(-pow(p, 3))));
                double y2 = -2 * sqrt(-p) * cos((1 / 3) * (acos(-q / sqrt(-pow(p, 3))) + M_PI));
                double y3 = -2 * sqrt(-p) * cos((1 / 3) * (acos(-q / sqrt(-pow(p, 3))) - M_PI));
                intersections.insert(intersections.end(), {y1, y2, y3});
            }

            return intersections;
        }

        template <typename T>
        std::pair<std::vector<Point3>, std::pair<size_t, size_t>> findDiamond(T tEdge)
        {
            std::pair<std::vector<Point3>, std::pair<size_t, size_t>> diamondElement;
            std::pair<size_t, size_t> sharedEdge = tEdge.first;
            std::vector<size_t> tetrIndices = tEdge.second;
            size_t kTetrCount = tetrIndices.size();

            if (kTetrCount >= 4)
            {
                std::vector<Point3> d;
                std::vector<std::vector<size_t>> tetrPointIndices;
                std::vector<size_t> diamondIndices;

                for (auto &tI : tetrIndices)
                {
                    Cell cell = grid->cell(tI);
                    std::vector<size_t> tempTetrI;
                    for (size_t i = 0; i < 4; i++)
                    {
                        if (cell.index(i) != sharedEdge.first && cell.index(i) != sharedEdge.second)
                            tempTetrI.push_back(cell.index(i));
                    }
                    tetrPointIndices.push_back(tempTetrI);
                }

                size_t k = 0;
                std::vector<size_t> visitedInd;
                for (size_t i = 0; i < kTetrCount; i++)
                {
                    for (size_t j = 0; j < kTetrCount; j++)
                    {
                        if (k == j || std::find(visitedInd.begin(), visitedInd.end(), j) != visitedInd.end())
                            continue;

                        if (tetrPointIndices[k][0] == tetrPointIndices[j][0] || tetrPointIndices[k][0] == tetrPointIndices[j][1])
                        {
                            if (k == 0)
                                diamondIndices.push_back(tetrPointIndices[k][1]);
                            diamondIndices.push_back(tetrPointIndices[k][0]);
                            visitedInd.push_back(k);
                            k = j;
                            break;
                        }
                        else if (tetrPointIndices[k][1] == tetrPointIndices[j][0] || tetrPointIndices[k][1] == tetrPointIndices[j][1])
                        {
                            if (k == 0)
                                diamondIndices.push_back(tetrPointIndices[k][0]);
                            diamondIndices.push_back(tetrPointIndices[k][1]);
                            visitedInd.push_back(k);
                            k = j;
                            break;
                        }
                    }
                }

                for (auto &dI : diamondIndices)
                {
                    d.push_back(gridPoints[dI]);
                }

                d.push_back(gridPoints[sharedEdge.first]);
                d.push_back(gridPoints[sharedEdge.second]);

                auto edgeE = std::make_pair(sharedEdge.first, sharedEdge.second);

                if (d.size() == kTetrCount + 2)
                    diamondElement = std::make_pair(d, edgeE);
            }
            return diamondElement;
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
            p += (a2 * diamondPoints[k] + a3 * diamondPoints[k + 1]);

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