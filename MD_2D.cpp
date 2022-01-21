#include <fantom/algorithm.hpp>
#include <fantom/datastructures/interfaces/Field.hpp>
#include <fantom/graphics.hpp>
#include <fantom/register.hpp>
#include <math.h>

#include <fantom-plugins/utils/Graphics/HelperFunctions.hpp>
#include <fantom-plugins/utils/Graphics/ObjectRenderer.hpp>

#include <stdexcept>
#include <vector>
#include <unordered_map>

using namespace fantom;

namespace
{

    class MarchingDiamondsAlgorithm : public VisAlgorithm
    {

    protected:
        std::shared_ptr<const Field<2, Scalar>> field;
        std::shared_ptr<const Grid<2>> grid;

        //aggregates intersection for each triangle by index
        std::unordered_map<size_t, std::vector<Point2>> trIntersectPoints;

    public:
        struct Options : public VisAlgorithm::Options
        {
            Options(fantom::Options::Control &control)
                : VisAlgorithm::Options(control)
            {
                add<Field<2, Scalar>>("Field", "A 2D scalar field", definedOn<Grid<2>>(Grid<2>::Points));
                add<double>("Isovalue", "The desired isovalue to show.", 1.05);
                add<Color>("Color", "The color of the graphics.", Color(1.0, 0.0, 0.0));
                add<int>("Max Splitting", "Max amount of splitting steps.", 2, &acceptNumber);
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
                //addGraphics("Diamonds");
                addGraphics("Iso-segments");
            }
        };

        MarchingDiamondsAlgorithm(InitData &data)
            : VisAlgorithm(data)
        {
        }

        virtual void execute(const Algorithm::Options &options, const volatile bool &abortFlag) override
        {
            field = options.get<Field<2, Scalar>>("Field");
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

            grid = std::dynamic_pointer_cast<const Grid<2>>(function->domain());

            //check if grid can be used or is correct
            if (!grid)
            {
                throw std::logic_error("Wrong type of grid!");
            }

            //get control points of the grid
            const ValueArray<Point2> &points = grid->points();

            //control points of grid cells and initial cell count
            std::unordered_map<size_t, std::vector<Point2>> trianglePoints;
            size_t trCellCount = grid->numCells();

            //points for final iso segment drawing
            std::vector<VectorF<2>> isoSegments;

            //load all cells into vector
            for (size_t i = 0; i < trCellCount; i++)
            {
                Cell cell = grid->cell(i);
                trianglePoints[i] = std::vector<Point2>({points[cell.index(0)], points[cell.index(1)], points[cell.index(2)]});
            }

            marchingDiamonds(trianglePoints, trCellCount, isovalue, maxSplit);

            //add triangle intersection points to final iso segment vector
            for (auto &tr : trIntersectPoints)
            {
                if (tr.second.size() == 2)
                {
                    isoSegments.insert(isoSegments.end(), {VectorF<2>(tr.second[0]), VectorF<2>(tr.second[1])});
                }
            }

            /************LOGGING***************/

            /*std::vector<VectorF<2>> vertGrid;
            std::vector<Point2> fD = diamondPoints[indexChoice];

            vertGrid.insert(vertGrid.end(), {VectorF<2>(fD[0]), VectorF<2>(fD[1])});
            vertGrid.insert(vertGrid.end(), {VectorF<2>(fD[0]), VectorF<2>(fD[2])});
            vertGrid.insert(vertGrid.end(), {VectorF<2>(fD[3]), VectorF<2>(fD[1])});
            vertGrid.insert(vertGrid.end(), {VectorF<2>(fD[3]), VectorF<2>(fD[2])});
*/
            if (abortFlag)
            {
                return;
            }

            auto const &system = graphics::GraphicsSystem::instance();
            std::string resourcePath = PluginRegistrationService::getInstance().getResourcePath("utils/Graphics");
            /*
            std::shared_ptr<graphics::Drawable> gridLines = system.makePrimitive(
                graphics::PrimitiveConfig{graphics::RenderPrimitives::LINES}
                    .vertexBuffer("in_vertex", system.makeBuffer(vertGrid))
                    .uniform("u_lineWidth", 1.0f)
                    .uniform("u_color", color),
                system.makeProgramFromFiles(resourcePath + "shader/line/noShading/singleColor/vertex.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/fragment.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/geometry.glsl"));
            setGraphics("Diamonds", gridLines);
*/
            std::shared_ptr<graphics::Drawable> isocontour = system.makePrimitive(
                graphics::PrimitiveConfig{graphics::RenderPrimitives::LINES}
                    .vertexBuffer("in_vertex", system.makeBuffer(isoSegments))
                    .uniform("u_lineWidth", 1.5f)
                    .uniform("u_color", color),
                system.makeProgramFromFiles(resourcePath + "shader/line/noShading/singleColor/vertex.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/fragment.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/geometry.glsl"));
            setGraphics("Iso-segments", isocontour);
        }

        bool marchingDiamonds(std::unordered_map<size_t, std::vector<Point2>> trianglePoints, size_t trCellCount, double isovalue, int splitCount)
        {
            if (splitCount == 0)
            {
                return true;
            }
            else
            {
                trIntersectPoints.clear();
            }

            auto evaluator = field->makeEvaluator();

            //check if splitting was neccesary
            bool splitting = false;

            //vector of point and indices vector for diamonds
            std::vector<std::vector<Point2>> diamondPoints;
            std::vector<std::vector<size_t>> diamondCells;

            //find all diamonds on grid
            for (size_t i = 0; i < trCellCount; i++)
            {
                if (trianglePoints[i].empty())
                    continue;
                for (size_t j = 1; j < trCellCount; j++)
                {
                    if (i >= j || trianglePoints[j].empty())
                        continue;

                    std::vector<Point2> diamP = findDiamond(trianglePoints[i], trianglePoints[j]); //new diamond vector element

                    if (!diamP.empty())
                    {
                        diamondPoints.push_back(diamP);
                        diamondCells.push_back({i, j});
                    }
                }
            }

            //find all intersections on grid
            for (size_t i = 0; i < diamondPoints.size(); i++)
            {
                std::vector<Point2> trIntPoints = findIntersections(diamondPoints[i], isovalue, evaluator);

                if (trIntPoints.size() == 2)
                {
                    splitting = true;

                    Point2 newVertex = (trIntPoints[0] + trIntPoints[1]) / 2;
                    std::vector<Point2> tempDiamond1 = {diamondPoints[i][0], diamondPoints[i][1], newVertex, diamondPoints[i][3]};
                    std::vector<Point2> tempDiamond2 = {diamondPoints[i][0], diamondPoints[i][2], newVertex, diamondPoints[i][3]};

                    std::vector<Point2> newTriangle1 = {diamondPoints[i][0], diamondPoints[i][1], newVertex};
                    std::vector<Point2> newTriangle2 = {diamondPoints[i][0], diamondPoints[i][2], newVertex};
                    std::vector<Point2> newTriangle3 = {diamondPoints[i][3], diamondPoints[i][1], newVertex};
                    std::vector<Point2> newTriangle4 = {diamondPoints[i][3], diamondPoints[i][2], newVertex};

                    trianglePoints.erase(diamondCells[i][0]);
                    trianglePoints.erase(diamondCells[i][1]);

                    trianglePoints[trCellCount + 1] = newTriangle1;
                    trianglePoints[trCellCount + 2] = newTriangle2;
                    trianglePoints[trCellCount + 3] = newTriangle3;
                    trianglePoints[trCellCount + 4] = newTriangle4;

                    trCellCount += 4;
                }
                if (trIntPoints.size() == 1)
                {
                    trIntersectPoints[diamondCells[i][0]].push_back(trIntPoints[0]);
                    trIntersectPoints[diamondCells[i][1]].push_back(trIntPoints[0]);
                }
            }

            if (splitting)
            {
                splitCount--;
                return marchingDiamonds(trianglePoints, trCellCount, isovalue, splitCount);
            }
            return splitting;
        }

    private:
        template <typename T>
        std::vector<Point2> findIntersections(std::vector<Point2> diamondPoint, double isovalue, T &evaluator)
        {
            std::vector<Point2> dP = diamondPoint;
            evaluator->reset(dP[0]);
            double s0 = norm(evaluator->value());
            evaluator->reset(dP[1]);
            double s1 = norm(evaluator->value());
            evaluator->reset(dP[2]);
            double s2 = norm(evaluator->value());
            evaluator->reset(dP[3]);
            double s3 = norm(evaluator->value());

            //the first two points have to be opposite of each other in relation to edge e
            std::vector<Point2> intersect = calcIntersection(isovalue, s0, s3, s1, s2);
            std::vector<Point2> trIntPoints;

            for (auto &i : intersect)
            {
                if (i[0] >= 0 && i[0] <= 1)
                {
                    //mapping has to be the same, as with the intersection calc
                    trIntPoints.push_back(bilinearTransform(dP[0], dP[3], dP[1], dP[2], i[0]));
                }
            }
            if (trIntPoints.size() == 2)
            {
                Point2 newVertex = (trIntPoints[0] + trIntPoints[1]) / 2;
                std::vector<Point2> tempDiamond1 = {dP[0], dP[1], newVertex, dP[3]};
                std::vector<Point2> tempDiamond2 = {dP[0], dP[2], newVertex, dP[3]};
                //return findIntersections(tempDiamond1, isovalue, evaluator);
            }

            return trIntPoints;
        }

        /*
        void splitDiamond(std::vector<Point2> diamondPoint, Point2 newVertex, double isovalue, T &evaluator)
        {
            std::vector<Point2> tempDiamond1 = {diamondPoints[i][0], diamondPoints[i][1], newVertex, diamondPoints[i][3]};
            std::vector<Point2> tempDiamond2 = {diamondPoints[i][0], diamondPoints[i][2], newVertex, diamondPoints[i][3]};
            Point2 newVertex;

            for (int i = 0; i < 3; i++)
            {
                std::vector<Point2> trIntPoints = findIntersections(tempDiamond1, isovalue, evaluator);

                if (trIntPoints.size() == 2)
                {
                    newVertex = (trIntPoints[0] + trIntPoints[1]) / 2;
                }
                else if (trIntPoints.size() == 1)
                {
                    trIntersectPoints[diamondCells[i][0]].push_back(trIntPoints[0]);
                    trIntersectPoints[diamondCells[i][1]].push_back(trIntPoints[0]);
                    break;
                }
            }
        }
        */

        /**
         * @brief calc intersections of edge e with bilinear plane (hyperbola for fixed iso-value c) in mapped unit square

                    s3                   
                   /  \                 s1 -> (0,1)         s1----s2
                  /    \                s2 -> (1,1)         |    / |
                s1     s0      --->     s0 -> (1,0)   --->  |   e  |
                  \    /                s3 -> (0,0)         | /    |
                   \  /                                     s3----s0
                    s2                  
        
         * @param c - desired iso-value
         * @param s0 - iso-value at point r0
         * @param s1 - iso-value at point r1
         * @param s2 - iso-value at point r2
         * @param s3 - iso-value at point r3
         * @return std::vector<Point2> - vector of intersection points
         */
        std::vector<Point2> calcIntersection(double c, double s0, double s1, double s2, double s3)
        {
            std::vector<Point2> intersections;
            double denominator = 2 * (s0 + s1 - s2 + s3);
            double result;

            //check if denominator is too close/equals zero, return empty as error
            if (abs(denominator) < 0.001)
            {
                return intersections;
            }

            double w = 4 * (s3 - c) * (s0 + s1 - s2 + s3) + ((s0 + s1) * (s0 + s1));
            //check if root is real, return empty as error (only real values)
            if (w < 0)
            {
                return intersections;
            }

            w = sqrt(w);
            if (w == 0)
            {
                result = (s0 + s1) / denominator;
                intersections.push_back({Point2({result, result})});
            }
            else
            {
                result = (s0 + s1 + w) / denominator;
                intersections.push_back({Point2({result, result})});

                result = (s0 + s1 - w) / denominator;
                intersections.push_back({Point2({result, result})});
            }

            return intersections;
        }

        /**
         * @brief Find triangles sharing an edge by comparing cell roations against one another
         * 
         * @param cell - main triangle to find shared edge triangle for
         * @param tempCell - triangle to compare with
         * @param points - point array for converting point indices to grid coords
         * @return std::vector<Point2> - diamond control points, first and last points are opposite of each other
         */
        std::vector<Point2> findDiamond(std::vector<Point2> trPoints, std::vector<Point2> trTempPoints)
        {
            Point2 d0, d1, d2, d3;

            Point2 t0 = trTempPoints[0];
            Point2 t1 = trTempPoints[1];
            Point2 t2 = trTempPoints[2];

            for (int p = 0; p < 3; p++)
            {
                d0 = trPoints[p];
                d1 = trPoints[(p + 1) % 3];
                d2 = trPoints[(p + 2) % 3];

                if (d0 == t0)
                {
                    if (d1 == t1)
                    { //i2 opposite of i3
                        d3 = t2;
                        return std::vector<Point2>({d2, d0, d1, d3});
                    }
                    else if (d2 == t1)
                    { //i1 opposite of i3
                        d3 = t2;
                        return std::vector<Point2>({d1, d0, d2, d3});
                    }
                    else if (d1 == t2)
                    { //i2 opposite of i3
                        d3 = t1;
                        return std::vector<Point2>({d2, d0, d1, d3});
                    }
                    else if (d2 == t2)
                    { //i1 opposite of i3
                        d3 = t1;
                        return std::vector<Point2>({d1, d0, d2, d3});
                    }
                }
            }
            return std::vector<Point2>();
        }

        /**
         * @brief maps coordinats inside the unit square to the reference diamond R via bilinear transformation
         * ONLY for specific intersection condition where x = y
         * 
         * @param y - coord of intersection of edge with hyperbola in unit square
         * @return Point2 - intersection point in reference space of dimaond R
         */
        Point2
        bilinearTransformUTR(double y)
        {
            Point2 v1 = {0, 2};
            Point2 v2 = {0, -2};

            return v1 + (v2 * y * y);
        }

        /**
         * @brief maps coordinats inside the unit square to diamond D via bilinear transformation
         * ONLY for specific intersection condition where x = y
         * 
                    d3                   
                   /  \           
                  /    \         
                d1      d0      
                  \    /      
                   \  /        
                    d2   
         * 
         * @param d0 - diamond control point in grid space
         * @param d1 - diamond control point in grid space
         * @param d2 - diamond control point in grid space
         * @param d3 - coord of intersection of edge with hyperbola in unit square
         * @param y - coord of intersection of edge with hyperbola in unit square
         * @return Point2 - intersection point in normal grid space of diamond D
         */
        Point2 bilinearTransform(Point2 d0, Point2 d1, Point2 d2, Point2 d3, double y)
        {
            return d3 + (d0 + d1 - 2 * d3) * y + (d3 - d0 + d2 - d1) * y * y;
        }
    };

    AlgorithmRegister<MarchingDiamondsAlgorithm> dummy("Iso-Surface/MarchingDiamonds2D",
                                                       "Show scalar values in 2D grid over certain threshold.");
} // namespace
