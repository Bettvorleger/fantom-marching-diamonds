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

using namespace fantom;

namespace
{

    class MarchingDiamondsAlgorithm : public VisAlgorithm
    {

    protected:
        std::shared_ptr<const Field<2, Scalar>> field;
        std::shared_ptr<const Grid<2>> grid;

    public:
        struct Options : public VisAlgorithm::Options
        {
            Options(fantom::Options::Control &control)
                : VisAlgorithm::Options(control)
            {
                add<Field<2, Scalar>>("Field", "A 2D scalar field", definedOn<Grid<2>>(Grid<2>::Points));
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
            std::map<size_t, std::vector<Point2>> trianglePoints;
            size_t trCellCount = grid->numCells();

            //vector of point and indices vector for diamonds
            std::map<std::vector<size_t>, std::vector<Point2>> diamondPoints;

            //aggregates intersection for each triangle by index
            std::unordered_map<size_t, std::vector<Point2>> trIntersectPoints;

            //points for final iso segment drawing
            std::vector<VectorF<2>> isoSegments;

            //load all cells into vector
            for (size_t i = 0; i < trCellCount; i++)
            {
                Cell cell = grid->cell(i);
                trianglePoints[i] = std::vector<Point2>({points[cell.index(0)], points[cell.index(1)], points[cell.index(2)]});
            }

            trIntersectPoints = marchingDiamonds(trianglePoints, trIntersectPoints, diamondPoints, trianglePoints, trCellCount, isovalue, maxSplit + 1);

            //add triangle intersection points to final iso segment vector
            for (auto &tr : trIntersectPoints)
            {
                if (tr.second.size() >= 2)
                {
                    isoSegments.insert(isoSegments.end(), {VectorF<2>(tr.second[0]), VectorF<2>(tr.second[1])});
                }
            }

            if (abortFlag)
            {
                return;
            }

            auto const &system = graphics::GraphicsSystem::instance();
            std::string resourcePath = PluginRegistrationService::getInstance().getResourcePath("utils/Graphics");

            std::shared_ptr<graphics::Drawable> isocontour = system.makePrimitive(
                graphics::PrimitiveConfig{graphics::RenderPrimitives::LINES}
                    .vertexBuffer("in_vertex", system.makeBuffer(isoSegments))
                    .uniform("u_lineWidth", 2.0f)
                    .uniform("u_color", color),
                system.makeProgramFromFiles(resourcePath + "shader/line/noShading/singleColor/vertex.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/fragment.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/geometry.glsl"));
            setGraphics("Iso-segments", isocontour);
        }

        /**
         * @brief Search for diamonds and find intersections inside them. If splitting is neccessary the function is called recursively with new triangles.
         * 
         * @param trianglePoints - all triangle cells as point vectors
         * @param trIntersectPoints - map for aggregating intersection points for segment connection
         * @param diamondPoints - point vector for all found diamonds, first call empty, then reused for splitting
         * @param newTriangles - new triangles found through splitting
         * @param trCellCount - triangle cell count changed through splitting triangles
         * @param isovalue - isovalue to be displayed
         * @param splitCount - recursion depth controling splitting steps
         * @return std::unordered_map<size_t, std::vector<Point2>> - map for aggregating intersection points for segment connection
         */
        template <typename T, typename U, typename V, typename W>
        std::unordered_map<size_t, std::vector<Point2>> marchingDiamonds(T trianglePoints, U trIntersectPoints, V diamondPoints, W newTriangles, size_t trCellCount, double isovalue, int splitCount)
        {
            if (splitCount == 0)
            {
                return trIntersectPoints;
            }
            trIntersectPoints.clear();

            auto evaluator = field->makeEvaluator();

            //check if splitting was neccesary
            bool splitting = false;

            //aggregate triangle points used in diamond over triangle indices
            std::map<size_t, std::vector<Point2>> trIndices;

            //find all diamonds on grid
            for (auto &i : trianglePoints)
            {
                for (auto &j : newTriangles)
                {
                    //avoid duplicate diamonds by only using the "upper" part of the symmetric matrix on first call
                    if (i.first >= j.first)
                        continue;
                    std::vector<Point2> diamP = findDiamond(i.second, j.second); //new diamond vector element

                    if (!diamP.empty())
                    {
                        //add new diamond point to map by triangle indices used
                        diamondPoints[{i.first, j.first}] = diamP;

                        //add for edge cases
                        trIndices[i.first].push_back(diamP[1]);
                        trIndices[i.first].push_back(diamP[2]);
                        trIndices[j.first].push_back(diamP[1]);
                        trIndices[j.first].push_back(diamP[2]);
                    }
                }
            }

            //add new triangle points now for later splitting
            trianglePoints = newTriangles;
            //reset for new triangles to be added from splitting
            newTriangles.clear();

            //find all intersections on grid-edge triangles without diamonds on all three edges
            for (auto &trI : trIndices)
            {
                Point2 intersectPoint = findEdgeIntersection(trI.second, isovalue, evaluator);
                if (intersectPoint[0] != 0 && intersectPoint[0] != 0)
                    trIntersectPoints[trI.first].push_back(intersectPoint);
            }

            std::map<std::vector<size_t>, std::vector<Point2>> newDiamondPoints;

            //find all intersections on grid, using constant iterator for manipulating the original structure
            for (auto dP = diamondPoints.cbegin(); dP != diamondPoints.cend();)
            {
                std::vector<Point2> trIntPoints = findIntersections(dP->second, isovalue, evaluator);

                //splitting for more than one intersection found
                if (trIntPoints.size() == 2)
                {
                    splitting = true;

                    //new point halfway between the two intersections
                    Point2 newVertex = (trIntPoints[0] + trIntPoints[1]) / 2;

                    //new triangles and diamonds created through splitting with new vertex, added by expanding cell indices
                    trianglePoints[trCellCount] = {dP->second[0], dP->second[1], newVertex};
                    trianglePoints[trCellCount + 1] = {dP->second[0], dP->second[2], newVertex};
                    trianglePoints[trCellCount + 2] = {dP->second[3], dP->second[1], newVertex};
                    trianglePoints[trCellCount + 3] = {dP->second[3], dP->second[2], newVertex};

                    newDiamondPoints[{trCellCount, trCellCount + 2}] = {dP->second[0], dP->second[1], newVertex, dP->second[3]};     //triangle 0 + 2
                    newDiamondPoints[{trCellCount + 1, trCellCount + 3}] = {dP->second[0], dP->second[2], newVertex, dP->second[3]}; //triangle 1 + 3
                    newDiamondPoints[{trCellCount, trCellCount + 1}] = {dP->second[1], dP->second[0], newVertex, dP->second[2]};     //triangle 0 + 1
                    newDiamondPoints[{trCellCount + 2, trCellCount + 3}] = {dP->second[1], dP->second[3], newVertex, dP->second[2]}; //triangle 2 + 3

                    trCellCount += 4;

                    //erase the two original triangles and found diamond
                    trianglePoints.erase(dP->first[0]);
                    trianglePoints.erase(dP->first[1]);
                    dP = diamondPoints.erase(dP);
                    continue;
                }
                else if (trIntPoints.size() == 1)
                {
                    //add the one intersection two both triangles
                    trIntersectPoints[dP->first[0]].push_back(trIntPoints[0]);
                    trIntersectPoints[dP->first[1]].push_back(trIntPoints[0]);
                }
                ++dP;
            }

            //merge old with new triangles because we want to check all non-splitted triangles against all non-splitted PLUS the splitted triangles
            newTriangles.insert(trianglePoints.begin(), trianglePoints.end());

            if (splitting)
            {
                return marchingDiamonds(trianglePoints, trIntersectPoints, newDiamondPoints, newTriangles, trCellCount, isovalue, splitCount - 1);
            }
            return trIntersectPoints;
        }

    private:
        /**
         * @brief Handle intersection search for given diamond point, isovalue and evaluator
         * 
         * @param diamondPoint
         * @param isovalue 
         * @param evaluator
         * @return std::vector<Point2> - vector of intersection points, two at most
         */
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
            std::vector<Point2> intersect = calcIntersection(isovalue, s0, s3, s2, s1);
            std::vector<Point2> trIntPoints;

            for (auto &i : intersect)
            {
                if (i[0] >= 0 && i[0] <= 1)
                {
                    //mapping has to be the same, as with the intersection calc
                    trIntPoints.push_back(bilinearTransform(dP[0], dP[3], dP[2], dP[1], i[0]));
                }
            }

            return trIntPoints;
        }

        /**
         * @brief calc intersections of edge e with bilinear plane (hyperbola for fixed iso-value c) in mapped unit square

                    s3                   
                   /  \                 s1 -> (0,1)         s1----s2
                  /    \                s2 -> (1,1)         |    / |
                s1     s0      --->     s0 -> (1,0)   --->  |   e  |
                  \    /                s3 -> (0,0)         | /    |
                   \  /                                     s3----s0
                    s2                  
        
            solve for y in c = (s0+s1)(y-y^2)+s3(1-y)^2+s2*y^2

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
            double denominator = 2 * (s0 + s1 - s2 - s3);
            double result;

            //check if denominator is too close/equals zero, return empty as error
            if (abs(denominator) < 0.000001)
            {
                return intersections;
            }

            double w = pow((s0 + s1), 2) - 4 * c * (s0 + s1 - s2 - s3) - (4 * s2 * s3);
            //check if root is real, return empty as error (only real values)
            if (w < 0)
            {
                return intersections;
            }

            w = sqrt(w);
            if (w == 0)
            {
                result = (s0 + s1 - 2 * s3) / denominator;
                intersections.push_back({Point2({result, result})});
            }
            else
            {
                result = (s0 + s1 - 2 * s3 + w) / denominator;
                intersections.push_back({Point2({result, result})});

                result = (s0 + s1 - 2 * s3 - w) / denominator;
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
            Point2 t0, t1, t2;

            for (int t = 0; t < 3; t++)
            {
                t0 = trTempPoints[t];
                t1 = trTempPoints[(t + 1) % 3];
                t2 = trTempPoints[(t + 2) % 3];

                for (int p = 0; p < 3; p++)
                {
                    d0 = trPoints[p];
                    d1 = trPoints[(p + 1) % 3];
                    d2 = trPoints[(p + 2) % 3];

                    if (d0 == t0)
                    {
                        if (d1 == t1)
                        { //d2 opposite of d3
                            d3 = t2;
                            return std::vector<Point2>({d2, d0, d1, d3});
                        }
                        else if (d2 == t1)
                        { //d1 opposite of d3
                            d3 = t2;
                            return std::vector<Point2>({d1, d0, d2, d3});
                        }
                        else if (d1 == t2)
                        { //d2 opposite of d3
                            d3 = t1;
                            return std::vector<Point2>({d2, d0, d1, d3});
                        }
                        else if (d2 == t2)
                        { //d1 opposite of d3
                            d3 = t1;
                            return std::vector<Point2>({d1, d0, d2, d3});
                        }
                    }
                }
            }
            return std::vector<Point2>();
        }

        /**
         * @brief Find intersection for edge cells via the fact, that edge-triangles have at most 2 diamonds, ergo points are only used at most twice
         * 
         * @param trPoints - trinagle point vector, one entry for each usage in a diamond
         * @param isovalue 
         * @param evaluator 
         * @return Point2 - point of intersection
         */
        template <typename T>
        Point2 findEdgeIntersection(std::vector<Point2> trPoints, double isovalue, T &evaluator)
        {
            Point2 p0, p1;
            if (trPoints.size() == 4)
            {
                if (trPoints[0] == trPoints[2])
                {
                    p0 = trPoints[1];
                    p1 = trPoints[3];
                }
                else if (trPoints[0] == trPoints[3])
                {
                    p0 = trPoints[1];
                    p1 = trPoints[2];
                }
                else if (trPoints[1] == trPoints[2])
                {
                    p0 = trPoints[0];
                    p1 = trPoints[3];
                }
                else if (trPoints[1] == trPoints[3])
                {
                    p0 = trPoints[0];
                    p1 = trPoints[2];
                }
            }
            return linearInterpolation(p0, p1, isovalue, evaluator);
        }

        template <typename T>
        Point2 linearInterpolation(Point2 p0, Point2 p1, double isovalue, T &evaluator)
        {
            Point2 intersectPoint;

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
            return d3 + (d0 + d1 - 2 * d3) * y + (d3 - d0 + d2 - d1) * pow(y, 2);
        }
    };

    AlgorithmRegister<MarchingDiamondsAlgorithm> dummy("Iso-Surface/MarchingDiamonds2D",
                                                       "Show iso-segment in 2D grid for given isovalue.");
} // namespace