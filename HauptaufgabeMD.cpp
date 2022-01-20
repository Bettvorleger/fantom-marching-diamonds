
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

    public:
        struct Options : public VisAlgorithm::Options
        {
            Options(fantom::Options::Control &control)
                : VisAlgorithm::Options(control)
            {
                add<Field<2, Scalar>>("Field", "A 2D scalar field", definedOn<Grid<2>>(Grid<2>::Points));
                add<double>("Isovalue", "The desired isovalue to show.", 1);
                add<Color>("Color", "The color of the graphics.", Color(1.0, 0.0, 0.0));
                add<int>("Index", "Index of Diamond.", 0);
            }
        };

        struct VisOutputs : public VisAlgorithm::VisOutputs
        {
            VisOutputs(fantom::VisOutputs::Control &control)
                : VisAlgorithm::VisOutputs(control)
            {
                addGraphics("Diamonds");
                addGraphics("Iso-segments");
            }
        };

        MarchingDiamondsAlgorithm(InitData &data)
            : VisAlgorithm(data)
        {
        }

        virtual void execute(const Algorithm::Options &options, const volatile bool &abortFlag) override
        {
            std::shared_ptr<const Field<2, Scalar>> field = options.get<Field<2, Scalar>>("Field");
            std::shared_ptr<const Function<Scalar>> function = options.get<Function<Scalar>>("Field");
            double isovalue = options.get<double>("Isovalue");
            Color color = options.get<Color>("Color");
            int indexChoice = options.get<int>("Index");

            //check if input field is set
            if (!field)
            {
                debugLog() << "Input field not set." << std::endl;
                return;
            }

            std::shared_ptr<const Grid<2>> grid = std::dynamic_pointer_cast<const Grid<2>>(function->domain());

            //check if grid can be used or is correct
            if (!grid)
            {
                throw std::logic_error("Wrong type of grid!");
            }

            //get control points of the grid
            const ValueArray<Point2> &points = grid->points();

            auto evaluator = field->makeEvaluator();

            std::vector<std::vector<Point2>> diamondPoints;
            std::vector<std::vector<size_t>> diamondCells;

            std::unordered_map<size_t, std::vector<Point2>> trIntersectPoints;
            
             std::vector<Point2> splitPoints;

            //find all diamonds on grid
            for (size_t i = 0; i < grid->numCells(); i++)
            {
                Cell cell = grid->cell(i);
				
                for (size_t j = 1; j < grid->numCells(); j++)
                {
                    if (i >= j)
                        continue;

                    Cell tempCell = grid->cell(j);
                    std::vector<Point2> diamP = findDiamond(cell, tempCell, points); //new diamond vector element

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
                std::vector<Point2> dP = diamondPoints[i];
                evaluator->reset(dP[0]);
                double s0 = norm(evaluator->value());
                evaluator->reset(dP[1]);
                double s1 = norm(evaluator->value());
                evaluator->reset(dP[2]);
                double s2 = norm(evaluator->value());
                evaluator->reset(dP[3]);
                double s3 = norm(evaluator->value());

                //the last two points have to be opposite of each other in relation to edge e
                std::vector<Point2> intersect = calcIntersection(isovalue, s0, s3, s1, s2);
                int validIntersects = 0;
                Point2 trIntPoint;
                std::vector<Point2> intersectPoints; // *Collect intersection points

                for (auto &i : intersect)
                {
                    if (i[0] >= 0 && i[0] <= 1)
                    {
                        validIntersects++;
                        trIntPoint = bilinearTransform(dP[0], dP[3], dP[1], dP[2], i[0]);
                        intersectPoints.push_back(trIntPoint);
                    }
                }

                if (validIntersects == 2)
                {
                    //SPLITTING-FUNCTION 
					Point2 i1, i2;
					i1 = intersectPoints[0];
					i2 = intersectPoints[1];
					// New Point v
					Point2 v;
					v = {((i1[0] + i2[0]) / 2), ((i1[1] + i2[1]) / 2)};
					// Call Splitting Function
					infoLog() << "Splitting Diamond "  << std::endl;
					splitPoints = splitting(dP[0], dP[3], dP[1], dP[2], v, isovalue, evaluator);                  
				}
                else if (validIntersects == 1)
                {
                    trIntersectPoints[diamondCells[i][0]].push_back(trIntPoint);
                    trIntersectPoints[diamondCells[i][1]].push_back(trIntPoint);
                }
            }

            std::vector<VectorF<2>> isoSegments;

            for (auto &tr : trIntersectPoints)
            {
                if (tr.second.size() == 2)
                {
                    isoSegments.insert(isoSegments.end(), {VectorF<2>(tr.second[0]), VectorF<2>(tr.second[1])});
                }
            }

            /************LOGGING***************/
            std::vector<VectorF<2>> vertGrid;
            std::vector<Point2> fD = diamondPoints[indexChoice];

            vertGrid.insert(vertGrid.end(), {VectorF<2>(fD[0]), VectorF<2>(fD[1])});
            vertGrid.insert(vertGrid.end(), {VectorF<2>(fD[0]), VectorF<2>(fD[2])});
            vertGrid.insert(vertGrid.end(), {VectorF<2>(fD[3]), VectorF<2>(fD[1])});
            vertGrid.insert(vertGrid.end(), {VectorF<2>(fD[3]), VectorF<2>(fD[2])});
            

            if (abortFlag)
            {
                return;
            }

            auto const &system = graphics::GraphicsSystem::instance();
            std::string resourcePath = PluginRegistrationService::getInstance().getResourcePath("utils/Graphics");

            std::shared_ptr<graphics::Drawable> gridLines = system.makePrimitive(
                graphics::PrimitiveConfig{graphics::RenderPrimitives::LINES}
                    .vertexBuffer("in_vertex", system.makeBuffer(vertGrid))
                    .uniform("u_lineWidth", 1.0f)
                    .uniform("u_color", color),
                system.makeProgramFromFiles(resourcePath + "shader/line/noShading/singleColor/vertex.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/fragment.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/geometry.glsl"));
            setGraphics("Diamonds", gridLines);

            std::shared_ptr<graphics::Drawable> isocontour = system.makePrimitive(
                graphics::PrimitiveConfig{graphics::RenderPrimitives::LINES}
                    .vertexBuffer("in_vertex", system.makeBuffer(isoSegments))
                    .uniform("u_lineWidth", 1.5f)
                    .uniform("u_color", color),
                system.makeProgramFromFiles(resourcePath + "shader/line/noShading/singleColor/vertex.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/fragment.glsl",
                                            resourcePath + "shader/line/noShading/singleColor/geometry.glsl"));
            setGraphics("Iso-segments", isocontour);

            evaluator->reset(fD[0]);
            double s0 = norm(evaluator->value());
            evaluator->reset(fD[1]);
            double s1 = norm(evaluator->value());
            evaluator->reset(fD[2]);
            double s2 = norm(evaluator->value());
            evaluator->reset(fD[3]);
            double s3 = norm(evaluator->value());

            std::vector<Point2> intersect = calcIntersection(isovalue, s0, s3, s1, s2);
            for (auto &i : intersect)
            {
                if (i[0] >= 0 && i[0] <= 1)
                {
                    infoLog() << "Intersection in D: " << bilinearTransform(fD[0], fD[3], fD[1], fD[2], i[0]) << std::endl;
                }
            }
        }

    private:
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
        std::vector<Point2> findDiamond(Cell cell, Cell tempCell, const ValueArray<Point2> &points)
        {
            size_t i0, i1, i2, i3;

            size_t t0 = tempCell.index(0);
            size_t t1 = tempCell.index(1);
            size_t t2 = tempCell.index(2);

            for (int p = 0; p < 3; p++)
            {
                i0 = cell.index(p);
                i1 = cell.index((p + 1) % 3);
                i2 = cell.index((p + 2) % 3);

                if (i0 == t0)
                {
                    if (i1 == t1)
                    { //i2 opposite of i3
                        i3 = t2;
                        return std::vector<Point2>({Point2(points[i2]), Point2(points[i0]), Point2(points[i1]), Point2(points[i3])});
                    }
                    else if (i2 == t1)
                    { //i1 opposite of i3
                        i3 = t2;
                        return std::vector<Point2>({Point2(points[i1]), Point2(points[i0]), Point2(points[i2]), Point2(points[i3])});
                    }
                    else if (i1 == t2)
                    { //i2 opposite of i3
                        i3 = t1;
                        return std::vector<Point2>({Point2(points[i2]), Point2(points[i0]), Point2(points[i1]), Point2(points[i3])});
                    }
                    else if (i2 == t2)
                    { //i1 opposite of i3
                        i3 = t1;
                        return std::vector<Point2>({Point2(points[i1]), Point2(points[i0]), Point2(points[i2]), Point2(points[i3])});
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
        Point2 bilinearTransformUTR(double y)
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
        
        void splitting(Point2 d0, Point2 d1, Point2 d2, Point2 d3, Point2 v, double isovalue, auto evaluator) {
        	
        	std::vector<Point2> splitPoints;
        	
        	evaluator->reset(d0);
               double s0 = norm(evaluator->value());
               evaluator->reset(d1);
               double s1 = norm(evaluator->value());
               evaluator->reset(d2);
               double s2 = norm(evaluator->value());
               evaluator->reset(d3);
               double s3 = norm(evaluator->value());
               evaluator->reset(v);
               double s4 = norm(evaluator->value());
               
		std::vector<Point2> intersect1 = calcIntersection(isovalue, s0, s3, s1, s4);
		std::vector<Point2> intersect2 = calcIntersection(isovalue, s0, s1, s4, s2);
		std::vector<Point2> intersect3 = calcIntersection(isovalue, s3, s1, s2, s4);
		std::vector<Point2> intersect4 = calcIntersection(isovalue, s3, s2, s4, s0);
		
		std::vector<std::vector<Point2>> intersect;
		intersect.push_back(intersect1);
		intersect.push_back(intersect2);
		intersect.push_back(intersect3);
		intersect.push_back(intersect4);
		
		std::vector<Point2> dp1 = {d0, d3, d1, v};
		std::vector<Point2> dp2 = {d0, d1, v, d2};
		std::vector<Point2> dp3 = {d3, d1, d2, v};
		std::vector<Point2> dp4 = {d3, d2, v, d0};
		
		std::vector<std::vector<Point2>> dp;
		dp.push_back(dp1);
		dp.push_back(dp2);
		dp.push_back(dp3);
		dp.push_back(dp4);
               
               for(int z = 0; z < 4; z++) {
               
               	std::vector<Point2> dp_tmp = dp[z];
               	std::vector<Point2> intersect_tmp = intersect[z];
               	
				int validIntersects = 0;
				Point2 trIntPoint;
				std::vector<Point2> intersectPoints; // *Collect intersection points

		        for (auto &i : intersect_tmp)
		        {
		            if (i[0] >= 0 && i[0] <= 1)
		            {
		                validIntersects++;
		                trIntPoint = bilinearTransform(dp_tmp[0], dp_tmp[1], dp_tmp[2], dp_tmp[3], i[0]);
		                intersectPoints.push_back(trIntPoint);
		            }
		        }
				if (validIntersects == 1)
                	{
						splitPoints.push_back(trIntPoint);
						splitPoints.push_back(trIntPoint);
                	}
						intersectPoints.clear();
		        else if (validIntersects == 2)
                	{
		            	//SPLITTING-FUNCTION 
						Point2 i1, i2;
						i1 = intersectPoints[0];
						i2 = intersectPoints[1];
						// New Point v
						Point2 v;
						v = {((i1[0] + i2[0]) / 2), ((i1[1] + i2[1]) / 2)};
						// Call Splitting Function
						splitting(dp_tmp[0], dp_tmp[1], dp_tmp[2], dp_tmp[3], v, isovalue, evaluator);                  
                	}

                }

        }
    };

    AlgorithmRegister<MarchingDiamondsAlgorithm> dummy("Iso-Surface/MarchingDiamonds2D",
                                                       "Show scalar values in 2D grid over certain threshold.");
} // namespace 


