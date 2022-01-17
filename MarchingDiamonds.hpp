#include<iostream>
#include <fantom/algorithm.hpp>
#include <fantom/datastructures/interfaces/Field.hpp>
#include <fantom/register.hpp>
#include <fantom/math.hpp>
#include <fantom/datastructures/Function.hpp>
#include <fantom/dataset.hpp>


#include <stdexcept>
#include <vector>

namespace fantom {

	class MarchingDiamonds {

		public:
		
			Point2 marchingDiamonds2D (double isovalue, std::shared_ptr< const Field< 2, Scalar > > field, const Function< Scalar > > function) 
			{
				
				std::shared_ptr< const Grid< 2 > > grid = std::dynamic_pointer_cast< const Grid< 2 > >( function->domain() );
					
				// DiscreteFunctionEvaluator
				auto discreteEvaluator = function->makeDiscreteEvaluator();
								
				// Iterate through each triangle cell	
				for( size_t i = 0; i < grid->numCells(); ++i ) {
					Cell cell = grid->cell( i );
					
					// The four points of a diamond, where d2 and d3 are the start and endpoint of edge e
					Point2 d0, d1, d2, d3;
					
					// Save intersection points from each edge of a triangle in this vector to connect them later.
					std::vector< Point< 2 > > connections;
					
					// Process each of the three edges in the triangle successively
					for(int count = 0; count < 3; ++count) { 
		 				int a, b, c;
		 				if(count = 0) a = 1; b = 2; c = 3;
		 				if(count = 1) a = 1; b = 3; c = 2;
		 				if(count = 2) a = 2; b = 3; c= 1;

						d2 = points[ cell.index(a) ];
						d3 = points[ cell.index(b) ];
						d0 = points[ cell.index(c) ];
						
						
						// Search for a matching triangle to complete the diamond
						for(size_t j = 0; j < grid->numCells(); ++j ) {
							if(i==j) continue;
							Cell cell2 = grid->cell( j );
							Point2 point1 = points[ cell2.index(0) ];
							Point2 point2 = points[ cell2.index(1) ];
							Point2 point3 = points[ cell2.index(2) ];
							

							if(d2 == point1 && d3 == point2) d1 = point3;
								else if(d2 == point1 && d3 == point3) d1 = point2;
								else if(d2 == point2 && d3 == point1) d1 = point3;
								else if(d2 == point2 && d3 == point3) d1 = point1;
								else if(d2 == point3 && d3 == point1) d1 = point2;
								else if(d2 == point3 && d3 == point2) d1 = point1;

							else continue;
							
							// Points of the reference diamond R, where r2 and r3 are the start and endpoint of edge e
							Point2 r0, r1, r2, r3;
							r0 = {1,1}; 
							r1 = {-1,1}; 
							r2 = {0,0}; 
							r3 = {0,2};
							
							// Get scalar values at r0,r1,r2,r3   
 							auto r0_eval = discreteEvaluator->value(r0);
 							double r0_value = norm(r0_eval);
 							auto r1_eval = discreteEvaluator->value(r1);
 							double r1_value = norm(r1_eval);
 							auto r2_eval = discreteEvaluator->value(r2);
 							double r2_value = norm(r2_eval);
 							auto r3_eval = discreteEvaluator->value(r3);
 							double r3_value = norm(r3_eval);
 							
 							//Bilinear interpolation
 							// f(α,β)==(1−α)(1−β)*r2_value+α(1−β)f10+(1−α)βf01+αβ*r0_value
							// Coordinates of intersections of the hyperbola
							
							// Calculate intersections with e in the unit square
							double a = sqrt(4*(r3_value - isovalue)*(r0_value + r1_value - r2_value + r3_value) + (pow(r0_value + r1_value, 2)))
							double b = 2* (r0_value + r1_value - r2_value + r3_value)
							
							double y1,y2;
							if(a >= 0) {
								if(b == 0) {
									y1 = (r0_value + r1_value + a) /b
									y2 = (r0_value + r1_value - a)/b
								}
								else {
									y1 = (r0_value + r1_value + a)
									y2 = (r0_value + r1_value - a)
								}
							}
								
							// Map back to R
							Point2D x1, x2; // intersection points in R
							
							x1 = {0,2} + {0,-2} * pow(y1,2)
							x2= {0,2} + {0,-2} * pow(y2,2)
							
							// Back to D
							// intersection Z = w*d0 + w*d1 + w*d2+ w*d3
							std::vector< Point< 2 > > intersections;
							
							// calculate intersection points
							
							if (intersections.size() == 0) {
								continue;
							}
							if (intersections.size() == 1) {
								connections.push_back(intersections.index(0));
							}
							if (intersections.size() == 2) {
								//splitting
								Point2D i1, i2;
								i1 = intersections.index(0);
								i2 = intersections.index(1);
								// New Point v
    								Point2D v;
    								v.x = (i1.x + i2.x) / 2;
    								v.y = (i1.y + i2.y) / 2;								
							}
							

 														

						
						}
					}

	 				
				}			
	
				return;
			}		

	};
	


} 



