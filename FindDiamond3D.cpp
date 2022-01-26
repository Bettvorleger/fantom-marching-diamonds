

        /**
         * @brief Find tetrahedrons sharing a triangle surface by comparing cell roations against one another
         * 
         * @param cell - main tetrahedron
         * @param tempCell - tetrahedron to compare with
         * @param points - point array for converting point indices to grid coords
         * @return std::vector<Point3> - diamond control points, first and last points are opposite of each other
         */
        std::vector<Point3> findDiamond(Cell cell, Cell tempCell, const ValueArray<Point3> &points)
        {
            size_t i0, i1, i2, i3, i4;

            size_t t0 = tempCell.index(0);
            size_t t1 = tempCell.index(1);
            size_t t2 = tempCell.index(2);
            size_t t3 = tempCell.index(3);

            for (int p = 0; p < 4; p++)
            {
                i0 = cell.index(p);
                i1 = cell.index((p + 1) % 4);
                i2 = cell.index((p + 2) % 4);
                i3 = cell.index((p + 3) % 4);

                if (i0 == t0)
                {
                	if (i2 == t1)
                	{
		            if (i3 == t2)
		            { //i1 opposite of i4
		                i4 = t3;
		                return std::vector<Point3>({Point3(points[i1]), Point3(points[i0]), Point3(points[i2]), Point3(points[i3]), Point3(points[i4])});
		            }
		            else if (i3 == t3)
		            { //i1 opposite of i4
		                i4 = t2;
		                return std::vector<Point3>({Point3(points[i1]), Point3(points[i0]), Point3(points[i2]), Point3(points[i3]), Point3(points[i4])});
		            }
		            else if (i1 == t2)
		            { //i3 opposite of i4
		                i4 = t3;
		                return std::vector<Point3>({Point3(points[i3]), Point3(points[i0]), Point3(points[i2]), Point3(points[i1]), Point3(points[i4])});
		            }
		            else if (i1 == t3)
		            { //i3 opposite of i4
		                i4 = t2;
		                return std::vector<Point3>({Point3(points[i3]), Point3(points[i0]), Point3(points[i2]), Point3(points[i1]), Point3(points[i4])});
		            }
	            	}
                	if (i2 == t2)
                	{
		            if (i3 == t1)
		            { //i1 opposite of i4
		                i4 = t3;
		                return std::vector<Point3>({Point3(points[i1]), Point3(points[i0]), Point3(points[i2]), Point3(points[i3]), Point3(points[i4])});
		            }
		            else if (i3 == t3)
		            { //i1 opposite of i4
		                i4 = t1;
		                return std::vector<Point3>({Point3(points[i1]), Point3(points[i0]), Point3(points[i2]), Point3(points[i3]), Point3(points[i4])});
		            }
		            else if (i1 == t1)
		            { //i3 opposite of i4
		                i4 = t3;
		                return std::vector<Point3>({Point3(points[i3]), Point3(points[i0]), Point3(points[i2]), Point3(points[i1]), Point3(points[i4])});
		            }
		            else if (i1 == t3)
		            { //i3 opposite of i4
		                i4 = t1;
		                return std::vector<Point3>({Point3(points[i3]), Point3(points[i0]), Point3(points[i2]), Point3(points[i1]), Point3(points[i4])});
		            }
	            	}
                	if (i2 == t3)
                	{
		            if (i3 == t1)
		            { //i1 opposite of i4
		                i4 = t2;
		                return std::vector<Point3>({Point3(points[i1]), Point3(points[i0]), Point3(points[i2]), Point3(points[i3]), Point3(points[i4])});
		            }
		            else if (i3 == t2)
		            { //i1 opposite of i4
		                i4 = t1;
		                return std::vector<Point3>({Point3(points[i1]), Point3(points[i0]), Point3(points[i2]), Point3(points[i3]), Point3(points[i4])});
		            }
		            else if (i1 == t1)
		            { //i3 opposite of i4
		                i4 = t2;
		                return std::vector<Point3>({Point3(points[i3]), Point3(points[i0]), Point3(points[i2]), Point3(points[i1]), Point3(points[i4])});
		            }
		            else if (i1 == t2)
		            { //i3 opposite of i4
		                i4 = t1;
		                return std::vector<Point3>({Point3(points[i3]), Point3(points[i0]), Point3(points[i2]), Point3(points[i1]), Point3(points[i4])});
		            }
	            	}
                	if (i3 == t1)
                	{
		            if (i1 == t2)
		            { //i2 opposite of i4
		                i4 = t3;
		                return std::vector<Point3>({Point3(points[i2]), Point3(points[i0]), Point3(points[i1]), Point3(points[i3]), Point3(points[i4])});
		            }
		            else if (i1 == t2)
		            { //i2 opposite of i4
		                i4 = t3;
		                return std::vector<Point3>({Point3(points[i2]), Point3(points[i0]), Point3(points[i1]), Point3(points[i3]), Point3(points[i4])});
		            }
	            	}
                	if (i3 == t2)
                	{
		            if (i1 == t1)
		            { //i2 opposite of i4
		                i4 = t3;
		                return std::vector<Point3>({Point3(points[i2]), Point3(points[i0]), Point3(points[i1]), Point3(points[i3]), Point3(points[i4])});
		            }
		            else if (i1 == t3)
		            { //i2 opposite of i4
		                i4 = t1;
		                return std::vector<Point3>({Point3(points[i2]), Point3(points[i0]), Point3(points[i1]), Point3(points[i3]), Point3(points[i4])});
		            }
	            	}
                	if (i3 == t3)
                	{
		            if (i1 == t1)
		            { //i2 opposite of i4
		                i4 = t2;
		                return std::vector<Point3>({Point3(points[i2]), Point3(points[i0]), Point3(points[i1]), Point3(points[i3]), Point3(points[i4])});
		            }
		            else if (i1 == t2)
		            { //i2 opposite of i4
		                i4 = t1;
		                return std::vector<Point3>({Point3(points[i2]), Point3(points[i0]), Point3(points[i1]), Point3(points[i3]), Point3(points[i4])});
		            }
	            	}	            		            		            		            		            	
		           
                }
            }
            return std::vector<Point3>();
        }

               



