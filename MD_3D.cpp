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
        std::shared_ptr<const Field<3, Scalar>> field;
        std::shared_ptr<const Grid<3>> grid;

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
                addGraphics("Iso-segments");
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

            //get control points of the grid
            const ValueArray<Point3> &points = grid->points();

            //control points of grid cells and initial cell count
            std::map<size_t, std::vector<Point3>> tetrahedraPoints;
            size_t trCellCount = grid->numCells();

            //vector of point and indices vector for diamonds
            std::map<std::vector<size_t>, std::vector<Point3>> diamondPoints;

            //aggregates intersection for each triangle by index
            std::unordered_map<size_t, std::vector<Point3>> tetrIntersectPoints;

            //points for final iso segment drawing
            std::vector<VectorF<3>> isoSurfaces;

            //load all cells into vector
            for (size_t i = 0; i < trCellCount; i++)
            {
                Cell cell = grid->cell(i);
                tetrahedraPoints[i] = std::vector<Point3>({points[cell.index(0)], points[cell.index(1)], points[cell.index(2)], points[cell.index(3)]});
            }

            tetrIntersectPoints = marchingDiamonds(tetrahedraPoints, tetrIntersectPoints, diamondPoints, tetrahedraPoints, trCellCount, isovalue, maxSplit + 1);

            //add triangle intersection points to final iso segment vector
            for (auto &tr : tetrIntersectPoints)
            {
                if (tr.second.size() >= 3)
                {
                    isoSurfaces.insert(isoSurfaces.end(), {VectorF<3>(tr.second[0]), VectorF<3>(tr.second[1]), VectorF<3>(tr.second[3])});
                }
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
        }

        /**
         * @brief Search for diamonds and find intersections inside them. If splitting is neccessary the function is called recursively with new Tetrahedra.
         * 
         * @param trianglePoints - all triangle cells as point vectors
         * @param trIntersectPoints - map for aggregating intersection points for segment connection
         * @param diamondPoints - point vector for all found diamonds, first call empty, then reused for splitting
         * @param newTetrahedra - new Tetrahedra found through splitting
         * @param trCellCount - triangle cell count changed through splitting Tetrahedra
         * @param isovalue - isovalue to be displayed
         * @param splitCount - recursion depth controling splitting steps
         * @return std::unordered_map<size_t, std::vector<Point3>> - map for aggregating intersection points for segment connection
         */
        template <typename T, typename U, typename V, typename W>
        std::unordered_map<size_t, std::vector<Point3>> marchingDiamonds(T tetrahedraPoints, U tetrIntersectPoints, V diamondPoints, W newTetrahedra, size_t trCellCount, double isovalue, int splitCount)
        {
            if (splitCount == 0)
            {
                return tetrIntersectPoints;
            }
            tetrIntersectPoints.clear();

            auto evaluator = field->makeEvaluator();

            //check if splitting was neccesary
            bool splitting = false;

            //find all intersections on grid, using constant iterator for manipulating the original structure
            for (auto dP = diamondPoints.cbegin(); dP != diamondPoints.cend();)
            {
                std::vector<Point3> tetrIntPoints = findIntersections(evaluator);

                //splitting for more than one intersection found
                if (tetrIntPoints.size() == 2)
                {
                    splitting = true;

                    continue;
                }
                else if (tetrIntPoints.size() == 1)
                {
                    //add the one intersection two both Tetrahedra
                    tetrIntersectPoints[dP->first[0]].push_back(tetrIntPoints[0]);
                    tetrIntersectPoints[dP->first[1]].push_back(tetrIntPoints[0]);
                }
                ++dP;
            }

            if (splitting)
            {
                return marchingDiamonds(tetrahedraPoints, tetrIntersectPoints, diamondPoints, newTetrahedra, trCellCount, isovalue, splitCount - 1);
            }
            return tetrIntersectPoints;
        }

    private:
        
        template <typename T>
        std::vector<Point3> findIntersections(T &evaluator)
        {
            std::vector<Point3> trIntPoints;

            return trIntPoints;
        }

        std::vector<Point3> calcIntersection()
        {
            std::vector<Point3> intersections;

            return intersections;
        }

        

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
