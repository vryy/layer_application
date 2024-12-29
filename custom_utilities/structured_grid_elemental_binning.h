//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 27 Dec 2024 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_STRUCTURED_GRID_ELEMENTAL_BINNING_H_INCLUDED )
#define KRATOS_LAYER_APP_STRUCTURED_GRID_ELEMENTAL_BINNING_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_utilities/spatial_grid_elemental_binning.h"


namespace Kratos
{
///@addtogroup LayerApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/***
 * This utility class supports for element binning into a structured grid, assuming that element mesh is also structured mesh.
 */
class StructuredGridElementalBinning : public SpatialGridElementalBinning
{
public:
    ///@name Type Definitions
    ///@{

    typedef SpatialGridElementalBinning BaseType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(StructuredGridElementalBinning);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StructuredGridElementalBinning(double X0, double Y0, double Z0, double tol)
    : BaseType(X0, Y0, Z0, 1.0, 1.0, 1.0, tol)
    {
    }

    /// Destructor.
    ~StructuredGridElementalBinning() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Initialize the elements binning
    void Initialize( const ElementsContainerType& pElements ) override
    {
        // determine the points on the axes
        std::vector<std::set<double> > axes_point_set(3);
        std::vector<PointType> corners;
        for ( ElementsContainerType::const_iterator it = pElements.begin(); it != pElements.end(); ++it )
        {
            this->GetCorners(corners, it->GetGeometry());
            for (std::size_t i = 0; i < corners.size(); ++i)
            {
                for (int dim = 0; dim < 3; ++dim)
                    axes_point_set[dim].insert(corners[i][dim]);
            }
        }

        const double TOL = this->Tolerance();
        std::vector<std::vector<double> > axes_points(3);
        for (int dim = 0; dim < 3; ++dim)
        {
            auto it = axes_point_set[dim].begin();
            axes_points[dim].push_back(*it);
            double prev = *it;
            ++it;
            for (; it != axes_point_set[dim].end(); ++it)
            {
                if (fabs(*it - prev) > TOL)
                {
                    axes_points[dim].push_back(*it);
                    prev = *it;
                }
            }

            std::cout << "StructuredGridElementalBinning: Coordinates array on " << dim << "-axis (size = " << axes_points[dim].size() << "):";
            for (std::size_t i = 0; i < axes_points[dim].size(); ++i)
                std::cout << " " << axes_points[dim][i];
            std::cout << std::endl;
        }

        // find the smallest increment in each axis
        std::vector<double> dmin(3);
        for (int dim = 0; dim < 3; ++dim)
        {
            if (axes_points.size() > 1)
            {
                dmin[dim] = 1e99;
                for (unsigned int i = 0; i < axes_points.size() - 1; ++i)
                {
                    double d = axes_points[dim][i+1] - axes_points[dim][i];
                    if (d < dmin[dim])
                        dmin[dim] = d;
                }
            }
            else
                dmin[dim] = 0.0;
        }

        if (dmin[0] > 0.0) mDx = dmin[0];
        if (dmin[1] > 0.0) mDy = dmin[1];
        if (dmin[2] > 0.0) mDz = dmin[2];

        std::cout << "StructuredGridElementalBinning: Dx = " << mDx << ", Dy = " << mDy << ", Dz = " << mDz << std::endl;

        // initialize the binning
        BaseType::Initialize(pElements);
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "StructuredGridElementalBinning";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StructuredGridElementalBinning";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {}

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /// Get corners of a geometry
    void GetCorners(std::vector<PointType>& Corners, const GeometryType& rGeometry) const;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class StructuredGridElementalBinning

///@}

void StructuredGridElementalBinning::GetCorners(std::vector<PointType>& Corners, const GeometryType& rGeometry) const
{
    if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear)
    {
        if (Corners.size() != 2)
            Corners.resize(2);

        for (unsigned int i = 0; i < 2; ++i)
            noalias(Corners[i]) = rGeometry[i];
    }
    else if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Triangle)
    {
        if (Corners.size() != 3)
            Corners.resize(3);

        for (unsigned int i = 0; i < 3; ++i)
            noalias(Corners[i]) = rGeometry[i];
    }
    else if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Quadrilateral)
    {
        if (Corners.size() != 4)
            Corners.resize(4);

        for (unsigned int i = 0; i < 4; ++i)
            noalias(Corners[i]) = rGeometry[i];
    }
    else if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Tetrahedra)
    {
        if (Corners.size() != 4)
            Corners.resize(4);

        for (unsigned int i = 0; i < 4; ++i)
            noalias(Corners[i]) = rGeometry[i];
    }
    else if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Hexahedra)
    {
        if (Corners.size() != 8)
            Corners.resize(8);

        for (unsigned int i = 0; i < 8; ++i)
            noalias(Corners[i]) = rGeometry[i];
    }
    else if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_NURBS)
    {
        if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Bezier1D)
        {
            CoordinatesArrayType p;
            noalias(p) = ZeroVector(3);

            if (Corners.size() != 2)
                Corners.resize(2);

            p[0] = 0.0;
            rGeometry.GlobalCoordinates(Corners[0], p);
            p[0] = 1.0;
            rGeometry.GlobalCoordinates(Corners[1], p);
        }
        else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Bezier2D
         || rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Bezier2D3)
        {
            CoordinatesArrayType p;
            noalias(p) = ZeroVector(3);

            if (Corners.size() != 4)
                Corners.resize(4);

            p[0] = 0.0; p[1] = 0.0;
            rGeometry.GlobalCoordinates(Corners[0], p);
            p[0] = 1.0; p[1] = 0.0;
            rGeometry.GlobalCoordinates(Corners[1], p);
            p[0] = 1.0; p[1] = 1.0;
            rGeometry.GlobalCoordinates(Corners[2], p);
            p[0] = 0.0; p[1] = 1.0;
            rGeometry.GlobalCoordinates(Corners[3], p);
        }
        else if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Bezier3D)
        {
            CoordinatesArrayType p;
            noalias(p) = ZeroVector(3);

            if (Corners.size() != 8)
                Corners.resize(8);

            p[0] = 0.0; p[1] = 0.0;
            rGeometry.GlobalCoordinates(Corners[0], p);
            p[0] = 1.0; p[1] = 0.0;
            rGeometry.GlobalCoordinates(Corners[1], p);
            p[0] = 1.0; p[1] = 1.0;
            rGeometry.GlobalCoordinates(Corners[2], p);
            p[0] = 0.0; p[1] = 1.0;
            rGeometry.GlobalCoordinates(Corners[3], p);

            p[0] = 0.0; p[1] = 0.0; p[2] = 1.0;
            rGeometry.GlobalCoordinates(Corners[4], p);
            p[0] = 1.0; p[1] = 0.0; p[2] = 1.0;
            rGeometry.GlobalCoordinates(Corners[5], p);
            p[0] = 1.0; p[1] = 1.0; p[2] = 1.0;
            rGeometry.GlobalCoordinates(Corners[6], p);
            p[0] = 0.0; p[1] = 1.0; p[2] = 1.0;
            rGeometry.GlobalCoordinates(Corners[7], p);
        }
        else
            KRATOS_ERROR << "Invalid geometry type " << rGeometry.GetGeometryType();
    }
    else
        KRATOS_ERROR << "Invalid geometry family " << rGeometry.GetGeometryFamily();
}


///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_LAYER_APP_STRUCTURED_GRID_ELEMENTAL_BINNING_H_INCLUDED defined
