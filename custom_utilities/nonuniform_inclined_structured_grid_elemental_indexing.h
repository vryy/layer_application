//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 24 Aug 2021 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_NONUNIFORM_INCLINED_STRUCTURED_GRID_ELEMENTAL_INDEXING_H_INCLUDED )
#define KRATOS_LAYER_APP_NONUNIFORM_INCLINED_STRUCTURED_GRID_ELEMENTAL_INDEXING_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/progress.h"
#include "custom_utilities/mesh_query_tool.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// Short class definition.
/*** Detail class definition.
 * This utility class supports for structured mesh indexing.
 * The structured grid is assumed to be non-uniform in every dimension.
 * There is no restriction on the element type, providing that the nodes are populated in structured manner.
 */
template<int TDim>
class NonuniformInclinedStructuredGridElementalIndexing : public MeshQueryTool<Element>
{
public:
    ///@name Type Definitions
    ///@{

    typedef MeshQueryTool<Element> BaseType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::PointType PointType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::EntitiesContainerType ElementsContainerType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NonuniformInclinedStructuredGridElementalIndexing);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NonuniformInclinedStructuredGridElementalIndexing(const double Tolerance)
    : BaseType()
    {
        this->SetTolerance(Tolerance);

        // default axis is Cartesian axis
        mAxes.resize(3);
        mAxes[0][0] = 1.0; mAxes[0][1] = 0.0; mAxes[0][2] = 0.0;
        mAxes[1][0] = 0.0; mAxes[1][1] = 1.0; mAxes[1][2] = 0.0;
        mAxes[2][0] = 0.0; mAxes[2][1] = 0.0; mAxes[2][2] = 1.0;
    }

    /// Destructor.
    ~NonuniformInclinedStructuredGridElementalIndexing() override
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Set the axis
    void SetAxis(const std::size_t i, const double x, const double y, const double z)
    {
        double s = sqrt(x*x + y*y + z*z);
        mAxes[i][0] = x/s; mAxes[i][1] = y/s; mAxes[i][2] = z/s;
    }

    /// Initialize the elements binning
    void Initialize( const ElementsContainerType& pElements ) final
    {
        std::cout << "Initialize the structured grid indexing, number of elements = " << pElements.size() << std::endl;

#ifdef _OPENMP
        double start_init = omp_get_wtime();
#endif

        const double TOL = this->Tolerance();

        // obtain all the coordinates
        mAxesPoints.clear();
        mAxesPoints.resize(TDim);

        std::vector<std::set<double> > axes_point_set(TDim);
        std::vector<PointType> corners;
        for ( ElementsContainerType::const_iterator it = pElements.begin(); it != pElements.end(); ++it )
        {
            this->GetCorners(corners, it->GetGeometry());
            for (std::size_t i = 0; i < corners.size(); ++i)
            {
                for (int dim = 0; dim < TDim; ++dim)
                {
                    const double ct = this->ComputeCoordinates(dim, corners[i][0], corners[i][1], corners[i][2]);
                    axes_point_set[dim].insert(ct);
                }
            }
        }

        for (int dim = 0; dim < TDim; ++dim)
        {
            auto it = axes_point_set[dim].begin();
            mAxesPoints[dim].push_back(*it);
            double prev = *it;
            ++it;
            for (; it != axes_point_set[dim].end(); ++it)
            {
                if (fabs(*it - prev) > TOL)
                {
                    mAxesPoints[dim].push_back(*it);
                    prev = *it;
                }
            }

            std::cout << "Coordinates array on " << dim << "-axis (size = " << mAxesPoints[dim].size() << "):";
            for (std::size_t i = 0; i < mAxesPoints[dim].size(); ++i)
                std::cout << " " << mAxesPoints[dim][i];
            std::cout << std::endl;
        }

        IndexType Nx=0, Ny=0;
        Nx = mAxesPoints[0].size();
        if constexpr (TDim > 1)
            Ny = mAxesPoints[1].size();

        Kratos::progress_display show_progress( pElements.size() );

        PointType C;
        std::size_t i1, i2=0, i3=0, I;
        for ( ElementsContainerType::const_iterator it = pElements.begin(); it != pElements.end(); ++it )
        {
            // compute center of the element
            noalias(C) = ZeroVector(3);
            this->GetCorners(corners, it->GetGeometry());
            for (std::size_t i = 0; i < corners.size(); ++i)
                noalias(C) += corners[i];
            C /= corners.size();

            // find the index number of the element
            const double ctx = this->ComputeCoordinates(0, C[0], C[1], C[2]);
            i1 = static_cast<std::size_t>(FindSpan(ctx, mAxesPoints[0], TOL));
            if constexpr (TDim > 1)
            {
                const double cty = this->ComputeCoordinates(1, C[0], C[1], C[2]);
                i2 = static_cast<std::size_t>(FindSpan(cty, mAxesPoints[1], TOL));
            }
            if constexpr (TDim > 2)
            {
                const double ctz = this->ComputeCoordinates(2, C[0], C[1], C[2]);
                i3 = static_cast<std::size_t>(FindSpan(ctz, mAxesPoints[2], TOL));
            }
            I = (i3*Ny + i2)*Nx + i1;

            mElemenIndexing[I] = it->Id();

            ++show_progress;
        }

        std::cout << "Initialize indexing completed";
#ifdef _OPENMP
        double stop_init = omp_get_wtime();
        std::cout << ", time = " << (stop_init-start_init) << "s";
#endif
        std::cout << std::endl;
    }

    /**
     * Get the list of neighbour elements of a point
     */
    void FindPotentialPartners( const PointType& rSourcePoint, std::vector<IndexType>& master_elements ) const final
    {
        // get the containing elements from the bin
        const double TOL = this->Tolerance();
        std::size_t i1, i2=0, i3=0;
        const double xt = this->ComputeCoordinates(0, rSourcePoint.X(), rSourcePoint.Y(), rSourcePoint.Z());
        i1 = static_cast<std::size_t>(FindSpan(xt, mAxesPoints[0], TOL));
        if constexpr (TDim > 1)
        {
            const double yt = this->ComputeCoordinates(1, rSourcePoint.X(), rSourcePoint.Y(), rSourcePoint.Z());
            i2 = static_cast<std::size_t>(FindSpan(yt, mAxesPoints[1], TOL));
        }
        if constexpr (TDim > 2)
        {
            const double zt = this->ComputeCoordinates(2, rSourcePoint.X(), rSourcePoint.Y(), rSourcePoint.Z());
            i3 = static_cast<std::size_t>(FindSpan(zt, mAxesPoints[2], TOL));
        }

        IndexType Nx=0, Ny=0;
        Nx = mAxesPoints[0].size();
        if constexpr (TDim > 1)
            Ny = mAxesPoints[1].size();

        const std::size_t I = (i3*Ny + i2)*Nx + i1;
        auto it = mElemenIndexing.find(I);
        if (it == mElemenIndexing.end())
        {
            KRATOS_WATCH(i1)
            KRATOS_WATCH(i2)
            KRATOS_WATCH(i3)
            KRATOS_WATCH(I)
            KRATOS_WATCH(Nx)
            KRATOS_WATCH(Ny)
            std::vector<double> vmin(TDim), vmax(TDim);
            for (int d = 0; d < TDim; ++d)
            {
                std::cout << "mAxesPoints[" << d << "]:" << std::endl;
                for (int i = 0; i < mAxesPoints[d].size(); ++i)
                    std::cout << " " << mAxesPoints[d][i];
                std::cout << std::endl;
                this->FindMinMax(d, vmin[d], vmax[d]);
            }
            // std::cout << "mElemenIndexing:" << std::endl;
            // for (auto it = mElemenIndexing.begin(); it != mElemenIndexing.end(); ++it)
            //     std::cout << it->first << " " << it->second << std::endl;
            std::cout << "Mesh region:";
            for (int d = 0; d < TDim; ++d)
                std::cout << " [" << vmin[d] << ", " << vmax[d] << "]";
            std::cout << std::endl;
            KRATOS_ERROR << "Point " << rSourcePoint << " does not locate in the mesh region";
            master_elements.resize(0);
            return;
        }

        master_elements.resize(1);
        master_elements[0] = it->second;
    }

    /// Find the containing element of a source point
    bool SearchPartner( const PointType& rSourcePoint, const ElementsContainerType& rAllElements,
        const std::vector<IndexType>& master_elements,
        Element::Pointer& pMatchedMaster, CoordinatesArrayType& local_coords) const final
    {
        auto it = rAllElements.find(master_elements[0]);
        if (it == rAllElements.end())
            KRATOS_ERROR << "The master element is not yet determined";
        pMatchedMaster = *(it.base());

        const GeometryType& rGeometry = pMatchedMaster->GetGeometry();
        // noalias(local_coords) = ZeroVector(3);

        this->PredictLocalPoint(rSourcePoint, local_coords, rGeometry);
        // CoordinatesArrayType predict_local_coords = local_coords;

        if (rGeometry.IsInside(rSourcePoint, local_coords))
        // if (IsInside(rGeometry, rSourcePoint, local_coords))
        {
            // if constexpr (TDim == 2)
            // {
            //     if (abs(local_coords[2]) > 1.0e-7)
            //     {
            //         KRATOS_WATCH(rSourcePoint)
            //         KRATOS_WATCH(pMatchedMaster->Id())
            //         for (std::size_t i = 0; i < rGeometry.size(); ++i)
            //             KRATOS_WATCH(rGeometry[i].GetSolutionStepValue(DISPLACEMENT))
            //         std::cout << "predict local point: " << predict_local_coords << std::endl;
            //         std::cout << "found point: " << local_coords << std::endl;
            //         KRATOS_ERROR << "Error computing local point";
            //     }
            // }
            return true;
        }

        return false;
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
        buffer << "NonuniformInclinedStructuredGridElementalIndexing" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
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

    std::vector<CoordinatesArrayType> mAxes;
    std::vector<std::vector<double> > mAxesPoints;
    std::map<std::size_t, std::size_t> mElemenIndexing;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /// Compute the transformed coordinates
    double ComputeCoordinates(const std::size_t i, const double x, const double y, const double z) const
    {
        return (x*mAxes[i][0] + y*mAxes[i][1] + z*mAxes[i][2]);
    }

    /// Find the span of the value in an array
    int FindSpan(const double v, const std::vector<double>& values, const double TOL) const
    {
        if (v < (values[0] - TOL))
        {
            return -1;
        }
        else
        {
            for (std::size_t i = 1; i < values.size(); ++i)
            {
                if (v < (values[i] + TOL))
                    return i-1;
            }
            return values.size();
        }
    }

    /// Find the minimum and maximum value on the axis
    void FindMinMax(const int dim, double& vmin, double& vmax) const
    {
        vmin = 1e99;
        vmax = -1e99;
        for (std::size_t i = 0; i < mAxesPoints[dim].size(); ++i)
        {
            if (mAxesPoints[dim][i] > vmax)
                vmax = mAxesPoints[dim][i];
            if (mAxesPoints[dim][i] < vmin)
                vmin = mAxesPoints[dim][i];
        }
    }

    /// Get corners of a geometry
    void GetCorners(std::vector<PointType>& Corners, const GeometryType& rGeometry) const;

    /// Compute a good estimation of the local coordinates of the point in the geometry
    void PredictLocalPoint(const CoordinatesArrayType& rSourcePoint, CoordinatesArrayType& local_coords, const GeometryType& rGeometry) const;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    NonuniformInclinedStructuredGridElementalIndexing& operator=(NonuniformInclinedStructuredGridElementalIndexing const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    NonuniformInclinedStructuredGridElementalIndexing(NonuniformInclinedStructuredGridElementalIndexing const& rOther)
    {
    }

    ///@}

}; // Class NonuniformInclinedStructuredGridElementalIndexing

///@}

template<>
void NonuniformInclinedStructuredGridElementalIndexing<1>::GetCorners(std::vector<PointType>& Corners, const GeometryType& rGeometry) const
{
    if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear)
    {
        if (Corners.size() != 2)
            Corners.resize(2);

        for (unsigned int i = 0; i < 2; ++i)
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
    }
    else
        KRATOS_ERROR << "Unsupported geometry family " << rGeometry.GetGeometryFamily();
}

template<>
void NonuniformInclinedStructuredGridElementalIndexing<2>::GetCorners(std::vector<PointType>& Corners, const GeometryType& rGeometry) const
{
    if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Triangle)
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
    else if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_NURBS)
    {
        if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Bezier2D
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
    }
    else
        KRATOS_ERROR << "Unsupported geometry family " << rGeometry.GetGeometryFamily();
}

template<>
void NonuniformInclinedStructuredGridElementalIndexing<3>::GetCorners(std::vector<PointType>& Corners, const GeometryType& rGeometry) const
{
    if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Tetrahedra)
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
        if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Bezier3D)
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
    }
    else
        KRATOS_ERROR << "Unsupported geometry family " << rGeometry.GetGeometryFamily();
}

template<>
void NonuniformInclinedStructuredGridElementalIndexing<1>::PredictLocalPoint(const CoordinatesArrayType& rSourcePoint, CoordinatesArrayType& local_coords, const GeometryType& rGeometry) const
{
    if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear)
    {
        std::vector<PointType> corners;
        this->GetCorners(corners, rGeometry);

        const double D1 = norm_2(corners[1] - corners[0]);

        noalias(local_coords) = ZeroVector(3);

        double xi = inner_prod(rSourcePoint - corners[0], corners[1] - corners[0]) / pow(D1, 2);
        local_coords[0] = 2.0*xi - 1.0;
    }
    else
    {
        noalias(local_coords) = ZeroVector(3);
    }
}

template<>
void NonuniformInclinedStructuredGridElementalIndexing<2>::PredictLocalPoint(const CoordinatesArrayType& rSourcePoint, CoordinatesArrayType& local_coords, const GeometryType& rGeometry) const
{
    if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Quadrilateral)
    {
        std::vector<PointType> corners;
        this->GetCorners(corners, rGeometry);

        const double D1 = norm_2(corners[1] - corners[0]);
        const double D2 = norm_2(corners[3] - corners[0]);

        noalias(local_coords) = ZeroVector(3);

        double xi = inner_prod(rSourcePoint - corners[0], corners[1] - corners[0]) / pow(D1, 2);
        local_coords[0] = 2.0*xi - 1.0;

        double eta = inner_prod(rSourcePoint - corners[0], corners[3] - corners[0]) / pow(D2, 2);
        local_coords[1] = 2.0*eta - 1.0;
    }
    else if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_NURBS)
    {
        if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Bezier2D
         || rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Bezier2D3)
        {
            std::vector<PointType> corners;
            this->GetCorners(corners, rGeometry);

            const double D1 = norm_2(corners[1] - corners[0]);
            const double D2 = norm_2(corners[3] - corners[0]);

            noalias(local_coords) = ZeroVector(3);
            local_coords[0] = inner_prod(rSourcePoint - corners[0], corners[1] - corners[0]) / pow(D1, 2);
            local_coords[1] = inner_prod(rSourcePoint - corners[0], corners[3] - corners[0]) / pow(D2, 2);
        }
    }
    else
    {
        noalias(local_coords) = ZeroVector(3);
    }
}

template<>
void NonuniformInclinedStructuredGridElementalIndexing<3>::PredictLocalPoint(const CoordinatesArrayType& rSourcePoint, CoordinatesArrayType& local_coords, const GeometryType& rGeometry) const
{
    if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Hexahedra)
    {
        std::vector<PointType> corners;
        this->GetCorners(corners, rGeometry);

        const double D1 = norm_2(corners[1] - corners[0]);
        const double D2 = norm_2(corners[3] - corners[0]);
        const double D3 = norm_2(corners[4] - corners[0]);

        double xi = inner_prod(rSourcePoint - corners[0], corners[1] - corners[0]) / pow(D1, 2);
        local_coords[0] = 2.0*xi - 1.0;

        double eta = inner_prod(rSourcePoint - corners[0], corners[3] - corners[0]) / pow(D2, 2);
        local_coords[1] = 2.0*eta - 1.0;

        double zeta = inner_prod(rSourcePoint - corners[0], corners[4] - corners[0]) / pow(D3, 2);
        local_coords[2] = 2.0*zeta - 1.0;
    }
    else if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_NURBS)
    {
        if (rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Bezier3D)
        {
            std::vector<PointType> corners;
            this->GetCorners(corners, rGeometry);

            const double D1 = norm_2(corners[1] - corners[0]);
            const double D2 = norm_2(corners[3] - corners[0]);
            const double D3 = norm_2(corners[4] - corners[0]);

            local_coords[0] = inner_prod(rSourcePoint - corners[0], corners[1] - corners[0]) / pow(D1, 2);
            local_coords[1] = inner_prod(rSourcePoint - corners[0], corners[3] - corners[0]) / pow(D2, 2);
            local_coords[2] = inner_prod(rSourcePoint - corners[0], corners[4] - corners[0]) / pow(D3, 2);
        }
    }
    else
    {
        noalias(local_coords) = ZeroVector(3);
    }
}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_NONUNIFORM_INCLINED_STRUCTURED_GRID_ELEMENTAL_INDEXING_H_INCLUDED defined

