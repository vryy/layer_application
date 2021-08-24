//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 27 Jul 2021 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_STRUCTURED_GRID_ELEMENTAL_INDEXING_H_INCLUDED )
#define KRATOS_LAYER_APP_STRUCTURED_GRID_ELEMENTAL_INDEXING_H_INCLUDED

// System includes

// External includes
#include "boost/progress.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
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
 * The structured grid is assumed to be uniform in every dimension.
 * The element connectivity should be aligned with the axis
 * Only quadrilateral/hexahedra element is supported
 */
template<int TDim>
class StructuredGridElementalIndexing : public MeshQueryTool<Element>
{
public:
    ///@name Type Definitions
    ///@{

    typedef MeshQueryTool<Element> BaseType;
    typedef typename BaseType::IndexType IndexType;
    typedef std::set<IndexType> SetType;
    typedef std::map<SpatialKey, std::set<IndexType> > BinType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::PointType PointType;
    typedef typename BaseType::EntitiesContainerType ElementsContainerType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(StructuredGridElementalIndexing);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StructuredGridElementalIndexing(double X0, double Dx, std::size_t Nx)
    : BaseType()
    , mX0(X0), mY0(0.0), mZ0(0.0), mDx(Dx), mDy(1.0), mDz(1.0)
    , mNx(Nx), mNy(0), mNz(0)
    {}

    StructuredGridElementalIndexing(double X0, double Y0,
        double Dx, double Dy,
        std::size_t Nx, std::size_t Ny)
    : BaseType()
    , mX0(X0), mY0(Y0), mZ0(0.0), mDx(Dx), mDy(Dy), mDz(1.0)
    , mNx(Nx), mNy(Ny), mNz(0)
    {}

    StructuredGridElementalIndexing(double X0, double Y0, double Z0,
        double Dx, double Dy, double Dz,
        std::size_t Nx, std::size_t Ny, std::size_t Nz)
    : BaseType()
    , mX0(X0), mY0(Y0), mZ0(Z0), mDx(Dx), mDy(Dy), mDz(Dz)
    , mNx(Nx), mNy(Ny), mNz(Nz)
    {}

    /// Destructor.
    virtual ~StructuredGridElementalIndexing()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Initialize the elements binning
    void Initialize( const ElementsContainerType& pElements ) final
    {
        std::cout << "Initialize the structured grid indexing" << std::endl;

#ifdef _OPENMP
        double start_init = omp_get_wtime();
#endif

        boost::progress_display show_progress( pElements.size() );

        PointType C;
        const double& TOL = this->Tolerance();
        std::size_t i1, i2=0, i3=0, I;
        for ( ElementsContainerType::const_iterator it = pElements.begin(); it != pElements.end(); ++it )
        {
            // compute center of the element
            noalias(C) = ZeroVector(3);
            for (std::size_t i = 0; i < it->GetGeometry().size(); ++i)
            {
                noalias(C) += it->GetGeometry()[i].GetInitialPosition();
            }
            C /= it->GetGeometry().size();

            // find the index number of the element
            i1 = static_cast<std::size_t>(floor((C[0] - mX0 - TOL) / mDx));
            if (TDim > 1)
                i2 = static_cast<std::size_t>(floor((C[1] - mY0 - TOL) / mDy));
            if (TDim > 2)
                i3 = static_cast<std::size_t>(floor((C[2] - mZ0 - TOL) / mDz));
            I = (i3*mNy + i2)*mNx + i1;

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
        const double& TOL = this->Tolerance();
        std::size_t i1, i2=0, i3=0;
        i1 = static_cast<std::size_t>(floor((rSourcePoint.X() - mX0 - TOL) / mDx));
        if (TDim > 1)
            i2 = static_cast<std::size_t>(floor((rSourcePoint.Y() - mY0 - TOL) / mDy));
        if (TDim > 2)
            i3 = static_cast<std::size_t>(floor((rSourcePoint.Z() - mZ0 - TOL) / mDz));

        const std::size_t I = (i3*mNy + i2)*mNx + i1;
        auto it = mElemenIndexing.find(I);
        if (it == mElemenIndexing.end())
        {
            KRATOS_WATCH(i1)
            KRATOS_WATCH(i2)
            KRATOS_WATCH(i3)
            KRATOS_WATCH(rSourcePoint)
            KRATOS_THROW_ERROR(std::logic_error, "The point does not locate in the mesh region", "")
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
            KRATOS_THROW_ERROR(std::logic_error, "The master element is not yet determined", "")
        pMatchedMaster = *(it.base());

        const GeometryType& rGeometry = pMatchedMaster->GetGeometry();
        noalias(local_coords) = ZeroVector(3);

        double xi = inner_prod(rSourcePoint - rGeometry[0], rGeometry[1] - rGeometry[0]) / pow(mDx, 2);
        local_coords[0] = 2.0*xi - 1.0;

        if (TDim > 1)
        {
            double eta = inner_prod(rSourcePoint - rGeometry[0], rGeometry[3] - rGeometry[0]) / pow(mDy, 2);
            local_coords[1] = 2.0*eta - 1.0;
        }

        if (TDim > 2)
        {
            double zeta = inner_prod(rSourcePoint - rGeometry[0], rGeometry[4] - rGeometry[0]) / pow(mDz, 2);
            local_coords[2] = 2.0*zeta - 1.0;
        }

        return true;
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "StructuredGridElementalIndexing" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    double mX0, mY0, mZ0, mDx, mDy, mDz;
    std::size_t mNx, mNy, mNz;
    std::map<std::size_t, std::size_t> mElemenIndexing;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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
    StructuredGridElementalIndexing& operator=(StructuredGridElementalIndexing const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    StructuredGridElementalIndexing(StructuredGridElementalIndexing const& rOther)
    {
    }

    ///@}

}; // Class StructuredGridElementalIndexing

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
//inline std::istream& operator >>(std::istream& rIStream, StructuredGridElementalIndexing& rThis)
//{
//    return rIStream;
//}

///// output stream function
//inline std::ostream& operator <<(std::ostream& rOStream,
//        const StructuredGridElementalIndexing& rThis)
//{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);

//    return rOStream;
//}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_STRUCTURED_GRID_ELEMENTAL_INDEXING_H_INCLUDED defined

