//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Jul 2021 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_SPATIAL_GRID_ELEMENTAL_BINNING_H_INCLUDED )
#define KRATOS_LAYER_APP_SPATIAL_GRID_ELEMENTAL_BINNING_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/progress.h"
#include "custom_utilities/spatial_key.h"
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

/***
 * This utility class supports for element binning into a structured grid.
 */
class SpatialGridElementalBinning : public MeshQueryTool<Element>
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
    KRATOS_CLASS_POINTER_DEFINITION(SpatialGridElementalBinning);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpatialGridElementalBinning(double X0, double Y0, double Z0, double Dx, double Dy, double Dz, double tol)
    : mX0(X0), mY0(Y0), mZ0(Z0), mDx(Dx), mDy(Dy), mDz(Dz)
    {
        this->SetTolerance(tol);
    }

    /// Destructor.
    ~SpatialGridElementalBinning() override
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
        std::cout << "Initialize the spatial binning, number of elements = " << pElements.size() << std::endl;

#ifdef _OPENMP
        double start_init = omp_get_wtime();
#endif

        Kratos::progress_display show_progress( pElements.size() );
        std::vector<double> vmin(3);
        std::vector<double> vmax(3);
        for ( ElementsContainerType::ptr_const_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it )
        {
            mAllElements.insert((*it)->Id());
            this->FindBoundingBox(vmin, vmax, (*it)->GetGeometry());

            int i_min = (int) std::floor((vmin[0] - mX0) / mDx), i_max = (int) std::floor((vmax[0] - mX0) / mDx);
            int j_min = (int) std::floor((vmin[1] - mY0) / mDy), j_max = (int) std::floor((vmax[1] - mY0) / mDy);
            int k_min = (int) std::floor((vmin[2] - mZ0) / mDz), k_max = (int) std::floor((vmax[2] - mZ0) / mDz);

            for (int i = i_min; i <= i_max; ++i)
            {
                for (int j = j_min; j <= j_max; ++j)
                {
                    for (int k = k_min; k <= k_max; ++k)
                    {
                        SpatialKey key(i, j, k);
                        mBinElements[key].insert((*it)->Id());
                    }
                }
            }

            ++show_progress;
        }

        KRATOS_WATCH(mBinElements.size())

#ifdef _OPENMP
        double stop_init = omp_get_wtime();
        std::cout << "Initialize binning completed, time = " << (stop_init-start_init) << "s" << std::endl;
#endif
    }

    /**
     * Get the list of neighbour elements of a point
     */
    void FindPotentialPartners( const PointType& rSourcePoint, std::vector<IndexType>& master_elements ) const final
    {
        // get the containing elements from the bin
        int ix = (int) std::floor((rSourcePoint.X() - mX0) / mDx);
        int iy = (int) std::floor((rSourcePoint.Y() - mY0) / mDy);
        int iz = (int) std::floor((rSourcePoint.Z() - mZ0) / mDz);

        SpatialKey key(ix, iy, iz);
        auto it_bin_elements = mBinElements.find(key);

        if (it_bin_elements != mBinElements.end())
        {
            for (auto it = it_bin_elements->second.begin(); it != it_bin_elements->second.end(); ++it)
            {
                auto it_elem = mAllElements.find(*it);
                if (it_elem != mAllElements.end())
                    master_elements.push_back(*it_elem);
            }
        }

        #ifdef ENABLE_ERROR_CHECK
        if (master_elements.size() == 0)
        {
            if (it_bin_elements == mBinElements.end())
                std::cout << "Binning does not contain point " << rSourcePoint << std::endl;
            KRATOS_WATCH(ix)
            KRATOS_WATCH(iy)
            KRATOS_WATCH(iz)
            KRATOS_ERROR << "Point " << rSourcePoint << " does not locate in the mesh region";
        }
        #endif
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
        buffer << "SpatialGridElementalBinning";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SpatialGridElementalBinning";
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

    double mDx, mDy, mDz;

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

    double mX0, mY0, mZ0;

    SetType mAllElements;
    BinType mBinElements;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void FindBoundingBox(std::vector<double>& vmin, std::vector<double>& vmax, const GeometryType& rGeometry) const
    {
        vmin[0] = rGeometry[0].X(); vmin[1] = rGeometry[0].Y(); vmin[2] = rGeometry[0].Z();
        vmax[0] = rGeometry[0].X(); vmax[1] = rGeometry[0].Y(); vmax[2] = rGeometry[0].Z();
        for (std::size_t i = 1; i < rGeometry.size(); ++i)
        {
            if (rGeometry[i].X() < vmin[0]) vmin[0] = rGeometry[i].X();
            if (rGeometry[i].X() > vmax[0]) vmax[0] = rGeometry[i].X();
            if (rGeometry[i].Y() < vmin[1]) vmin[1] = rGeometry[i].Y();
            if (rGeometry[i].Y() > vmax[1]) vmax[1] = rGeometry[i].Y();
            if (rGeometry[i].Z() < vmin[2]) vmin[2] = rGeometry[i].Z();
            if (rGeometry[i].Z() > vmax[2]) vmax[2] = rGeometry[i].Z();
        }
    }

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
    SpatialGridElementalBinning& operator=(SpatialGridElementalBinning const& rOther)
    {
        mX0 = rOther.mX0;
        mY0 = rOther.mY0;
        mZ0 = rOther.mZ0;
        mDx = rOther.mDx;
        mDy = rOther.mDy;
        mDz = rOther.mDz;
        BaseType::operator=(rOther);
        return *this;
    }

    /// Copy constructor.
    SpatialGridElementalBinning(SpatialGridElementalBinning const& rOther)
    : mX0(rOther.mX0), mY0(rOther.mY0), mZ0(rOther.mZ0)
    , mDx(rOther.mDx), mDy(rOther.mDy), mDz(rOther.mDz)
    , BaseType(rOther)
    {
    }

    ///@}

}; // Class SpatialGridElementalBinning

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos.

#undef ENABLE_ERROR_CHECK

#endif // KRATOS_LAYER_APP_SPATIAL_GRID_ELEMENTAL_BINNING_H_INCLUDED defined
