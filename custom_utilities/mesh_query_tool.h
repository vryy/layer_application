//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Jul 2021 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_MESH_QUERY_TOOL_H_INCLUDED )
#define KRATOS_LAYER_APP_MESH_QUERY_TOOL_H_INCLUDED

// System includes
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <map>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/pointer_vector_set.h"
#include "utilities/indexed_object.h"


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
 * Abstract class for mesh query tool of element/condition
 */
template<class TEntityType>
class MeshQueryTool
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef typename TEntityType::GeometryType GeometryType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename GeometryType::PointType NodeType;
    typedef typename NodeType::PointType PointType;
    typedef PointerVectorSet<TEntityType, IndexedObject> EntitiesContainerType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MeshQueryTool);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MeshQueryTool() : mTolerance(0.0)
    {
    }

    /// Destructor.
    virtual ~MeshQueryTool()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Get and Set the tolerance
    double Tolerance() const
    {
        return mTolerance;
    }

    void SetTolerance(const double& v)
    {
        mTolerance = v;
    }

    /// Initialize the elements binning
    virtual void Initialize( const EntitiesContainerType& pElements )
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    /**
     * Get the list of neighbour elements of a point
     */
    virtual void FindPotentialPartners( const PointType& rSourcePoint, std::vector<IndexType>& master_elements ) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    /// Find the containing element of a source point
    virtual bool SearchPartner( const PointType& rSourcePoint, const EntitiesContainerType& rAllElements,
        const std::vector<IndexType>& master_elements,
        typename TEntityType::Pointer& pMatchedMaster, CoordinatesArrayType& local_coords) const
    {
        for (auto elem_id : master_elements)
        {
            auto it = rAllElements.find(elem_id);
            if (it == rAllElements.end())
                continue; // that should not happen but we do it for safety

            if (it->GetGeometry().IsInside(rSourcePoint, local_coords))
            // if (IsInside(it->GetGeometry(), rSourcePoint, local_coords))
            {
                pMatchedMaster = *(it.base());
                return true;
            }
        }

        pMatchedMaster = NULL;
        local_coords.clear();
        return false;
    }

    // /// Compute the local coordinates and check if the point is inside the geometry
    // static bool IsInside(const GeometryType& rGeometry, const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
    // {
    //     PointLocalCoordinates( rGeometry, rResult, rPoint );
    //     return rGeometry.IsInside(rResult);
    // }

    // /// Better version of PointLocalCoordinates in the geometry. This version uses stringer tolerance and allows for external
    // /// prediction of local coordinates
    // static CoordinatesArrayType& PointLocalCoordinates( const GeometryType& rGeometry, CoordinatesArrayType& rResult,
    //         const CoordinatesArrayType& rPoint )
    // {
    //     const unsigned int working_dim = rGeometry.WorkingSpaceDimension();
    //     const unsigned int local_dim = rGeometry.LocalSpaceDimension();

    //     if (working_dim != local_dim)
    //         KRATOS_THROW_ERROR(std::logic_error, "Attention, the Point Local Coordinates must be specialized for the current geometry", "");

    //     Matrix J = ZeroMatrix( local_dim, local_dim );

    //     Vector DeltaXi = ZeroVector( local_dim );

    //     CoordinatesArrayType CurrentGlobalCoords( ZeroVector( 3 ) );

    //     //reset the prediction if it is outside the valid domain
    //     if (norm_2(rResult) >= working_dim)
    //         rResult.clear();

    //     //Newton iteration:
    //     constexpr double tol_xi = 1.0e-8;
    //     constexpr double tol_g = 1.0e-20;

    //     constexpr int maxiter = 1000;

    //     constexpr double max_norm_xi = 30.0;

    //     // CoordinatesArrayType PredictResult = rResult;

    //     int k;
    //     double norm_dxi;
    //     for ( k = 0; k < maxiter; k++ )
    //     {
    //         CurrentGlobalCoords.clear();
    //         rGeometry.GlobalCoordinates( CurrentGlobalCoords, rResult );
    //         noalias( CurrentGlobalCoords ) = rPoint - CurrentGlobalCoords;

    //         if (norm_2(CurrentGlobalCoords) < tol_g)
    //             break;

    //         DeltaXi.clear();
    //         rGeometry.InverseOfJacobian( J, rResult );
    //         for(unsigned int i = 0; i < local_dim; i++)
    //         {
    //             for(unsigned int j = 0; j < local_dim; j++)
    //             {
    //                 DeltaXi[i] += J(i,j)*CurrentGlobalCoords[j];
    //             }
    //             rResult[i] += DeltaXi[i];
    //         }

    //         norm_dxi = norm_2(DeltaXi);

    //         if ( norm_dxi > max_norm_xi )
    //         {
    //             KRATOS_THROW_ERROR(std::logic_error, "Computation of point local coordinates fails at step", k)
    //             break;
    //         }

    //         if ( norm_dxi < tol_xi )
    //         {
    //             break;
    //         }
    //     }

    //     if (k > 100)
    //     {
    //         // KRATOS_WATCH(PredictResult)
    //         KRATOS_WATCH(rResult)
    //         KRATOS_WATCH(rPoint)
    //         KRATOS_WATCH(CurrentGlobalCoords)
    //         KRATOS_WATCH(DeltaXi)
    //         KRATOS_WATCH(MathUtils<double>::Norm3( DeltaXi ))
    //         KRATOS_WATCH(norm_2( DeltaXi ))
    //         std::cout << "PointLocalCoordinates iterations: " << k << std::endl;
    //         KRATOS_THROW_ERROR(std::logic_error, "too much iterations", "")
    //     }

    //     return( rResult );
    // }

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
        buffer << "MeshQueryTool";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MeshQueryTool";
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

    double mTolerance;

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
    MeshQueryTool& operator=(MeshQueryTool const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    MeshQueryTool(MeshQueryTool const& rOther)
    {
    }

    ///@}

}; // Class MeshQueryTool

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
//inline std::istream& operator >>(std::istream& rIStream, MeshQueryTool& rThis)
//{
//    return rIStream;
//}

///// output stream function
//inline std::ostream& operator <<(std::ostream& rOStream,
//        const MeshQueryTool& rThis)
//{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);

//    return rOStream;
//}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_MESH_QUERY_TOOL_H_INCLUDED defined
