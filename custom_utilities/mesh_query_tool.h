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
#ifdef SD_APP_FORWARD_COMPATIBILITY
#include "includes/indexed_object.h"
#else
#include "utilities/indexed_object.h"
#endif


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
template<class TEntityType, class TEntitiesContainerType = PointerVectorSet<TEntityType, IndexedObject> >
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
    typedef TEntitiesContainerType EntitiesContainerType;

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

    void SetTolerance(const double v)
    {
        mTolerance = v;
    }

    /// Initialize the elements binning
    virtual void Initialize( const EntitiesContainerType& pElements )
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /**
     * Get the list of neighbour elements of a point
     */
    virtual void FindPotentialPartners( const PointType& rSourcePoint, std::vector<IndexType>& master_elements ) const
    {
        KRATOS_ERROR << "Error calling base class function";
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

    /// Assignment operator.
    MeshQueryTool& operator=(MeshQueryTool const& rOther)
    {
        mTolerance = rOther.mTolerance;
        return *this;
    }

    /// Copy constructor.
    MeshQueryTool(MeshQueryTool const& rOther)
    : mTolerance(rOther.mTolerance)
    {
    }

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

    ///@}

}; // Class MeshQueryTool

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TEntityType, class TEntitiesContainerType>
inline std::istream& operator >>(std::istream& rIStream, MeshQueryTool<TEntityType, TEntitiesContainerType>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TEntityType, class TEntitiesContainerType>
inline std::ostream& operator <<(std::ostream& rOStream, const MeshQueryTool<TEntityType, TEntitiesContainerType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_MESH_QUERY_TOOL_H_INCLUDED defined
