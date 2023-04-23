//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Jul 2019 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_MODEL_PART_UTILITY_H_INCLUDED )
#define  KRATOS_LAYER_APP_MODEL_PART_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>
#include <set>
#include <ctime>
#include <iomanip>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


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
 */
class ModelPartUtilities
{
public:

    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ModelPartUtilities);


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ModelPartUtilities()
    {
    }

    /// Destructor.
    virtual ~ModelPartUtilities()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Export the nodal coordinates to output stream
    static void ExportNodalCoordinates(std::ostream& rOStream, ModelPart& r_model_part,
        const bool& with_id, const bool& with_x, const bool& with_y, const bool& with_z,
        const std::string& separator, const int& precision)
    {
        rOStream << std::setprecision(precision);

        for(ModelPart::NodeIterator i_node = r_model_part.NodesBegin();
                i_node != r_model_part.NodesEnd(); ++i_node)
        {
            bool init = false;

            if (with_id)
            {
                rOStream << i_node->Id();
                init = true;
            }

            if (with_x)
            {
                if (init)
                    rOStream << separator;
                rOStream << i_node->X0();
                init = true;
            }

            if (with_y)
            {
                if (init)
                    rOStream << separator;
                rOStream << i_node->Y0();
            }

            if (with_z)
            {
                if (init)
                    rOStream << separator;
                rOStream << i_node->Z0();
            }

            rOStream << std::endl;
        }
    }

    /// Export the nodal coordinates to output stream
    static void ExportNodalCoordinatesToGiD(std::ostream& rOStream, ModelPart& r_model_part,
        const int& precision)
    {
        rOStream << std::setprecision(precision);

        for(ModelPart::NodeIterator i_node = r_model_part.NodesBegin();
                i_node != r_model_part.NodesEnd(); ++i_node)
        {
            rOStream << "Mescape Geometry Create Point "
                     << i_node->X0() << " " << i_node->Y0() << " " << i_node->Z0()
                     << std::endl;
        }
    }

    /// Export the edge information to output stream
    /// REMARK: for second order element, only the 2-node edge on the side is returned
    static void ExportEdgeInformation(std::ostream& rOStream, ModelPart::ElementsContainerType& rpElements,
        const std::string& separator)
    {
        typedef std::pair<std::size_t, std::size_t> edge_t;
        typedef std::set<edge_t> edge_container_t;

        edge_container_t edges;
        ExtractEdgeInformation(edges, rpElements);

        for (edge_container_t::iterator it = edges.begin(); it != edges.end(); ++it)
        {
            rOStream << it->first << separator << it->second << std::endl;
        }
    }

    /// Export the edge information to output stream
    /// REMARK: for second order element, only the 2-node edge on the side is returned
    static void ExportEdgeInformationToGiD(std::ostream& rOStream, ModelPart::ElementsContainerType& rpElements)
    {
        typedef std::pair<std::size_t, std::size_t> edge_t;
        typedef std::set<edge_t> edge_container_t;

        edge_container_t edges;
        ExtractEdgeInformation(edges, rpElements);

        for (edge_container_t::iterator it = edges.begin(); it != edges.end(); ++it)
        {
            rOStream << "Mescape Geometry Create Line Join " << it->first << " " << it->second << std::endl;
        }
    }

    /// Create a new entity out from a geometry
    template<class TEntityType>
    static typename TEntityType::Pointer CreateEntity(const std::string& sample_entity_name,
        const std::size_t& Id, Properties::Pointer pProperties, typename TEntityType::GeometryType& rGeometry)
    {
        if(!KratosComponents<TEntityType>::Has(sample_entity_name))
            KRATOS_THROW_ERROR(std::logic_error, sample_entity_name, "is not registered to the KRATOS kernel")
        TEntityType const& r_clone_entity = KratosComponents<TEntityType>::Get(sample_entity_name);

        typename TEntityType::Pointer pNewEntity = r_clone_entity.Create(Id, rGeometry, pProperties);

        return pNewEntity;
    }

    /// Create a new entity out from a list of nodes
    template<class TEntityType>
    static typename TEntityType::Pointer CreateEntity(ModelPart& r_model_part, const std::string& sample_entity_name,
        const std::size_t& Id, Properties::Pointer pProperties, const std::vector<std::size_t>& node_ids)
    {
        if(!KratosComponents<TEntityType>::Has(sample_entity_name))
            KRATOS_THROW_ERROR(std::logic_error, sample_entity_name, "is not registered to the KRATOS kernel")
        TEntityType const& r_clone_entity = KratosComponents<TEntityType>::Get(sample_entity_name);

        // create the points array
        typename TEntityType::GeometryType::PointsArrayType Points;
        for(std::size_t i = 0; i < node_ids.size(); ++i)
            Points.push_back(r_model_part.pGetNode(node_ids[i]));

        // create new entity
        typename TEntityType::Pointer pNewEntity = r_clone_entity.Create(Id, Points, pProperties);

        return pNewEntity;
    }

    /// Invoke the CalculateLocalSystem of an Element/Condition. It is useful when debugging them.
    template<class TEntityType>
    static void CalculateLocalSystem(TEntityType& rEntity, const ProcessInfo& rCurrentProcessInfo, const int& echo_level)
    {
        Matrix LHS;
        Vector RHS;

        rEntity.CalculateLocalSystem(LHS, RHS, rCurrentProcessInfo);

        if (echo_level > 0)
        {
            rEntity.PrintInfo(std::cout);
            std::cout << ":" << std::endl;
            KRATOS_WATCH(LHS)
            KRATOS_WATCH(RHS)
        }
    }

    /// Invoke the CalculateLocalSystem of an Element/Condition. It is useful when debugging them.
    template<class TEntityType>
    static void CalculateMassMatrix(TEntityType& rEntity, const ProcessInfo& rCurrentProcessInfo, const int& echo_level)
    {
        Matrix MassMatrix;

        rEntity.CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

        if (echo_level > 0)
        {
            rEntity.PrintInfo(std::cout);
            std::cout << ":" << std::endl;
            KRATOS_WATCH(MassMatrix)
        }
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
        buffer << "ModelPartUtilities";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    void Print()
    {
        PrintInfo(std::cout);
        std::cout << std::endl;
        PrintData(std::cout);
    }

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

    // decode the 1st order edges of the geometry
    static const int msT3Edges[][2];
    static const int msQ4Edges[][2];
    static const int msT4Edges[][2];
    static const int msH8Edges[][2];

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    template<typename edge_container_t>
    static void ExtractEdgeInformation(edge_container_t& edges, ModelPart::ElementsContainerType& rpElements)
    {
        typedef typename edge_container_t::value_type edge_t;
        std::size_t n1, n2, n;
        for (typename ModelPart::ElementsContainerType::ptr_iterator it = rpElements.ptr_begin();
                it != rpElements.ptr_end(); ++it)
        {
            if ( (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3
              || (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D6
              || (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3
              || (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D6 )
            {
                for (unsigned int i = 0; i < 3; ++i)
                {
                    n1 = (*it)->GetGeometry()[msT3Edges[i][0]].Id();
                    n2 = (*it)->GetGeometry()[msT3Edges[i][1]].Id();
                    if (n2 > n1)
                    {
                        n = n1;
                        n1 = n2;
                        n2 = n;
                    }
                    edge_t e(n1, n2);
                    edges.insert(e);
                }
            }
            else if ( (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4
                   || (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8
                   || (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9
                   || (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4
                   || (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8
                   || (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9 )
            {
                for (unsigned int i = 0; i < 4; ++i)
                {
                    n1 = (*it)->GetGeometry()[msQ4Edges[i][0]].Id();
                    n2 = (*it)->GetGeometry()[msQ4Edges[i][1]].Id();
                    if (n2 > n1)
                    {
                        n = n1;
                        n1 = n2;
                        n2 = n;
                    }
                    edge_t e(n1, n2);
                    edges.insert(e);
                }
            }
            else if ( (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4
                   || (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10 )
            {
                for (unsigned int i = 0; i < 6; ++i)
                {
                    n1 = (*it)->GetGeometry()[msT4Edges[i][0]].Id();
                    n2 = (*it)->GetGeometry()[msT4Edges[i][1]].Id();
                    if (n2 > n1)
                    {
                        n = n1;
                        n1 = n2;
                        n2 = n;
                    }
                    edge_t e(n1, n2);
                    edges.insert(e);
                }
            }
            else if ( (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8
                   || (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20
                   || (*it)->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D27 )
            {
                for (unsigned int i = 0; i < 12; ++i)
                {
                    n1 = (*it)->GetGeometry()[msH8Edges[i][0]].Id();
                    n2 = (*it)->GetGeometry()[msH8Edges[i][1]].Id();
                    if (n2 > n1)
                    {
                        n = n1;
                        n1 = n2;
                        n2 = n;
                    }
                    edge_t e(n1, n2);
                    edges.insert(e);
                }
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "Invalid geometry", static_cast<int>((*it)->GetGeometry().GetGeometryType()))
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
    ModelPartUtilities& operator=(ModelPartUtilities const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    ModelPartUtilities(ModelPartUtilities const& rOther)
    {
    }

    ///@}

}; // Class ModelPartUtilities

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, ModelPartUtilities& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const ModelPartUtilities& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_MODEL_PART_UTILITY_H_INCLUDED

