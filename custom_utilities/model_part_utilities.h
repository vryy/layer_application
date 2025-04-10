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

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "custom_utilities/gidpost_reader.h"


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

/*** Utilities for various operations on Kratos::ModelPart
 */
class ModelPartUtilities
{
public:

    ///@name Type Definitions
    ///@{

    typedef ModelPart::IndexType IndexType;
    typedef ModelPart::SizeType SizeType;
    typedef ModelPart::ElementType::GeometryType GeometryType;

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
        const bool with_id, const bool with_x, const bool with_y, const bool with_z,
        const std::string& separator, const int precision);

    /// Export the nodal coordinates to output stream
    static void ExportNodalCoordinatesToGiD(std::ostream& rOStream, ModelPart& r_model_part,
        const int precision);

    /// Export the edge information to output stream
    /// REMARK: for second order element, only the 2-node edge on the side is returned
    static void ExportEdgeInformation(std::ostream& rOStream, ModelPart::ElementsContainerType& rpElements,
        const std::string& separator);

    /// Export the edge information to output stream
    /// REMARK: for second order element, only the 2-node edge on the side is returned
    static void ExportEdgeInformationToGiD(std::ostream& rOStream, ModelPart::ElementsContainerType& rpElements);

    /// Create a new entity out from a geometry
    template<class TEntityType>
    static typename TEntityType::Pointer CreateEntity(const std::string& sample_entity_name,
        const std::size_t& Id, Properties::Pointer pProperties, typename TEntityType::GeometryType& rGeometry)
    {
        if(!KratosComponents<TEntityType>::Has(sample_entity_name))
            KRATOS_ERROR << sample_entity_name << " is not registered to the KRATOS kernel";
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
            KRATOS_ERROR << sample_entity_name << " is not registered to the KRATOS kernel";
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
    static void CalculateLocalSystem(TEntityType& rEntity, const ProcessInfo& rCurrentProcessInfo, const int echo_level)
    {
        Matrix LHS;
        Vector RHS;

        rEntity.CalculateLocalSystem(LHS, RHS, rCurrentProcessInfo);

        if (echo_level > 0)
        {
            rEntity.PrintInfo(std::cout);
            std::cout << std::endl;
        }

        if (echo_level > 1)
        {
            KRATOS_WATCH(LHS)
            KRATOS_WATCH(RHS)
        }
    }

    /// Invoke the CalculateLocalSystem of an Element/Condition. It is useful when debugging them.
    template<class TEntityType>
    static void CalculateMassMatrix(TEntityType& rEntity, const ProcessInfo& rCurrentProcessInfo, const int echo_level)
    {
        Matrix MassMatrix;

        rEntity.CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

        if (echo_level > 0)
        {
            rEntity.PrintInfo(std::cout);
            std::cout << std::endl;
        }

        if (echo_level > 1)
        {
            KRATOS_WATCH(MassMatrix)
        }
    }

    /// Guess the integration order, providing the number of integration points
    static int GetIntegrationOrder(const GeometryType& rGeometry, const unsigned int npoints);

    /// Clear everything from the model_part, make it fresh again
    static void ClearModelPart(ModelPart& r_model_part);

    /// Fill the model_part based on GiD post result data reader
    /// A json parameters can be optionally provided to assist with search and read process. An example
    /// of the parameters looks like
    // """
    // {
    //     "Kratos_Hexahedra3D27_Mesh_1":
    //     {
    //         "type": "element",
    //         "name": "KinematicLinear3D27N",
    //         "prop_id": 1
    //     },
    //     "echo_level": 1
    // }
    // """
    static void GiDPost2ModelPart(GiDPostReader& reader, ModelPart& r_model_part, const Parameters& mesh_info, VariablesList* pElementalVariablesList = nullptr);

    /// Fill the model_part based on GiD binary data
    /// A json parameters can be optionally provided to assist with search and read process. An example
    /// of the parameters looks like
    // """
    // {
    //     "Kratos_Hexahedra3D27_Mesh_1":
    //     {
    //         "type": "element",
    //         "name": "KinematicLinear3D27N",
    //         "prop_id": 1
    //     },
    //     "echo_level": 1
    // }
    // """
    static void GiDPostBin2ModelPart(const std::string& fileName, ModelPart& r_model_part, const Parameters& mesh_info, VariablesList* pElementalVariablesList = nullptr);

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
                KRATOS_ERROR << "Invalid geometry " << (*it)->GetGeometry().GetGeometryType();
        }
    }

    static int ExtractPropertiesId(const std::string& Name);

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
