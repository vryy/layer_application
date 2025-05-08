//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Hoang-Giang Bui $
//   Date:                $Date: 23 Jun 2023 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_GID_MFEM_MESH_CONTAINER_H_INCLUDED)
#define  KRATOS_GID_MFEM_MESH_CONTAINER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_io/gid_sd_mesh_container.h"


namespace Kratos
{

/**
 * Auxiliary class to store meshes of different element types and to
 * write these meshes to an output file
 */
template<class TModelPartType>
class GidMfemMeshContainer : public GidSDMeshContainer<TModelPartType>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(GidMfemMeshContainer);

    typedef GidSDMeshContainer<TModelPartType> BaseType;

    ///Constructor
    GidMfemMeshContainer ( GeometryData::KratosGeometryType geometryType, std::string mesh_title )
    : BaseType(geometryType, mesh_title)
    {
    }

    void WriteMesh ( GiD_FILE MeshFile, bool deformed ) override
    {
        KRATOS_TRY

        std::vector<int> ns = GetNodalSequence(BaseType::mGeometryType);

        bool nodes_written = false;
        if ( BaseType::mMeshElements.size() != 0 )
        {
            //compute number of layers
            int max_id = 0;
            for ( typename TModelPartType::ElementsContainerType::iterator it = BaseType::mMeshElements.begin();
                    it != BaseType::mMeshElements.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                if (max_id < prop_id) max_id = prop_id;
            }
            if (max_id > 10000)
                std::cout<< "a property Id > 10000 found. Are u sure you need so many properties?" << std::endl;
            std::vector<int> elements_per_layer (max_id+1,0);

            //fill layer list
            for ( typename TModelPartType::ElementsContainerType::iterator it = BaseType::mMeshElements.begin();
                    it != BaseType::mMeshElements.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                elements_per_layer[prop_id] += 1;
            }
            //std::cout << "start printing elements" <<std::endl;
            for (unsigned int current_layer = 0; current_layer < elements_per_layer.size(); current_layer++)
            {
                if (elements_per_layer[current_layer] > 0)
                {
                    //create an appropiate name
                    std::stringstream current_layer_name (std::stringstream::in | std::stringstream::out);
                    current_layer_name << BaseType::mMeshTitle << "_" << current_layer;
                    for ( typename TModelPartType::ElementsContainerType::iterator it = BaseType::mMeshElements.begin();
                            it != BaseType::mMeshElements.end(); ++it )
                    {
                        if ( it->GetProperties().Id() == current_layer )
                        {
                            if (it->GetProperties().Has(LAYER_NAME))
                                current_layer_name << "_" << it->GetProperties().GetValue(LAYER_NAME);
                            break;
                        }
                    }

                    //begin mesh
                    if ( BaseType::mMeshElements.begin()->GetGeometry().WorkingSpaceDimension() == 2 )
                    {
                        //std::cout << " -print element 2D mesh: layer ["<<current_layer<<"]-"<<std::endl;
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_2D, BaseType::mGidElementType, BaseType::mMeshElements.begin()->GetGeometry().size() );
                    }
                    else if ( BaseType::mMeshElements.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        //std::cout << " -print element 3D mesh: layer ["<<current_layer<<"]-"<<std::endl;
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_3D, BaseType::mGidElementType, BaseType::mMeshElements.begin()->GetGeometry().size() );
                    }
                    else
                        KRATOS_THROW_ERROR (std::logic_error,"check working space dimension of model","");

                    //printing nodes
                    if(nodes_written == false)
                    {
                        GiD_fBeginCoordinates(MeshFile);
                        for ( typename TModelPartType::NodesContainerType::iterator it = BaseType::mMeshNodes.begin();
                                it != BaseType::mMeshNodes.end(); ++it )
                        {
                            if ( deformed )
                                GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X(),
                                                       (it)->Y(), (it)->Z() );
                            else
                                GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X0(),
                                                       (it)->Y0(), (it)->Z0() );
                            /// DEBUGGING
                            // std::cout << "node " << (it)->Id() << " is written to GiD" << std::endl;
                            /// END DEBUGGING
                        }
                        GiD_fEndCoordinates(MeshFile);

                        nodes_written = true;
                    }

                    //printing elements
                    GiD_fBeginElements(MeshFile);
                    int* nodes_id = new int[BaseType::mMeshElements.begin()->GetGeometry().size() + 1];
                    for ( typename TModelPartType::ElementsContainerType::iterator it = BaseType::mMeshElements.begin();
                            it != BaseType::mMeshElements.end(); ++it )
                    {
                        for ( unsigned int i=0; i< (it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry()[ns[i]].Id();

                        if ( BaseType::mGeometryType == GeometryData::KratosGeometryType::Kratos_Line2D3
                                || BaseType::mGeometryType == GeometryData::KratosGeometryType::Kratos_Line3D3 )
                        {
                            nodes_id[0] = (it)->GetGeometry()[0].Id();
                            nodes_id[1] = (it)->GetGeometry()[2].Id();
                            nodes_id[2] = (it)->GetGeometry()[1].Id();
                        }
                        // nodes_id[(it)->GetGeometry().size()] = (it)->GetProperties().Id()+1;
                        nodes_id[(it)->GetGeometry().size()] = (it)->GetProperties().Id();

                        bool element_is_active = true;
                        if( it->IsDefined( ACTIVE ) )
                            element_is_active = it->Is(ACTIVE);
                        if ( element_is_active )
                        {
                            if ((it)->GetProperties().Id()==current_layer)
                            {
                                GiD_fWriteElementMat ( MeshFile, (it)->Id(), nodes_id );
                                /// DEBUGGING
                                // std::cout << "element " << (it)->Id() << " is written to GiD" << std::endl;
                                /// END DEBUGGING
                            }
                        }
                    }

                    delete [] nodes_id;
                    GiD_fEndElements(MeshFile);
                    GiD_fEndMesh(MeshFile);
                }
            }
            //std::cout << "end printing elements" <<std::endl;
        }

        if ( BaseType::mMeshConditions.size() != 0 )
        {
            //compute number of layers
            int max_id = 0;
            for ( typename TModelPartType::ConditionsContainerType::iterator it = BaseType::mMeshConditions.begin();
                    it != BaseType::mMeshConditions.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                if (max_id < prop_id) max_id = prop_id;
            }
            if (max_id > 10000)
                std::cout<< "a property Id > 10000 found. Are u sure you need so many properties?" << std::endl;
            std::vector<int> conditions_per_layer (max_id+1,0);
            //fill layer list
            for ( typename TModelPartType::ConditionsContainerType::iterator it = BaseType::mMeshConditions.begin();
                    it != BaseType::mMeshConditions.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                conditions_per_layer[prop_id] += 1;
            }
            //std::cout << "start printing conditions" <<std::endl;
            for (unsigned int current_layer = 0; current_layer < conditions_per_layer.size(); current_layer++)
            {
                if (conditions_per_layer[current_layer] > 0)
                {
                    // determine mesh name
                    std::stringstream current_layer_name (std::stringstream::in | std::stringstream::out);
                    current_layer_name << BaseType::mMeshTitle << "_" << current_layer;

                    for ( typename TModelPartType::ConditionsContainerType::iterator it = BaseType::mMeshConditions.begin(  );
                            it != BaseType::mMeshConditions.end(); ++it )
                    {
                        if ( it->GetProperties().Id() == current_layer )
                        {
                            if (it->GetProperties().Has(LAYER_NAME))
                                current_layer_name << "_" << it->GetProperties().GetValue(LAYER_NAME);
                            break;
                        }
                    }

                    // begin mesh
                    if ( BaseType::mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() == 2 )
                    {
                        //std::cout << " -print condition 2D mesh: layer ["<<current_layer<<"]-"<<std::endl;
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_2D, BaseType::mGidElementType,
                                        BaseType::mMeshConditions.begin()->GetGeometry().size() );
                    }
                    else if ( BaseType::mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        //std::cout << " -print condition 3D mesh: layer ["<<current_layer<<"]-"<<std::endl;
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_3D, BaseType::mGidElementType,
                                        BaseType::mMeshConditions.begin()->GetGeometry().size() );
                    }
                    else
                        KRATOS_THROW_ERROR (std::logic_error,"check working space dimension of model","");

                    //printing nodes
                    if(nodes_written == false)
                    {
                        GiD_fBeginCoordinates(MeshFile);
                        for ( typename TModelPartType::NodesContainerType::iterator it = BaseType::mMeshNodes.begin();
                                it != BaseType::mMeshNodes.end(); ++it )
                        {
                            if ( deformed )
                                GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X(),
                                                       (it)->Y(), (it)->Z() );
                            else
                                GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X0(),
                                                       (it)->Y0(), (it)->Z0() );
                        }
                        GiD_fEndCoordinates(MeshFile);
                        nodes_written = true;
                    }

                    //printing elements
                    GiD_fBeginElements(MeshFile);
                    int* nodes_id = new int[BaseType::mMeshConditions.begin()->GetGeometry().size() + 1];
                    for ( typename TModelPartType::ConditionsContainerType::iterator it = BaseType::mMeshConditions.begin(  );
                            it != BaseType::mMeshConditions.end(); ++it )
                    {
                        for ( unsigned int i=0; i< (it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry()[ns[i]].Id();
                        // nodes_id[ (it)->GetGeometry().size()]= (it)->GetProperties().Id()+1;
                        nodes_id[ (it)->GetGeometry().size()] = (it)->GetProperties().Id();

                        bool condition_is_active = true;
                        if( it->IsDefined( ACTIVE ) )
                            condition_is_active = it->Is(ACTIVE);
                        if ( condition_is_active )
                        {
                            if ((it)->GetProperties().Id()==current_layer)
                                GiD_fWriteElementMat ( MeshFile, (it)->Id(), nodes_id );
                        }
                    }

                    delete [] nodes_id;
                    GiD_fEndElements(MeshFile);
                    GiD_fEndMesh(MeshFile);
                }
            }
            //std::cout << "end printing conditions" <<std::endl;
        }

        KRATOS_CATCH ("")
    }

private:

    std::vector<int> GetNodalSequence(const GeometryData::KratosGeometryType& geometryType) const
    {
        switch (geometryType)
        {
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                return std::vector<int>{0, 1, 2};
            case GeometryData::KratosGeometryType::Kratos_Triangle2D6:
                return std::vector<int>{0, 1, 2, 3, 4, 5};
            case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4:
                return std::vector<int>{0, 1, 2, 3};
            case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8:
                return std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7};
            case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9:
                return std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8};
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                return std::vector<int>{0, 1, 2, 3};
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10:
                return std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
            case GeometryData::KratosGeometryType::Kratos_Hexahedra3D8:
                return std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7};
            case GeometryData::KratosGeometryType::Kratos_Hexahedra3D20:
                return std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                    12, 13, 14, 15, 16, 17, 18, 19};
            // mfem nodal sequence for Hex27 can be found in constructor of RefinedTriLinear3DFiniteElement, fe_fixed_order.cpp
            case GeometryData::KratosGeometryType::Kratos_Hexahedra3D27:
                return std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                    16, 17, 18, 19,
                    12, 13, 14, 15,
                    20, 21, 22, 23, 24, 25, 26};
            default:
                return std::vector<int>{};
        }
    }

}; // class GidMfemMeshContainer

} // namespace Kratos.

#endif // KRATOS_GID_MFEM_MESH_CONTAINER_H_INCLUDED defined
