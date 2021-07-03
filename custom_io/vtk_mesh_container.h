//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 layer_application/LICENSE.txt
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//

#if !defined(KRATOS_VTK_MESH_CONTAINER_H_INCLUDED)
#define  KRATOS_VTK_MESH_CONTAINER_H_INCLUDED
// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstddef>

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry_data.h"
#include "includes/deprecated_variables.h"
#include "custom_utilities/vtk.h"

namespace Kratos
{

/**
 * Auxiliary class to store meshes of different element types and to
 * write these meshes to an output file
 */
class VtkMeshContainer
{
public:
    /**
     * Type definitions
     */
    typedef GeometryData::IntegrationMethod IntegrationMethodType;
    typedef GeometryData::KratosGeometryFamily KratosGeometryFamily;
    typedef Element::GeometryType GeometryType;
    typedef std::map<std::size_t, ModelPart::ElementsContainerType> MeshElementsContainerType;
    typedef std::map<std::size_t, ModelPart::ConditionsContainerType> MeshConditionsContainerType;
    typedef std::map<std::size_t, ModelPart::NodesContainerType> MeshNodesContainerType;

    ///Constructor
    VtkMeshContainer ( GeometryData::KratosGeometryType geometryType,
                       VTK_ElementType elementType, const char* mesh_title )
        : mGeometryType (geometryType), mVtkElementType (elementType), mMeshTitle (mesh_title) {}

    const std::string MeshTitle() const {return mMeshTitle;}

    bool AddElement ( const ModelPart::ElementsContainerType::iterator pElemIt )
    {
        KRATOS_TRY

        bool element_is_active;

        if ( pElemIt->GetGeometry().GetGeometryType() == mGeometryType )
        {
            element_is_active = true;
            if( pElemIt->IsDefined( ACTIVE ) )
                element_is_active = pElemIt->Is(ACTIVE);

            if (element_is_active)
            {
                mMeshElements[pElemIt->GetProperties().Id()].push_back ( * (pElemIt.base() ) );
                GeometryType& geom = pElemIt->GetGeometry();
                for ( Element::GeometryType::iterator it = geom.begin(); it != geom.end(); ++it)
                {
                    mMeshElementNodes[pElemIt->GetProperties().Id()].push_back ( * (it.base() ) );
                }
                return true;
            }
            else
                return false;
        }
        else
            return false;

        KRATOS_CATCH ("")
    }

    bool AddCondition (const ModelPart::ConditionsContainerType::iterator pCondIt)
    {
        KRATOS_TRY

        bool condition_is_active;

        if ( pCondIt->GetGeometry().GetGeometryType() == mGeometryType )
        {
            condition_is_active = true;
            if( pCondIt->IsDefined( ACTIVE ) )
                condition_is_active = pCondIt->Is(ACTIVE);

            if (condition_is_active)
            {
                mMeshConditions[pCondIt->GetProperties().Id()].push_back ( * (pCondIt.base() ) );
                GeometryType& geom = pCondIt->GetGeometry();
                for ( Condition::GeometryType::iterator it = geom.begin(); it != geom.end(); ++it)
                {
                    mMeshConditionNodes[pCondIt->GetProperties().Id()].push_back ( * (it.base() ) );
                }
                return true;
            }
            else
                return false;
        }
        else
            return false;

        KRATOS_CATCH ("")
    }

    void FinalizeMeshCreation()
    {
        for(std::map<std::size_t, ModelPart::ElementsContainerType>::iterator it = mMeshElements.begin();
                it != mMeshElements.end(); ++it)
        {
            if(it->second.size() != 0)
                mMeshElementNodes[it->first].Unique();
        }

        for(std::map<std::size_t, ModelPart::ConditionsContainerType>::iterator it = mMeshConditions.begin();
                it != mMeshConditions.end(); ++it)
        {
            if(it->second.size() != 0)
                mMeshConditionNodes[it->first].Unique();
        }
    }

    template<class TNodesContainerType, class TElementsContainerType>
    void WriteMesh(FILE* MeshFile, TNodesContainerType& MeshNodes, TElementsContainerType& MeshElements, bool deformed, VTK_PostMode mode)
    {
        // printing nodes
        std::map<std::size_t, std::size_t> NodeIdMap;
        std::size_t cnt = 0;

        VTK_fBeginCoordinates(MeshFile, mode);
        if (mode == VTK_PostAscii)
        {
            for ( typename TNodesContainerType::iterator it = MeshNodes.begin();
                    it != MeshNodes.end(); ++it )
            {
                if ( deformed )
                    VTK_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X(),
                                           (it)->Y(), (it)->Z() );
                else
                    VTK_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X0(),
                                           (it)->Y0(), (it)->Z0() );
                NodeIdMap[(it)->Id()] = cnt++;
            }
        }
        else if (mode == VTK_PostBinary)
        {
            std::vector<float> data_list;
            data_list.reserve(3*(MeshNodes.end() - MeshNodes.begin()));

            for ( typename TNodesContainerType::iterator it = MeshNodes.begin();
                    it != MeshNodes.end(); ++it )
            {
                if ( deformed )
                {
                    data_list.push_back((it)->X());
                    data_list.push_back((it)->Y());
                    data_list.push_back((it)->Z());
                }
                else
                {
                    data_list.push_back((it)->X0());
                    data_list.push_back((it)->Y0());
                    data_list.push_back((it)->Z0());
                }
                NodeIdMap[(it)->Id()] = cnt++;
            }

            float* tmp = (float*)(&data_list[0]);
            vtk_write_compressed ( MeshFile, (char*)tmp, sizeof(*tmp)*data_list.size());
        }
        VTK_fEndCoordinates(MeshFile);

        // printing elements
        VTK_fBeginElements(MeshFile);
        VTK_fBeginElementsConnectivity(MeshFile, mode);
        unsigned int nodes_size = MeshElements.begin()->GetGeometry().size();

        bool is_active;
        if (mode == VTK_PostAscii)
        {
            int* nodes_id = new int[nodes_size];
            for ( typename TElementsContainerType::iterator it = MeshElements.begin();
                    it != MeshElements.end(); ++it )
            {
                for ( unsigned int i=0; i< (it)->GetGeometry().size(); i++ )
                    nodes_id[i] = NodeIdMap[(it)->GetGeometry()[i].Id()];

                is_active = true;
                if(it->IsDefined ( ACTIVE ))
                {
                    is_active = it->Is( ACTIVE );
                #ifdef SD_APP_FORWARD_COMPATIBILITY
                }
                #else
                    if(it->Has ( IS_INACTIVE ))
                        is_active = is_active && (!it->GetValue( IS_INACTIVE ));
                }
                else
                {
                    if(it->Has ( IS_INACTIVE ))
                        is_active = !it->GetValue( IS_INACTIVE );
                }
                #endif

                if ( is_active )
                {
                    VTK_fWriteElementConnectivity ( MeshFile, (it)->Id(), nodes_id, nodes_size );
                }
            }
            delete [] nodes_id;
        }
        else if (mode == VTK_PostBinary)
        {
            std::vector<int> data_list;
            data_list.reserve(nodes_size*(MeshElements.end() - MeshElements.begin()));
            for ( typename TElementsContainerType::iterator it = MeshElements.begin();
                    it != MeshElements.end(); ++it )
            {
                is_active = true;
                if(it->IsDefined ( ACTIVE ))
                {
                    is_active = it->Is( ACTIVE );
                #ifdef SD_APP_FORWARD_COMPATIBILITY
                }
                #else
                    if(it->Has ( IS_INACTIVE ))
                        is_active = is_active && (!it->GetValue( IS_INACTIVE ));
                }
                else
                {
                    if(it->Has ( IS_INACTIVE ))
                        is_active = !it->GetValue( IS_INACTIVE );
                }
                #endif

                if ( is_active )
                {
                    for ( unsigned int i=0; i< (it)->GetGeometry().size(); i++ )
                        data_list.push_back(NodeIdMap[(it)->GetGeometry()[i].Id()]);
                }
            }
            int* tmp = (int*)(&data_list[0]);
            vtk_write_compressed(MeshFile, (char*)tmp, sizeof(*tmp)*data_list.size());
        }
        VTK_fEndElementsConnectivity(MeshFile);

        VTK_fBeginElementsOffsets(MeshFile, mode);
        std::size_t offset = 0;

        if (mode == VTK_PostAscii)
        {
            for ( typename TElementsContainerType::iterator it = MeshElements.begin();
                    it != MeshElements.end(); ++it )
            {
                is_active = true;
                if(it->IsDefined ( ACTIVE ))
                {
                    is_active = it->Is( ACTIVE );
                #ifdef SD_APP_FORWARD_COMPATIBILITY
                }
                #else
                    if(it->Has ( IS_INACTIVE ))
                        is_active = is_active && (!it->GetValue( IS_INACTIVE ));
                }
                else
                {
                    if(it->Has ( IS_INACTIVE ))
                        is_active = !it->GetValue( IS_INACTIVE );
                }
                #endif

                if ( is_active )
                {
                    offset += nodes_size;
                    VTK_fWriteElementOffset ( MeshFile, offset );
                }
            }
        }
        else if (mode == VTK_PostBinary)
        {
            std::vector<int> data_list;
            data_list.reserve(MeshElements.end() - MeshElements.begin());
            for ( typename TElementsContainerType::iterator it = MeshElements.begin();
                    it != MeshElements.end(); ++it )
            {
                is_active = true;
                if(it->IsDefined ( ACTIVE ))
                {
                    is_active = it->Is( ACTIVE );
                #ifdef SD_APP_FORWARD_COMPATIBILITY
                }
                #else
                    if(it->Has ( IS_INACTIVE ))
                        is_active = is_active && (!it->GetValue( IS_INACTIVE ));
                }
                else
                {
                    if(it->Has ( IS_INACTIVE ))
                        is_active = !it->GetValue( IS_INACTIVE );
                }
                #endif

                if ( is_active )
                {
                    offset += nodes_size;
                    data_list.push_back(offset);
                }
            }
            int* tmp = (int*)(&data_list[0]);
            vtk_write_compressed(MeshFile, (char*)tmp, sizeof(*tmp)*data_list.size());
        }
        VTK_fEndElementsOffsets(MeshFile);

        VTK_fBeginElementsTypes(MeshFile, mode);
        if (mode == VTK_PostAscii)
        {
            for ( typename TElementsContainerType::iterator it = MeshElements.begin();
                    it != MeshElements.end(); ++it )
            {
                is_active = true;
                if(it->IsDefined ( ACTIVE ))
                {
                    is_active = it->Is( ACTIVE );
                #ifdef SD_APP_FORWARD_COMPATIBILITY
                }
                #else
                    if(it->Has ( IS_INACTIVE ))
                        is_active = is_active && (!it->GetValue( IS_INACTIVE ));
                }
                else
                {
                    if(it->Has ( IS_INACTIVE ))
                        is_active = !it->GetValue( IS_INACTIVE );
                }
                #endif

                if ( is_active )
                {
                    VTK_fWriteElementType ( MeshFile, mVtkElementType );
                }
            }
        }
        else if (mode == VTK_PostBinary)
        {
            std::vector<uint8_t> data_list;
            data_list.reserve(MeshElements.end() - MeshElements.begin());
            for ( typename TElementsContainerType::iterator it = MeshElements.begin();
                    it != MeshElements.end(); ++it )
            {
                is_active = true;
                if(it->IsDefined ( ACTIVE ))
                {
                    is_active = it->Is( ACTIVE );
                #ifdef SD_APP_FORWARD_COMPATIBILITY
                }
                #else
                    if(it->Has ( IS_INACTIVE ))
                        is_active = is_active && (!it->GetValue( IS_INACTIVE ));
                }
                else
                {
                    if(it->Has ( IS_INACTIVE ))
                        is_active = !it->GetValue( IS_INACTIVE );
                }
                #endif

                if ( is_active )
                {
                    data_list.push_back(mVtkElementType);
                }
            }
            uint8_t* tmp = (uint8_t*)(&data_list[0]);
            vtk_write_compressed(MeshFile, (char*)tmp, sizeof(*tmp)*data_list.size());
        }
        VTK_fEndElementsTypes(MeshFile);

        VTK_fEndElements(MeshFile);
    }

    void Reset()
    {
        mMeshElementNodes.clear();
        mMeshConditionNodes.clear();
        mMeshElements.clear();
        mMeshConditions.clear();
    }

    const std::string GetMeshElementsName(const std::size_t& prop_id)
    {
        std::stringstream ss;
        ss << mMeshTitle << "_element_" << prop_id;
        return ss.str();
    }

    const std::string GetMeshConditionsName(const std::size_t& prop_id)
    {
        std::stringstream ss;
        ss << mMeshTitle << "_condition_" << prop_id;
        return ss.str();
    }

    ModelPart::NodesContainerType& GetMeshElementNodes(const std::size_t& prop_id)
    {
        return mMeshElementNodes[prop_id];
    }

    ModelPart::NodesContainerType& GetMeshConditionNodes(const std::size_t& prop_id)
    {
        return mMeshConditionNodes[prop_id];
    }

    MeshElementsContainerType& GetMeshElements()
    {
        return mMeshElements;
    }

    MeshConditionsContainerType& GetMeshConditions()
    {
        return mMeshConditions;
    }

protected:
    ///member variables
    GeometryData::KratosGeometryType mGeometryType;
    VTK_ElementType mVtkElementType;

    // here we sort out the mesh based on Properties Id and element/condition. In sort, different element/condition type with
    // different Properties Id will be contained as a separate mesh, and will be export as a Piece in Paraview's Vtk format.
    // each mesh will have respective nodes container.

    MeshNodesContainerType mMeshElementNodes;

    MeshNodesContainerType mMeshConditionNodes;

    MeshElementsContainerType mMeshElements;

    MeshConditionsContainerType mMeshConditions;

    const char* mMeshTitle;
};//class VtkMeshContainer
}// namespace Kratos.
#endif // KRATOS_VTK_MESH_CONTAINER_H_INCLUDED defined

