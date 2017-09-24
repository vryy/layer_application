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

// reference: vtkCellType.h
enum VTK_ElementType {
    VTK_Empty_Cell = 0,
    VTK_Vertex = 1,
    VTK_Poly_Vertex = 2,
    VTK_Line = 3, 
    VTK_Poly_Line = 4,
    VTK_Triangle = 5,
    VTK_Triangle_Strip = 6,
    VTK_Polygon = 7, 
    VTK_Pixel = 8,
    VTK_Quad = 9,
    VTK_Tetra = 10,
    VTK_Voxel = 11, 
    VTK_Hexahedron = 12,
    VTK_Wedge = 13,
    VTK_Pyramid = 14,
    VTK_Pentagonal_Prism = 15, 
    VTK_Hexagonal_Prism = 16,
    VTK_Quadratic_Edge = 21,
    VTK_Quadratic_Triangle = 22,
    VTK_Quadratic_Quad = 23, 
    VTK_Quadratic_Polygon = 36,
    VTK_Quadratic_Tetra = 24,
    VTK_Quadratic_Hexahedron = 25,
    VTK_Quadratic_Wedge = 26, 
    VTK_Quadratic_Pyramid = 27,
    VTK_Biquadratic_Quad = 28,
    VTK_Triquadratic_Hexahedron = 29,
    VTK_Quadratic_Linear_Quad = 30, 
    VTK_Quadratic_Linear_Wedge = 31,
    VTK_Biquadratic_Quadratic_Wedge = 32,
    VTK_Biquadratic_Quadratic_Hexahedron = 33,
    VTK_Biquadratic_Triangle = 34, 
    VTK_Cubix_Line = 35,
    VTK_Convex_Point_Set = 41,
    VTK_Polyhedron = 42,
    VTK_Parametric_Curve = 51, 
    VTK_Parametric_Surface = 52,
    VTK_Parametric_Tri_Surface = 53,
    VTK_Parametric_Quad_Surface = 54,
    VTK_Parametric_Tetra_Region = 55, 
    VTK_Parametric_Hex_Region = 56,
    VTK_Higher_Order_Edge = 60,
    VTK_Higher_Order_Triangle = 61,
    VTK_Higher_Order_Quad = 62, 
    VTK_Higher_Order_Polygon = 63,
    VTK_Higher_Order_Tetrahedron = 64,
    VTK_Higher_Order_Wedge = 65,
    VTK_Higher_Order_Pyramid = 66, 
    VTK_Higher_Order_Hexahedron = 67,
    VTK_Number_Of_Cell_Types 
};

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
        if ( pElemIt->GetGeometry().GetGeometryType() == mGeometryType )
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
        KRATOS_CATCH ("")
    }

    bool AddCondition (const ModelPart::ConditionsContainerType::iterator pCondIt)
    {
        KRATOS_TRY
        if ( pCondIt->GetGeometry().GetGeometryType() == mGeometryType )
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

        if (mode == VTK_PostAscii)
        {
            int* nodes_id = new int[nodes_size];
            for ( typename TElementsContainerType::iterator it = MeshElements.begin();
                    it != MeshElements.end(); ++it )
            {
                for ( unsigned int i=0; i< (it)->GetGeometry().size(); i++ )
                    nodes_id[i] = NodeIdMap[(it)->GetGeometry()[i].Id()];

                if ( it->Has ( IS_INACTIVE ) )
                {
                    if ( ! it->GetValue ( IS_INACTIVE ) )
                    {
                        VTK_fWriteElementConnectivity ( MeshFile, (it)->Id(), nodes_id, nodes_size );
                    }
                }
                else
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
                if ( it->Has ( IS_INACTIVE ) )
                {
                    if ( ! it->GetValue ( IS_INACTIVE ) )
                    {
                        for ( unsigned int i=0; i< (it)->GetGeometry().size(); i++ )
                            data_list.push_back(NodeIdMap[(it)->GetGeometry()[i].Id()]);
                    }
                }
                else
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
                if ( it->Has ( IS_INACTIVE ) )
                {
                    if ( ! it->GetValue ( IS_INACTIVE ) )
                    {
                        offset += nodes_size;
                        VTK_fWriteElementOffset ( MeshFile, offset );
                    }
                }
                else
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
                if ( it->Has ( IS_INACTIVE ) )
                {
                    if ( ! it->GetValue ( IS_INACTIVE ) )
                    {
                        offset += nodes_size;
                        data_list.push_back(offset);
                    }
                }
                else
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
                if ( it->Has ( IS_INACTIVE ) )
                {
                    if ( ! it->GetValue ( IS_INACTIVE ) )
                    {
                        VTK_fWriteElementType ( MeshFile, mVtkElementType );
                    }
                }
                else
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
                if ( it->Has ( IS_INACTIVE ) )
                {
                    if ( ! it->GetValue ( IS_INACTIVE ) )
                    {
                        data_list.push_back(mVtkElementType);
                    }
                }
                else
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

