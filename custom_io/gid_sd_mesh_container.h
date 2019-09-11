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
//   Date:                $Date: 14 Jun 2016 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_GID_SD_MESH_CONTAINER_H_INCLUDED)
#define  KRATOS_GID_SD_MESH_CONTAINER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstddef>

// External includes
#include "gidpost/source/gidpost.h"

// Project includes
#include "includes/define.h"
#include "geometries/geometry_data.h"
#include "includes/deprecated_variables.h"
#include "layer_application/layer_application.h"


namespace Kratos
{

/**
 * Type definitions
 */
typedef ModelPart::ElementsContainerType ElementsArrayType;
typedef ModelPart::NodesContainerType NodesArrayType;
typedef ModelPart::ConditionsContainerType ConditionsArrayType;
typedef GeometryData::IntegrationMethod IntegrationMethodType;
typedef GeometryData::KratosGeometryFamily KratosGeometryFamily;

/**
 * Auxiliary class to store meshes of different element types and to
 * write these meshes to an output file
 */
class GidSDMeshContainer
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(GidSDMeshContainer);

    ///Constructor
    GidSDMeshContainer ( GeometryData::KratosGeometryType geometryType, std::string mesh_title )
    : mMeshTitle(mesh_title)
    {
        mGeometryType = geometryType;

        if(     mGeometryType == GeometryData::Kratos_Hexahedra3D8
            ||  mGeometryType == GeometryData::Kratos_Hexahedra3D20
            ||  mGeometryType == GeometryData::Kratos_Hexahedra3D27
        )
        {
            mGidElementType = GiD_Hexahedra;
        }
        else if(mGeometryType == GeometryData::Kratos_Tetrahedra3D4
            ||  mGeometryType == GeometryData::Kratos_Tetrahedra3D10
        )
        {
            mGidElementType = GiD_Tetrahedra;
        }
        else if(mGeometryType == GeometryData::Kratos_Prism3D6
            ||  mGeometryType == GeometryData::Kratos_Prism3D15
        )
        {
            mGidElementType = GiD_Prism;
        }
        else if(mGeometryType == GeometryData::Kratos_Quadrilateral3D4
            ||  mGeometryType == GeometryData::Kratos_Quadrilateral3D8
            ||  mGeometryType == GeometryData::Kratos_Quadrilateral3D9
            ||  mGeometryType == GeometryData::Kratos_Quadrilateral2D4
            ||  mGeometryType == GeometryData::Kratos_Quadrilateral2D8
            ||  mGeometryType == GeometryData::Kratos_Quadrilateral2D9
        )
        {
            mGidElementType = GiD_Quadrilateral;
        }
        else if(mGeometryType == GeometryData::Kratos_Triangle3D3
            ||  mGeometryType == GeometryData::Kratos_Triangle3D6
            ||  mGeometryType == GeometryData::Kratos_Triangle2D3
            ||  mGeometryType == GeometryData::Kratos_Triangle2D6
        )
        {
            mGidElementType = GiD_Triangle;
        }
        else if(mGeometryType == GeometryData::Kratos_Line3D3
            ||  mGeometryType == GeometryData::Kratos_Line3D2
            ||  mGeometryType == GeometryData::Kratos_Line2D3
            ||  mGeometryType == GeometryData::Kratos_Line2D2
        )
        {
            mGidElementType = GiD_Linear;
        }
        else if(mGeometryType == GeometryData::Kratos_Point3D
            ||  mGeometryType == GeometryData::Kratos_Point2D
        )
        {
            mGidElementType = GiD_Point;
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error, "Unknown geometry type", mGeometryType)
        }
    }

    bool AddElement ( const ModelPart::ElementsContainerType::iterator pElemIt )
    {
        KRATOS_TRY

        if ( pElemIt->GetGeometry().GetGeometryType() == mGeometryType )
        {
//            std::cout << "element " << pElemIt->Id() << " of geometryType "
//                      << mGeometryType << " is added to " << mMeshTitle << std::endl;
            mMeshElements.push_back ( * (pElemIt.base() ) );
            Geometry<Node<3> >&geom = pElemIt->GetGeometry();
            for ( Element::GeometryType::iterator it = geom.begin(); it != geom.end(); it++)
            {
                mMeshNodes.push_back ( * (it.base() ) );
            }
            return true;
        }
        else
            return false;

        KRATOS_CATCH ("")
    }

    bool AddCondition ( const ModelPart::ConditionsContainerType::iterator pCondIt )
    {
        KRATOS_TRY

        if ( pCondIt->GetGeometry().GetGeometryType() == mGeometryType )
        {
//            std::cout << "condition " << pCondIt->Id() << " of geometryType "
//                      << mGeometryType << " is added to " << mMeshTitle << std::endl;
            mMeshConditions.push_back ( * (pCondIt.base() ) );
            Geometry<Node<3> >&geom = pCondIt->GetGeometry();
            for ( Condition::GeometryType::iterator it = geom.begin(); it != geom.end(); it++)
            {
                mMeshNodes.push_back ( * (it.base() ) );
            }
            return true;
        }
        else
            return false;

        KRATOS_CATCH ("")
    }

    void FinalizeMeshCreation()
    {
        if ( mMeshElements.size() != 0 )
        {
            mMeshNodes.Unique();
        }
        if ( mMeshConditions.size() != 0 )
        {
            mMeshNodes.Unique();
        }
    }

    void WriteMesh ( GiD_FILE MeshFile, bool deformed )
    {
        KRATOS_TRY

        bool nodes_written = false;
        if ( mMeshElements.size() != 0 )
        {
            //compute number of layers
            int max_id = 0;
            for ( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                    it != mMeshElements.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                if (max_id < prop_id) max_id = prop_id;
            }
            if (max_id > 10000)
                std::cout<< "a property Id > 10000 found. Are u sure you need so many properties?" << std::endl;
            std::vector<int> elements_per_layer (max_id+1,0);

            //fill layer list
            for ( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                    it != mMeshElements.end(); ++it )
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
                    current_layer_name << mMeshTitle << "_" << current_layer;
                    for ( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                            it != mMeshElements.end(); ++it )
                    {
                        if ( it->GetProperties().Id() == current_layer )
                        {
                            if (it->GetProperties().Has(LAYER_NAME))
                                current_layer_name << "_" << it->GetProperties().GetValue(LAYER_NAME);
                            break;
                        }
                    }

                    //begin mesh
                    if ( mMeshElements.begin()->GetGeometry().WorkingSpaceDimension() == 2 )
                    {
                        //std::cout << " -print element 2D mesh: layer ["<<current_layer<<"]-"<<std::endl;
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_2D, mGidElementType,mMeshElements.begin()->GetGeometry().size() );
                    }
                    else if ( mMeshElements.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        //std::cout << " -print element 3D mesh: layer ["<<current_layer<<"]-"<<std::endl;
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_3D, mGidElementType,mMeshElements.begin()->GetGeometry().size() );
                    }
                    else
                        KRATOS_THROW_ERROR (std::logic_error,"check working space dimension of model","");

                    //printing nodes
                    if(nodes_written == false)
                    {
                        GiD_fBeginCoordinates(MeshFile);
                        for ( ModelPart::NodesContainerType::iterator it = mMeshNodes.begin();
                                it != mMeshNodes.end(); ++it )
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
                    int* nodes_id = new int[mMeshElements.begin()->GetGeometry().size() + 1];
                    for ( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                            it != mMeshElements.end(); ++it )
                    {
                        for ( unsigned int i=0; i< (it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry() [i].Id();

                        if ( mGeometryType == GeometryData::Kratos_Line2D3
                                || mGeometryType == GeometryData::Kratos_Line3D3 )
                        {
                            nodes_id[0] = (it)->GetGeometry() [0].Id();
                            nodes_id[1] = (it)->GetGeometry() [2].Id();
                            nodes_id[2] = (it)->GetGeometry() [1].Id();
                        }
                        nodes_id[(it)->GetGeometry().size()] = (it)->GetProperties().Id()+1;

                        bool element_is_active = true;
                        if( it->IsDefined( ACTIVE ) )
                            element_is_active = it->Is(ACTIVE);
                        if ( element_is_active )
                        {
                            if ((it)->GetProperties().Id()==current_layer)
                                GiD_fWriteElementMat ( MeshFile, (it)->Id(), nodes_id);
                        }
                    }

                    delete [] nodes_id;
                    GiD_fEndElements(MeshFile);
                    GiD_fEndMesh(MeshFile);
                }
            }
            //std::cout << "end printing elements" <<std::endl;
        }

        if ( mMeshConditions.size() != 0 )
        {
            //compute number of layers
            int max_id = 0;
            for ( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                    it != mMeshConditions.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                if (max_id < prop_id) max_id = prop_id;
            }
            if (max_id > 10000)
                std::cout<< "a property Id > 10000 found. Are u sure you need so many properties?" << std::endl;
            std::vector<int> conditions_per_layer (max_id+1,0);
            //fill layer list
            for ( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                    it != mMeshConditions.end(); ++it )
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
                    current_layer_name << mMeshTitle << "_" << current_layer;

                    for ( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin(  );
                            it != mMeshConditions.end(); ++it )
                    {
                        if ( it->GetProperties().Id() == current_layer )
                        {
                            if (it->GetProperties().Has(LAYER_NAME))
                                current_layer_name << "_" << it->GetProperties().GetValue(LAYER_NAME);
                            break;
                        }
                    }

                    // begin mesh
                    if ( mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() == 2 )
                    {
                        //std::cout << " -print condition 2D mesh: layer ["<<current_layer<<"]-"<<std::endl;
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_2D, mGidElementType,
                                        mMeshConditions.begin()->GetGeometry().size() );
                    }
                    else if ( mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        //std::cout << " -print condition 3D mesh: layer ["<<current_layer<<"]-"<<std::endl;
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_3D, mGidElementType,
                                        mMeshConditions.begin()->GetGeometry().size() );
                    }
                    else
                        KRATOS_THROW_ERROR (std::logic_error,"check working space dimension of model","");

                    //printing nodes
                    if(nodes_written == false)
                    {
                        GiD_fBeginCoordinates(MeshFile);
                        for ( ModelPart::NodesContainerType::iterator it = mMeshNodes.begin();
                                it != mMeshNodes.end(); ++it )
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
                    int* nodes_id = new int[mMeshConditions.begin()->GetGeometry().size() + 1];
                    for ( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin(  );
                            it != mMeshConditions.end(); ++it )
                    {
                        for ( unsigned int i=0; i< (it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry() [i].Id();
                        nodes_id[ (it)->GetGeometry().size()]= (it)->GetProperties().Id()+1;

                        bool condition_is_active = true;
                        if( it->IsDefined( ACTIVE ) )
                            condition_is_active = it->Is(ACTIVE);
                        if ( condition_is_active )
                        {
                            if ((it)->GetProperties().Id()==current_layer)
                                GiD_fWriteElementMat ( MeshFile, (it)->Id(), nodes_id);
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

    void Reset()
    {
        mMeshNodes.clear();
        mMeshElements.clear();
        mMeshConditions.clear();
    }

    ModelPart::NodesContainerType GetMeshNodes()
    {
        return mMeshNodes;
    }

protected:
    ///member variables
    GeometryData::KratosGeometryType mGeometryType;
    GiD_ElementType mGidElementType;
    ModelPart::NodesContainerType mMeshNodes;
    ModelPart::ElementsContainerType mMeshElements;
    ModelPart::ConditionsContainerType mMeshConditions;
    std::string mMeshTitle;
};//class GidSDMeshContainer

}// namespace Kratos.

#endif // KRATOS_GID_SD_MESH_CONTAINER_H_INCLUDED defined

