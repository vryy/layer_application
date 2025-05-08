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
#ifndef SD_APP_FORWARD_COMPATIBILITY
#include "includes/deprecated_variables.h"
#endif
#include "layer_application_variables.h"


namespace Kratos
{

/**
 * Auxiliary class to store meshes of different element types and to
 * write these meshes to an output file
 * In case the model part is ComplexModelPart, the second option will determine real (1) or imaginary (2)
 * component used for the coordinate
 */
template<class TModelPartType, int type = 1>
class GidSDMeshContainer
{
public:

    /// Type definitions
    typedef TModelPartType ModelPartType;
    typedef typename TModelPartType::CoordinateType CoordinateType;
    typedef typename TModelPartType::ElementType::GeometryType GeometryType;

    static constexpr int complex_coordinate_type = type;

    KRATOS_CLASS_POINTER_DEFINITION(GidSDMeshContainer);

    ///Constructor
    GidSDMeshContainer ( GeometryData::KratosGeometryType geometryType, std::string mesh_title )
    : mMeshTitle(mesh_title)
    {
        mGeometryType = geometryType;

        if(     mGeometryType == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Hexahedra3D27
        )
        {
            mGidElementType = GiD_Hexahedra;
        }
        else if(mGeometryType == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10
        )
        {
            mGidElementType = GiD_Tetrahedra;
        }
        else if(mGeometryType == GeometryData::KratosGeometryType::Kratos_Prism3D6
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Prism3D15
        )
        {
            mGidElementType = GiD_Prism;
        }
        else if(mGeometryType == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9
        )
        {
            mGidElementType = GiD_Quadrilateral;
        }
        else if(mGeometryType == GeometryData::KratosGeometryType::Kratos_Triangle3D3
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Triangle3D6
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Triangle2D3
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Triangle2D6
        )
        {
            mGidElementType = GiD_Triangle;
        }
        else if(mGeometryType == GeometryData::KratosGeometryType::Kratos_Line3D3
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Line3D2
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Line2D3
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Line2D2
        )
        {
            mGidElementType = GiD_Linear;
        }
        else if(mGeometryType == GeometryData::KratosGeometryType::Kratos_Point3D
            ||  mGeometryType == GeometryData::KratosGeometryType::Kratos_Point2D
        )
        {
            mGidElementType = GiD_Point;
        }
        else
        {
            KRATOS_ERROR << "Unknown geometry type " << static_cast<int>(mGeometryType);
        }
    }

    bool AddElement ( const typename TModelPartType::ElementsContainerType::const_iterator pElemIt )
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
                mMeshElements.push_back ( * (pElemIt.base() ) );
                GeometryType& geom = pElemIt->GetGeometry();
                for ( typename GeometryType::iterator it = geom.begin(); it != geom.end(); it++)
                {
                    mMeshNodes.push_back ( * (it.base() ) );
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

    bool AddCondition ( const typename TModelPartType::ConditionsContainerType::const_iterator pCondIt )
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
                mMeshConditions.push_back ( * (pCondIt.base() ) );
                GeometryType& geom = pCondIt->GetGeometry();
                for ( typename GeometryType::iterator it = geom.begin(); it != geom.end(); it++)
                {
                    mMeshNodes.push_back ( * (it.base() ) );
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
        if ( mMeshElements.size() != 0 )
        {
            mMeshNodes.Unique();
        }
        if ( mMeshConditions.size() != 0 )
        {
            mMeshNodes.Unique();
        }
    }

    virtual void WriteMesh ( GiD_FILE MeshFile, bool deformed )
    {
        KRATOS_TRY

        if ( mMeshElements.size() != 0 )
        {
            //compute number of layers
            int max_id = 0;
            for ( typename TModelPartType::ElementsContainerType::iterator it = mMeshElements.begin();
                    it != mMeshElements.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                if (max_id < prop_id) max_id = prop_id;
            }
            if (max_id > 10000)
                std::cout<< "a property Id > 10000 found. Are u sure you need so many properties?" << std::endl;
            std::vector<int> elements_per_layer(max_id+1, 0);

            //fill layer list
            for ( typename TModelPartType::ElementsContainerType::iterator it = mMeshElements.begin();
                    it != mMeshElements.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                elements_per_layer[prop_id] += 1;
            }
            //std::cout << "start printing elements" <<std::endl;
            for (unsigned int current_layer = 0; current_layer < elements_per_layer.size(); current_layer++)
            {
                bool nodes_written = false;
                if (elements_per_layer[current_layer] > 0)
                {
                    //create an appropriate name
                    std::stringstream current_layer_name (std::stringstream::in | std::stringstream::out);
                    current_layer_name << mMeshTitle << "_" << current_layer;
                    for ( typename TModelPartType::ElementsContainerType::iterator it = mMeshElements.begin();
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
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_2D, mGidElementType, mMeshElements.begin()->GetGeometry().size() );
                    }
                    else if ( mMeshElements.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_3D, mGidElementType, mMeshElements.begin()->GetGeometry().size() );
                    }
                    else
                        KRATOS_ERROR << "check working space dimension of model";

                    //printing nodes
                    if(nodes_written == false)
                    {
                        GiD_fBeginCoordinates(MeshFile);
                        for ( typename TModelPartType::NodesContainerType::iterator it = mMeshNodes.begin();
                                it != mMeshNodes.end(); ++it )
                        {
                            if constexpr (std::is_arithmetic<CoordinateType>::value)
                            {
                                if ( deformed )
                                    GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X(),
                                                           (it)->Y(), (it)->Z() );
                                else
                                    GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X0(),
                                                           (it)->Y0(), (it)->Z0() );
                            }
                            else if constexpr (std::is_same<CoordinateType, KRATOS_COMPLEX_TYPE>::value)
                            {
                                if constexpr (complex_coordinate_type == 1)
                                {
                                    if ( deformed )
                                        GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X().real(),
                                                               (it)->Y().real(), (it)->Z().real() );
                                    else
                                        GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X0().real(),
                                                               (it)->Y0().real(), (it)->Z0().real() );
                                }
                                else if constexpr (complex_coordinate_type == 2)
                                {
                                    if ( deformed )
                                        GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X().imag(),
                                                               (it)->Y().imag(), (it)->Z().imag() );
                                    else
                                        GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X0().imag(),
                                                               (it)->Y0().imag(), (it)->Z0().imag() );
                                }
                            }
                        }
                        GiD_fEndCoordinates(MeshFile);

                        nodes_written = true;
                    }

                    //printing elements
                    GiD_fBeginElements(MeshFile);
                    std::vector<int> nodes_id(mMeshElements.begin()->GetGeometry().size() + 1);
                    for ( typename TModelPartType::ElementsContainerType::iterator it = mMeshElements.begin();
                            it != mMeshElements.end(); ++it )
                    {
                        for ( unsigned int i=0; i< (it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry()[i].Id();

                        if ( mGeometryType == GeometryData::KratosGeometryType::Kratos_Line2D3
                                || mGeometryType == GeometryData::KratosGeometryType::Kratos_Line3D3 )
                        {
                            nodes_id[0] = (it)->GetGeometry()[0].Id();
                            nodes_id[1] = (it)->GetGeometry()[2].Id();
                            nodes_id[2] = (it)->GetGeometry()[1].Id();
                        }
                        nodes_id[(it)->GetGeometry().size()] = (it)->GetProperties().Id();

                        bool element_is_active = true;
                        if( it->IsDefined( ACTIVE ) )
                            element_is_active = it->Is(ACTIVE);
                        if ( element_is_active )
                        {
                            if ((it)->GetProperties().Id()==current_layer)
                            {
                                GiD_fWriteElementMat ( MeshFile, (it)->Id(), nodes_id.data() );
                            }
                        }
                    }

                    GiD_fEndElements(MeshFile);
                    GiD_fEndMesh(MeshFile);
                }
            }
        }

        if ( mMeshConditions.size() != 0 )
        {
            //compute number of layers
            int max_id = 0;
            for ( typename TModelPartType::ConditionsContainerType::iterator it = mMeshConditions.begin();
                    it != mMeshConditions.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                if (max_id < prop_id) max_id = prop_id;
            }
            if (max_id > 10000)
                std::cout<< "a property Id > 10000 found. Are u sure you need so many properties?" << std::endl;
            std::vector<int> conditions_per_layer (max_id+1,0);
            //fill layer list
            for ( typename TModelPartType::ConditionsContainerType::iterator it = mMeshConditions.begin();
                    it != mMeshConditions.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                conditions_per_layer[prop_id] += 1;
            }
            for (unsigned int current_layer = 0; current_layer < conditions_per_layer.size(); current_layer++)
            {
                bool nodes_written = false;
                if (conditions_per_layer[current_layer] > 0)
                {
                    // determine mesh name
                    std::stringstream current_layer_name (std::stringstream::in | std::stringstream::out);
                    current_layer_name << mMeshTitle << "_" << current_layer;

                    for ( typename TModelPartType::ConditionsContainerType::iterator it = mMeshConditions.begin(  );
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
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_2D, mGidElementType,
                                        mMeshConditions.begin()->GetGeometry().size() );
                    }
                    else if ( mMeshConditions.begin()->GetGeometry().WorkingSpaceDimension() == 3 )
                    {
                        GiD_fBeginMesh ( MeshFile, (char *) (current_layer_name.str() ).c_str(), GiD_3D, mGidElementType,
                                        mMeshConditions.begin()->GetGeometry().size() );
                    }
                    else
                        KRATOS_ERROR << "check working space dimension of model";

                    //printing nodes
                    if(nodes_written == false)
                    {
                        GiD_fBeginCoordinates(MeshFile);
                        for ( typename TModelPartType::NodesContainerType::iterator it = mMeshNodes.begin();
                                it != mMeshNodes.end(); ++it )
                        {
                            if constexpr (std::is_arithmetic<CoordinateType>::value)
                            {
                                if ( deformed )
                                    GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X(),
                                                           (it)->Y(), (it)->Z() );
                                else
                                    GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X0(),
                                                           (it)->Y0(), (it)->Z0() );
                            }
                            else if constexpr (std::is_same<CoordinateType, KRATOS_COMPLEX_TYPE>::value)
                            {
                                if constexpr (complex_coordinate_type == 1)
                                {
                                    if ( deformed )
                                        GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X().real(),
                                                               (it)->Y().real(), (it)->Z().real() );
                                    else
                                        GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X0().real(),
                                                               (it)->Y0().real(), (it)->Z0().real() );
                                }
                                else if constexpr (complex_coordinate_type == 2)
                                {
                                    if ( deformed )
                                        GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X().imag(),
                                                               (it)->Y().imag(), (it)->Z().imag() );
                                    else
                                        GiD_fWriteCoordinates ( MeshFile, (it)->Id(), (it)->X0().imag(),
                                                               (it)->Y0().imag(), (it)->Z0().imag() );
                                }
                            }
                        }
                        GiD_fEndCoordinates(MeshFile);
                        nodes_written = true;
                    }

                    //printing elements
                    GiD_fBeginElements(MeshFile);
                    std::vector<int> nodes_id(mMeshConditions.begin()->GetGeometry().size() + 1);
                    for ( typename TModelPartType::ConditionsContainerType::iterator it = mMeshConditions.begin(  );
                            it != mMeshConditions.end(); ++it )
                    {
                        for ( unsigned int i=0; i< (it)->GetGeometry().size(); i++ )
                            nodes_id[i] = (it)->GetGeometry()[i].Id();
                        nodes_id[(it)->GetGeometry().size()] = (it)->GetProperties().Id();

                        bool condition_is_active = true;
                        if( it->IsDefined( ACTIVE ) )
                            condition_is_active = it->Is(ACTIVE);
                        if ( condition_is_active )
                        {
                            if ((it)->GetProperties().Id()==current_layer)
                                GiD_fWriteElementMat ( MeshFile, (it)->Id(), nodes_id.data() );
                        }
                    }

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

    const typename TModelPartType::NodesContainerType& GetMeshNodes() const
    {
        return mMeshNodes;
    }

protected:
    ///member variables
    GeometryData::KratosGeometryType mGeometryType;
    GiD_ElementType mGidElementType;
    typename TModelPartType::NodesContainerType mMeshNodes;
    typename TModelPartType::ElementsContainerType mMeshElements;
    typename TModelPartType::ConditionsContainerType mMeshConditions;
    std::string mMeshTitle;
}; //class GidSDMeshContainer

} // namespace Kratos.

#endif // KRATOS_GID_SD_MESH_CONTAINER_H_INCLUDED defined
