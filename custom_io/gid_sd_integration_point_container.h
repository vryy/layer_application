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





#if !defined(KRATOS_GID_SD_INTEGRATION_POINT_CONTAINER_H_INCLUDED)
#define  KRATOS_GID_SD_INTEGRATION_POINT_CONTAINER_H_INCLUDED

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
#include "includes/kratos_flags.h"
#include "includes/deprecated_variables.h"
#include "geometries/geometry_data.h"

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
 * Auxiliary class to store integration point containers and perform result printing
 * on integration points
 */
class GidSDIntegrationPointsContainer
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(GidSDIntegrationPointsContainer);

    ///Constructor
    GidSDIntegrationPointsContainer( const char * gp_title, KratosGeometryFamily geometryFamily,
                             GeometryData::IntegrationMethod integration_method )
    : mGPTitle(gp_title), mIntegrationMethod(integration_method)
    {
        mKratosElementFamily = geometryFamily;

        if(     mKratosElementFamily == GeometryData::Kratos_Hexahedra )
        {
            mGidElementFamily = GiD_Hexahedra;
        }
        else if(mKratosElementFamily == GeometryData::Kratos_Tetrahedra )
        {
            mGidElementFamily = GiD_Tetrahedra;
        }
        else if(mKratosElementFamily == GeometryData::Kratos_Prism )
        {
            mGidElementFamily = GiD_Prism;
        }
        else if(mKratosElementFamily == GeometryData::Kratos_Quadrilateral)
        {
            mGidElementFamily = GiD_Quadrilateral;
        }
        else if(mKratosElementFamily == GeometryData::Kratos_Triangle)
        {
            mGidElementFamily = GiD_Triangle;
        }
        else if(mKratosElementFamily == GeometryData::Kratos_Linear)
        {
            mGidElementFamily = GiD_Linear;
        }
        else if(mKratosElementFamily == GeometryData::Kratos_Point)
        {
            mGidElementFamily = GiD_Point;
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error, "Unknown geometry family type", mKratosElementFamily)
        }
    }

    bool AddElement( const ModelPart::ElementsContainerType::iterator pElemIt )
    {
        KRATOS_TRY
        if( pElemIt->GetGeometry().GetGeometryFamily() == mKratosElementFamily
                && pElemIt->GetIntegrationMethod() == mIntegrationMethod )
        {
            bool element_is_active = true;
            if( pElemIt->IsDefined( ACTIVE ) )
                element_is_active = pElemIt->Is(ACTIVE);

            if (element_is_active)
            {
                mMeshElements.push_back( *(pElemIt.base() ) );
                return true;
            }
        }
        else return false;
        KRATOS_CATCH("")
    }

    bool AddCondition( const ModelPart::ConditionsContainerType::iterator pCondIt )
    {
        KRATOS_TRY
        if( pCondIt->GetGeometry().GetGeometryFamily() == mKratosElementFamily
                && pCondIt->GetIntegrationMethod() == mIntegrationMethod )
        {
            bool condition_is_active = true;
            if( pCondIt->IsDefined( ACTIVE ) )
                condition_is_active = pCondIt->Is(ACTIVE);

            if (condition_is_active)
            {
                mMeshConditions.push_back( *(pCondIt.base() ) );
                return true;
            }
        }
        else return false;
        KRATOS_CATCH("")
    }


    void PrintPartitionIndex( GiD_FILE ResultFile, const Variable<double>& rVariable, ModelPart& r_model_part,
                               double SolutionTag, unsigned int value_index, int rank )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile, mGPTitle);

            GiD_fBeginResult(ResultFile,  (char *)(rVariable.Name()).c_str(), (char *)("Kratos"), SolutionTag,
                             GiD_Scalar, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                        it != mMeshElements.end(); ++it )
                {
                    for(unsigned int i=0; i<mSize; i++)
                    {
                        GiD_fWriteScalar( ResultFile, it->Id(), (double) rank );
                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                        it != mMeshConditions.end(); ++it )
                {
                    for(unsigned int i=0; i<mSize; i++)
                    {
                        GiD_fWriteScalar( ResultFile, it->Id(), (double) rank );
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }


    virtual void PrintResults( GiD_FILE ResultFile, const Variable<double>& rVariable, ModelPart& r_model_part,
                               double SolutionTag, unsigned int value_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile, mGPTitle);

            GiD_fBeginResult(ResultFile,  (char *)(rVariable.Name()).c_str(), (char *)("Kratos"), SolutionTag,
                             GiD_Scalar, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<double> ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                        it != mMeshElements.end(); ++it )
                {
                    it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                     r_model_part.GetProcessInfo() );
                    for(unsigned int i=0; i<mSize; i++)
                    {
                        GiD_fWriteScalar( ResultFile, it->Id(), ValuesOnIntPoint[i] );
                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                        it != mMeshConditions.end(); ++it )
                {
                    it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                     r_model_part.GetProcessInfo() );
                    for(unsigned int i=0; i<mSize; i++)
                    {
                        GiD_fWriteScalar( ResultFile, it->Id(), ValuesOnIntPoint[i] );
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }


    virtual void PrintResults( GiD_FILE ResultFile, const Variable<array_1d<double,3> >& rVariable, ModelPart& r_model_part,
                               double SolutionTag, unsigned int value_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile, mGPTitle);

            GiD_fBeginResult( ResultFile,  (char *)(rVariable.Name()).c_str(), (char *)("Kratos"), SolutionTag,
                             GiD_Vector, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<array_1d<double,3> > ValuesOnIntPoint(mSize,ZeroVector(3));
            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                        it != mMeshElements.end(); ++it )
                {
                    it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                     r_model_part.GetProcessInfo() );
                    for(unsigned int i=0; i<mSize; i++)
                    {
                        if( ValuesOnIntPoint[0].size() == 3 )
                            GiD_fWriteVector( ResultFile, it->Id(), ValuesOnIntPoint[i][0],
                                             ValuesOnIntPoint[i][1], ValuesOnIntPoint[i][2] );
                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                        it != mMeshConditions.end(); ++it )
                {
                    it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                     r_model_part.GetProcessInfo() );
                    for(unsigned int i=0; i<mSize; i++)
                    {
                        GiD_fWriteVector( ResultFile, it->Id(), ValuesOnIntPoint[i][0],
                                         ValuesOnIntPoint[i][1], ValuesOnIntPoint[i][2] );
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, const Variable<array_1d<double,6> >& rVariable, ModelPart& r_model_part,
                               double SolutionTag, unsigned int value_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile, mGPTitle);

            GiD_fBeginResult( ResultFile, (char *)(rVariable.Name()).c_str(), ( char*)("Kratos"),
                             SolutionTag, GiD_Matrix, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<array_1d<double, 6> > ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                        it != mMeshElements.end(); ++it )
                {
                    it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                     r_model_part.GetProcessInfo() );
                    for(unsigned int i=0; i<mSize; i++)
                    {
                        GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[i][0],
                                           ValuesOnIntPoint[i][1], ValuesOnIntPoint[i][2],
                                           ValuesOnIntPoint[i][3], ValuesOnIntPoint[i][4],
                                           ValuesOnIntPoint[i][5] );
                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                        it != mMeshConditions.end(); ++it )
                {
                    it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                     r_model_part.GetProcessInfo() );
                    for(unsigned int i=0; i<mSize; i++)
                    {
                        GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[i][0],
                                           ValuesOnIntPoint[i][1], ValuesOnIntPoint[i][2],
                                           ValuesOnIntPoint[i][3], ValuesOnIntPoint[i][4],
                                           ValuesOnIntPoint[i][5] );
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }


    virtual void PrintResults( GiD_FILE ResultFile, const Variable<Vector>& rVariable, ModelPart& r_model_part,
                               double SolutionTag, unsigned int value_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile, mGPTitle);

            if( rVariable.Name() == std::string("INSITU_STRESS")
             || rVariable.Name() == std::string("PRESTRESS")
             || rVariable.Name() == std::string("STRESSES")
             || rVariable.Name() == std::string("PLASTIC_STRAIN_VECTOR") )
            {
                GiD_fBeginResult( ResultFile, (char *)(rVariable.Name()).c_str(), "Kratos", SolutionTag,
                                 GiD_Matrix, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            }
            else if( rVariable.Name() == std::string("MATERIAL_PARAMETERS")
                  || rVariable.Name() == std::string("INTERNAL_VARIABLES") )
            {
                std::stringstream param_index;
                param_index << value_index;
                GiD_fBeginResult( ResultFile, (char *)(rVariable.Name() + param_index.str() ).c_str(),
                                 "Kratos", SolutionTag, GiD_Scalar, GiD_OnGaussPoints,
                                 mGPTitle, NULL, 0, NULL );
            }
            else
            {
                GiD_fBeginResult( ResultFile, (char *)(rVariable.Name()).c_str(), "Kratos", SolutionTag,
                                 GiD_Vector, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            }

            std::vector<Vector> ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                        it != mMeshElements.end(); ++it )
                {
                    it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                     r_model_part.GetProcessInfo() );
                    for(unsigned int i=0; i<mSize; i++)
                    {
                        if( rVariable.Name() == std::string("INSITU_STRESS")
                         || rVariable.Name() == std::string("PRESTRESS")
                         || rVariable.Name() == std::string("STRESSES")
                         || rVariable.Name() == std::string("PLASTIC_STRAIN_VECTOR") )
                        {
                            if(ValuesOnIntPoint[i].size() ==6 )
                                GiD_fWrite3DMatrix( ResultFile, it->Id(),
                                                   ValuesOnIntPoint[i](0),
                                                   ValuesOnIntPoint[i](1),
                                                   ValuesOnIntPoint[i](2),
                                                   ValuesOnIntPoint[i](3),
                                                   ValuesOnIntPoint[i](4),
                                                   ValuesOnIntPoint[i](5) );
                            if(ValuesOnIntPoint[i].size() ==3 )
                                GiD_fWrite3DMatrix( ResultFile, it->Id(),
                                                   ValuesOnIntPoint[i](0),
                                                   ValuesOnIntPoint[i](1),
                                                   0.0,
                                                   ValuesOnIntPoint[i](2),
                                                   0.0,
                                                   0.0 );
                            if(ValuesOnIntPoint[i].size() ==4 )
                                GiD_fWrite3DMatrix( ResultFile, it->Id(),
                                                   ValuesOnIntPoint[i](0), //o_xx
                                                   ValuesOnIntPoint[i](1), //o_yy
                                                   ValuesOnIntPoint[i](3), //o_zz
                                                   ValuesOnIntPoint[i](2), //o_xy
                                                   0.0,
                                                   0.0 );
                        }
                        else if( rVariable.Name() == std::string("MATERIAL_PARAMETERS")
                              || rVariable.Name() == std::string("INTERNAL_VARIABLES") )
                        {
                            double value = 0.0;
                            if( ValuesOnIntPoint[i].size() > value_index )
                                value = ValuesOnIntPoint[i][value_index];
                            GiD_fWriteScalar( ResultFile, it->Id(), value );
                        }
                        else if( ValuesOnIntPoint[0].size() == 3 )
                            GiD_fWriteVector( ResultFile, it->Id(), ValuesOnIntPoint[i][0],
                                             ValuesOnIntPoint[i][1], ValuesOnIntPoint[i][2] );
                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                        it != mMeshConditions.end(); ++it )
                {
                    it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                     r_model_part.GetProcessInfo() );
                    for(unsigned int i=0; i<mSize; i++)
                    {
                        if( rVariable.Name() == std::string("INSITU_STRESS")
                         || rVariable.Name() == std::string("PRESTRESS")
                         || rVariable.Name() == std::string("STRESSES")
                         || rVariable.Name() == std::string("PLASTIC_STRAIN_VECTOR") )
                        {
                            if(ValuesOnIntPoint[i].size() ==6 )
                                GiD_fWrite3DMatrix( ResultFile, it->Id(),
                                                   ValuesOnIntPoint[i](0),
                                                   ValuesOnIntPoint[i](1),
                                                   ValuesOnIntPoint[i](2),
                                                   ValuesOnIntPoint[i](3),
                                                   ValuesOnIntPoint[i](4),
                                                   ValuesOnIntPoint[i](5) );
                            if(ValuesOnIntPoint[i].size() ==3 )
                                GiD_fWrite3DMatrix( ResultFile, it->Id(),
                                                   ValuesOnIntPoint[i](0),
                                                   ValuesOnIntPoint[i](1),
                                                   0.0,
                                                   ValuesOnIntPoint[i](2),
                                                   0.0,
                                                   0.0 );
                            if(ValuesOnIntPoint[i].size() ==4 )
                                GiD_fWrite3DMatrix( ResultFile, it->Id(),
                                                   ValuesOnIntPoint[i](0), //o_xx
                                                   ValuesOnIntPoint[i](1), //o_yy
                                                   ValuesOnIntPoint[i](3), //o_zz
                                                   ValuesOnIntPoint[i](2), //o_xy
                                                   0.0,
                                                   0.0 );
                        }
                        else if( rVariable.Name() == std::string("MATERIAL_PARAMETERS")
                              || rVariable.Name() == std::string("INTERNAL_VARIABLES") )
                        {
                            double value = 0.0;
                            if( ValuesOnIntPoint[i].size() > value_index )
                                value = ValuesOnIntPoint[i][value_index];
                            GiD_fWriteScalar( ResultFile, it->Id(), value );
                        }
                        else if( ValuesOnIntPoint[0].size() == 3 )
                            GiD_fWriteVector( ResultFile, it->Id(), ValuesOnIntPoint[i][0],
                                             ValuesOnIntPoint[i][1], ValuesOnIntPoint[i][2] );
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    virtual void PrintResults( GiD_FILE ResultFile, const Variable<Matrix>& rVariable, ModelPart& r_model_part,
                               double SolutionTag, int value_index )
    {
        if( mMeshElements.size() != 0 || mMeshConditions.size() != 0 )
        {
            WriteGaussPoints(ResultFile, mGPTitle);

            GiD_fBeginResult(ResultFile,  (char *)(rVariable.Name()).c_str(), (char *)("Kratos"),
                             SolutionTag, GiD_Matrix, GiD_OnGaussPoints, mGPTitle, NULL, 0, NULL );
            std::vector<Matrix> ValuesOnIntPoint(mSize);
            if( mMeshElements.size() != 0 )
            {
                for( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                        it != mMeshElements.end(); ++it )
                {
                    it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                     r_model_part.GetProcessInfo() );
                    for(unsigned int i=0; i<mSize; i++)
                    {
                        if(ValuesOnIntPoint[i].size1() ==3
                                && ValuesOnIntPoint[i].size2() ==3)
                            GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[i](0,0),
                                               ValuesOnIntPoint[i](1,1), ValuesOnIntPoint[i](2,2),
                                               ValuesOnIntPoint[i](0,1), ValuesOnIntPoint[i](1,2),
                                               ValuesOnIntPoint[i](0,2) );
                        else if(ValuesOnIntPoint[i].size1() ==2
                                && ValuesOnIntPoint[i].size2() ==2)
                            GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[i](0,0),
                                               ValuesOnIntPoint[i](1,1), 0.0,
                                               ValuesOnIntPoint[i](0,1), 0.0, 0.0);
                        else if(ValuesOnIntPoint[i].size1() ==1
                                && ValuesOnIntPoint[i].size2() ==3)
                            GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[i](0,0),
                                               ValuesOnIntPoint[i](0,1), 0.0,
                                               ValuesOnIntPoint[i](0,2), 0.0, 0.0);
                        else if(ValuesOnIntPoint[i].size1() ==1
                                && ValuesOnIntPoint[i].size2() ==4)
                            GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[i](0,0),
                                               ValuesOnIntPoint[i](0,1), ValuesOnIntPoint[i](0,2),
                                               ValuesOnIntPoint[i](0,3), 0.0, 0.0);
                        else if(ValuesOnIntPoint[i].size1() ==1
                                && ValuesOnIntPoint[i].size2() ==6)
                            GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[i](0,0),
                                               ValuesOnIntPoint[i](0,1), ValuesOnIntPoint[i](0,2),
                                               ValuesOnIntPoint[i](0,3), ValuesOnIntPoint[i](0,4),
                                               ValuesOnIntPoint[i](0,5) );

                    }
                }
            }
            if( mMeshConditions.size() != 0 )
            {
                for( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                        it != mMeshConditions.end(); ++it )
                {
                    it->GetValueOnIntegrationPoints( rVariable, ValuesOnIntPoint,
                                                     r_model_part.GetProcessInfo() );
                    for(unsigned int i=0; i<mSize; i++)
                    {
                        if(ValuesOnIntPoint[i].size1() ==3
                                && ValuesOnIntPoint[i].size2() ==3)
                            GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[i](0,0),
                                               ValuesOnIntPoint[i](1,1), ValuesOnIntPoint[i](2,2),
                                               ValuesOnIntPoint[i](0,1), ValuesOnIntPoint[i](1,2),
                                               ValuesOnIntPoint[i](0,2) );
                        else if(ValuesOnIntPoint[i].size1() ==1
                                && ValuesOnIntPoint[i].size2() ==6)
                            GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[i](0,0),
                                               ValuesOnIntPoint[i](0,1), ValuesOnIntPoint[i](0,2),
                                               ValuesOnIntPoint[i](0,3), ValuesOnIntPoint[i](0,4),
                                               ValuesOnIntPoint[i](0,5) );
                        else if(ValuesOnIntPoint[i].size1() ==1
                                && ValuesOnIntPoint[i].size2() ==3)
                            GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[i](0,0),
                                               ValuesOnIntPoint[i](0,1), 0.0,
                                               ValuesOnIntPoint[i](0,2), 0.0, 0.0);
                        else if(ValuesOnIntPoint[i].size1() ==1
                                && ValuesOnIntPoint[i].size2() ==4)
                            GiD_fWrite3DMatrix( ResultFile, it->Id(), ValuesOnIntPoint[i](0,0),
                                               ValuesOnIntPoint[i](0,1), ValuesOnIntPoint[i](0,2),
                                               ValuesOnIntPoint[i](0,3), 0.0, 0.0);
                    }
                }
            }
            GiD_fEndResult(ResultFile);
        }
    }

    void Reset()
    {
        mMeshElements.clear();
        mMeshConditions.clear();
    }

public:

    void WriteGaussPoints(GiD_FILE MeshFile, const char* GPTitle)
    {
        //setting up gauss points
        bool gp_written = false;
        if( !gp_written && mMeshElements.size() != 0 )
        {
            const Element::GeometryType::IntegrationPointsArrayType& integration_points = mMeshElements.begin()->GetGeometry().IntegrationPoints( mIntegrationMethod );
            mSize = integration_points.size();

            GiD_fBeginGaussPoint( MeshFile, GPTitle, mGidElementFamily, NULL, integration_points.size(), 0, 0 );
            if(mGidElementFamily == GiD_Tetrahedra || mGidElementFamily == GiD_Hexahedra)
            {
                for(std::size_t i = 0; i < integration_points.size(); ++i)
                    GiD_fWriteGaussPoint3D( MeshFile, integration_points[i].X(), integration_points[i].Y(), integration_points[i].Z());
            }
            else if(mGidElementFamily == GiD_Triangle || mGidElementFamily == GiD_Quadrilateral)
            {
                for(std::size_t i = 0; i < integration_points.size(); ++i)
                    GiD_fWriteGaussPoint2D( MeshFile, integration_points[i].X(), integration_points[i].Y());
            }
            GiD_fEndGaussPoint(MeshFile);

            gp_written = true;
        }

        if( !gp_written && mMeshConditions.size() != 0 )
        {
            const Element::GeometryType::IntegrationPointsArrayType& integration_points = mMeshConditions.begin()->GetGeometry().IntegrationPoints( mIntegrationMethod );
            mSize = integration_points.size();

            if(mGidElementFamily != GiD_Linear)
                GiD_fBeginGaussPoint( MeshFile, GPTitle, mGidElementFamily, NULL, integration_points.size(), 0, 0 );
            else
                GiD_fBeginGaussPoint( MeshFile, GPTitle, mGidElementFamily, NULL, integration_points.size(), 0, 1 );
            // Note: By now, Natural Coordinates for line elements cannot be "Given"

            if(mGidElementFamily == GiD_Tetrahedra || mGidElementFamily == GiD_Hexahedra)
            {
                for(std::size_t i = 0; i < integration_points.size(); ++i)
                    GiD_fWriteGaussPoint3D( MeshFile, integration_points[i].X(), integration_points[i].Y(), integration_points[i].Z());
            }
            else if(mGidElementFamily == GiD_Triangle || mGidElementFamily == GiD_Quadrilateral)
            {
                for(std::size_t i = 0; i < integration_points.size(); ++i)
                    GiD_fWriteGaussPoint2D( MeshFile, integration_points[i].X(), integration_points[i].Y());
            }
            GiD_fEndGaussPoint(MeshFile);
        }
    }


protected:


    ///member variables
    const char* mGPTitle;
    KratosGeometryFamily mKratosElementFamily;
    GiD_ElementType mGidElementFamily;
    GeometryData::IntegrationMethod mIntegrationMethod;
    std::size_t mSize;
    ModelPart::ElementsContainerType mMeshElements;
    ModelPart::ConditionsContainerType mMeshConditions;
};//class GidSDIntegrationPointsContainer
}// namespace Kratos.

#endif // KRATOS_GID_SD_INTEGRATION_POINT_CONTAINER_H_INCLUDED defined

