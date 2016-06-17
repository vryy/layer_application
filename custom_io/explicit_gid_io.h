//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Jun 2016 $
//   Revision:            $Revision: 1.0 $
//
//







#if !defined(KRATOS_EXPLICIT_GID_IO_H_INCLUDED)
#define  KRATOS_EXPLICIT_GID_IO_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <iomanip>

// External includes
#define USE_CONST
#include "gidpost/source/gidpost.h"

// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "geometries/geometry_data.h"

#include "utilities/timer.h"

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

///Flags for mesh writing
enum WriteDeformedMeshFlag {WriteDeformed, WriteUndeformed};
enum WriteConditionsFlag {WriteConditions, WriteElementsOnly, WriteConditionsOnly};
enum MultiFileFlag {SingleFile, MultipleFiles};


/**
 * This class defines an interface to the GiDPost library
 * in order to provide GiD compliant I/O functionality
 */
template<class TGaussPointContainer, class TMeshContainer>
class ExplicitGidIO : public IO
{
public:
    ///pointer definition of ExplicitGidIO
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitGidIO);

    ///typedefs
	typedef IO BaseType;

    ///Flags for mesh writing
//             enum WriteDeformedMeshFlag{WriteDeformed, WriteUndeformed};
//             enum WriteConditionsFlag{WriteConditions, WriteElementsOnly};
//             enum MultiFileFlag{SingleFile, MultipleFiles};

    ///Constructor
    ///single stream IO constructor
    ExplicitGidIO( const std::string& rDatafilename,
           GiD_PostMode Mode,
           MultiFileFlag use_multiple_files_flag,
           WriteDeformedMeshFlag write_deformed_flag,
           WriteConditionsFlag write_conditions_flag
         )
    {
        mMode = Mode;
        mResultFileOpen = false;
        mMeshFileOpen = false;
        mWriteDeformed = write_deformed_flag;
        mWriteConditions = write_conditions_flag;
        mUseMultiFile = use_multiple_files_flag;
        mResultFileName = rDatafilename;
        InitializeResultFile(mResultFileName);
        mMeshFileName = rDatafilename;
//                 mMeshFileName += ".post.msh";
        SetUpMeshContainers();
        SetUpGaussPointContainers();
    }

    ///Destructor.
    virtual ~ExplicitGidIO()
    {
        Timer::PrintTimingInformation();

        if ( mResultFileOpen )
        {
            GiD_fClosePostResultFile( mResultFile );
            mResultFileOpen = false;
        }
    }

    ///initialization functions
    /**
     * creates the mesh containers for all different element types.
     * Note that the containers are not filled yet in here!
     */
    void SetUpMeshContainers()
    {
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Hexahedra3D20,
                                          GiD_Hexahedra, "Kratos_Hexahedra3D20_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Hexahedra3D27,
                                          GiD_Hexahedra, "Kratos_Hexahedra3D27_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Hexahedra3D8,
                                          GiD_Hexahedra, "Kratos_Hexahedra3D8_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Prism3D15,
                                          GiD_Prism, "Kratos_Prism3D15_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Prism3D6,
                                          GiD_Prism, "Kratos_Prism3D6_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Quadrilateral2D4,
                                          GiD_Quadrilateral, "Kratos_Quadrilateral2D4_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Quadrilateral2D8,
                                          GiD_Quadrilateral, "Kratos_Quadrilateral2D8_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Quadrilateral2D9,
                                          GiD_Quadrilateral, "Kratos_Quadrilateral2D9_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Quadrilateral3D4,
                                          GiD_Quadrilateral, "Kratos_Quadrilateral3D4_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Quadrilateral3D8,
                                          GiD_Quadrilateral, "Kratos_Quadrilateral3D8_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Quadrilateral3D9,
                                          GiD_Quadrilateral, "Kratos_Quadrilateral3D9_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Tetrahedra3D10,
                                          GiD_Tetrahedra, "Kratos_Tetrahedra3D10_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Tetrahedra3D4,
                                          GiD_Tetrahedra, "Kratos_Tetrahedra3D4_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Triangle2D3,
                                          GiD_Triangle, "Kratos_Triangle2D3_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Triangle2D6,
                                          GiD_Triangle, "Kratos_Triangle2D6_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Triangle3D3,
                                          GiD_Triangle, "Kratos_Triangle3D3_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Triangle3D6,
                                          GiD_Triangle, "Kratos_Triangle3D6_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Line2D2,
                                          GiD_Linear, "Kratos_Line2D2_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Line3D2,
                                          GiD_Linear, "Kratos_Line3D2_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Line2D3,
                                          GiD_Linear, "Kratos_Line2D3_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Line3D3,
                                          GiD_Linear, "Kratos_Line3D3_Mesh" ) );
        mGidMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Point3D,
                                          GiD_Point, "Kratos_Point3D_Mesh" ) );


    }//SetUpMeshContainers

    /**
     * creates the gauss point containers for all different element types.
     * Note that the containers are not filled yet in here!
     */
    virtual void SetUpGaussPointContainers()
    {
        //case Triangle
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tri_gauss_legendre_1_element_gp",
                                            GeometryData::Kratos_Triangle, GiD_Triangle, GeometryData::GI_GAUSS_1 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tri_gauss_legendre_2_element_gp",
                                            GeometryData::Kratos_Triangle, GiD_Triangle, GeometryData::GI_GAUSS_2 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tri_gauss_legendre_3_element_gp",
                                            GeometryData::Kratos_Triangle, GiD_Triangle, GeometryData::GI_GAUSS_3 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tri_gauss_legendre_4_element_gp",
                                            GeometryData::Kratos_Triangle, GiD_Triangle, GeometryData::GI_GAUSS_4 ) );

        //case Quadrilateral
        mGidGaussPointContainers.push_back( TGaussPointContainer( "quad_gauss_legendre_1_element_gp",
                                            GeometryData::Kratos_Quadrilateral, GiD_Quadrilateral, GeometryData::GI_GAUSS_1 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "quad_gauss_legendre_2_element_gp",
                                            GeometryData::Kratos_Quadrilateral, GiD_Quadrilateral, GeometryData::GI_GAUSS_2 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "quad_gauss_legendre_3_element_gp",
                                            GeometryData::Kratos_Quadrilateral, GiD_Quadrilateral, GeometryData::GI_GAUSS_3 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "quad_gauss_legendre_4_element_gp",
                                            GeometryData::Kratos_Quadrilateral, GiD_Quadrilateral, GeometryData::GI_GAUSS_4 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "quad_gauss_legendre_5_element_gp",
                                            GeometryData::Kratos_Quadrilateral, GiD_Quadrilateral, GeometryData::GI_GAUSS_5 ) );

        //case Tetrahedra
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tet_gauss_legendre_1_element_gp",
                                            GeometryData::Kratos_Tetrahedra, GiD_Tetrahedra, GeometryData::GI_GAUSS_1 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tet_gauss_legendre_2_element_gp",
                                            GeometryData::Kratos_Tetrahedra, GiD_Tetrahedra, GeometryData::GI_GAUSS_2 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tet_gauss_legendre_3_element_gp",
                                            GeometryData::Kratos_Tetrahedra, GiD_Tetrahedra, GeometryData::GI_GAUSS_3 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tet_gauss_legendre_4_element_gp",
                                            GeometryData::Kratos_Tetrahedra, GiD_Tetrahedra, GeometryData::GI_GAUSS_4 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "tet_gauss_legendre_5_element_gp",
                                            GeometryData::Kratos_Tetrahedra, GiD_Tetrahedra, GeometryData::GI_GAUSS_5 ) );

        //case Hexahedra
        mGidGaussPointContainers.push_back( TGaussPointContainer( "hex_gauss_legendre_1_element_gp",
                                            GeometryData::Kratos_Hexahedra, GiD_Hexahedra, GeometryData::GI_GAUSS_1 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "hex_gauss_legendre_2_element_gp",
                                            GeometryData::Kratos_Hexahedra, GiD_Hexahedra, GeometryData::GI_GAUSS_2 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "hex_gauss_legendre_3_element_gp",
                                            GeometryData::Kratos_Hexahedra, GiD_Hexahedra, GeometryData::GI_GAUSS_3 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "hex_gauss_legendre_4_element_gp",
                                            GeometryData::Kratos_Hexahedra, GiD_Hexahedra, GeometryData::GI_GAUSS_4 ) );
        mGidGaussPointContainers.push_back( TGaussPointContainer( "hex_gauss_legendre_5_element_gp",
                                            GeometryData::Kratos_Hexahedra, GiD_Hexahedra, GeometryData::GI_GAUSS_5 ) );
    }//SetUpGaussPointContainers


    ///general ExplicitGidIO related functions
    /**
     * TODO: to be removed
     */
    void ChangeOutputName(const std::string& rDatafilename )
    {
        KRATOS_TRY
        mMeshFileName = rDatafilename;
        mResultFileName = rDatafilename;
        KRATOS_CATCH("")
    }

    /**
     * sets up the file names and opens the result file in case there
     * is ASCII mode and only one file written
     */
    void InitializeResultFile( std::string const& rResultFileName )
    {
        //std::cout << "initializing result files" << std::endl;
        mResultFileName = rResultFileName;
    }

    /**
     * TODO: check whether this is still necessary!
     */
    void  CloseResultFile()
    {
        if ( mResultFileOpen )
            GiD_fClosePostResultFile( mResultFile );
    }

    /**
     * TODO: check whether this is still necessary!
     */
    void Flush()
    {
        GiD_fFlushPostFile( mResultFile );
    }

    /**
     * Turn back information as a string.
     */
    virtual std::string Info() const
    {
        return "gid io";
    }

    /**
     * Print information about this object.
     */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /**
     * Print object's data.
     */
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///result functions
    /**
     * This has to be called for each solution step BEFORE any results
     * (on nodes and on gauss points) is written
     * @param SolutionTag the current solution step (i.e. time)
     * @param conditions_flag states whether results should also be written on conditions
     */
    virtual void InitializeResults( double name, MeshType rThisMesh )
    {
        if ( mMode == GiD_PostAscii && ! mResultFileOpen )
        {
            std::stringstream file_name;
            file_name << mResultFileName << std::setprecision(12) << "_" << name << ".post.res";
            mResultFile = GiD_fOpenPostResultFile((char*)(file_name.str()).c_str(), mMode);
            mResultFileOpen = true;
            
        }
        //initializing gauss points containers
        if ( mWriteConditions != WriteConditionsOnly )
        {
            int i=0;
            for ( MeshType::ElementIterator element_iterator = rThisMesh.ElementsBegin();
                    element_iterator != rThisMesh.ElementsEnd(); ++element_iterator )
            {            
                for ( typename std::vector<TGaussPointContainer>::iterator it =
                            mGidGaussPointContainers.begin();
                        it != mGidGaussPointContainers.end(); it++ )
                {                      
                    i++;
                    if ( it->AddElement( element_iterator ) )
                        break;
                }
                
            }
        }

        if ( mWriteConditions == WriteConditions || mWriteConditions == WriteConditionsOnly )
        {
            for ( MeshType::ConditionsContainerType::iterator conditions_iterator =
                        rThisMesh.ConditionsBegin(); conditions_iterator
                    != rThisMesh.ConditionsEnd(); conditions_iterator++ )
            {
                for ( typename std::vector<TGaussPointContainer>::iterator it =
                            mGidGaussPointContainers.begin();
                        it != mGidGaussPointContainers.end(); it++ )
                {                      

                    if ( it->AddCondition( conditions_iterator ) )
                        break;
                }
            }
        }
    }

    /**
     * This has to be called for each solution step AFTER all the results
     * have been written
     */
    void FinalizeResults()
    {
        if ( mUseMultiFile == MultipleFiles || mMode == GiD_PostAscii )
        {
            GiD_fClosePostResultFile( mResultFile );
            mResultFileOpen = false;
        }
        //resetting gauss point containers
        for ( typename std::vector<TGaussPointContainer>::iterator it =
                    mGidGaussPointContainers.begin();
                it != mGidGaussPointContainers.end(); it++ )
        {
            it->Reset();
        }
    }



    ///functions for writing nodal results

	///////////////////////////////////////////////////////////////////////
	//////                  HISTORICAL DATABASE BLOCK                 /////
	///////////////////////////////////////////////////////////////////////
	 /**
     * writes nodal results for variables of type bool
     */
    void WriteNodalResults( Variable<bool> const& rVariable,
                            NodesContainerType& rNodes, double SolutionTag,
                            std::size_t SolutionStepNumber)
    {

        Timer::Start("Writing Results");
        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Scalar,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
            GiD_fWriteScalar( mResultFile, i_node->Id(), i_node->GetSolutionStepValue(rVariable,
                             SolutionStepNumber) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }


    ///functions for writing nodal results
    /**
     * writes nodal results for variables of type double
     */
    void WriteNodalResults( Variable<double> const& rVariable,
                            NodesContainerType& rNodes, double SolutionTag,
                            std::size_t SolutionStepNumber)
    {

        Timer::Start("Writing Results");
        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Scalar,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
            GiD_fWriteScalar( mResultFile, i_node->Id(), i_node->GetSolutionStepValue(rVariable,
                             SolutionStepNumber) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }

    /**
     * writes nodal results for variables of type array_1d<double, 3>
     * (e.g. DISPLACEMENT)
     */
    void WriteNodalResults( Variable<array_1d<double, 3> > const& rVariable,
                            NodesContainerType& rNodes,
                            double SolutionTag, std::size_t SolutionStepNumber)
    {

        Timer::Start("Writing Results");

        GiD_fBeginResult(mResultFile,(char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Vector,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for (NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
        {
            array_1d<double, 3>& temp = i_node->GetSolutionStepValue( rVariable,
                                        SolutionStepNumber );
            GiD_fWriteVector( mResultFile, i_node->Id(), temp[0], temp[1], temp[2] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }


    /**
     * writes nodal results for variables of type Vector
     * (note that only vectors with 3 components can be printed)
     */
    void WriteNodalResults( Variable<Vector> const& rVariable,
                            NodesContainerType& rNodes,
                            double SolutionTag, std::size_t SolutionStepNumber)
    {

        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Matrix,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for (NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
        {
            Vector& tempVector = i_node->GetSolutionStepValue(rVariable,
                                 SolutionStepNumber);
            if (tempVector.size() ==3 )
                GiD_fWriteVector(mResultFile, i_node->Id(), tempVector(0), tempVector(1), tempVector(2) );
            else if (tempVector.size() == 6 )
                GiD_fWrite3DMatrix( mResultFile, i_node->Id(), tempVector(0), tempVector(1), tempVector(2),
                                    tempVector(3), tempVector(4), tempVector(5) );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }

    /**
     * writes nodal results for variables of type Matrix
     */
    void WriteNodalResults( Variable<Matrix> const& rVariable,
                            NodesContainerType& rNodes,
                            double SolutionTag, std::size_t SolutionStepNumber)
    {

        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Matrix,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for (NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
        {
            Matrix& tempMatrix = i_node->GetSolutionStepValue(rVariable,
                    SolutionStepNumber);
            //Matrix& tempMatrix = i_node->GetValue(rVariable);
            if (tempMatrix.size1() ==3 && tempMatrix.size2() ==3)
            {
                GiD_fWrite3DMatrix( mResultFile,  i_node->Id(), tempMatrix(0,0), tempMatrix(1,1),
                                    tempMatrix(2,2), tempMatrix(0,1), tempMatrix(1,2),
                                    tempMatrix(0,2) );
            }
            else if (tempMatrix.size1() ==2 && tempMatrix.size2() ==2)
            {
                GiD_fWrite2DMatrix( mResultFile, i_node->Id(), tempMatrix(0,0), tempMatrix(1,1), tempMatrix(0,1));
            }

            else if (tempMatrix.size1() ==1 && tempMatrix.size2() ==3)
            {

                GiD_fWrite3DMatrix( mResultFile, i_node->Id(), tempMatrix(0,0), tempMatrix(0,1), 0.00,
                                   tempMatrix(0,2), 0.00, 0.00);
            }
            else if (tempMatrix.size1() ==1 && tempMatrix.size2() ==6)
            {
                GiD_fWrite3DMatrix( mResultFile, i_node->Id(), tempMatrix(0,0), tempMatrix(0,1), tempMatrix(0,2),
                                   tempMatrix(0,3), tempMatrix(0,4), tempMatrix(0,5) );
            }
            //i_node->GetValue(rVariable) = tempMatrix;

        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }

   void WriteLocalAxesOnNodes( Variable<array_1d<double, 3> > const& rVariable,
                            NodesContainerType& rNodes,
                            double SolutionTag, std::size_t SolutionStepNumber)
    {

        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_LocalAxes,
                         GiD_OnNodes, NULL, NULL, 0, NULL );

        for (NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
        {
            array_1d<double, 3>& temp = i_node->GetSolutionStepValue( rVariable,
                                        SolutionStepNumber );
            GiD_fWriteLocalAxes( mResultFile, i_node->Id(), temp[0], temp[1], temp[2] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }
   	///////////////////////////////////////////////////////////////////////
	//////                 NON- HISTORICAL DATABASE BLOCK                 /////
	///////////////////////////////////////////////////////////////////////
	 /**
     * writes nodal results for variables of type bool
     */
    void WriteNodalResultsNonHistorical( Variable<bool> const& rVariable, NodesContainerType& rNodes, double SolutionTag)
    {

        Timer::Start("Writing Results");
        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Scalar,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
            GiD_fWriteScalar( mResultFile, i_node->Id(), i_node->GetValue(rVariable) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }


    ///functions for writing nodal results
    /**
     * writes nodal results for variables of type double
     */
    void WriteNodalResultsNonHistorical( Variable<double> const& rVariable, NodesContainerType& rNodes, double SolutionTag)
    {

        Timer::Start("Writing Results");
        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Scalar,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
            GiD_fWriteScalar( mResultFile, i_node->Id(), i_node->GetValue(rVariable) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }

    /**
     * writes nodal results for variables of type array_1d<double, 3>
     * (e.g. DISPLACEMENT)
     */
    void WriteNodalResultsNonHistorical( Variable<array_1d<double, 3> > const& rVariable, NodesContainerType& rNodes, double SolutionTag)
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult(mResultFile,(char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Vector,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for (NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
        {
            array_1d<double, 3>& temp = i_node->GetValue( rVariable);
            GiD_fWriteVector( mResultFile, i_node->Id(), temp[0], temp[1], temp[2] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }


    /**
     * writes nodal results for variables of type Vector
     * (note that only vectors with 3 components can be printed)
     */
    void WriteNodalResultsNonHistorical( Variable<Vector> const& rVariable, NodesContainerType& rNodes, double SolutionTag )
    {

        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Matrix,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for (NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
        {
            Vector& tempVector = i_node->GetSolutionStepValue(rVariable);
            if (tempVector.size() ==3 )
                GiD_fWriteVector(mResultFile, i_node->Id(), tempVector(0), tempVector(1), tempVector(2) );
            else if (tempVector.size() == 6 )
                GiD_fWrite3DMatrix( mResultFile, i_node->Id(), tempVector(0), tempVector(1), tempVector(2),
                                    tempVector(3), tempVector(4), tempVector(5) );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }

    /**
     * writes nodal results for variables of type Matrix
     */
    void WriteNodalResultsNonHistorical( Variable<Matrix> const& rVariable, NodesContainerType& rNodes, double SolutionTag)
    {

        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Matrix,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for (NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
        {
            Matrix& tempMatrix = i_node->GetSolutionStepValue(rVariable);
            //Matrix& tempMatrix = i_node->GetValue(rVariable);
            if (tempMatrix.size1() ==3 && tempMatrix.size2() ==3)
            {
                GiD_fWrite3DMatrix( mResultFile,  i_node->Id(), tempMatrix(0,0), tempMatrix(1,1),
                                    tempMatrix(2,2), tempMatrix(0,1), tempMatrix(1,2),
                                    tempMatrix(0,2) );
            }
            else if (tempMatrix.size1() ==2 && tempMatrix.size2() ==2)
            {
                GiD_fWrite2DMatrix( mResultFile, i_node->Id(), tempMatrix(0,0), tempMatrix(1,1), tempMatrix(0,1));
            }

            else if (tempMatrix.size1() ==1 && tempMatrix.size2() ==3)
            {

                GiD_fWrite3DMatrix( mResultFile, i_node->Id(), tempMatrix(0,0), tempMatrix(0,1), 0.00,
                                   tempMatrix(0,2), 0.00, 0.00);
            }
            else if (tempMatrix.size1() ==1 && tempMatrix.size2() ==6)
            {
                GiD_fWrite3DMatrix( mResultFile, i_node->Id(), tempMatrix(0,0), tempMatrix(0,1), tempMatrix(0,2),
                                   tempMatrix(0,3), tempMatrix(0,4), tempMatrix(0,5) );
            }
            //i_node->GetValue(rVariable) = tempMatrix;

        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }

   void WriteLocalAxesOnNodesNonHistorical( Variable<array_1d<double, 3> > const& rVariable, NodesContainerType& rNodes, double SolutionTag)
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_LocalAxes,
                         GiD_OnNodes, NULL, NULL, 0, NULL );

        for (NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node)
        {
            array_1d<double, 3>& temp = i_node->GetSolutionStepValue( rVariable);
            GiD_fWriteLocalAxes( mResultFile, i_node->Id(), temp[0], temp[1], temp[2] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");

    }








    ///mesh writing functions
    /**
     * opens a new mesh group
     */
    void InitializeMesh( double name )
    {
        if ( mUseMultiFile == MultipleFiles )
        {
            if ( mMode == GiD_PostAscii && ! mMeshFileOpen )
            {
                std::stringstream file_name;
                file_name << std::setprecision(12) << mMeshFileName << "_" << name << ".post.msh";
                mMeshFile = GiD_fOpenPostMeshFile( (char *)(file_name.str()).c_str(), mMode);
                mMeshFileOpen = true;
            }
            if ( (mMode == GiD_PostBinary || mMode == GiD_PostHDF5) && ! mResultFileOpen )
            {
                std::stringstream file_name;
                file_name << std::setprecision(12) << mResultFileName << "_" << name << ".post.bin";
                if ( ! mResultFileOpen )
                {
                    mResultFile = GiD_fOpenPostResultFile((char*)(file_name.str()).c_str(), mMode);
                    mResultFileOpen = true;
                }
				mMeshFile = mResultFile;
            }
        }
        if ( mUseMultiFile == SingleFile )
        {
            if ( (mMode == GiD_PostBinary || mMode == GiD_PostHDF5) && ! mResultFileOpen )
            {
                std::stringstream file_name;
                file_name << mResultFileName << ".post.bin";
		         //KRATOS_WATCH(file_name.str())
		        mResultFile = GiD_fOpenPostResultFile((char*)(file_name.str()).c_str(), mMode);
                if ( mResultFile == 0) //error handler can not be zero
                {
                    std::stringstream buffer;
                    buffer << "error opening results file:" << "/" <<  file_name.str()   << "/";
                    KRATOS_THROW_ERROR(std::runtime_error, buffer.str(), "");
                }
                mResultFileOpen = true;
				mMeshFile = mResultFile;
            }
            if ( mMode == GiD_PostAscii && ! mMeshFileOpen )
            {
                std::stringstream file_name;
                file_name << mMeshFileName << "_" << name << ".post.msh";
                mMeshFile = GiD_fOpenPostMeshFile( (char *)(file_name.str()).c_str(), mMode);
                mMeshFileOpen = true;
            }

        }
    }

    /**
     * closes a mesh group
     */
    void FinalizeMesh()
    {
        if ( mUseMultiFile == MultipleFiles && mMode == GiD_PostAscii )
        {
            GiD_fClosePostMeshFile(mMeshFile);
            mMeshFileOpen = false;
        }
        if ( mUseMultiFile == SingleFile && mMode == GiD_PostAscii )
        {
            GiD_fClosePostMeshFile(mMeshFile);
            mMeshFileOpen = false;
        }
    }

    /**
     * Writes a node mesh.
     * @param rThisMesh the given mesh to be written to the output file
     * @param solution_step the current solution step
     * @param deformed_flag indicates whether the mesh shall be written in deformed
     * or undeformed state
     * @param Mode either GiD_PostAscii (default) or GiD_PostBinary
     */
    void WriteNodeMesh( MeshType& rThisMesh )
    {
        KRATOS_TRY

        Timer::Start("Writing Mesh");

        GiD_fBeginMesh(mMeshFile,  "Kratos Mesh",GiD_3D,GiD_Point,1);
        GiD_fBeginCoordinates(mMeshFile);
        for ( MeshType::NodeIterator node_iterator = rThisMesh.NodesBegin();
                node_iterator != rThisMesh.NodesEnd();
                ++node_iterator)
        {
            if ( mWriteDeformed == WriteUndeformed )
                GiD_fWriteCoordinates(mMeshFile, node_iterator->Id(), node_iterator->X0(),
                                      node_iterator->Y0(), node_iterator->Z0() );
            else if ( mWriteDeformed == WriteDeformed )
                GiD_fWriteCoordinates(mMeshFile, node_iterator->Id(), node_iterator->X(),
                                      node_iterator->Y(), node_iterator->Z() );
            else
                KRATOS_THROW_ERROR( std::logic_error,"undefined WriteDeformedMeshFlag","" );
        }
        GiD_fEndCoordinates(mMeshFile);
        int nodes_id[1];
        GiD_fBeginElements(mMeshFile);

//         mNodeList.clear();

        for ( MeshType::NodeIterator node_iterator = rThisMesh.NodesBegin();
                node_iterator != rThisMesh.NodesEnd();
                ++node_iterator)
        {
            nodes_id[0] = node_iterator->Id();
            GiD_fWriteElement(mMeshFile,node_iterator->Id(), nodes_id);
//             mNodeList.push_back(*node_iterator);
        }
        GiD_fEndElements(mMeshFile);
        GiD_fEndMesh(mMeshFile);

//         mNodeList.Unique();

        Timer::Stop("Writing Mesh");

        KRATOS_CATCH("")
    }//WriteNodeMesh



    /**
     * This is a multi-purpose function that writes arbitrary meshes of elements
     * and conditions in either deformed or undeformed state
     * @param rThisMesh the current mesh to be written
     * @param deformed_flag states whether the mesh should be written in deformed configuration
     * @param conditions_flag states whether conditions should also be written
     */
    void WriteMesh( MeshType& rThisMesh )
    {
        KRATOS_TRY

        Timer::Start("Writing Mesh");
        
        if ( mWriteConditions != WriteConditionsOnly )
        {
            for ( MeshType::ElementIterator element_iterator = rThisMesh.ElementsBegin();
                    element_iterator != rThisMesh.ElementsEnd(); ++element_iterator)
                for ( typename std::vector<TMeshContainer>::iterator it = mGidMeshContainers.begin();
                        it != mGidMeshContainers.end(); it++ )
                    if ( it->AddElement( element_iterator ) )
                        break;
        }
        if ( mWriteConditions == WriteConditions || mWriteConditions == WriteConditionsOnly )
		{
            for ( MeshType::ConditionsContainerType::iterator conditions_iterator =
                        rThisMesh.ConditionsBegin();
                    conditions_iterator != rThisMesh.ConditionsEnd(); conditions_iterator++ )
                for ( typename std::vector<TMeshContainer>::iterator it = mGidMeshContainers.begin();
                        it != mGidMeshContainers.end(); it++ )
                    if ( it->AddCondition( conditions_iterator ) )
                        break;
		}
//         mNodeList.clear();

        for ( typename std::vector<TMeshContainer>::iterator it = mGidMeshContainers.begin();
                it != mGidMeshContainers.end(); it++ )
        {
            it->FinalizeMeshCreation();
            if ( mWriteDeformed == WriteDeformed )
                it->WriteMesh(mMeshFile,true);
            else if ( mWriteDeformed == WriteUndeformed )
                it->WriteMesh(mMeshFile,false);
            else
                KRATOS_THROW_ERROR( std::logic_error, "undefined WriteDeformedMeshFlag" , "" );

            ModelPart::NodesContainerType tempNodes = it->GetMeshNodes();
            for( ModelPart::NodesContainerType::iterator iter = tempNodes.begin(); iter != tempNodes.end(); iter++ )
            {
//                 mNodeList.push_back(*iter);

            }

            it->Reset();
        }

//         mNodeList.Unique();



        Timer::Stop("Writing Mesh");

        KRATOS_CATCH("")
    }//WriteMesh


    ///functions for printing results on gauss points
    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param r_model_part the current model part
     */
    virtual void PrintOnGaussPoints( const Variable<double>& rVariable, ModelPart& r_model_part,
                                     double SolutionTag, int value_index = 0 )
    {
        KRATOS_TRY;

        Timer::Start("Writing Results");

        for ( typename std::vector<TGaussPointContainer>::iterator it =
                    mGidGaussPointContainers.begin();
                it != mGidGaussPointContainers.end(); it++ )
        {

            it->PrintResults( mResultFile, rVariable, r_model_part, SolutionTag, value_index );
        }

        Timer::Stop("Writing Results");

        KRATOS_CATCH("");
    }

    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param r_model_part the current model part
     */
    virtual void PrintOnGaussPoints( const Variable<array_1d<double,3> >& rVariable, ModelPart& r_model_part, double SolutionTag, int value_index = 0 )
    {
        KRATOS_TRY;

        Timer::Start("Writing Results");

        for ( typename std::vector<TGaussPointContainer>::iterator it =
                    mGidGaussPointContainers.begin();
                it != mGidGaussPointContainers.end(); it++ )
        {
            it->PrintResults(  mResultFile, rVariable, r_model_part, SolutionTag, value_index );
        }

        Timer::Stop("Writing Results");

        KRATOS_CATCH("");
    }

    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param r_model_part the current model part
     */
    virtual void PrintOnGaussPoints( const Variable<Vector>& rVariable, ModelPart& r_model_part,
                                     double SolutionTag, int value_index = 0 )
    {
        KRATOS_TRY;
        Timer::Start("Writing Results");

        for ( typename std::vector<TGaussPointContainer>::iterator it =
                    mGidGaussPointContainers.begin();
                it != mGidGaussPointContainers.end(); it++ )
        {
            it->PrintResults(  mResultFile, rVariable, r_model_part, SolutionTag, value_index );

        }

        Timer::Stop("Writing Results");

        KRATOS_CATCH("");
    }

    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param r_model_part the current model part
     */
    virtual void PrintOnGaussPoints( const Variable<Matrix>& rVariable, ModelPart& r_model_part,
                                     double SolutionTag, int value_index = 0 )
    {
        KRATOS_TRY;
        Timer::Start("Writing Results");
        for ( typename std::vector<TGaussPointContainer>::iterator it =
                    mGidGaussPointContainers.begin();
                it != mGidGaussPointContainers.end(); it++ )
        {
         
            it->PrintResults(  mResultFile, rVariable, r_model_part, SolutionTag, value_index );
        }

        Timer::Stop("Writing Results");

        KRATOS_CATCH("");
    }

protected:
    /**
     * File names
     */
    std::string mResultFileName;
    std::string mMeshFileName;
    
    GiD_FILE mMeshFile;
    GiD_FILE mResultFile;

    /**
     * Flags
     */
    WriteDeformedMeshFlag mWriteDeformed;
    WriteConditionsFlag mWriteConditions;
    MultiFileFlag mUseMultiFile;
    GiD_PostMode mMode;

    /**
     * member variables
     */
    std::vector<TMeshContainer> mGidMeshContainers;
    std::vector<TGaussPointContainer> mGidGaussPointContainers;
    bool mMeshFileOpen;
    bool mResultFileOpen;
//     ModelPart::NodesContainerType mNodeList;

private:
    /**
     * assignment operator
     */
    ExplicitGidIO& operator=(ExplicitGidIO const& rOther);

    /**
     * Copy constructor
     */
    ExplicitGidIO(ExplicitGidIO const& rOther);
}; // Class ExplicitGidIO


///**
// * Input and output
// */
///*    ExplicitGidIO& operator >> (ExplicitGidIO& rInput, IO::NodeType& rNode)
//    {
//        rInput.ReadNode(rNode);
//        return rInput;
//    }

//    ExplicitGidIO& operator >> (ExplicitGidIO& rInput, IO::NodesContainerType& rNodes)
//    {
//        rInput.ReadNodes(rNodes);
//        return rInput;
//    }

//    ExplicitGidIO& operator >> (ExplicitGidIO& rInput, IO::PropertiesContainerType& rProperties)
//    {
//        rInput.ReadProperties(rProperties);
//        return rInput;
//    }

//    ExplicitGidIO& operator >> (ExplicitGidIO& rInput, IO::MeshType& rMesh)
//    {
//        rInput.ReadMesh(rMesh);
//        return rInput;
//    }

//    ExplicitGidIO& operator << (ExplicitGidIO& rOutput, IO::NodesContainerType& rNodes)
//    {
//        rOutput.WriteNodes(rNodes);
//        return rOutput;
//    }

//    ExplicitGidIO& operator << (ExplicitGidIO& rOutput, IO::ElementsContainerType& rElements)
//    {
//        rOutput.WriteElements(rElements);
//        return rOutput;
//    }*/

///**
// * output stream function
// */
//inline std::ostream& operator << (std::ostream& rOStream, const ExplicitGidIO<>& rThis)
//{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);
//    return rOStream;
//}

}// namespace Kratos.

#undef KRATOS_INDEX_PARSER
#undef KRATOS_COORDINATES_PARSER
#undef KRATOS_NODE_PARSER
#undef KRATOS_NODE_DATA_PARSER
#undef KRATOS_PROPERTIES_LHS_PARSER
#undef KRATOS_PROPERTIES_TEMPORARY_VARIABLES
#undef KRATOS_ARRAY_1D_3_PARSER
#undef KRATOS_VECTOR_PARSER
#undef KRATOS_MATRIX_PARSER
#undef KRATOS_CONDITIONS_TEMPORARY_VARIABLES
#undef KRATOS_CONDITIONS_FIX_PARSER

#endif // KRATOS_GID_IO_BASE_H_INCLUDED  defined 

