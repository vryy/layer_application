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







#if !defined(KRATOS_SD_GID_POST_IO_H_INCLUDED)
#define  KRATOS_SD_GID_POST_IO_H_INCLUDED

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
#include "includes/mesh.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"
#include "utilities/timer.h"

#include "custom_io/sd_post_io.h"

namespace Kratos
{


/**
 */
template<class TGaussPointContainer, class TMeshContainer>
class SDGidPostIO : public SDPostIO<TGaussPointContainer, TMeshContainer>
{
public:

    /// Type definitions
    KRATOS_CLASS_POINTER_DEFINITION(SDGidPostIO);

    typedef SDPostIO<TGaussPointContainer, TMeshContainer> BaseType;

    typedef Node<3> NodeType;

    typedef Mesh<NodeType, Properties, Element, Condition> MeshType;

    typedef MeshType::NodesContainerType NodesContainerType;

    typedef MeshType::PropertiesContainerType PropertiesContainerType;

    typedef MeshType::ElementsContainerType ElementsContainerType;

    typedef MeshType::ConditionsContainerType ConditionsContainerType;

    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef GeometryData::IntegrationMethod IntegrationMethodType;
    typedef GeometryData::KratosGeometryFamily KratosGeometryFamily;

    typedef typename BaseType::MeshContainerVectorType MeshContainerVectorType;
    typedef typename BaseType::GaussPointContainerVectorType GaussPointContainerVectorType;

    ///Constructor
    ///single stream IO constructor
    SDGidPostIO( const std::string& rDatafilename,
           GiD_PostMode Mode,
           MultiFileFlag use_multiple_files_flag,
           WriteDeformedMeshFlag write_deformed_flag,
           WriteConditionsFlag write_conditions_flag )
    : BaseType(rDatafilename)
    {
        mMode = Mode;
        mResultFileOpen = false;
        mMeshFileOpen = false;
        mWriteDeformed = write_deformed_flag;
        mWriteConditions = write_conditions_flag;
        mUseMultiFile = use_multiple_files_flag;
        InitializeResultFile(BaseType::mResultFileName);
    }

    ///Destructor.
    virtual ~SDGidPostIO()
    {
        Timer::PrintTimingInformation();

        if ( mResultFileOpen )
        {
            GiD_fClosePostResultFile( mResultFile );
            mResultFileOpen = false;
        }
    }

    /***************************************************************************************************/
    /***************************************************************************************************/
    /***************************************************************************************************/


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
                file_name << std::setprecision(12) << BaseType::mMeshFileName << "_" << name << ".post.msh";
                mMeshFile = GiD_fOpenPostMeshFile( (char *)(file_name.str()).c_str(), mMode);
                mMeshFileOpen = true;
            }
            if ( (mMode == GiD_PostBinary || mMode == GiD_PostHDF5) && ! mResultFileOpen )
            {
                std::stringstream file_name;
                file_name << std::setprecision(12) << BaseType::mResultFileName << "_" << name << ".post.bin";
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
                file_name << BaseType::mResultFileName << ".post.bin";
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
                file_name << BaseType::mMeshFileName << "_" << name << ".post.msh";
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
    void WriteNodeMesh( const MeshType& rThisMesh )
    {
        KRATOS_TRY

        Timer::Start("Writing Mesh");

        GiD_fBeginMesh(mMeshFile, "Kratos Mesh", GiD_3D, GiD_Point, 1);

        GiD_fBeginCoordinates(mMeshFile);
        for ( MeshType::NodeConstantIterator node_iterator = rThisMesh.NodesBegin();
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
        for ( MeshType::NodeConstantIterator node_iterator = rThisMesh.NodesBegin();
                node_iterator != rThisMesh.NodesEnd();
                ++node_iterator)
        {
            nodes_id[0] = node_iterator->Id();
            GiD_fWriteElement(mMeshFile, node_iterator->Id(), nodes_id);
        }
        GiD_fEndElements(mMeshFile);

        GiD_fEndMesh(mMeshFile);

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
    void WriteMesh( const MeshType& rThisMesh )
    {
        KRATOS_TRY

        Timer::Start("Writing Mesh");

        if ( mWriteConditions != WriteConditionsOnly )
        {
            for ( MeshType::ElementConstantIterator element_iterator = rThisMesh.ElementsBegin();
                    element_iterator != rThisMesh.ElementsEnd(); ++element_iterator)
                for ( typename MeshContainerVectorType::iterator it = BaseType::mMeshContainers.begin();
                        it != BaseType::mMeshContainers.end(); ++it )
                    if ( it->AddElement( element_iterator ) )
                        break;
        }
        if ( mWriteConditions == WriteConditions || mWriteConditions == WriteConditionsOnly )
        {
            for ( MeshType::ConditionConstantIterator conditions_iterator =
                        rThisMesh.ConditionsBegin();
                    conditions_iterator != rThisMesh.ConditionsEnd(); ++conditions_iterator )
                for ( typename MeshContainerVectorType::iterator it = BaseType::mMeshContainers.begin();
                        it != BaseType::mMeshContainers.end(); ++it )
                    if ( it->AddCondition( conditions_iterator ) )
                        break;
        }

        for ( typename MeshContainerVectorType::iterator it = BaseType::mMeshContainers.begin();
                it != BaseType::mMeshContainers.end(); ++it )
        {
            it->FinalizeMeshCreation();
            if ( mWriteDeformed == WriteDeformed )
                it->WriteMesh(mMeshFile, true);
            else if ( mWriteDeformed == WriteUndeformed )
                it->WriteMesh(mMeshFile, false);
            else
                KRATOS_THROW_ERROR( std::logic_error, "undefined WriteDeformedMeshFlag" , "" );
        }

        Timer::Stop("Writing Mesh");

        KRATOS_CATCH("")
    }//WriteMesh


    /***************************************************************************************************/
    /***************************************************************************************************/
    /***************************************************************************************************/


    /**
     * sets up the file names and opens the result file in case there
     * is ASCII mode and only one file written
     */
    void InitializeResultFile( std::string const& rResultFileName )
    {
        //std::cout << "initializing result files" << std::endl;
        BaseType::mResultFileName = rResultFileName;
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
        return "sd post io";
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
            file_name << BaseType::mResultFileName << std::setprecision(12) << "_" << name << ".post.res";
            mResultFile = GiD_fOpenPostResultFile((char*)(file_name.str()).c_str(), mMode);
            mResultFileOpen = true;

        }

        // initializing gauss points containers
        if ( mWriteConditions != WriteConditionsOnly )
        {
            for ( MeshType::ElementIterator element_iterator = rThisMesh.ElementsBegin();
                    element_iterator != rThisMesh.ElementsEnd(); ++element_iterator )
            {
                for ( typename GaussPointContainerVectorType::iterator it =
                            BaseType::mGaussPointContainers.begin();
                        it != BaseType::mGaussPointContainers.end(); it++ )
                {
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
                for ( typename GaussPointContainerVectorType::iterator it =
                            BaseType::mGaussPointContainers.begin();
                        it != BaseType::mGaussPointContainers.end(); it++ )
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
    }

    /**
     * writes nodal results for flags using the nodes in the mesh container
     */
    void WriteNodalResults( const char* FlagName, Flags const& rFlag,
                            double SolutionTag)
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, FlagName, "Kratos",
                         SolutionTag, GiD_Scalar,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( typename MeshContainerVectorType::iterator it = BaseType::mMeshContainers.begin();
                        it != BaseType::mMeshContainers.end(); ++it )
        {
            for ( NodesContainerType::const_iterator i_node = it->GetMeshNodes().begin();
                    i_node != it->GetMeshNodes().end(); ++i_node )
                GiD_fWriteScalar( mResultFile, i_node->Id(), i_node->Is(rFlag) );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for flags
     */
    void WriteNodalResults( const char* FlagName, Flags const& rFlag,
                            NodesContainerType& rNodes, double SolutionTag)
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, FlagName, "Kratos",
                         SolutionTag, GiD_Scalar,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node )
            GiD_fWriteScalar( mResultFile, i_node->Id(), i_node->Is(rFlag) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type bool using nodes from mesh container
     */
    void WriteNodalResults( Variable<bool> const& rVariable,
                            double SolutionTag,
                            std::size_t SolutionStepNumber)
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Scalar,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( typename MeshContainerVectorType::iterator it = BaseType::mMeshContainers.begin();
                        it != BaseType::mMeshContainers.end(); ++it )
        {
            for ( NodesContainerType::const_iterator i_node = it->GetMeshNodes().begin();
                    i_node != it->GetMeshNodes().end() ; ++i_node )
                GiD_fWriteScalar( mResultFile, i_node->Id(), i_node->GetSolutionStepValue(rVariable,
                                 SolutionStepNumber) );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

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
                i_node != rNodes.end() ; ++i_node )
            GiD_fWriteScalar( mResultFile, i_node->Id(), i_node->GetSolutionStepValue(rVariable,
                             SolutionStepNumber) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type double using nodes in the mesh container
     */
    void WriteNodalResults( Variable<double> const& rVariable,
                            double SolutionTag,
                            std::size_t SolutionStepNumber)
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Scalar,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( typename MeshContainerVectorType::iterator it = BaseType::mMeshContainers.begin();
                        it != BaseType::mMeshContainers.end(); ++it )
        {
            for ( NodesContainerType::const_iterator i_node = it->GetMeshNodes().begin();
                    i_node != it->GetMeshNodes().end() ; ++i_node )
                GiD_fWriteScalar( mResultFile, i_node->Id(), i_node->GetSolutionStepValue(rVariable,
                                 SolutionStepNumber) );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

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
                i_node != rNodes.end() ; ++i_node )
            GiD_fWriteScalar( mResultFile, i_node->Id(), i_node->GetSolutionStepValue(rVariable,
                             SolutionStepNumber) );
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type array_1d<double, 3> using nodes from the mesh container
     * (e.g. DISPLACEMENT)
     */
    void WriteNodalResults( Variable<array_1d<double, 3> > const& rVariable,
                            double SolutionTag, std::size_t SolutionStepNumber)
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult(mResultFile,(char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Vector,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( typename MeshContainerVectorType::iterator it = BaseType::mMeshContainers.begin();
                        it != BaseType::mMeshContainers.end(); ++it )
        {
            for ( NodesContainerType::const_iterator i_node = it->GetMeshNodes().begin();
                    i_node != it->GetMeshNodes().end() ; ++i_node )
            {
                array_1d<double, 3>& temp = i_node->GetSolutionStepValue( rVariable,
                                            SolutionStepNumber );
                GiD_fWriteVector( mResultFile, i_node->Id(), temp[0], temp[1], temp[2] );
            }
        }
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
        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node )
        {
            array_1d<double, 3>& temp = i_node->GetSolutionStepValue( rVariable,
                                        SolutionStepNumber );
            GiD_fWriteVector( mResultFile, i_node->Id(), temp[0], temp[1], temp[2] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type Vector using nodes from the mesh container
     * (note that only vectors with 3 or 6 components can be printed)
     */
    void WriteNodalResults( Variable<Vector> const& rVariable,
                            double SolutionTag, std::size_t SolutionStepNumber)
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Matrix,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( typename MeshContainerVectorType::iterator it = BaseType::mMeshContainers.begin();
                        it != BaseType::mMeshContainers.end(); ++it )
        {
            for ( NodesContainerType::const_iterator i_node = it->GetMeshNodes().begin();
                    i_node != it->GetMeshNodes().end() ; ++i_node )
            {
                Vector& tempVector = i_node->GetSolutionStepValue(rVariable,
                                     SolutionStepNumber);
                if (tempVector.size() == 3 )
                    GiD_fWriteVector(mResultFile, i_node->Id(), tempVector(0), tempVector(1), tempVector(2) );
                else if (tempVector.size() == 6 )
                    GiD_fWrite3DMatrix( mResultFile, i_node->Id(), tempVector(0), tempVector(1), tempVector(2),
                                        tempVector(3), tempVector(4), tempVector(5) );
            }
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type Vector
     * (note that only vectors with 3 or 6 components can be printed)
     */
    void WriteNodalResults( Variable<Vector> const& rVariable,
                            NodesContainerType& rNodes,
                            double SolutionTag, std::size_t SolutionStepNumber)
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Matrix,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node )
        {
            Vector& tempVector = i_node->GetSolutionStepValue(rVariable,
                                 SolutionStepNumber);
            if (tempVector.size() == 3 )
                GiD_fWriteVector(mResultFile, i_node->Id(), tempVector(0), tempVector(1), tempVector(2) );
            else if (tempVector.size() == 6 )
                GiD_fWrite3DMatrix( mResultFile, i_node->Id(), tempVector(0), tempVector(1), tempVector(2),
                                    tempVector(3), tempVector(4), tempVector(5) );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type Matrix using nodes from the mesh container
     */
    void WriteNodalResults( Variable<Matrix> const& rVariable,
                            double SolutionTag, std::size_t SolutionStepNumber)
    {
        Timer::Start("Writing Results");

        GiD_fBeginResult( mResultFile, (char*)(rVariable.Name().c_str()), "Kratos",
                         SolutionTag, GiD_Matrix,
                         GiD_OnNodes, NULL, NULL, 0, NULL );
        for ( typename MeshContainerVectorType::iterator it = BaseType::mMeshContainers.begin();
                        it != BaseType::mMeshContainers.end(); ++it )
        {
            for ( NodesContainerType::const_iterator i_node = it->GetMeshNodes().begin();
                    i_node != it->GetMeshNodes().end() ; ++i_node )
            {
                Matrix& tempMatrix = i_node->GetSolutionStepValue(rVariable, SolutionStepNumber);
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
            }
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
        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node )
        {
            Matrix& tempMatrix = i_node->GetSolutionStepValue(rVariable, SolutionStepNumber);
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

        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node )
        {
            array_1d<double, 3>& temp = i_node->GetSolutionStepValue( rVariable,
                                        SolutionStepNumber );
            GiD_fWriteLocalAxes( mResultFile, i_node->Id(), temp[0], temp[1], temp[2] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }

    ///////////////////////////////////////////////////////////////////////
    //////                 NON- HISTORICAL DATABASE BLOCK             /////
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
                i_node != rNodes.end() ; ++i_node )
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
                i_node != rNodes.end() ; ++i_node )
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
        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node )
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
        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node )
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
        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node )
        {
            Matrix& tempMatrix = i_node->GetSolutionStepValue(rVariable);
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
        for ( NodesContainerType::iterator i_node = rNodes.begin();
                i_node != rNodes.end() ; ++i_node )
        {
            array_1d<double, 3>& temp = i_node->GetSolutionStepValue( rVariable);
            GiD_fWriteLocalAxes( mResultFile, i_node->Id(), temp[0], temp[1], temp[2] );
        }
        GiD_fEndResult(mResultFile);

        Timer::Stop("Writing Results");
    }


    /***************************************************************************************************/
    /***************************************************************************************************/
    /***************************************************************************************************/


    ///functions for printing results on gauss points

    /**
     * Prints variables of type double on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param r_model_part the current model part
     */
    void PrintElementalPartitionIndex( const Variable<double>& rVariable, ModelPart& r_model_part,
                                     double SolutionTag, int value_index, int rank )
    {
        KRATOS_TRY;

        if(rVariable.Name() == "PARTITION_INDEX")
        {
            Timer::Start("Writing Results");

            for ( typename GaussPointContainerVectorType::iterator it =
                        BaseType::mGaussPointContainers.begin();
                    it != BaseType::mGaussPointContainers.end(); it++ )
            {
                it->PrintPartitionIndex( mResultFile, rVariable, r_model_part, SolutionTag, value_index, rank );
            }

            Timer::Stop("Writing Results");
        }

        KRATOS_CATCH("");
    }

    /**
     * Prints variables of type int on gauss points of the complete mesh
     * @param rVariable the given variable name
     * @param r_model_part the current model part
     */
    virtual void PrintOnGaussPoints( const Variable<int>& rVariable, ModelPart& r_model_part,
                                     double SolutionTag, int value_index = 0 )
    {
        KRATOS_TRY;

        Timer::Start("Writing Results");

        for ( typename GaussPointContainerVectorType::iterator it =
                    BaseType::mGaussPointContainers.begin();
                it != BaseType::mGaussPointContainers.end(); it++ )
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
    virtual void PrintOnGaussPoints( const Variable<double>& rVariable, ModelPart& r_model_part,
                                     double SolutionTag, int value_index = 0 )
    {
        KRATOS_TRY;

        Timer::Start("Writing Results");

        for ( typename GaussPointContainerVectorType::iterator it =
                    BaseType::mGaussPointContainers.begin();
                it != BaseType::mGaussPointContainers.end(); it++ )
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

        for ( typename GaussPointContainerVectorType::iterator it =
                    BaseType::mGaussPointContainers.begin();
                it != BaseType::mGaussPointContainers.end(); it++ )
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

        for ( typename GaussPointContainerVectorType::iterator it =
                    BaseType::mGaussPointContainers.begin();
                it != BaseType::mGaussPointContainers.end(); it++ )
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
        for ( typename GaussPointContainerVectorType::iterator it =
                    BaseType::mGaussPointContainers.begin();
                it != BaseType::mGaussPointContainers.end(); it++ )
        {
            it->PrintResults(  mResultFile, rVariable, r_model_part, SolutionTag, value_index );
        }

        Timer::Stop("Writing Results");

        KRATOS_CATCH("");
    }

protected:
    /**
     * File pointers
     */
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
    bool mMeshFileOpen;
    bool mResultFileOpen;

private:
    /**
     * assignment operator
     */
    SDGidPostIO& operator=(SDGidPostIO const& rOther);

    /**
     * Copy constructor
     */
    SDGidPostIO(SDGidPostIO const& rOther);
}; // Class SDGidPostIO


///**
// * Input and output
// */
//SDGidPostIO& operator >> (SDGidPostIO& rInput, IO::NodeType& rNode)
//{
//    rInput.ReadNode(rNode);
//    return rInput;
//}

//SDGidPostIO& operator >> (SDGidPostIO& rInput, IO::NodesContainerType& rNodes)
//{
//    rInput.ReadNodes(rNodes);
//    return rInput;
//}

//SDGidPostIO& operator >> (SDGidPostIO& rInput, IO::PropertiesContainerType& rProperties)
//{
//    rInput.ReadProperties(rProperties);
//    return rInput;
//}

//SDGidPostIO& operator >> (SDGidPostIO& rInput, IO::MeshType& rMesh)
//{
//    rInput.ReadMesh(rMesh);
//    return rInput;
//}

//SDGidPostIO& operator << (SDGidPostIO& rOutput, IO::NodesContainerType& rNodes)
//{
//    rOutput.WriteNodes(rNodes);
//    return rOutput;
//}

//SDGidPostIO& operator << (SDGidPostIO& rOutput, IO::ElementsContainerType& rElements)
//{
//    rOutput.WriteElements(rElements);
//    return rOutput;
//}

///**
// * output stream function
// */
//inline std::ostream& operator << (std::ostream& rOStream, const SDGidPostIO<>& rThis)
//{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);
//    return rOStream;
//}

}// namespace Kratos.

#endif // KRATOS_SD_GID_POST_IO_H_INCLUDED  defined

