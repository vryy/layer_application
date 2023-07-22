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
//   Date:                $Date: 8 Sep 2016 $
//   Revision:            $Revision: 1.0 $
//
//







#if !defined(KRATOS_SD_TIKZ_POST_IO_H_INCLUDED)
#define  KRATOS_SD_TIKZ_POST_IO_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

// External includes

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
class SDTikzPostIO : public SDPostIO<TGaussPointContainer, TMeshContainer>
{
public:

    /// Type definitions
    KRATOS_CLASS_POINTER_DEFINITION(SDTikzPostIO);

    typedef SDPostIO<TGaussPointContainer, TMeshContainer> BaseType;

    typedef typename BaseType::MeshType MeshType;

    typedef typename BaseType::NodesContainerType NodesContainerType;

    typedef typename BaseType::PropertiesContainerType PropertiesContainerType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    typedef typename BaseType::ConditionsContainerType ConditionsContainerType;


    ///Constructor
    ///single stream IO constructor
    SDTikzPostIO( const std::string& rDatafilename,
        WriteDeformedMeshFlag write_deformed_flag,
        WriteConditionsFlag write_conditions_flag )
    : BaseType(rDatafilename)
    {
        mWriteDeformed = write_deformed_flag;
        mWriteConditions = write_conditions_flag;
        mWriteId = true;
        mAlpha = 60.0;
        mBeta = 110.0;
        mGamma = 0.0;
    }

    ///Destructor.
    virtual ~SDTikzPostIO()
    {
    }

    /***************************************************************************************************/
    /***************************************************************************************************/
    /***************************************************************************************************/

    void SetWriteId(bool WriteId)
    {
        mWriteId = WriteId;
    }

    void SetCamera(double Alpha, double Beta, double Gamma)
    {
        mAlpha = Alpha;
        mBeta = Beta;
        mGamma = Gamma;
    }

    void SetStyle(std::string name, std::string style)
    {
        mStyles[name] = style;
    }

    void SetCurrentStyle(std::string NodeStyle, double NodeSize, std::string ElementStyle, std::string ConditionStyle)
    {
        mCurrentNodeStyle = NodeStyle;
        mCurrentNodeSize = NodeSize;
        mCurrentElementStyle = ElementStyle;
        mCurrentConditionStyle = ConditionStyle;
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
        std::stringstream ss;
        ss << BaseType::mMeshFileName << "_" << name << ".tex";

        mMeshFile.open(ss.str().c_str());

        mMeshFile << "\\documentclass[]{article}" << std::endl;
        mMeshFile << "\\usepackage{tikz,tikz-3dplot}" << std::endl;
        mMeshFile << "%%%<" << std::endl;
        mMeshFile << "\\usepackage[active,tightpage]{preview}" << std::endl;
        mMeshFile << "\\PreviewEnvironment{tikzpicture}" << std::endl;
        mMeshFile << "\\setlength\\PreviewBorder{5pt}%" << std::endl;
        mMeshFile << "%%%>" << std::endl;

        mMeshFile << "\\usetikzlibrary{calc,3d}" << std::endl;

        mMeshFile << "\\begin{document}" << std::endl;
        mMeshFile << "\\tdplotsetmaincoords{" << mAlpha << "}{" << mBeta << "}{" << mGamma << "}" << std::endl;
        mMeshFile << "\\thispagestyle{empty}" << std::endl;

        mMeshFile << "\\begin{center}" << std::endl;
        mMeshFile << "\\begin{tikzpicture}[tdplot_main_coords]" << std::endl;

        mMeshFile << "\\tikzset{" << std::endl;
        for(std::map<std::string, std::string>::iterator it = mStyles.begin(); it != mStyles.end(); ++it)
        {
            mMeshFile << "," << it->first << "/.style={" << it->second << "}" << std::endl;
        }
        mMeshFile << "}" << std::endl;

        std::cout << "initialized tikz mesh on " << ss.str() << std::endl;
    }

    /**
     * closes a mesh group
     */
    void FinalizeMesh()
    {
        mMeshFile << "\\end{tikzpicture}" << std::endl;
        mMeshFile << "\\end{center}" << std::endl;
        mMeshFile << "\\end{document}" << std::endl;
        mMeshFile.close();
        std::cout << "finalized tikz mesh" << std::endl;
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
    }



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
            for ( typename MeshType::ElementIterator element_iterator = rThisMesh.ElementsBegin();
                    element_iterator != rThisMesh.ElementsEnd(); ++element_iterator)
                for ( typename std::vector<TMeshContainer>::iterator it = BaseType::mMeshContainers.begin();
                        it != BaseType::mMeshContainers.end(); it++ )
                    if ( it->AddElement( element_iterator ) )
                        break;
        }
        if ( mWriteConditions == WriteConditions || mWriteConditions == WriteConditionsOnly )
        {
            for ( typename MeshType::ConditionsContainerType::iterator conditions_iterator =
                        rThisMesh.ConditionsBegin();
                    conditions_iterator != rThisMesh.ConditionsEnd(); ++conditions_iterator )
                for ( typename std::vector<TMeshContainer>::iterator it = BaseType::mMeshContainers.begin();
                        it != BaseType::mMeshContainers.end(); it++ )
                    if ( it->AddCondition( conditions_iterator ) )
                        break;
        }

        for ( typename std::vector<TMeshContainer>::iterator it = BaseType::mMeshContainers.begin();
                it != BaseType::mMeshContainers.end(); ++it )
        {
            it->FinalizeMeshCreation();
            if ( mWriteDeformed == WriteDeformed )
                it->WriteMesh(mMeshFile, true, mWriteId, mCurrentNodeStyle, mCurrentNodeSize, mCurrentElementStyle, mCurrentConditionStyle);
            else if ( mWriteDeformed == WriteUndeformed )
                it->WriteMesh(mMeshFile, false, mWriteId, mCurrentNodeStyle, mCurrentNodeSize, mCurrentElementStyle, mCurrentConditionStyle);
            else
                KRATOS_THROW_ERROR( std::logic_error, "undefined WriteDeformedMeshFlag" , "" );

            it->Reset();
        }

        Timer::Stop("Writing Mesh");

        KRATOS_CATCH("")
    }


protected:

    /// File handles
    std::ofstream mMeshFile;

    /// member variables
    WriteDeformedMeshFlag mWriteDeformed;
    WriteConditionsFlag mWriteConditions;
    bool mWriteId;

    /// camera angle
    double mAlpha;
    double mBeta;
    double mGamma;

    /// list of styles
    std::map<std::string, std::string> mStyles;

    /// current style
    std::string mCurrentNodeStyle;
    double mCurrentNodeSize;
    std::string mCurrentElementStyle;
    std::string mCurrentConditionStyle;

private:
    /**
     * assignment operator
     */
    SDTikzPostIO& operator=(SDTikzPostIO const& rOther);

    /**
     * Copy constructor
     */
    SDTikzPostIO(SDTikzPostIO const& rOther);
}; // Class SDTikzPostIO


///**
// * Input and output
// */
//SDTikzPostIO& operator >> (SDTikzPostIO& rInput, IO::NodeType& rNode)
//{
//    rInput.ReadNode(rNode);
//    return rInput;
//}

//SDTikzPostIO& operator >> (SDTikzPostIO& rInput, IO::NodesContainerType& rNodes)
//{
//    rInput.ReadNodes(rNodes);
//    return rInput;
//}

//SDTikzPostIO& operator >> (SDTikzPostIO& rInput, IO::PropertiesContainerType& rProperties)
//{
//    rInput.ReadProperties(rProperties);
//    return rInput;
//}

//SDTikzPostIO& operator >> (SDTikzPostIO& rInput, IO::MeshType& rMesh)
//{
//    rInput.ReadMesh(rMesh);
//    return rInput;
//}

//SDTikzPostIO& operator << (SDTikzPostIO& rOutput, IO::NodesContainerType& rNodes)
//{
//    rOutput.WriteNodes(rNodes);
//    return rOutput;
//}

//SDTikzPostIO& operator << (SDTikzPostIO& rOutput, IO::ElementsContainerType& rElements)
//{
//    rOutput.WriteElements(rElements);
//    return rOutput;
//}

///**
// * output stream function
// */
//inline std::ostream& operator << (std::ostream& rOStream, const SDTikzPostIO<>& rThis)
//{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);
//    return rOStream;
//}

}// namespace Kratos.

#endif // KRATOS_SD_TIKZ_POST_IO_H_INCLUDED  defined

