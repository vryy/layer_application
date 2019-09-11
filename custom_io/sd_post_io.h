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







#if !defined(KRATOS_SD_POST_IO_H_INCLUDED)
#define  KRATOS_SD_POST_IO_H_INCLUDED

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
//#include "includes/io.h"
#include "includes/mesh.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

#include "utilities/timer.h"

namespace Kratos
{

/**
 * Type definitions
 */



///Flags for mesh writing
enum WriteDeformedMeshFlag {WriteDeformed, WriteUndeformed};
enum WriteConditionsFlag {WriteConditions, WriteElementsOnly, WriteConditionsOnly};
enum MultiFileFlag {SingleFile, MultipleFiles};



/**
 * This class defines an interface to the post processing library
 * In particular, it can be used to interface with GiDPost, in order to provide GiD compliant I/O functionality
 */
template<class TGaussPointContainer, class TMeshContainer>
class SDPostIO// : public IO
{
public:

    /// Type definitions
    KRATOS_CLASS_POINTER_DEFINITION(SDPostIO);

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

    typedef std::vector<TMeshContainer> MeshContainerVectorType;
    typedef std::vector<TGaussPointContainer> GaussPointContainerVectorType;

    ///Constructor
    ///single stream IO constructor
    SDPostIO( const std::string& rDatafilename)
    {
        mResultFileName = rDatafilename;
        mMeshFileName = rDatafilename;
        SetUpMeshContainers(); // setup default mesh containers
        SetUpGaussPointContainers(); // setup default integration point containers
    }

    ///Destructor.
    virtual ~SDPostIO()
    {
    }

    ///initialization functions
    /**
     * creates the mesh containers for all different element types.
     * Note that the containers are not filled yet in here!
     */
    virtual void SetUpMeshContainers()
    {
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Hexahedra3D20, "Kratos_Hexahedra3D20_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Hexahedra3D27, "Kratos_Hexahedra3D27_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Hexahedra3D8, "Kratos_Hexahedra3D8_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Prism3D15, "Kratos_Prism3D15_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Prism3D6, "Kratos_Prism3D6_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Quadrilateral2D4, "Kratos_Quadrilateral2D4_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Quadrilateral2D8, "Kratos_Quadrilateral2D8_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Quadrilateral2D9, "Kratos_Quadrilateral2D9_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Quadrilateral3D4, "Kratos_Quadrilateral3D4_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Quadrilateral3D8, "Kratos_Quadrilateral3D8_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Quadrilateral3D9, "Kratos_Quadrilateral3D9_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Tetrahedra3D10, "Kratos_Tetrahedra3D10_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Tetrahedra3D4, "Kratos_Tetrahedra3D4_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Triangle2D3, "Kratos_Triangle2D3_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Triangle2D6, "Kratos_Triangle2D6_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Triangle3D3, "Kratos_Triangle3D3_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Triangle3D6, "Kratos_Triangle3D6_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Line2D2, "Kratos_Line2D2_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Line3D2, "Kratos_Line3D2_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Line2D3, "Kratos_Line2D3_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Line3D3, "Kratos_Line3D3_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Point3D, "Kratos_Point3D_Mesh" ) );
        mMeshContainers.push_back( TMeshContainer( GeometryData::Kratos_Point2D, "Kratos_Point2D_Mesh" ) );
    }//SetUpMeshContainers

    /**
     * creates the gauss point containers for all different element types.
     * Note that the containers are not filled yet in here!
     */
    virtual void SetUpGaussPointContainers()
    {
        //case Line
        mGaussPointContainers.push_back( TGaussPointContainer( "lin_gauss_legendre_1_element_gp", GeometryData::Kratos_Linear, GeometryData::GI_GAUSS_1 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "lin_gauss_legendre_2_element_gp", GeometryData::Kratos_Linear, GeometryData::GI_GAUSS_2 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "lin_gauss_legendre_3_element_gp", GeometryData::Kratos_Linear, GeometryData::GI_GAUSS_3 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "lin_gauss_legendre_4_element_gp", GeometryData::Kratos_Linear, GeometryData::GI_GAUSS_4 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "lin_gauss_legendre_5_element_gp", GeometryData::Kratos_Linear, GeometryData::GI_GAUSS_5 ) );

        //case Triangle
        mGaussPointContainers.push_back( TGaussPointContainer( "tri_gauss_legendre_1_element_gp", GeometryData::Kratos_Triangle, GeometryData::GI_GAUSS_1 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "tri_gauss_legendre_2_element_gp", GeometryData::Kratos_Triangle, GeometryData::GI_GAUSS_2 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "tri_gauss_legendre_3_element_gp", GeometryData::Kratos_Triangle, GeometryData::GI_GAUSS_3 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "tri_gauss_legendre_4_element_gp", GeometryData::Kratos_Triangle, GeometryData::GI_GAUSS_4 ) );

        //case Quadrilateral
        mGaussPointContainers.push_back( TGaussPointContainer( "quad_gauss_legendre_1_element_gp", GeometryData::Kratos_Quadrilateral, GeometryData::GI_GAUSS_1 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "quad_gauss_legendre_2_element_gp", GeometryData::Kratos_Quadrilateral, GeometryData::GI_GAUSS_2 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "quad_gauss_legendre_3_element_gp", GeometryData::Kratos_Quadrilateral, GeometryData::GI_GAUSS_3 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "quad_gauss_legendre_4_element_gp", GeometryData::Kratos_Quadrilateral, GeometryData::GI_GAUSS_4 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "quad_gauss_legendre_5_element_gp", GeometryData::Kratos_Quadrilateral, GeometryData::GI_GAUSS_5 ) );

        //case Tetrahedra
        mGaussPointContainers.push_back( TGaussPointContainer( "tet_gauss_legendre_1_element_gp", GeometryData::Kratos_Tetrahedra, GeometryData::GI_GAUSS_1 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "tet_gauss_legendre_2_element_gp", GeometryData::Kratos_Tetrahedra, GeometryData::GI_GAUSS_2 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "tet_gauss_legendre_3_element_gp", GeometryData::Kratos_Tetrahedra, GeometryData::GI_GAUSS_3 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "tet_gauss_legendre_4_element_gp", GeometryData::Kratos_Tetrahedra, GeometryData::GI_GAUSS_4 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "tet_gauss_legendre_5_element_gp", GeometryData::Kratos_Tetrahedra, GeometryData::GI_GAUSS_5 ) );

        //case Hexahedra
        mGaussPointContainers.push_back( TGaussPointContainer( "hex_gauss_legendre_1_element_gp", GeometryData::Kratos_Hexahedra, GeometryData::GI_GAUSS_1 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "hex_gauss_legendre_2_element_gp", GeometryData::Kratos_Hexahedra, GeometryData::GI_GAUSS_2 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "hex_gauss_legendre_3_element_gp", GeometryData::Kratos_Hexahedra, GeometryData::GI_GAUSS_3 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "hex_gauss_legendre_4_element_gp", GeometryData::Kratos_Hexahedra, GeometryData::GI_GAUSS_4 ) );
        mGaussPointContainers.push_back( TGaussPointContainer( "hex_gauss_legendre_5_element_gp", GeometryData::Kratos_Hexahedra, GeometryData::GI_GAUSS_5 ) );
    }//SetUpGaussPointContainers

    ///general SDPostIO related functions
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

protected:
    /**
     * File names
     */
    std::string mResultFileName;
    std::string mMeshFileName;

    /**
     * member variables
     */
    MeshContainerVectorType mMeshContainers;
    GaussPointContainerVectorType mGaussPointContainers;

private:
    /**
     * assignment operator
     */
    SDPostIO& operator=(SDPostIO const& rOther);

    /**
     * Copy constructor
     */
    SDPostIO(SDPostIO const& rOther);
}; // Class SDPostIO


///**
// * Input and output
// */
//SDPostIO& operator >> (SDPostIO& rInput, IO::NodeType& rNode)
//{
//    rInput.ReadNode(rNode);
//    return rInput;
//}

//SDPostIO& operator >> (SDPostIO& rInput, IO::NodesContainerType& rNodes)
//{
//    rInput.ReadNodes(rNodes);
//    return rInput;
//}

//SDPostIO& operator >> (SDPostIO& rInput, IO::PropertiesContainerType& rProperties)
//{
//    rInput.ReadProperties(rProperties);
//    return rInput;
//}

//SDPostIO& operator >> (SDPostIO& rInput, IO::MeshType& rMesh)
//{
//    rInput.ReadMesh(rMesh);
//    return rInput;
//}

//SDPostIO& operator << (SDPostIO& rOutput, IO::NodesContainerType& rNodes)
//{
//    rOutput.WriteNodes(rNodes);
//    return rOutput;
//}

//SDPostIO& operator << (SDPostIO& rOutput, IO::ElementsContainerType& rElements)
//{
//    rOutput.WriteElements(rElements);
//    return rOutput;
//}

///**
// * output stream function
// */
//inline std::ostream& operator << (std::ostream& rOStream, const SDPostIO<>& rThis)
//{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);
//    return rOStream;
//}

}// namespace Kratos.

#endif // KRATOS_SD_POST_IO_H_INCLUDED  defined

