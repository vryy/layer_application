//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Oct 30, 2014$
//   Revision:            $Revision: 1.0 $
//
//
//Change log:
//  +   30/10/2014: create layer_application.cpp


// System includes


// External includes


// Project includes
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "layer_application.h"
#include "layer_application_variables.h"

#ifdef SD_APP_FORWARD_COMPATIBILITY
#define LAYER_APP_CREATE_ELEMENT(element_type, geometry_type, number_of_nodes) \
    element_type( 0, Element::GeometryType::Pointer( new geometry_type <Node>( Element::GeometryType::PointsArrayType( number_of_nodes ) ) ) )
#define LAYER_APP_CREATE_CONDITION(condition_type, geometry_type, number_of_nodes) \
    condition_type( 0, Condition::GeometryType::Pointer( new geometry_type <Node>( Condition::GeometryType::PointsArrayType( number_of_nodes ) ) ) )
#else
#define LAYER_APP_CREATE_ELEMENT(element_type, geometry_type, number_of_nodes) \
    element_type( 0, Element::GeometryType::Pointer( new geometry_type <Node<3> >( Element::GeometryType::PointsArrayType( number_of_nodes, Node<3>() ) ) ) )
#define LAYER_APP_CREATE_CONDITION(condition_type, geometry_type, number_of_nodes) \
    condition_type( 0, Condition::GeometryType::Pointer( new geometry_type <Node<3> >( Condition::GeometryType::PointsArrayType( number_of_nodes, Node<3>() ) ) ) )
#endif

namespace Kratos
{

    KratosLayerApplication::KratosLayerApplication()
    #ifdef SD_APP_FORWARD_COMPATIBILITY
    : KratosApplication("KratosLayerApplication")
    #else
    : KratosApplication()
    #endif
    , LAYER_APP_CREATE_ELEMENT( mPostElement2D3N, Triangle2D3, 3 )
    , LAYER_APP_CREATE_ELEMENT( mPostElement2D4N, Quadrilateral2D4, 4 )
    , LAYER_APP_CREATE_ELEMENT( mPostElement2D6N, Triangle2D6, 6 )
    , LAYER_APP_CREATE_ELEMENT( mPostElement2D8N, Quadrilateral2D8, 8 )
    , LAYER_APP_CREATE_ELEMENT( mPostElement2D9N, Quadrilateral2D9, 9 )
    , LAYER_APP_CREATE_ELEMENT( mPostFaceElement3D3N, Triangle3D3, 3 )
    , LAYER_APP_CREATE_ELEMENT( mPostFaceElement3D4N, Quadrilateral3D4, 4 )
    , LAYER_APP_CREATE_ELEMENT( mPostFaceElement3D6N, Triangle3D6, 6 )
    , LAYER_APP_CREATE_ELEMENT( mPostFaceElement3D8N, Quadrilateral3D8, 8 )
    , LAYER_APP_CREATE_ELEMENT( mPostFaceElement3D9N, Quadrilateral3D9, 9 )
    , LAYER_APP_CREATE_ELEMENT( mPostElement3D4N, Tetrahedra3D4, 4 )
    , LAYER_APP_CREATE_ELEMENT( mPostElement3D10N, Tetrahedra3D10, 10 )
    , LAYER_APP_CREATE_ELEMENT( mPostElement3D8N, Hexahedra3D8, 8 )
    , LAYER_APP_CREATE_ELEMENT( mPostElement3D20N, Hexahedra3D20, 20 )
    , LAYER_APP_CREATE_ELEMENT( mPostElement3D27N, Hexahedra3D27, 27 )
    , LAYER_APP_CREATE_ELEMENT( mPostElement3D6N, Prism3D6, 6 )
    , LAYER_APP_CREATE_ELEMENT( mPostElement3D15N, Prism3D15, 15 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSElement2D3N, Triangle2D3, 3 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSElement2D4N, Quadrilateral2D4, 4 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSElement2D6N, Triangle2D6, 6 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSElement2D8N, Quadrilateral2D8, 8 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSElement2D9N, Quadrilateral2D9, 9 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSFaceElement3D3N, Triangle3D3, 3 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSFaceElement3D4N, Quadrilateral3D4, 4 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSFaceElement3D6N, Triangle3D6, 6 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSFaceElement3D8N, Quadrilateral3D8, 8 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSFaceElement3D9N, Quadrilateral3D9, 9 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSElement3D4N, Tetrahedra3D4, 4 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSElement3D10N, Tetrahedra3D10, 10 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSElement3D8N, Hexahedra3D8, 8 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSElement3D20N, Hexahedra3D20, 20 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSElement3D27N, Hexahedra3D27, 27 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSElement3D6N, Prism3D6, 6 )
    , LAYER_APP_CREATE_ELEMENT( mPostUPSElement3D15N, Prism3D15, 15 )
    , LAYER_APP_CREATE_CONDITION( mPostSurfaceCondition3D3N, Triangle3D3, 3 )
    , LAYER_APP_CREATE_CONDITION( mPostSurfaceCondition3D6N, Triangle3D6, 6 )
    , LAYER_APP_CREATE_CONDITION( mPostSurfaceCondition3D4N, Quadrilateral3D4, 4 )
    , LAYER_APP_CREATE_CONDITION( mPostSurfaceCondition3D8N, Quadrilateral3D8, 8 )
    , LAYER_APP_CREATE_CONDITION( mPostSurfaceCondition3D9N, Quadrilateral3D9, 9 )
    {}

    void KratosLayerApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosLayerApplication... " << std::endl;

        // register variables to Kratos kernel
        KRATOS_REGISTER_VARIABLE(LAYER_ENTITY_TYPE)
        KRATOS_REGISTER_VARIABLE(LAYER_ENTITY_NAME)
        KRATOS_REGISTER_VARIABLE(LAYER_NAME)
        KRATOS_REGISTER_VARIABLE(LAYER_PROP_ID)
        #ifdef LAYER_APP_USE_MMG
        KRATOS_REGISTER_VARIABLE(NODAL_MMG_SCALAR_METRIC)
        KRATOS_REGISTER_VARIABLE(NODAL_MMG_VECTOR_METRIC)
        KRATOS_REGISTER_VARIABLE(NODAL_MMG_TENSOR_METRIC)
        KRATOS_REGISTER_VARIABLE(NODAL_MMG_LEVEL_SET)
        KRATOS_REGISTER_VARIABLE(MMG_GRADATION)
        KRATOS_REGISTER_VARIABLE(MMG_HAUSDORFF_DISTANCE)
        KRATOS_REGISTER_VARIABLE(MMG_MINIMAL_MESH_SIZE)
        KRATOS_REGISTER_VARIABLE(MMG_MAXIMAL_MESH_SIZE)
        KRATOS_REGISTER_VARIABLE(MMG_CONSTANT_MESH_SIZE)
        KRATOS_REGISTER_VARIABLE(MMG_RMC_VOLUME_FRACTION)
        #endif

        // register element to the kernel
        KRATOS_REGISTER_ELEMENT( "PostElement2D3N", mPostElement2D3N )
        KRATOS_REGISTER_ELEMENT( "PostElement2D4N", mPostElement2D4N )
        KRATOS_REGISTER_ELEMENT( "PostElement2D6N", mPostElement2D6N )
        KRATOS_REGISTER_ELEMENT( "PostElement2D8N", mPostElement2D8N )
        KRATOS_REGISTER_ELEMENT( "PostElement2D9N", mPostElement2D9N )
        KRATOS_REGISTER_ELEMENT( "PostFaceElement3D3N", mPostFaceElement3D3N )
        KRATOS_REGISTER_ELEMENT( "PostFaceElement3D4N", mPostFaceElement3D4N )
        KRATOS_REGISTER_ELEMENT( "PostFaceElement3D6N", mPostFaceElement3D6N )
        KRATOS_REGISTER_ELEMENT( "PostFaceElement3D8N", mPostFaceElement3D8N )
        KRATOS_REGISTER_ELEMENT( "PostFaceElement3D9N", mPostFaceElement3D9N )
        KRATOS_REGISTER_ELEMENT( "PostElement3D4N", mPostElement3D4N )
        KRATOS_REGISTER_ELEMENT( "PostElement3D10N", mPostElement3D10N )
        KRATOS_REGISTER_ELEMENT( "PostElement3D8N", mPostElement3D8N )
        KRATOS_REGISTER_ELEMENT( "PostElement3D20N", mPostElement3D20N )
        KRATOS_REGISTER_ELEMENT( "PostElement3D27N", mPostElement3D27N )
        KRATOS_REGISTER_ELEMENT( "PostElement3D6N", mPostElement3D6N )
        KRATOS_REGISTER_ELEMENT( "PostElement3D15N", mPostElement3D15N )

        KRATOS_REGISTER_ELEMENT( "PostUPSElement2D3N", mPostUPSElement2D3N )
        KRATOS_REGISTER_ELEMENT( "PostUPSElement2D4N", mPostUPSElement2D4N )
        KRATOS_REGISTER_ELEMENT( "PostUPSElement2D6N", mPostUPSElement2D6N )
        KRATOS_REGISTER_ELEMENT( "PostUPSElement2D8N", mPostUPSElement2D8N )
        KRATOS_REGISTER_ELEMENT( "PostUPSElement2D9N", mPostUPSElement2D9N )
        KRATOS_REGISTER_ELEMENT( "PostUPSFaceElement3D3N", mPostUPSFaceElement3D3N )
        KRATOS_REGISTER_ELEMENT( "PostUPSFaceElement3D4N", mPostUPSFaceElement3D4N )
        KRATOS_REGISTER_ELEMENT( "PostUPSFaceElement3D6N", mPostUPSFaceElement3D6N )
        KRATOS_REGISTER_ELEMENT( "PostUPSFaceElement3D8N", mPostUPSFaceElement3D8N )
        KRATOS_REGISTER_ELEMENT( "PostUPSFaceElement3D9N", mPostUPSFaceElement3D9N )
        KRATOS_REGISTER_ELEMENT( "PostUPSElement3D4N", mPostUPSElement3D4N )
        KRATOS_REGISTER_ELEMENT( "PostUPSElement3D10N", mPostUPSElement3D10N )
        KRATOS_REGISTER_ELEMENT( "PostUPSElement3D8N", mPostUPSElement3D8N )
        KRATOS_REGISTER_ELEMENT( "PostUPSElement3D20N", mPostUPSElement3D20N )
        KRATOS_REGISTER_ELEMENT( "PostUPSElement3D27N", mPostUPSElement3D27N )
        KRATOS_REGISTER_ELEMENT( "PostUPSElement3D6N", mPostUPSElement3D6N )
        KRATOS_REGISTER_ELEMENT( "PostUPSElement3D15N", mPostUPSElement3D15N )

        KRATOS_REGISTER_CONDITION( "PostSurfaceCondition3D3N", mPostSurfaceCondition3D3N )
        KRATOS_REGISTER_CONDITION( "PostSurfaceCondition3D4N", mPostSurfaceCondition3D4N )
        KRATOS_REGISTER_CONDITION( "PostSurfaceCondition3D6N", mPostSurfaceCondition3D6N )
        KRATOS_REGISTER_CONDITION( "PostSurfaceCondition3D8N", mPostSurfaceCondition3D8N )
        KRATOS_REGISTER_CONDITION( "PostSurfaceCondition3D9N", mPostSurfaceCondition3D9N )
    }

} // namespace Kratos
