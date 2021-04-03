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


namespace Kratos
{
    KRATOS_CREATE_VARIABLE(int, LAYER_ENTITY_TYPE)
    KRATOS_CREATE_VARIABLE(int, LAYER_PROP_ID)
    KRATOS_CREATE_VARIABLE(std::string, LAYER_ENTITY_NAME)
    KRATOS_CREATE_VARIABLE(std::string, LAYER_NAME)

    KratosLayerApplication::KratosLayerApplication()
    : KratosApplication()
    , mPostElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) )
    , mPostElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) )
    , mPostElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) )
    , mPostElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) )
    , mPostElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9 ) ) ) )
    , mPostFaceElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) )
    , mPostFaceElement3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) )
    , mPostFaceElement3D6N( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) )
    , mPostFaceElement3D8N( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) )
    , mPostFaceElement3D9N( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9 ) ) ) )
    , mPostElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) )
    , mPostElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) )
    , mPostElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) )
    , mPostElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) )
    , mPostElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) )
    , mPostElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) )
    , mPostElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) )
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
    }
} // namespace Kratos

