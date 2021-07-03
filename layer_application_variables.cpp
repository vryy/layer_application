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
#include "layer_application_variables.h"


namespace Kratos
{
    KRATOS_CREATE_VARIABLE(int, LAYER_ENTITY_TYPE)
    KRATOS_CREATE_VARIABLE(int, LAYER_PROP_ID)
    KRATOS_CREATE_VARIABLE(std::string, LAYER_ENTITY_NAME)
    KRATOS_CREATE_VARIABLE(std::string, LAYER_NAME)
    #ifdef LAYER_APP_USE_MMG
    KRATOS_CREATE_VARIABLE(double, NODAL_MMG_SCALAR_METRIC)
    KRATOS_CREATE_VARIABLE(Vector, NODAL_MMG_VECTOR_METRIC)
    KRATOS_CREATE_VARIABLE(Matrix, NODAL_MMG_TENSOR_METRIC)
    KRATOS_CREATE_VARIABLE(double, NODAL_MMG_LEVEL_SET)
    KRATOS_CREATE_VARIABLE(double, MMG_GRADATION)
    KRATOS_CREATE_VARIABLE(double, MMG_HAUSDORFF_DISTANCE)
    KRATOS_CREATE_VARIABLE(double, MMG_MINIMAL_MESH_SIZE)
    KRATOS_CREATE_VARIABLE(double, MMG_MAXIMAL_MESH_SIZE)
    KRATOS_CREATE_VARIABLE(double, MMG_CONSTANT_MESH_SIZE)
    KRATOS_CREATE_VARIABLE(double, MMG_RMC_VOLUME_FRACTION)
    #endif

}
