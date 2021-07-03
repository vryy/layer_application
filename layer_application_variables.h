//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Oct 30, 2014 $
//   Revision:            $Revision: 1.0 $
//
//
//Change log:
//  +   30/10/2014: create layer_application.h

#if !defined(KRATOS_LAYER_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_LAYER_APPLICATION_VARIABLES_H_INCLUDED


// System includes
#include <string>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"


namespace Kratos
{

    ///@name Kratos Globals
    ///@{
    // Variables definition
    KRATOS_DEFINE_VARIABLE(int, LAYER_ENTITY_TYPE)
    KRATOS_DEFINE_VARIABLE(int, LAYER_PROP_ID)
    KRATOS_DEFINE_VARIABLE(std::string, LAYER_ENTITY_NAME)
    KRATOS_DEFINE_VARIABLE(std::string, LAYER_NAME)
    #ifdef LAYER_APP_USE_MMG
    KRATOS_DEFINE_VARIABLE(double, NODAL_MMG_SCALAR_METRIC)
    KRATOS_DEFINE_VARIABLE(Vector, NODAL_MMG_VECTOR_METRIC)
    KRATOS_DEFINE_VARIABLE(Matrix, NODAL_MMG_TENSOR_METRIC)
    KRATOS_DEFINE_VARIABLE(double, NODAL_MMG_LEVEL_SET)
    KRATOS_DEFINE_VARIABLE(double, MMG_GRADATION)
    KRATOS_DEFINE_VARIABLE(double, MMG_HAUSDORFF_DISTANCE)
    KRATOS_DEFINE_VARIABLE(double, MMG_MINIMAL_MESH_SIZE)
    KRATOS_DEFINE_VARIABLE(double, MMG_MAXIMAL_MESH_SIZE)
    KRATOS_DEFINE_VARIABLE(double, MMG_CONSTANT_MESH_SIZE)
    KRATOS_DEFINE_VARIABLE(double, MMG_RMC_VOLUME_FRACTION) // associate with -rmc option
    #endif

}

#endif // KRATOS_LAYER_APPLICATION_VARIABLES_H_INCLUDED