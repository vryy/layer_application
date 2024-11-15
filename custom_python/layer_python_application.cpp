//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Oct 30, 2014 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes


// External includes
#if defined(KRATOS_PYTHON)


// Project includes
#include "includes/define_python.h"
#include "layer_application.h"
#include "layer_application_variables.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_python/add_io_to_python.h"
#include "custom_python/add_parameter_list_to_python.h"

namespace Kratos
{

namespace Python
{

    BOOST_PYTHON_MODULE(KratosLayerApplication)
    {
        using namespace boost::python;

        class_<KratosLayerApplication, KratosLayerApplication::Pointer,
               bases<KratosApplication>, boost::noncopyable>
               ("KratosLayerApplication");

        AddParameterListToPython();
        LayerApp_AddCustomUtilitiesToPython();
        LayerApplication_AddIOToPython();

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAYER_ENTITY_TYPE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAYER_NAME )
        #ifdef LAYER_APP_USE_MMG
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_MMG_SCALAR_METRIC )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_MMG_VECTOR_METRIC )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_MMG_TENSOR_METRIC )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_MMG_LEVEL_SET )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( MMG_GRADATION )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( MMG_HAUSDORFF_DISTANCE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( MMG_MINIMAL_MESH_SIZE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( MMG_MAXIMAL_MESH_SIZE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( MMG_CONSTANT_MESH_SIZE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( MMG_RMC_VOLUME_FRACTION )
        #endif
    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
