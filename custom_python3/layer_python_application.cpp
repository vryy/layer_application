//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: June 4, 2021 $
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
// #include "custom_python/add_utilities_to_python.h"
// #include "custom_python/add_io_to_python.h"
// #include "custom_python/add_parameter_list_to_python.h"

namespace Kratos
{

namespace Python
{

    PYBIND11_MODULE(KratosLayerApplication, m)
    {
        namespace py = pybind11;

        py::class_<KratosLayerApplication,
                KratosLayerApplication::Pointer,
                KratosApplication >(m, "KratosLayerApplication")
                .def(py::init<>())
                ;

        // AddParameterListToPython();
        // LayerApp_AddCustomUtilitiesToPython();
        // LayerApplication_AddIOToPython();

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, LAYER_ENTITY_TYPE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, LAYER_NAME )
        #ifdef LAYER_APP_USE_MMG
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NODAL_MMG_SCALAR_METRIC )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NODAL_MMG_VECTOR_METRIC )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NODAL_MMG_TENSOR_METRIC )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NODAL_MMG_LEVEL_SET )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, MMG_GRADATION )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, MMG_HAUSDORFF_DISTANCE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, MMG_MINIMAL_MESH_SIZE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, MMG_MAXIMAL_MESH_SIZE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, MMG_CONSTANT_MESH_SIZE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, MMG_RMC_VOLUME_FRACTION )
        #endif
    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON

