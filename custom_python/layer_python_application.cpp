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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "layer_application.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_python/add_io_to_python.h"
#include "custom_python/add_parameter_list_to_python.h"

namespace Kratos
{

namespace Python
{

    using namespace boost::python;
    
    BOOST_PYTHON_MODULE(KratosLayerApplication)
    {
        class_<KratosLayerApplication, KratosLayerApplication::Pointer,
               bases<KratosApplication>, boost::noncopyable>
               ("KratosLayerApplication");

        AddParameterListToPython();
        LayerApp_AddCustomUtilitiesToPython();
        LayerApplication_AddIOToPython();
        
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAYER_ENTITY_TYPE )
    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON

