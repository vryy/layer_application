//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: June 4, 2021 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_LAYER_APP_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_LAYER_APP_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/define_python.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

void  LayerApp_AddCustomUtilitiesToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_LAYER_APP_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED  defined

