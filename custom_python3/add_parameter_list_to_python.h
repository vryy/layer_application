/*
see layer_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: June 4, 2021 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_ADD_PARAMETER_LIST_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_PARAMETER_LIST_TO_PYTHON_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/define_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

void LayerApp_AddParameterListToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_PARAMETER_LIST_TO_PYTHON_H_INCLUDED  defined
