/*
see layer_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: June 4, 2021 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_LAYER_APPLICATION_ADD_IO_TO_PYTHON_H_INCLUDED )
#define  KRATOS_LAYER_APPLICATION_ADD_IO_TO_PYTHON_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/define_python.h"

namespace Kratos
{

namespace Python
{

void  LayerApplication_AddIOToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_LAYER_APPLICATION_ADD_IO_TO_PYTHON_H_INCLUDED  defined
