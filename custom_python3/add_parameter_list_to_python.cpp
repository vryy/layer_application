/*
see layer_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: June 4, 2021 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "add_parameter_list_to_python.h"
#include "custom_utilities/parameter_list.h"


namespace Kratos
{
namespace Python
{

typedef Kratos::ParameterList<std::string> ParameterListType;

void SetDoubleValue(ParameterListType& dummy, const std::string &name, double value)
{
    dummy.set(name, value);
}

void SetIntValue(ParameterListType& dummy, const std::string &name, int value)
{
    dummy.set(name, value);
}

// void SetUnsignedIntValue(ParameterListType& dummy, const std::string &name, unsigned int value)
// {
//    dummy.set(name, value);
// }

//void SetCharValue(ParameterListType& dummy, const std::string &name, const char value[])
//{
//    dummy.set(name, value);
//}

void SetStringValue(ParameterListType& dummy, const std::string &name, const std::string value)
{
    dummy.set(name, value);
}

void SetBoolValue(ParameterListType& dummy, const std::string &name, bool value)
{
    dummy.set(name, value);
}

void SetVectorValue(ParameterListType& dummy, const std::string &name, Vector value)
{
    dummy.set(name, value);
}

void SetMatrixValue(ParameterListType& dummy, const std::string &name, Matrix value)
{
    dummy.set(name, value);
}

ParameterListType& SubList(ParameterListType& dummy, const std::string &name)
{
    return dummy.sublist(name);
}

void LayerApp_AddParameterListToPython(pybind11::module& m)
{
    class_<ParameterListType, ParameterListType::Pointer>
    (m, "ParameterList")
    .def(init<>())
    .def("set", SetDoubleValue)
    .def("set", SetIntValue)
//    .def("set", SetUnsignedIntValue)
    .def("set", SetStringValue)
//    .def("set", SetBoolValue)
    .def("set", SetVectorValue)
    .def("set", SetMatrixValue)
    .def("sublist", SubList, pybind11::return_value_policy::reference_internal)
    ;
}

}  // namespace Python.

} // Namespace Kratos

