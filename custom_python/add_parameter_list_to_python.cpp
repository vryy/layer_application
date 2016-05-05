/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Jan 9, 2013 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes
#include <string>

// External includes
#include <boost/python.hpp>

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

void AddParameterListToPython()
{
    using namespace boost::python;

    class_<ParameterListType, ParameterListType::Pointer, boost::noncopyable>("ParameterList", init<>())
    .def("set", SetDoubleValue)
    .def("set", SetIntValue)
//    .def("set", SetUnsignedIntValue)
    .def("set", SetStringValue)
//    .def("set", SetBoolValue)
    .def("set", SetVectorValue)
    .def("set", SetMatrixValue)
    .def("sublist", SubList, return_internal_reference<>())
    ;
}

}  // namespace Python.

} // Namespace Kratos

