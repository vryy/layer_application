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
#include "layer_application.h"


namespace Kratos
{
    KRATOS_CREATE_VARIABLE(int, LAYER_ENTITY_TYPE)
    KRATOS_CREATE_VARIABLE(int, LAYER_PROP_ID)
    KRATOS_CREATE_VARIABLE(std::string, LAYER_ENTITY_NAME)
    
    void KratosLayerApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosLayerApplication... " << std::endl;
        
        // register variables to Kratos kernel
        KRATOS_REGISTER_VARIABLE(LAYER_ENTITY_TYPE)
        KRATOS_REGISTER_VARIABLE(LAYER_ENTITY_NAME)
        KRATOS_REGISTER_VARIABLE(LAYER_PROP_ID)
    }
} // namespace Kratos

