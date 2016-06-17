//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Oct 25, 2014 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_NODAL_MONITORING_NodalMonitoringProcess_H_INCLUDED )
#define  KRATOS_NODAL_MONITORING_NodalMonitoringProcess_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "NodalMonitoringProcesses/NodalMonitoringProcess.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for all NodalMonitoringProcesses in Kratos.
/** The NodalMonitoringProcess is the base class for all NodalMonitoringProcesses and defines a simple interface for them.
    Execute method is used to execute the NodalMonitoringProcess algorithms. While the parameters of this method
  can be very different from one NodalMonitoringProcess to other there is no way to create enough overridden
  versions of it. For this reason this method takes no argument and all NodalMonitoringProcess parameters must
  be passed at construction time. The reason is that each constructor can take different set of
  argument without any dependency to other NodalMonitoringProcesses or the base NodalMonitoringProcess class.
*/
class NodalMonitoringProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    typedef Process BaseType;

    /// Pointer definition of NodalMonitoringProcess
    KRATOS_CLASS_POINTER_DEFINITION(NodalMonitoringProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NodalMonitoringProcess() : BaseType()() {}
    NodalMonitoringProcess(Flags options) : BaseType( options ) {}

    /// Destructor.
    virtual ~NodalMonitoringProcess() {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the NodalMonitoringProcess as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Execute method is used to execute the NodalMonitoringProcess algorithms.
    virtual void Execute() {}

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    virtual void ExecuteInitialize()
    {
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    virtual void ExecuteBeforeSolutionLoop()
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    virtual void ExecuteInitializeSolutionStep()
    {
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    virtual void ExecuteFinalizeSolutionStep()
    {
    }


    /// this function will be executed at every time step BEFORE  writing the output
    virtual void ExecuteBeforeOutputStep()
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    virtual void ExecuteAfterOutputStep()
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
    {
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "NodalMonitoringProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NodalMonitoringProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{




    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    NodalMonitoringProcess& operator=(NodalMonitoringProcess const& rOther);

    /// Copy constructor.
    //NodalMonitoringProcess(NodalMonitoringProcess const& rOther);


    ///@}

}; // Class NodalMonitoringProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  NodalMonitoringProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const NodalMonitoringProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_NodalMonitoringProcess_H_INCLUDED  defined 


