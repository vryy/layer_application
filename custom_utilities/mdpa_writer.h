//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Nov 2014 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_MDPA_WRITER_H_INCLUDED )
#define  KRATOS_LAYER_APP_MDPA_WRITER_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>
#include <set>
#include <ctime>

// External includes
#include <boost/foreach.hpp>
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"


namespace Kratos
{

///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/*** Detail class definition.
Class defines interface for writing to MDPA data file
 */
class MDPAWriter
{
public:

    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MDPAWriter);
    
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MDPAWriter()
    {
    }

    /// Destructor.
    virtual ~MDPAWriter()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    void WriteMDPA(std::string fn)
    {
        std::ofstream fid;
        std::string new_fn = fn + std::string(".mdpa");
        fid.open(new_fn.c_str());
        fid.precision(15);
//        fid.setf(std::ios_base::scientific);
//        fid.setf(std::ios_base::floatfield);
//        fid.setf(std::ios_base::scientific, std::ios_base::floatfield);

        std::cout << "MDPA_Header starts" << std::endl;
        MDPA_Header(fid);
        std::cout << "MDPA_Header completed" << std::endl;

        std::cout << "MDPA_Data starts" << std::endl;
        MDPA_Data(fid);
        std::cout << "MDPA_Data completed" << std::endl;

        std::cout << "MDPA_Properties starts" << std::endl;
        MDPA_Properties(fid);
        std::cout << "MDPA_Properties completed" << std::endl;

        std::cout << "MDPA_Nodes starts" << std::endl;
        MDPA_Nodes(fid);
        std::cout << "MDPA_Nodes completed" << std::endl;

        std::cout << "MDPA_Elements starts" << std::endl;
        MDPA_Elements(fid);
        std::cout << "MDPA_Elements completed" << std::endl;

        std::cout << "MDPA_Conditions starts" << std::endl;
        MDPA_Conditions(fid);
        std::cout << "MDPA_Conditions completed" << std::endl;

        fid.close();
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MDPAWriter";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    void Print()
    {
        PrintInfo(std::cout);
        std::cout << std::endl;
        PrintData(std::cout);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void MDPA_Header(std::ostream& rOStream)
    {
        std::time_t curTime = std::time(NULL);
        std::tm* timePtr = localtime(&curTime);
        rOStream << "//KRATOS isogeometric application data file\n";
        rOStream << "//(c) " << (timePtr->tm_year + 1900) << " Hoang-Giang Bui, Ruhr-University Bochum\n";
        rOStream << "//This file is created at " << timePtr->tm_mday << "/" << timePtr->tm_mon << "/" << (timePtr->tm_year + 1900) % 100;
        rOStream << " " << timePtr->tm_hour << ":" << timePtr->tm_min << ":" << timePtr->tm_sec << "\n\n";
    }

    virtual void MDPA_Data(std::ostream& rOStream)
    {
        Begin(rOStream, "ModelPartData");
        End(rOStream, "ModelPartData");
    }
    
    virtual void MDPA_Properties(std::ostream& rOStream)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Call the virtual class function of MDPAWriter: ", __FUNCTION__)
    }

    virtual void MDPA_Nodes(std::ostream& rOStream)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Call the virtual class function of MDPAWriter: ", __FUNCTION__)
    }
    
    virtual void MDPA_Elements(std::ostream& rOStream)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Call the virtual class function of MDPAWriter: ", __FUNCTION__)
    }

    virtual void MDPA_Conditions(std::ostream& rOStream)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Call the virtual class function of MDPAWriter: ", __FUNCTION__)
    }
    
    void Begin(std::ostream& rOStream, const std::string& section)
    {
        rOStream << "Begin " << section << std::endl;
    }

    void Begin(std::ostream& rOStream, const std::string& section, int id)
    {
        rOStream << "Begin " << section << " " << id << std::endl;
    }

    void Begin(std::ostream& rOStream, const std::string& section, const std::string& header)
    {
        rOStream << "Begin " << section << " " << header << std::endl;
    }

    void End(std::ostream& rOStream, const std::string& section)
    {
        rOStream << "End " << section << std::endl << std::endl;
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{
    
    ///@}
    ///@name Member Variables
    ///@{
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
        
    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{
    
    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MDPAWriter& operator=(MDPAWriter const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    MDPAWriter(MDPAWriter const& rOther)
    {
    }

    ///@}

}; // Class MDPAWriter

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, MDPAWriter& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const MDPAWriter& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_MDPA_WRITER_H_INCLUDED

