//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Nov 2014 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_MDPA_MODEL_PART_WRITER_H_INCLUDED )
#define  KRATOS_LAYER_APP_MDPA_MODEL_PART_WRITER_H_INCLUDED

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
#include "includes/model_part.h"


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
class MDPAModelPartWriter : public MDPAWriter
{
public:

    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MDPAModelPartWriter);
    
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MDPAModelPartWriter(ModelPart::Pointer p_model_part) : mp_model_part(p_model_part)
    {
    }

    /// Destructor.
    virtual ~MDPAModelPartWriter()
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

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MDPAModelPartWriter";
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
        rOStream << "//KRATOS analysis data file\n";
        rOStream << "//(c) " << (timePtr->tm_year + 1900) << " Hoang Giang Bui, Ruhr-University Bochum\n";
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
        for(typename ModelPart::PropertiesContainerType::ContainerType::iterator it = mp_model_part->PropertiesArray().begin();
            it != mp_model_part->PropertiesArray().end(); ++it)
        {
            Begin(rOStream, "Properties", (*it)->Id());
            End(rOStream, "Properties");
        }
    }

    virtual void MDPA_Nodes(std::ostream& rOStream)
    {
        Begin(rOStream, "Nodes");
        for(ModelPart::NodeIterator it = mp_model_part->NodesBegin(); it != mp_model_part->NodesEnd(); ++it)
        {
            rOStream << it->Id()
                     << " " << it->X0()
                     << " " << it->Y0()
                     << " " << it->Z0()
                     << std::endl;
        }
        End(rOStream, "Nodes");
    }
    
    virtual void MDPA_Elements(std::ostream& rOStream)
    {
        //TODO
        // extract the list of element types in the model part and put into the container
//        std::set<std::string> ElementNames;
//        
//        ElementsArrayType& pElements = r_model_part.Elements();
//        for (typename ModelPart::ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
//        {
//            
//        }
    }

    virtual void MDPA_Conditions(std::ostream& rOStream)
    {
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
    ModelPart::Pointer mp_model_part;
    
    
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
    MDPAModelPartWriter& operator=(MDPAModelPartWriter const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    MDPAModelPartWriter(MDPAModelPartWriter const& rOther)
    {
    }

    ///@}

}; // Class MDPAModelPartWriter

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, MDPAModelPartWriter& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const MDPAModelPartWriter& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_MDPA_MODEL_PART_WRITER_H_INCLUDED

