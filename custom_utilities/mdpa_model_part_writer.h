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
    typedef ModelPart::ElementsContainerType ElementsContainerType;
    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MDPAModelPartWriter);


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MDPAModelPartWriter(const ModelPart& r_model_part)
    : mr_model_part(r_model_part), mnode_id_offset(0)
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

    void SetNodeIndexOffset(const std::size_t& Number)
    {
        mnode_id_offset = Number;
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
    std::string Info() const final
    {
        std::stringstream buffer;
        buffer << "MDPAModelPartWriter";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
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

    void MDPA_Header(std::ostream& rOStream) const final
    {
        std::time_t curTime = std::time(NULL);
        std::tm* timePtr = localtime(&curTime);
        rOStream << "//KRATOS analysis data file\n";
        rOStream << "//(c) " << (timePtr->tm_year + 1900) << " Hoang-Giang Bui, Ruhr-University Bochum\n";
        rOStream << "//This file is created at " << timePtr->tm_mday << "/" << (timePtr->tm_mon + 1) << "/" << (timePtr->tm_year + 1900) % 100;
        rOStream << " " << timePtr->tm_hour << ":" << timePtr->tm_min << ":" << timePtr->tm_sec << "\n\n";
    }

    void MDPA_Data(std::ostream& rOStream) const final
    {
        Begin(rOStream, "ModelPartData");
        End(rOStream, "ModelPartData");
    }

    void MDPA_Properties(std::ostream& rOStream) const final
    {
        std::cout << " Writing " << mr_model_part.NumberOfProperties() << " properties" << std::endl;
        for(typename ModelPart::PropertiesConstantIterator it = mr_model_part.PropertiesBegin();
            it != mr_model_part.PropertiesEnd(); ++it)
        {
            Begin(rOStream, "Properties", it->Id());
            End(rOStream, "Properties");
        }
    }

    void MDPA_Nodes(std::ostream& rOStream) const final
    {
        std::cout << " Writing " << mr_model_part.NumberOfNodes() << " nodes" << std::endl;
        Begin(rOStream, "Nodes");
        for(typename ModelPart::NodeConstantIterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); ++it)
        {
            rOStream << it->Id() + mnode_id_offset
                     << " " << it->X0()
                     << " " << it->Y0()
                     << " " << it->Z0()
                     << std::endl;
        }
        End(rOStream, "Nodes");
    }

    void MDPA_Elements(std::ostream& rOStream) const final
    {
        // extract the list of element types in the model part and put into the container
        std::map<std::string, std::vector<std::size_t> > ElementList;

        for (typename ModelPart::ElementConstantIterator it = mr_model_part.ElementsBegin(); it != mr_model_part.ElementsEnd(); ++it)
        {
            std::stringstream ElementName;
            ElementName << it->Info() << it->GetGeometry().WorkingSpaceDimension() << "D"
                        << it->GetGeometry().size() << "N";

            ElementList[ElementName.str()].push_back(it->Id());
        }

        // write the connectivities
        for (auto it = ElementList.begin(); it != ElementList.end(); ++it)
        {
            this->Begin(rOStream, "Elements", it->first);
            std::cout << " Writing element " << it->first << std::endl;

            for (std::size_t i = 0; i < it->second.size(); ++i)
            {
                const Element& rElement = mr_model_part.GetElement(it->second[i]);
                rOStream << rElement.Id() << " " << rElement.GetProperties().Id();
                for (std::size_t j = 0; j < rElement.GetGeometry().size(); ++j)
                    rOStream << " " << rElement.GetGeometry()[j].Id() + mnode_id_offset;
                rOStream << std::endl;
            }

            this->End(rOStream, "Elements");
        }
    }

    void MDPA_Conditions(std::ostream& rOStream) const final
    {
        // extract the list of Condition types in the model part and put into the container
        std::map<std::string, std::vector<std::size_t> > ConditionList;

        for (typename ModelPart::ConditionConstantIterator it = mr_model_part.ConditionsBegin(); it != mr_model_part.ConditionsEnd(); ++it)
        {
            std::stringstream ConditionName;
            ConditionName << it->Info() << it->GetGeometry().WorkingSpaceDimension() << "D"
                        << it->GetGeometry().size() << "N";

            ConditionList[ConditionName.str()].push_back(it->Id());
        }

        // write the connectivities
        for (auto it = ConditionList.begin(); it != ConditionList.end(); ++it)
        {
            this->Begin(rOStream, "Conditions", it->first);
            std::cout << " Writing condition " << it->first << std::endl;

            for (std::size_t i = 0; i < it->second.size(); ++i)
            {
                const Condition& rCondition = mr_model_part.GetCondition(it->second[i]);
                rOStream << rCondition.Id() << " " << rCondition.GetProperties().Id();
                for (std::size_t j = 0; j < rCondition.GetGeometry().size(); ++j)
                    rOStream << " " << rCondition.GetGeometry()[j].Id() + mnode_id_offset;
                rOStream << std::endl;
            }

            this->End(rOStream, "Conditions");
        }
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

    const ModelPart& mr_model_part;
    std::size_t mnode_id_offset;

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

    // /// Assignment operator.
    // MDPAModelPartWriter& operator=(MDPAModelPartWriter const& rOther)
    // {
    //     this->mr_model_part = rOther.mr_model_part;
    //     return *this;
    // }

    /// Copy constructor.
    MDPAModelPartWriter(MDPAModelPartWriter const& rOther)
    : mr_model_part(rOther.mr_model_part), mnode_id_offset(rOther.mnode_id_offset)
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

