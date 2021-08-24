//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2 Nov 2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_MODEL_STATE_H_INCLUDED )
#define  KRATOS_MODEL_STATE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "processes/process.h"
#include "utilities/math_utils.h"


namespace Kratos
{

class DataState
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( DataState );
    DataState() {}
    ~DataState() {}
};

template<typename TVariableType>
class NodalDataState : public DataState
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( NodalDataState );
    typedef DataState BaseType;
    typedef typename TVariableType::Type DataType;

    NodalDataState(const TVariableType& rVariable)
    : BaseType(), mrVariable(rVariable)
    {}

    DataType& operator[](const std::size_t& key) { return mData[key]; }
    const DataType& operator[](const std::size_t& key) const { return mData[key]; }

    const TVariableType& mrVariable;
    std::map<std::size_t, DataType> mData;
};

template<typename TVariableType>
class ElementalDataState : public DataState
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( ElementalDataState );
    typedef DataState BaseType;
    typedef typename TVariableType::Type DataType;

    ElementalDataState(const TVariableType& rVariable)
    : BaseType(), mrVariable(rVariable)
    {}

    std::vector<DataType>& operator[](const std::size_t& key) { return mData[key]; }
    const std::vector<DataType>& operator[](const std::size_t& key) const { return mData[key]; }

    const TVariableType& mrVariable;
    std::map<std::size_t, std::vector<DataType> > mData;
};

/**
 * Class to support storing the state of the model (nodal displacements, stresses at Gauss points, ...)
 */
class ModelState
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( ModelState );
    typedef std::map<std::size_t, DataState::Pointer> DataContainerType;

    ModelState() {}

    virtual ~ModelState()
    {
    }

    /// Helper function to create a new instance
    static ModelState::Pointer Create()
    {
        return ModelState::Pointer(new ModelState());
    }

    /// Get the nodal state
    template<typename TVariableType>
    typename NodalDataState<TVariableType>::Pointer GetNodalState(const TVariableType& rVariable)
    {
        typedef NodalDataState<TVariableType> NodalDataType;
        typename NodalDataType::Pointer pNodalData;
        DataContainerType::iterator it_data = mNodalData.find(rVariable.Key());

        if (it_data != mNodalData.end())
        {
            #ifdef SD_APP_FORWARD_COMPATIBILITY
            pNodalData = std::static_pointer_cast<NodalDataType>(it_data->second);
            #else
            pNodalData = boost::static_pointer_cast<NodalDataType>(it_data->second);
            #endif
        }
        else
        {
            std::stringstream ss;
            ss << "Nodal state " << rVariable.Name() << " does not exist in the database";
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }

        return pNodalData;
    }

    /// Get the elemental state
    template<typename TVariableType>
    typename ElementalDataState<TVariableType>::Pointer GetElementalState(const TVariableType& rVariable)
    {
        typedef typename TVariableType::Type DataType;
        typedef ElementalDataState<TVariableType> ElementalDataType;
        typename ElementalDataType::Pointer pElementalData;
        DataContainerType::iterator it_data = mElementalData.find(rVariable.Key());

        if (it_data != mElementalData.end())
        {
            #ifdef SD_APP_FORWARD_COMPATIBILITY
            pElementalData = std::static_pointer_cast<ElementalDataType>(it_data->second);
            #else
            pElementalData = boost::static_pointer_cast<ElementalDataType>(it_data->second);
            #endif
        }
        else
        {
            std::stringstream ss;
            ss << "Elemental state " << rVariable.Name() << " does not exist in the database";
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }

        return pElementalData;
    }

    /// Save the nodal state
    template<typename TVariableType>
    void SaveNodalState(ModelPart& r_model_part, const TVariableType& rVariable, const ProcessInfo& rCurrentProcessInfo)
    {
        this->StoreNodalState(r_model_part, rVariable, rVariable, rCurrentProcessInfo);
    }

    /// Store the nodal state from a source variable to a target variable
    template<typename TVariableType>
    inline void StoreNodalState(ModelPart& r_model_part, const TVariableType& rSourceVariable, const TVariableType& rTargetVariable, const ProcessInfo& rCurrentProcessInfo)
    {
        typedef NodalDataState<TVariableType> NodalDataType;
        typename NodalDataType::Pointer pNodalData;
        pNodalData = typename NodalDataType::Pointer(new NodalDataType(rTargetVariable));

        for(ModelPart::NodeIterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd() ; ++it)
        {
            (*pNodalData)[it->Id()] = it->GetSolutionStepValue(rSourceVariable);
        }

        mNodalData[rTargetVariable.Key()] = pNodalData;

        if (rSourceVariable == rTargetVariable)
            std::cout << "Nodal state " << rSourceVariable.Name() << " is saved" << std::endl;
        else
            std::cout << "Nodal state " << rSourceVariable.Name() << " is stored to " << rTargetVariable.Name() << std::endl;
    }

    /// Save the element state
    template<typename TVariableType>
    void SaveElementalState(ModelPart& r_model_part, const TVariableType& rVariable, const ProcessInfo& rCurrentProcessInfo)
    {
        typedef typename TVariableType::Type DataType;
        typedef ElementalDataState<TVariableType> ElementalDataType;
        typename ElementalDataType::Pointer pElementalData;
        pElementalData = typename ElementalDataType::Pointer(new ElementalDataType(rVariable));

        std::vector<ConstitutiveLaw::Pointer> pConstitutiveLaws;
        std::size_t npoints;

        for (ModelPart::ElementsContainerType::ptr_iterator it = r_model_part.Elements().ptr_begin();
            it != r_model_part.Elements().ptr_end(); ++it)
        {
            (*it)->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, pConstitutiveLaws, rCurrentProcessInfo);
            npoints = pConstitutiveLaws.size();
            (*pElementalData)[(*it)->Id()].resize(npoints);

            for (std::size_t i = 0; i < npoints; ++i)
            {
                std::vector<DataType>& elemental_data = (*pElementalData)[(*it)->Id()];
                pConstitutiveLaws[i]->GetValue(rVariable, elemental_data[i]);
            }
        }

        mElementalData[rVariable.Key()] = pElementalData;

        std::cout << "Elemental state " << rVariable.Name() << " is saved" << std::endl;
    }

    /// Restore the nodal state to the model_part
    template<typename TVariableType>
    void LoadNodalState(ModelPart& r_model_part, const TVariableType& rVariable, const ProcessInfo& rCurrentProcessInfo)
    {
        typedef NodalDataState<TVariableType> NodalDataType;
        typename NodalDataType::Pointer pNodalData = this->GetNodalState<TVariableType>(rVariable);

        for(ModelPart::NodeIterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd() ; ++it)
        {
            it->GetSolutionStepValue(rVariable) = (*pNodalData)[it->Id()];
        }

        std::cout << "Nodal state " << rVariable.Name() << " is loaded" << std::endl;
    }

    /// Restore the elemental state to the model_part
    template<typename TVariableType>
    void LoadElementalState(ModelPart& r_model_part, const TVariableType& rVariable, const ProcessInfo& rCurrentProcessInfo)
    {
        typedef typename TVariableType::Type DataType;
        typedef ElementalDataState<TVariableType> ElementalDataType;
        typename ElementalDataType::Pointer pElementalData = this->GetElementalState<TVariableType>(rVariable);

        std::vector<ConstitutiveLaw::Pointer> pConstitutiveLaws;
        std::size_t npoints;

        for (ModelPart::ElementsContainerType::ptr_iterator it = r_model_part.Elements().ptr_begin();
            it != r_model_part.Elements().ptr_end(); ++it)
        {
            (*it)->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, pConstitutiveLaws, rCurrentProcessInfo);
            npoints = pConstitutiveLaws.size();

            for (std::size_t i = 0; i < npoints; ++i)
            {
                const DataType& value = (*pElementalData)[(*it)->Id()][i];
                pConstitutiveLaws[i]->SetValue(rVariable, value, rCurrentProcessInfo);
            }
        }

        std::cout << "Elemental state " << rVariable.Name() << " is loaded" << std::endl;
    }

private:

    DataContainerType mNodalData;
    DataContainerType mElementalData;

};

}

#endif // KRATOS_MODEL_STATE_H_INCLUDED

