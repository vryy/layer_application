
/* *********************************************************
 *
 *   Last modified by:    $Author: hbui $
 *   Date:                $Date: Nov 2, 2014 $
 *   Revision:            $Revision: 1.2 $
 *
 * ***********************************************************/


#if !defined(KRATOS_LAYER_APP_PARAMETER_LIST_H_INCLUDED )
#define  KRATOS_LAYER_APP_PARAMETER_LIST_H_INCLUDED

// External includes
#include <boost/variant.hpp>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#ifndef SD_APP_FORWARD_COMPATIBILITY
#include "containers/vector_map.h"
#else
#include "custom_utilities/containers/vector_map.h"
#endif


namespace Kratos
{

typedef boost::variant<
        bool,
        int,
        unsigned int,
        long,
        long long,
        std::size_t,
        float,
        double,
        std::string,
        Vector,
        Matrix
        > KratosParameterListAcceptedType;

/**
A ParameterList class offers an alternative to Teuchos::ParameterList
*/
template<class TKeyType>
class ParameterList
{

    class parameter_list_visitor : public boost::static_visitor<int>
    {
    public:
        int operator()(bool b) const { return 0;}
        int operator()(int i) const { return 1;}
        int operator()(unsigned int ui) const { return 2;}
        int operator()(long l) const { return 3;}
        int operator()(long long ll) const { return 4;}
        int operator()(float f) const { return 5;}
        int operator()(double d) const { return 6;}
        int operator()(std::string s) const { return 7;}
        int operator()(Vector& v) const { return 8;}
        int operator()(Matrix& m) const { return 9;}
    };

public:

    KRATOS_CLASS_POINTER_DEFINITION(ParameterList);

    typedef VectorMap<TKeyType, KratosParameterListAcceptedType> DataContainerType;

    typedef typename DataContainerType::iterator iterator;

    typedef typename DataContainerType::const_iterator const_iterator;

    typedef typename DataContainerType::pair_iterator pair_iterator;

    typedef typename DataContainerType::pair_const_iterator pair_const_iterator;

    typedef typename DataContainerType::key_type KeyType; //TKeyType

    typedef typename DataContainerType::data_type DataType; //KratosParameterListAcceptedType

    typedef VectorMap<KeyType, ParameterList<TKeyType> > SubListType;

    inline ParameterList() : mData(), mSubList()
    {}

    inline ParameterList(const DataContainerType& rOther) : mData(rOther), mSubList()
    {}

    inline ParameterList(const ParameterList& rOther) : mData(rOther.mData), mSubList(rOther.mSubList)
    {}

    inline ~ParameterList()
    {}

    inline ParameterList& operator=(const ParameterList& rOther)
    {
        mData = rOther.mData;
        mSubList = rOther.mSubList;
        return *this;
    }

    inline DataType& operator[](const TKeyType& rKey)
    {
        return mData[rKey];
    }

    template<class TDataType>
    inline void set(const TKeyType& rKey, const TDataType& rData)
    {
        mData[rKey] = rData;
    }

    template<class TDataType>
    inline TDataType& get(const TKeyType& rKey)
    {
        return boost::get<TDataType&>(mData[rKey]);
    }

    template<class TDataType>
    inline TDataType& get(const TKeyType& rKey, const TDataType& rInitialData)
    {
        iterator i = mData.find(rKey);

        if(i == mData.end())
            set(rKey, rInitialData);

        return boost::get<TDataType&>(mData[rKey]);
    }

    inline std::string& get(const TKeyType& rKey, const char* rInitialData)
    {
        iterator i = mData.find(rKey);

        if(i == mData.end())
            set(rKey, std::string(rInitialData));

        return boost::get<std::string&>(mData[rKey]);
    }

    inline int type(const TKeyType& rKey)
    {
        return boost::apply_visitor(parameter_list_visitor(), mData[rKey]);
    }

    inline ParameterList& sublist(const TKeyType& rKey)
    {
        iterator i = mData.find(rKey);
        if(i != mData.end())
        {
            KRATOS_THROW_ERROR(std::logic_error, "Existing key has been associated with value", "");
        }

        typename SubListType::iterator ii = mSubList.find(rKey);
        if(ii == mSubList.end())
        {
            ParameterList pl;
            mSubList[rKey] = pl;
        }

        return mSubList[rKey];
    }

    iterator begin() { return mData.begin(); }
    const_iterator begin() const { return mData.begin(); }
    iterator end() { return mData.end(); }
    const_iterator end() const { return mData.end(); }

    pair_iterator pair_begin() { return mData.pair_begin(); }
    pair_const_iterator pair_begin() const { return mData.pair_begin(); }
    pair_iterator pair_end() { return mData.pair_end(); }
    pair_const_iterator pair_end() const { return mData.pair_end(); }

    iterator find(const KeyType& rKey) { return mData.find(rKey); }
    const_iterator find(const KeyType& rKey) const { return mData.find(rKey); }

    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "parameter list (size = " << mData.size() + mSubList.size() << ") : ";

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
        mData.PrintData(rOStream);
        for(typename SubListType::pair_const_iterator i = mSubList.pair_begin() ; i != mSubList.pair_end() ; ++i)
        {
            rOStream << "(" << (i->first) << " , " << (i->second) << ")" << std::endl;
        }
    }

    void Print(std::ostream& rOStream)
    {
        for(pair_iterator i = mData.pair_begin(); i != mData.pair_end(); ++i)
        {
            rOStream << "(" << (i->first) << ", " << (i->second) << ", type = " << type(i->first) << ")" << std::endl;
        }
        for(typename SubListType::pair_const_iterator i = mSubList.pair_begin() ; i != mSubList.pair_end() ; ++i)
            rOStream << "(" << (i->first) << " , " << (i->second) << ")" << std::endl;
    }

private:

    DataContainerType mData;
    SubListType mSubList;

    friend class Serializer;

    // TODO enable serializer when serialization is added to vector map
    virtual void save(Serializer& rSerializer) const
    {
        // rSerializer.save("Data", mData);
        // rSerializer.save("SubList", mSubList);
    }

    virtual void load(Serializer& rSerializer)
    {
        // rSerializer.load("Data", mData);
        // rSerializer.load("SubList", mSubList);
    }

};

template<class TKeyType>
inline std::istream& operator >> (std::istream& is, ParameterList<TKeyType>& rThis)
{
    return is;
}

template<class TKeyType>
inline std::ostream& operator << (std::ostream& os, const ParameterList<TKeyType>& rThis)
{
    rThis.PrintInfo(os);
    os << std::endl;
    rThis.PrintData(os);
    return os;
}

} // end namespace Kratos

#endif /* KRATOS_LAYER_APP_PARAMETER_LIST_H_INCLUDED */
