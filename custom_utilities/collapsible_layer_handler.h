//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 4 May 2016 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_COLLAPSIBLE_LAYER_HANDLER_H_INCLUDED )
#define  KRATOS_LAYER_APP_COLLAPSIBLE_LAYER_HANDLER_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>
#include <set>
#include <cmath>

// External includes
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

// Project includes
#include "includes/define.h"
#include "custom_utilities/parameter_list.h"
#include "custom_utilities/group.h"
#include "custom_utilities/spatial_key.h"
#include "custom_utilities/layer_handler.h"


namespace Kratos
{

///@addtogroup LayerApplication
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
    Layer handler keeps a global node information and handle layer collapse properly
 */
class CollapsibleLayerHandler : public LayerHandler
{
public:

    const double PI = std::atan(1.0)*4;

    ///@name Type Definitions
    ///@{

    typedef LayerHandler BaseType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::PointsContainerType PointsContainerType;
    typedef typename BaseType::LayersContainerType LayersContainerType;
    typedef std::map<std::string, Group::Pointer> GroupsContainerType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(CollapsibleLayerHandler);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CollapsibleLayerHandler() : BaseType()
    {
        mDx = 0.1 * PI;
        mDy = 0.1 * PI;
        mDz = 0.1 * PI;
    }

    /// Destructor.
    virtual ~CollapsibleLayerHandler()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetSpacing(double Dx, double Dy, double Dz)
    {
        mDx = Dx * PI; // this is to make sure that we have a non-trivial spacecing
        mDy = Dy * PI;
        mDz = Dz * PI;
    }

    virtual void AddLayer(std::string name, boost::python::dict& pyDictNodalSet, boost::python::dict& pyDictEntitySet, boost::python::dict& pyDictEntityInfoSet)
    {
        BaseType::AddLayer(name, pyDictNodalSet, pyDictEntitySet, pyDictEntityInfoSet);
        mLayerIsCollapsed.insert(std::pair<std::string, bool>(name, false));
    }

    void AddGroup(std::string name, boost::python::list& pyListLayerNames)
    {
        Group::Pointer pGroup = Group::Pointer(new Group(name));

        typedef boost::python::stl_input_iterator<std::string> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& layer_name,
                      std::make_pair(iterator_value_type(pyListLayerNames), // begin
                      iterator_value_type() ) ) // end
        {
            if(mpLayers.find(layer_name) != mpLayers.end())
                // check if layer layer_name has been existed in mpLayers
                pGroup->AddLayer(mpLayers[layer_name]);
            else
            {
                std::stringstream ss;
                ss << "The layer " << layer_name << " does not exist in the layer list";
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }
        }

        mpGroups.insert(GroupsContainerType::value_type(name, pGroup));
        mGroupIsCollapsed.insert(std::pair<std::string, bool>(name, false));
    }

    void CollapseLayer(std::string name, const double TOL)
    {
        if(mpLayers.find(name) == mpLayers.end())
        {
            std::cout << "Layer " << name << " does not exist" << std::endl;
            return;
        }

        Layer& thisLayer = *mpLayers[name];

        // firstly put nodes in layer into bin
        typedef std::map<SpatialKey, Layer::NodesContainerType> BinType;
        BinType Binning;
        for(Layer::NodesIteratorType it = thisLayer.NodesBegin(); it != thisLayer.NodesEnd(); ++it)
        {
            double X = (*it)->X();
            double Y = (*it)->Y();
            double Z = (*it)->Z();

            int kx = X / mDx;
            int ky = Y / mDy;
            int kz = Z / mDz;

            SpatialKey key(kx, ky, kz);
            Binning[key].push_back(*it);
        }

        // for each bin, collapse locally
        this->CollapseBin(Binning, TOL);

        mLayerIsCollapsed[name] = true;
    }

    void CollapseGroup(std::string name, const double TOL)
    {
        if(mpGroups.find(name) == mpGroups.end())
        {
            std::cout << "Group " << name << " does not exist" << std::endl;
            return;
        }

        // check if each layer in group is collapsed
        for(Group::LayerContainerIteratorType it = mpGroups[name]->LayerBegin(); it != mpGroups[name]->LayerEnd(); ++it)
        {
            if(mLayerIsCollapsed[*it] == false)
                CollapseLayer(*it, TOL);
        }

        // firstly put nodes in layer into bin
        typedef std::map<SpatialKey, Layer::NodesContainerType> BinType;
        BinType Binning;
        for(Group::LayerContainerIteratorType itl = mpGroups[name]->LayerBegin(); itl != mpGroups[name]->LayerEnd(); ++itl)
        {
            Layer& thisLayer = *mpLayers[*itl];
            for(Layer::NodesIteratorType it = thisLayer.NodesBegin(); it != thisLayer.NodesEnd(); ++it)
            {
                double X = (*it)->X();
                double Y = (*it)->Y();
                double Z = (*it)->Z();

                int kx = X / mDx;
                int ky = Y / mDy;
                int kz = Z / mDz;

                SpatialKey key(kx, ky, kz);
                Binning[key].push_back(*it);
            }
        }

        // for each bin, collapse locally
        this->CollapseBin(Binning, TOL);

        mGroupIsCollapsed[name] = true;
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
        std::stringstream buffer;
        buffer << "CollapsibleLayerHandler";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        BaseType::PrintData(rOStream);
        for(GroupsContainerType::const_iterator it = mpGroups.begin(); it != mpGroups.end(); ++it)
            std::cout << "+" << *(it->second);
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

    double mDx, mDy, mDz;

    GroupsContainerType mpGroups;
    std::map<std::string, bool> mLayerIsCollapsed;
    std::map<std::string, bool> mGroupIsCollapsed;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CollapseBin(std::map<SpatialKey, Layer::NodesContainerType>& rBin, const double TOL)
    {
        for(std::map<SpatialKey, Layer::NodesContainerType>::iterator it = rBin.begin(); it != rBin.end(); ++it)
        {
            std::size_t cnt1 = 0;
            for(Layer::NodesIteratorType it2 = it->second.begin(); it2 != it->second.end(); ++it2)
            {
                double X1 = (*it2)->X();
                double Y1 = (*it2)->Y();
                double Z1 = (*it2)->Z();

                std::size_t cnt2 = 0;
                for(Layer::NodesIteratorType it3 = it->second.begin(); it3 != it->second.end(); ++it3)
                {
                    double X2 = (*it3)->X();
                    double Y2 = (*it3)->Y();
                    double Z2 = (*it3)->Z();

                    if(cnt2 != cnt1)
                    {
                        double diff = sqrt(pow(X2 - X1, 2) + pow(Y2 - Y1, 2) + pow(Z2 - Z1, 2));
                        if(diff < TOL)
                            (*it3)->SetId((*it2)->Id());
                    }
                    ++cnt2;
                }
                ++cnt1;
            }
        }
    }

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
    CollapsibleLayerHandler& operator=(CollapsibleLayerHandler const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    CollapsibleLayerHandler(CollapsibleLayerHandler const& rOther)
    {
    }

    ///@}

}; // Class CollapsibleLayerHandler

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, CollapsibleLayerHandler& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const CollapsibleLayerHandler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_COLLAPSIBLE_LAYER_HANDLER_H_INCLUDED
