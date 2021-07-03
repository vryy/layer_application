//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 7 May 2021 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_MMG_MESHER_H_INCLUDED )
#define  KRATOS_LAYER_APP_MMG_MESHER_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/mmg_adapter_2d.h"
#include "custom_utilities/mmg_adapter_3d.h"


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

template<int TDim>
struct MMG;

template<>
struct MMG<2> : MMG2DAdapter {};

template<>
struct MMG<3> : MMG3DAdapter {};

///
/**
 * Simple utility to refine the mesh using MMG
 * The node shall be preset with NODAL_MMG_SCALAR_METRIC, NODAL_MMG_VECTOR_METRIC and NODAL_MMG_TENSOR_METRIC
 */
template<int TDim>
class MMGMesher
{
public:

    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MMGMesher);


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MMGMesher();

    /// Specialized Constructor.
    MMGMesher(const bool& use_level_set, const bool& use_metric);

    /// Destructor.
    virtual ~MMGMesher();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Set the interger value for controlling the remeshing
    void SetValue(const Variable<int>& rVariable, const int& Value);

    /// Set the interger value for controlling the remeshing
    void SetValue(const int& param, const int& Value);

    /// Set the double value for controlling the remeshing
    void SetValue(const Variable<double>& rVariable, const double& Value);

    /// Set the interger value for controlling the remeshing
    void SetValue(const int& param, const double& Value);

    /// Initialize the mesh from the model part
    /// The model_part shall contain only one type of element (triangle or tetrahedra)
    /// The nodal index shall be consecutive and starts from 1
    void Initialize(ModelPart& r_model_part);

    /// Export the mesh to a model part
    void Export(ModelPart& r_model_part, const std::size_t& lastNodeId, const std::size_t& lastElementId,
        const std::string& SampleElementName, const std::size_t& lastPropertiesIndex) const;

    /// Save the mesh
    void SaveMesh(const std::string& Name) const;

    // Save the level set
    void SaveLevelSet(const std::string& Name) const;

    // Save the metric
    void SaveMetric(const std::string& Name) const;

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
        buffer << "MMG" << TDim << "DMesher";
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

    typename MMG<TDim>::mesh_t      mmmgMesh;
    typename MMG<TDim>::met_t       mmmgMet;
    typename MMG<TDim>::ls_t        mmmgLs;

    bool muse_level_set;
    bool muse_metric;

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
    MMGMesher& operator=(MMGMesher const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    MMGMesher(MMGMesher const& rOther)
    {
    }

    ///@}

}; // Class MMGMesher

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<int TDim>
inline std::istream& operator >>(std::istream& rIStream, MMGMesher<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const MMGMesher<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_MMG_MESHER_H_INCLUDED
