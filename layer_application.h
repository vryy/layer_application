//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Oct 30, 2014 $
//   Revision:            $Revision: 1.0 $
//
//
//Change log:
//  +   30/10/2014: create layer_application.h

#if !defined(KRATOS_LAYER_APPLICATION_H_INCLUDED)
#define KRATOS_LAYER_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"
#include "custom_utilities/spatial_point.h"
#include "custom_elements/post_element.h"


namespace Kratos
{

    ///@name Kratos Globals
    ///@{

    // Variables definition
    KRATOS_DEFINE_VARIABLE(int, LAYER_ENTITY_TYPE)
    KRATOS_DEFINE_VARIABLE(int, LAYER_PROP_ID)
    KRATOS_DEFINE_VARIABLE(std::string, LAYER_ENTITY_NAME)
    KRATOS_DEFINE_VARIABLE(std::string, LAYER_NAME)

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Enum's
    ///@{

    ///@}
    ///@name Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /// Short class definition.
    /** Detail class definition.
    */
    class KratosLayerApplication : public KratosApplication
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KratosMultiphaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(KratosLayerApplication);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        KratosLayerApplication();

        /// Destructor.
        virtual ~KratosLayerApplication(){}

        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        virtual void Register();

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
            return "Application for layer handling for arbitrary mesh";
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << Info();
            PrintData(rOStream);
        }

        ///// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const
        {
            rOStream << "in KratosLayerApplication:";
            KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
            rOStream << "Variables:" << std::endl;
            KratosComponents<VariableData>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Elements:" << std::endl;
            KratosComponents<Element>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Conditions:" << std::endl;
            KratosComponents<Condition>().PrintData(rOStream);
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

        const PostElement mPostElement2D3N;
        const PostElement mPostElement2D4N;
        const PostElement mPostElement2D6N;
        const PostElement mPostElement2D8N;
        const PostElement mPostElement2D9N;
        const PostElement mPostFaceElement3D3N;
        const PostElement mPostFaceElement3D4N;
        const PostElement mPostFaceElement3D6N;
        const PostElement mPostFaceElement3D8N;
        const PostElement mPostFaceElement3D9N;
        const PostElement mPostElement3D4N;
        const PostElement mPostElement3D10N;
        const PostElement mPostElement3D8N;
        const PostElement mPostElement3D20N;
        const PostElement mPostElement3D27N;
        const PostElement mPostElement3D6N;
        const PostElement mPostElement3D15N;

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
        KratosLayerApplication& operator=(KratosLayerApplication const& rOther);

        /// Copy constructor.
        KratosLayerApplication(KratosLayerApplication const& rOther);


        ///@}

    }; // Class KratosLayerApplication

    ///@}


    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}


} // namespace Kratos

#endif // KRATOS_LAYER_APPLICATION_H_INCLUDED defined

