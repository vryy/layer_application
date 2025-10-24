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
#include "custom_conditions/post_condition.h"
#include "custom_elements/post_element.h"
#include "custom_elements/post_ups_element.h"


namespace Kratos
{

    ///@name Kratos Globals
    ///@{

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
    class KRATOS_API(LAYER_APPLICATION) KratosLayerApplication : public KratosApplication
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
        ~KratosLayerApplication() override {}

        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        void Register() override;

        void RegisterVariables() override;

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
        std::string Info() const override
        {
            return "Application for layer handling for arbitrary mesh";
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << Info();
            PrintData(rOStream);
        }

        ///// Print object's data.
        void PrintData(std::ostream& rOStream) const override
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

        const PostUPSElement mPostUPSElement2D3N;
        const PostUPSElement mPostUPSElement2D4N;
        const PostUPSElement mPostUPSElement2D6N;
        const PostUPSElement mPostUPSElement2D8N;
        const PostUPSElement mPostUPSElement2D9N;
        const PostUPSElement mPostUPSFaceElement3D3N;
        const PostUPSElement mPostUPSFaceElement3D4N;
        const PostUPSElement mPostUPSFaceElement3D6N;
        const PostUPSElement mPostUPSFaceElement3D8N;
        const PostUPSElement mPostUPSFaceElement3D9N;
        const PostUPSElement mPostUPSElement3D4N;
        const PostUPSElement mPostUPSElement3D10N;
        const PostUPSElement mPostUPSElement3D8N;
        const PostUPSElement mPostUPSElement3D20N;
        const PostUPSElement mPostUPSElement3D27N;
        const PostUPSElement mPostUPSElement3D6N;
        const PostUPSElement mPostUPSElement3D15N;

        const PostCondition mPostSurfaceCondition3D3N;
        const PostCondition mPostSurfaceCondition3D4N;
        const PostCondition mPostSurfaceCondition3D6N;
        const PostCondition mPostSurfaceCondition3D8N;
        const PostCondition mPostSurfaceCondition3D9N;

        const PostCondition mPostLineCondition2D2N;
        const PostCondition mPostLineCondition2D3N;
        const PostCondition mPostLineCondition3D2N;
        const PostCondition mPostLineCondition3D3N;

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
