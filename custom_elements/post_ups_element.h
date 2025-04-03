//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 10 Apr 21 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_POST_UPS_ELEMENT_H_INCLUDED )
#define  KRATOS_POST_UPS_ELEMENT_H_INCLUDED


// External includes
#include <unordered_map>

// Project includes
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{

/**
 * An element used to store the displacement, stress and pressure for post-processing
 */
class PostUPSElement : public Element
{
    public:
        // Counted pointer of PostUPSElement
        KRATOS_CLASS_POINTER_DEFINITION(PostUPSElement);

        /**
         * Default constructor.
         */
        PostUPSElement();
        PostUPSElement( IndexType NewId, GeometryType::Pointer pGeometry);
        PostUPSElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        /**
         * Destructor.
         */
        ~PostUPSElement() override;

        /**
         * Operations.
         */

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const final;

        Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const final;

        /**
         * Calculates the local system contributions for this contact element
         */
        void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                   VectorType& rRightHandSideVector,
                                   const ProcessInfo& rCurrentProcessInfo) final;

        void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                     const ProcessInfo& rCurrentProcessInfo) final;

        void EquationIdVector( EquationIdVectorType& rResult,
                               const ProcessInfo& rCurrentProcessInfo) const final;

        void GetDofList( DofsVectorType& ElementalDofList,
                         const ProcessInfo& CurrentProcessInfo) const final;

        void Initialize(const ProcessInfo& rCurrentProcessInfo) final;

        /**
         * Get and Set values
         */

        void SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
                             const std::vector<double>& rValues,
                             const ProcessInfo& rCurrentProcessInfo) final;

        void SetValuesOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                             const std::vector<array_1d<double, 3 > >& rValues,
                             const ProcessInfo& rCurrentProcessInfo) final;

        void SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
                             const std::vector<Vector>& rValues,
                             const ProcessInfo& rCurrentProcessInfo) final;

        void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                             std::vector<double>& rValues,
                             const ProcessInfo& rCurrentProcessInfo) final;

        void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                             std::vector<array_1d<double, 3 > >& rValues,
                             const ProcessInfo& rCurrentProcessInfo) final;

        void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                             std::vector<Vector>& rValues,
                             const ProcessInfo& rCurrentProcessInfo) final;

        /**
         * Turn back information as a string.
         * (DEACTIVATED)
         */
        std::string Info() const final
        {
            return "PostUPSElement";
        }

        /**
         * Print information about this object.
         * (DEACTIVATED)
         */
        void PrintInfo(std::ostream& rOStream) const final
        {
            rOStream << this->Info() << " #" << this->Id();
        }

        /**
         * Print object's data.
         * (DEACTIVATED)
         */
        void PrintData(std::ostream& rOStream) const final
        {
        }

    private:

        friend class Serializer;

        std::vector<double> mAirPressure;
        std::vector<double> mWaterPressure;
        std::vector<array_1d<double, 3> > mDisplacement;
        std::vector<Vector> mStress;

        void save ( Serializer& rSerializer ) const final
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Element )
        }

        void load ( Serializer& rSerializer ) final
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Element )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix,
                           VectorType& rRightHandSideVector,
                           const ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag);

}; // Class PostUPSElement

}  // namespace Kratos.

#endif // KRATOS_POST_UPS_ELEMENT_H_INCLUDED defined
