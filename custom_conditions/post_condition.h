//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 27 Mar 25 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_POST_CONDITION_H_INCLUDED )
#define  KRATOS_POST_CONDITION_H_INCLUDED


// External includes
#include <unordered_map>

// Project includes
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{

/**
 * A condition used to store the values for post-processing
 */
class PostCondition : public Condition
{
    public:
        // Counted pointer of PostCondition
        KRATOS_CLASS_POINTER_DEFINITION(PostCondition);

        /**
         * Default constructor.
         */
        PostCondition();
        PostCondition( IndexType NewId, GeometryType::Pointer pGeometry);
        PostCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        /**
         * Destructor.
         */
        ~PostCondition() override;

        /**
         * Operations.
         */

        Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const final;

        Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const final;

        /**
         * Calculates the local system contributions for this condition. In fact it does nothing.
         */
        void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                   VectorType& rRightHandSideVector,
                                   const ProcessInfo& rCurrentProcessInfo) final;

        void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                     const ProcessInfo& rCurrentProcessInfo) final;

        void EquationIdVector( EquationIdVectorType& rResult,
                               const ProcessInfo& rCurrentProcessInfo) const final;

        void GetDofList( DofsVectorType& ConditionalDofList,
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
         */
        std::string Info() const final
        {
            return "PostCondition";
        }

        /**
         * Print information about this object.
         */
        void PrintInfo(std::ostream& rOStream) const final
        {
            rOStream << this->Info() << " #" << this->Id();
        }

        /**
         * Print object's data.
         */
        void PrintData(std::ostream& rOStream) const final
        {
        }

    private:

        friend class Serializer;

        std::unordered_map<std::size_t, std::vector<double> > mDoubleValuesContainer;   // container for double variables, the key is variable key
        std::unordered_map<std::size_t, std::vector<Vector> > mVectorValuesContainer;
        std::unordered_map<std::size_t, std::vector<array_1d<double, 3> > > mArray1DValuesContainer;

        void save ( Serializer& rSerializer ) const final
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
        }

        void load ( Serializer& rSerializer ) final
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix,
                           VectorType& rRightHandSideVector,
                           const ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag );

}; // Class PostCondition

}  // namespace Kratos.

#endif // KRATOS_POST_CONDITION_H_INCLUDED defined
