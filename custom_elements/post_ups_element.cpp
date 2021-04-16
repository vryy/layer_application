//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 10 Apr 2021 $
//   Revision:            $Revision: 1.0 $
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/post_ups_element.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
PostUPSElement::PostUPSElement()
{
}

PostUPSElement::PostUPSElement( IndexType NewId,
                              GeometryType::Pointer pGeometry)
: Element( NewId, pGeometry )
{
}

PostUPSElement::PostUPSElement( IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: Element( NewId, pGeometry, pProperties )
{
}

/**
 * Destructor. Never to be called manually
 */
PostUPSElement::~PostUPSElement()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Element::Pointer PostUPSElement::Create( IndexType NewId, NodesArrayType const& ThisNodes,
                                      PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer(new PostUPSElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer PostUPSElement::Create( IndexType NewId, GeometryType::Pointer pGeom,
                                      PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer(new PostUPSElement(NewId, pGeom, pProperties));
}

void PostUPSElement::Initialize( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************
/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void PostUPSElement::CalculateRightHandSide( VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType matrix = Matrix();
    CalculateAll( matrix, rRightHandSideVector,
                  rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag,
                  CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************
/**
 * calculates this contact element's local contributions
 */
void PostUPSElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                        VectorType& rRightHandSideVector,
                                        const ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}
    //************************************************************************************
//************************************************************************************    /
/**
 * This function calculates all system contributions due to the contact problem
 * with regard to the current master and slave partners.
 * All Conditions are assumed to be defined in 2D/3D space and having 2/3 DOFs per node
 */
void PostUPSElement::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo,
                                bool CalculateStiffnessMatrixFlag,
                                bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    rLeftHandSideMatrix.resize(0, 0, false);
    rRightHandSideVector.resize(0, false);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
/**
* Setting up the EquationIdVector for the current partners.
* All conditions are assumed to be defined in 2D/3D space with 2/3 DOFs per node.
* All Equation IDs are given Master first, Slave second
*/
void PostUPSElement::EquationIdVector( EquationIdVectorType& rResult,
                                    const ProcessInfo& CurrentProcessInfo) const
{
    rResult.resize(0);
}

//************************************************************************************
//************************************************************************************
/**
 * Setting up the DOF list for the current partners.
 * All conditions are assumed to be defined in 2D/3D space with 2/3 DOFs per Node.
 * All DOF are given Master first, Slave second
 */
//************************************************************************************
//************************************************************************************
void PostUPSElement::GetDofList( DofsVectorType& ElementalDofList,
                              const ProcessInfo& CurrentProcessInfo) const
{
    ElementalDofList.resize(0);
}

//************************************************************************************
//************************************************************************************
void PostUPSElement::SetValuesOnIntegrationPoints( const Variable<double>& rVariable,
        const std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if (rVariable == AIR_PRESSURE)
    {
        if (mAirPressure.size() != rValues.size())
            mAirPressure.resize(rValues.size());
        std::copy(rValues.begin(), rValues.end(), mAirPressure.begin());
    }
    else if (rVariable == WATER_PRESSURE)
    {
        if (mWaterPressure.size() != rValues.size())
            mWaterPressure.resize(rValues.size());
        std::copy(rValues.begin(), rValues.end(), mWaterPressure.begin());
    }
}

void PostUPSElement::SetValuesOnIntegrationPoints( const Variable<array_1d<double, 3 > >& rVariable,
        const std::vector<array_1d<double, 3 > >& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if (rVariable == DISPLACEMENT)
    {
        if (mDisplacement.size() != rValues.size())
            mDisplacement.resize(rValues.size());
        std::copy(rValues.begin(), rValues.end(), mDisplacement.begin());
    }
}

void PostUPSElement::SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
        const std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == STRESSES)
    {
        if (mStress.size() != rValues.size())
            mStress.resize(rValues.size());
        std::copy(rValues.begin(), rValues.end(), mStress.begin());
    }
}

void PostUPSElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == AIR_PRESSURE)
    {
        if (rValues.size() != mAirPressure.size())
            rValues.resize(mAirPressure.size());
        std::copy(mAirPressure.begin(), mAirPressure.end(), rValues.begin());
    }
    else if (rVariable == WATER_PRESSURE)
    {
        if (rValues.size() != mWaterPressure.size())
            rValues.resize(mWaterPressure.size());
        std::copy(mWaterPressure.begin(), mWaterPressure.end(), rValues.begin());
    }
}

void PostUPSElement::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
       std::vector<array_1d<double, 3 > >& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == DISPLACEMENT)
    {
        if (rValues.size() != mDisplacement.size())
            rValues.resize(mDisplacement.size());
        std::copy(mDisplacement.begin(), mDisplacement.end(), rValues.begin());
    }
}

void PostUPSElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == STRESSES)
    {
        if (rValues.size() != mStress.size())
            rValues.resize(mStress.size());
        std::copy(mStress.begin(), mStress.end(), rValues.begin());
    }
}

} // Namespace Kratos

