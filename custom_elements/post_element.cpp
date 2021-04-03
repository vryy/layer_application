//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 20 Mar 2021 $
//   Revision:            $Revision: 1.0 $
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/post_element.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
PostElement::PostElement()
{
}

PostElement::PostElement( IndexType NewId,
                              GeometryType::Pointer pGeometry)
: Element( NewId, pGeometry )
{
}

PostElement::PostElement( IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: Element( NewId, pGeometry, pProperties )
{
}

/**
 * Destructor. Never to be called manually
 */
PostElement::~PostElement()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Element::Pointer PostElement::Create( IndexType NewId, NodesArrayType const& ThisNodes,
                                      PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer(new PostElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer PostElement::Create( IndexType NewId, GeometryType::Pointer pGeom,
                                      PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer(new PostElement(NewId, pGeom, pProperties));
}

void PostElement::Initialize( const ProcessInfo& rCurrentProcessInfo )
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
void PostElement::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
void PostElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
void PostElement::CalculateAll( MatrixType& rLeftHandSideMatrix,
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
void PostElement::EquationIdVector( EquationIdVectorType& rResult,
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
void PostElement::GetDofList( DofsVectorType& ElementalDofList,
                              const ProcessInfo& CurrentProcessInfo) const
{
    ElementalDofList.resize(0);
}

//************************************************************************************
//************************************************************************************
void PostElement::SetValuesOnIntegrationPoints( const Variable<double>& rVariable,
        const std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    auto it = mDoubleValuesContainer.find(rVariable.Key());
    if (it == mDoubleValuesContainer.end())
    {
        // it = mDoubleValuesContainer.insert(rVariable.Key(), std::vector<double>(rValues.size())); // only compile with C++17
        mDoubleValuesContainer[rVariable.Key()] = std::vector<double>(rValues.size());
        it = mDoubleValuesContainer.find(rVariable.Key());
    }
    else
    {
        if (it->second.size() != rValues.size())
            it->second.resize(rValues.size());
    }

    std::copy(rValues.begin(), rValues.end(), it->second.begin());
}

void PostElement::SetValuesOnIntegrationPoints( const Variable<array_1d<double, 3 > >& rVariable,
        const std::vector<array_1d<double, 3 > >& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    auto it = mArray1DValuesContainer.find(rVariable.Key());
    if (it == mArray1DValuesContainer.end())
    {
        // it = mArray1DValuesContainer.insert(rVariable.Key(), std::vector<double>(rValues.size())); // only compile with C++17
        mArray1DValuesContainer[rVariable.Key()] = std::vector<array_1d<double, 3> >(rValues.size());
        it = mArray1DValuesContainer.find(rVariable.Key());
    }
    else
    {
        if (it->second.size() != rValues.size())
            it->second.resize(rValues.size());
    }

    std::copy(rValues.begin(), rValues.end(), it->second.begin());
}

void PostElement::SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
        const std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    auto it = mVectorValuesContainer.find(rVariable.Key());
    if (it == mVectorValuesContainer.end())
    {
        // it = mVectorValuesContainer.insert(rVariable.Key(), std::vector<double>(rValues.size())); // only compile with C++17
        mVectorValuesContainer[rVariable.Key()] = std::vector<Vector>(rValues.size());
        it = mVectorValuesContainer.find(rVariable.Key());
    }
    else
    {
        if (it->second.size() != rValues.size())
            it->second.resize(rValues.size());
    }

    std::copy(rValues.begin(), rValues.end(), it->second.begin());
}

void PostElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    auto it = mDoubleValuesContainer.find(rVariable.Key());
    if (it == mDoubleValuesContainer.end())
    {
        return;
    }
    else
    {
        std::copy(it->second.begin(), it->second.end(), rValues.begin());
    }
}

void PostElement::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
       std::vector<array_1d<double, 3 > >& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    auto it = mArray1DValuesContainer.find(rVariable.Key());
    if (it == mArray1DValuesContainer.end())
    {
        return;
    }
    else
    {
        std::copy(it->second.begin(), it->second.end(), rValues.begin());
    }
}

void PostElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    auto it = mVectorValuesContainer.find(rVariable.Key());
    if (it == mVectorValuesContainer.end())
    {
        return;
    }
    else
    {
        std::copy(it->second.begin(), it->second.end(), rValues.begin());
    }
}

} // Namespace Kratos

