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
#include "includes/legacy_structural_app_vars.h"
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
void PostElement::EquationIdVector( EquationIdVectorType& rResult,
                                    const ProcessInfo& CurrentProcessInfo) const
{
    rResult.resize(0);
}

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
        if (rValues.size() != it->second.size())
            rValues.resize(it->second.size());
        std::copy(it->second.begin(), it->second.end(), rValues.begin());
    }
}

void PostElement::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
       std::vector<array_1d<double, 3 > >& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    // special handling to obtain the coordinates of the integration points
    if( rVariable == INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION )
    {
        auto it = mArray1DValuesContainer.find(INTEGRATION_POINT_LOCAL.Key());
        if (it != mArray1DValuesContainer.end())
        // if the integration points are provided for the element via INTEGRATION_POINT_LOCAL,
        // we will use it
        {
            if (rValues.size() != it->second.size())
                rValues.resize(it->second.size());

            GeometryType::IntegrationPointsArrayType integration_points(it->second.size());
            for(std::size_t point = 0; point < it->second.size(); ++point)
                noalias(integration_points[point]) = it->second[point];

            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            GetGeometry().Initialize(integration_points);
            #endif

            VectorType N( GetGeometry().size() );
            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                GetGeometry().ShapeFunctionsValues( N, integration_points[point] );

                noalias( rValues[point] ) = ZeroVector(3);
                for(std::size_t i = 0 ; i < GetGeometry().size() ; ++i)
                    noalias( rValues[point] ) += N[i] * GetGeometry()[i].GetInitialPosition();
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Clean();
            #endif
        }
        else
        // otherwise we use the default integration method, which can be controlled via INTEGRATION_ORDER
        {
            const IntegrationMethod ThisIntegrationMethod = this->GetIntegrationMethod();

            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            GetGeometry().Initialize(ThisIntegrationMethod);
            #endif

            //reading integration points and local gradients
            const GeometryType::IntegrationPointsArrayType& integration_points =
                GetGeometry().IntegrationPoints( ThisIntegrationMethod );

            if (rValues.size() != integration_points.size())
                rValues.resize(integration_points.size());

            const MatrixType& Ncontainer = GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

            VectorType N( GetGeometry().size() );
            for(std::size_t point = 0; point < integration_points.size(); ++point)
            {
                noalias(N) = row( Ncontainer, point );

                noalias( rValues[point] ) = ZeroVector(3);
                for(std::size_t i = 0 ; i < GetGeometry().size() ; ++i)
                    noalias( rValues[point] ) += N[i] * GetGeometry()[i].GetInitialPosition();
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            //clean the internal data of the geometry
            GetGeometry().Clean();
            #endif
        }

        return;
    }

    auto it = mArray1DValuesContainer.find(rVariable.Key());
    if (it == mArray1DValuesContainer.end())
    {
        return;
    }
    else
    {
        if (rValues.size() != it->second.size())
            rValues.resize(it->second.size());
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
        if (rValues.size() != it->second.size())
            rValues.resize(it->second.size());
        std::copy(it->second.begin(), it->second.end(), rValues.begin());
    }
}

PostElement::IntegrationMethod PostElement::GetIntegrationMethod() const
{
    if(this->Has( INTEGRATION_ORDER ))
    {
        if(this->GetValue(INTEGRATION_ORDER) == 1)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_1;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 2)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 3)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_3;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 4)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_4;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 5)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_5;
        }
        else
            KRATOS_ERROR << Info() << " does not support for integration order " << this->GetValue(INTEGRATION_ORDER);
    }
    else if(GetProperties().Has( INTEGRATION_ORDER ))
    {
        if(GetProperties()[INTEGRATION_ORDER] == 1)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_1;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 2)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 3)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_3;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 4)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_4;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 5)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_5;
        }
        else
            KRATOS_ERROR << Info() << " does not support for integration order " << GetProperties()[INTEGRATION_ORDER];
    }
    else
        return GetGeometry().GetDefaultIntegrationMethod(); // default method
}

} // Namespace Kratos
