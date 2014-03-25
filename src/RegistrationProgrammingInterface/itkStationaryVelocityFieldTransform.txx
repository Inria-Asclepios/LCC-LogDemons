#ifndef _itkStationaryVelocityFieldTransform_cxx_
#define _itkStationaryVelocityFieldTransform_cxx_

#include "itkStationaryVelocityFieldTransform.h"

#include <itkNeighborhoodAlgorithm.h>
#include <itkImageRegionIterator.h>
#include <vnl/vnl_det.h>
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>
#include <itkMultiplyByConstantImageFilter.h>
#include <itkWarpVectorImageFilter.h>
#include <itkImageRegionIterator.h>  
#include <itkDivideByConstantImageFilter.h>

namespace itk
{



template <class TScalarType, unsigned int NDimensions>
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
StationaryVelocityFieldTransform():Superclass(NDimensions)
{
    this->m_InterpolateFunction = InterpolateFunctionType::New();
}



template <class TScalarType, unsigned int NDimensions>
void
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
SetIdentity(void)
{
    TScalarType            value = 0;
    VectorFieldPointerType field = VectorFieldType::New();
    field->Allocate();
    field->FillBuffer(value);
    field->Register();
    this->m_VectorField = field;
}



template<class TScalarType, unsigned int NDimensions>
void
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
PrintSelf(std::ostream &os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
}



template<class TScalarType, unsigned int NDimensions>
const typename StationaryVelocityFieldTransform<TScalarType, NDimensions>::VectorFieldType *
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
GetParametersAsVectorField(void) const
{
    return this->m_VectorField.GetPointer();
}



template<class TScalarType, unsigned int NDimensions>
void
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
SetParametersAsVectorField(const VectorFieldType * field)
{
    this->m_VectorField = field;
    this->m_InterpolateFunction->SetInputImage(this->m_VectorField);
    this->Modified();
}



template<class TScalarType, unsigned int NDimensions>
bool
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
GetInverse( Self* inverse ) const
{
    // Initial field
    VectorFieldConstPointerType initial_field = this->m_VectorField;

    // Initialize the field inverter
    typedef itk::MultiplyByConstantImageFilter<VectorFieldType, int, VectorFieldType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(    initial_field );
    filter->SetConstant( -1 );
    filter->SetInPlace(  false);

    // Update the filter
    filter->UpdateLargestPossibleRegion();

    // Get the inverted field and set the orientation that has been lost
    VectorFieldPointerType inverted_field = filter->GetOutput();
    inverted_field->SetDirection( initial_field->GetDirection() );

    // Set the displacement field to the current objet
    inverse->SetParametersAsVectorField( inverted_field );

    return true;
}



template<class TScalarType, unsigned int NDimensions>
typename StationaryVelocityFieldTransform<TScalarType, NDimensions>::InverseTransformBasePointer
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
GetInverseTransform() const
{
   Pointer inv = New();
   return GetInverse(inv) ? inv.GetPointer() : NULL;
}



template<class TScalarType, unsigned int NDimensions>
typename StationaryVelocityFieldTransform<TScalarType, NDimensions>::OutputPointType
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
TransformPoint(const InputPointType & point) const
{
    double minpixelspacing = m_VectorField->GetSpacing()[0];
    for (unsigned int i = 0; i<3; ++i)
    {
        if ( m_VectorField->GetSpacing()[i] < minpixelspacing/2 )
        {
            minpixelspacing = m_VectorField->GetSpacing()[i];
        }
    }

    typedef typename itk::ImageRegionIterator<VectorFieldType> IteratorType;

    IteratorType InputIt(const_cast<VectorFieldType*>(m_VectorField.GetPointer()), m_VectorField->GetRequestedRegion());

    float norm2,maxnorm2=0;
    for( InputIt.GoToBegin(); !InputIt.IsAtEnd(); ++InputIt )
    {
        norm2 = InputIt.Get().GetSquaredNorm();
        if (maxnorm2<norm2) maxnorm2=norm2;
    }

    maxnorm2 /= vnl_math_sqr(minpixelspacing);
    
    float numiterfloat = 2.0 +
            0.5 * vcl_log(maxnorm2)/vnl_math::ln2;
    std::cout<<"Num iter float: "<<numiterfloat<<std::endl;

    unsigned int constant=static_cast<unsigned int>(1<<static_cast<unsigned int>(abs(numiterfloat)));

    typedef typename itk::DivideByConstantImageFilter<VectorFieldType,float,VectorFieldType> DividerType;
    typename DividerType::Pointer Divider=DividerType::New();
    Divider->SetInput(m_VectorField);
    Divider->SetConstant( constant );
    std::cout<<"divider "<<  static_cast<unsigned int>(1<<static_cast<unsigned int>(abs(numiterfloat))) <<std::endl;

    Divider->Update();

    VectorFieldPointerType UpdatedVector = Divider->GetOutput();

    typedef typename itk::WarpVectorImageFilter<VectorFieldType,VectorFieldType,VectorFieldType> VectorWarperType;
    typename VectorWarperType::Pointer VectorWarper=VectorWarperType::New();

    VectorWarper->SetOutputOrigin(UpdatedVector->GetOrigin());
    VectorWarper->SetOutputSpacing(UpdatedVector->GetSpacing());
    VectorWarper->SetOutputDirection(UpdatedVector->GetDirection());


    OutputPointType output = point;

    for (unsigned int i =0; i < constant;++i)
    { 
        this->m_InterpolateFunction->SetInputImage(UpdatedVector);
        typename InterpolateFunctionType::OutputType vector = m_InterpolateFunction->Evaluate(output);

        for (unsigned int i=0; i<NDimensions; i++)
            output[i] += vector[i];

        VectorWarper->SetInput(UpdatedVector);
        VectorWarper->SetDisplacementField(UpdatedVector);
        VectorWarper->Update();

        UpdatedVector=VectorWarper->GetOutput();
        UpdatedVector->DisconnectPipeline();
    }

    return output;
}



template<class TScalarType, unsigned int NDimensions>
typename StationaryVelocityFieldTransform<TScalarType, NDimensions>::OutputVectorType
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
TransformVector(const InputVectorType & vector) const
{
    // Convert vector into point
    InputPointType point_0;
    for (unsigned int i=0; i<NDimensions; i++)
        point_0[i] = vector[i];

    // Transform point
    InputPointType point_1 = TransformPoint(point_0);

    // Convert point into vector
    OutputVectorType vector_1;
    for (unsigned int i=0; i<NDimensions; i++)
        vector_1[i] = point_1[i];
    return vector_1;
}



template<class TScalarType, unsigned int NDimensions>
typename StationaryVelocityFieldTransform<TScalarType, NDimensions>::OutputVnlVectorType
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
TransformVector(const InputVnlVectorType & vector) const
{
    // Convert vector into point
    InputPointType point_0;
    for (unsigned int i=0; i<NDimensions; i++)
        point_0[i] = vector[i];

    // Transform point
    InputPointType point_1 = TransformPoint(point_0);
    return point_1.GetVnlVector();
}



template<class TScalarType, unsigned int NDimensions>
typename StationaryVelocityFieldTransform<TScalarType, NDimensions>::VectorFieldPointerType
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
GetDisplacementFieldAsVectorField(typename SVFExponentialType::NumericalScheme scheme) const
{
    typename SVFExponentialType::Pointer exp = SVFExponentialType::New();
    exp->SetInput( this->m_VectorField );
    exp->SetIterativeScheme( scheme );
    exp->UpdateLargestPossibleRegion();
    return exp->GetOutput();
}



template<class TScalarType, unsigned int NDimensions>
typename StationaryVelocityFieldTransform<TScalarType, NDimensions>::OriginType
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
GetOrigin(void) const
{
    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");
    return this->m_VectorField->GetOrigin();
}



template<class TScalarType, unsigned int NDimensions>
typename StationaryVelocityFieldTransform<TScalarType, NDimensions>::SpacingType
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
GetSpacing(void) const
{
    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");
    return this->m_VectorField->GetSpacing();
}



template<class TScalarType, unsigned int NDimensions>
typename StationaryVelocityFieldTransform<TScalarType, NDimensions>::DirectionType
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
GetDirection(void) const
{
    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");
    return this->m_VectorField->GetDirection();
}



template<class TScalarType, unsigned int NDimensions>
typename StationaryVelocityFieldTransform<TScalarType, NDimensions>::RegionType
StationaryVelocityFieldTransform<TScalarType, NDimensions>::
GetLargestPossibleRegion(void) const
{
    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");
    return this->m_VectorField->GetLargestPossibleRegion();
}


} // namespace

#endif // _itkStationaryVelocityFieldTransform_cxx_
