#ifndef _itkDisplacementFieldTransform_cxx_
#define _itkDisplacementFieldTransform_cxx_

#include "itkDisplacementFieldTransform.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkImageRegionIterator.h"
#include "vnl/vnl_det.h"

#include <itkInverseDeformationFieldImageFilter.h>
#include "itkFixedPointInverseDeformationFieldImageFilter.h"


namespace itk
{



template <class TScalarType, unsigned int NDimensions>
DisplacementFieldTransform<TScalarType, NDimensions>::
DisplacementFieldTransform() : Superclass( SpaceDimension, ParametersDimension )
{
    this->m_InterpolateFunction = InterpolateFunctionType::New();
}



template <class TScalarType, unsigned int NDimensions>
void
DisplacementFieldTransform<TScalarType, NDimensions>::
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
DisplacementFieldTransform<TScalarType, NDimensions>::
PrintSelf(std::ostream &os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
}



template<class TScalarType, unsigned int NDimensions>
const typename DisplacementFieldTransform<TScalarType, NDimensions>::VectorFieldType *
DisplacementFieldTransform<TScalarType, NDimensions>::
GetParametersAsVectorField(void) const
{
    return this->m_VectorField.GetPointer();
}



template<class TScalarType, unsigned int NDimensions>
void
DisplacementFieldTransform<TScalarType, NDimensions>::
SetParametersAsVectorField(const VectorFieldType * field)
{
    // Set displacement field and affect it to the interpolate object
    this->m_VectorField = field;
    this->m_InterpolateFunction->SetInputImage(this->m_VectorField);

    // Compute the derivative weights
    itk::Vector<double, NDimensions> spacing = GetSpacing();
    for (unsigned int i=0; i<NDimensions; i++)
        m_DerivativeWeights[i] = (double)(1.0/spacing[i]);

    this->Modified();
}



template<class TScalarType, unsigned int NDimensions>
bool
DisplacementFieldTransform<TScalarType, NDimensions>::
GetInverse( Self* inverse ) const
{
    // Initial field
    VectorFieldConstPointerType initial_field = this->m_VectorField;//GetParametersAsVectorField();

    // Initialize the field inverter
    typedef itk::FixedPointInverseDeformationFieldImageFilter<VectorFieldType, VectorFieldType> FPInverseType;
    typename FPInverseType::Pointer filter = FPInverseType::New();
    filter->SetInput(         initial_field );
    filter->SetOutputOrigin(  initial_field->GetOrigin() );
    filter->SetSize(          initial_field->GetLargestPossibleRegion().GetSize() );
    filter->SetOutputSpacing( initial_field->GetSpacing() );
    filter->SetNumberOfIterations( 20 );

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
typename DisplacementFieldTransform<TScalarType, NDimensions>::InverseTransformBasePointer
DisplacementFieldTransform<TScalarType, NDimensions>::
GetInverseTransform() const
{
    Pointer inv = New();
    return GetInverse(inv) ? inv.GetPointer() : NULL;
}



template<class TScalarType, unsigned int NDimensions>
typename DisplacementFieldTransform<TScalarType, NDimensions>::OutputPointType
DisplacementFieldTransform<TScalarType, NDimensions>::
TransformPoint(const InputPointType & point) const
{

    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");

    OutputPointType output = point;
    typename InterpolateFunctionType::OutputType vector = m_InterpolateFunction->Evaluate(output);
    for (unsigned int i=0; i<NDimensions; i++)
        output[i] += vector[i];

    return output;
}



template<class TScalarType, unsigned int NDimensions>
typename DisplacementFieldTransform<TScalarType, NDimensions>::OutputVectorType
DisplacementFieldTransform<TScalarType, NDimensions>::
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
typename DisplacementFieldTransform<TScalarType, NDimensions>::OutputVnlVectorType
DisplacementFieldTransform<TScalarType, NDimensions>::
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
typename DisplacementFieldTransform<TScalarType, NDimensions>::OriginType
DisplacementFieldTransform<TScalarType, NDimensions>::
GetOrigin(void)
{
    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");
    return this->m_VectorField->GetOrigin();
}



template<class TScalarType, unsigned int NDimensions>
typename DisplacementFieldTransform<TScalarType, NDimensions>::SpacingType
DisplacementFieldTransform<TScalarType, NDimensions>::
GetSpacing(void)
{
    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");
    return this->m_VectorField->GetSpacing();
}



template<class TScalarType, unsigned int NDimensions>
typename DisplacementFieldTransform<TScalarType, NDimensions>::DirectionType
DisplacementFieldTransform<TScalarType, NDimensions>::
GetDirection(void)
{
    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");
    return this->m_VectorField->GetDirection();
}



template<class TScalarType, unsigned int NDimensions>
typename DisplacementFieldTransform<TScalarType, NDimensions>::RegionType
DisplacementFieldTransform<TScalarType, NDimensions>::
GetLargestPossibleRegion(void)
{
    if (this->m_VectorField.IsNull())
        itkExceptionMacro("No field has been set.");
    return this->m_VectorField->GetLargestPossibleRegion();
}



template<class TScalarType, unsigned int NDimensions>
vnl_matrix_fixed<double,NDimensions,NDimensions>
DisplacementFieldTransform<TScalarType, NDimensions>::
GetSpatialJacobian( const InputPointType & point) const
{

    unsigned int i, j;
    vnl_matrix_fixed<double,NDimensions,NDimensions> J;

    if (!this->m_VectorField)
        itkExceptionMacro("No field has been set, cannot compute Jacobian !");

    double weight;
    itk::Vector<double, NDimensions> spacing = m_VectorField->GetSpacing();

    InputPointType point_prev, point_next;

    for (i = 0; i < NDimensions; ++i)
    {
        point_prev[i] = static_cast<double>(point[i]) - spacing[i];
        point_next[i] = static_cast<double>(point[i]) + spacing[i];
    }

    typename InterpolateFunctionType::OutputType vector_prev = m_InterpolateFunction->Evaluate(point_prev);
    typename InterpolateFunctionType::OutputType vector_next = m_InterpolateFunction->Evaluate(point_next);

    for (i = 0; i < NDimensions; ++i)
    {
        weight = 0.5 * m_DerivativeWeights[i];

        for (j = 0; j < NDimensions; ++j)
            J[i][j] = weight * ( static_cast<double>(vector_next[j]) - static_cast<double>(vector_prev[j]) );

        // Add one on the diagonal to consider the warp and not only the deformation field
        J[i][i] += 1.0;
    }

    return J;
}



template<class TScalarType, unsigned int NDimensions>
typename DisplacementFieldTransform<TScalarType, NDimensions>::ScalarType
DisplacementFieldTransform< TScalarType,  NDimensions >::
GetSpatialJacobianDeterminant( const InputPointType & point) const
{
    return static_cast<ScalarType>(vnl_det(this->GetSpatialJacobian(point)));
}


} // namespace

#endif // _itkDisplacementFieldTransform_cxx_
