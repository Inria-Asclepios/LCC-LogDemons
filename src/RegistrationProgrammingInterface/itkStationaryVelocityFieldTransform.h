#ifndef _itkStationaryVelocityFieldTransform_h_
#define _itkStationaryVelocityFieldTransform_h_

#include <itkTransform.h>
#include <itkProcessObject.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkImage.h>
#include <itkStationaryVelocityFieldExponential.h>


namespace itk
{


/**
 *
 * This class defines the non-linear stationary velocity field transformation.
 *
 * There are two templates for this class:
 *
 *   TScalarType  Type of the transformation parameters. Can be "float" or "double".
 *
 *   NDimensions  Dimension of the transformation space.
 *
 * The parameters of a stationary velocity field is a vector field. Since a vector field cannot
 * be simply represented by an simple linear array, we store this vector field into a VectorFieldType
 * which is simply an itk::Image where the type of voxels is itk::Vector. The methods
 * SetParametersAsVectorField and GetParametersAsVectorField respectively allow to set and get the
 * vector field.
 *
 * There is no I/O support as there is no parameters, the user might save the
 * vector field instead.
 *
 * @brief   Stationary velocity field transformation
 *
 * @class   StationaryVelocityFieldTransform (itk)
 *
 * @author  Vincent Garcia, Marco Lorenzi, Asclepios team, INRIA
 *
 */
template <class TScalarType=double, unsigned int NDimensions=3>
class ITK_EXPORT StationaryVelocityFieldTransform : public Transform<TScalarType, NDimensions, NDimensions>
{


public:

    typedef StationaryVelocityFieldTransform                  Self;
    typedef Transform<TScalarType, NDimensions, NDimensions>  Superclass;
    typedef SmartPointer<Self>                                Pointer;
    typedef SmartPointer<const Self>                          ConstPointer;

    typedef Superclass                                        TransformType;
    typedef typename TransformType::ConstPointer              TransformPointerType;

    /**
     * Base inverse transform type. This type should not be changed to the concrete inverse
     * transform type or inheritance would be lost.
     */
    typedef typename Superclass::InverseTransformBaseType     InverseTransformBaseType;
    typedef typename InverseTransformBaseType::Pointer        InverseTransformBasePointer;

    /**
     * Type of the scalar representing coordinate and vector elements
     */
    typedef TScalarType                                       ScalarType;

    /**
     * Type of the input parameters
     */
    typedef typename Superclass::ParametersType               ParametersType;

    /**
     * Type of the Jacobian matrix
     */
    typedef typename Superclass::JacobianType                 JacobianType;

    /**
     * Standard vector type for this class
     */
    typedef typename Superclass::InputVectorType              InputVectorType;
    typedef typename Superclass::OutputVectorType             OutputVectorType;

    /**
     * Standard covariant vector type for this class
     */
    typedef typename Superclass::InputCovariantVectorType     InputCovariantVectorType;
    typedef typename Superclass::OutputCovariantVectorType    OutputCovariantVectorType;

    /**
     * Standard vnl_vector type for this class
     */
    typedef typename Superclass::InputVnlVectorType           InputVnlVectorType;
    typedef typename Superclass::OutputVnlVectorType          OutputVnlVectorType;

    /**
     * Standard coordinate point type for this class
     */
    typedef typename Superclass::InputPointType               InputPointType;
    typedef typename Superclass::OutputPointType              OutputPointType;

    /**
     * Type related to the container of field
     */
    typedef itk::Vector<ScalarType, NDimensions>              VectorType;
    typedef itk::Image<VectorType, NDimensions>               VectorFieldType;
    typedef typename VectorFieldType::Pointer                 VectorFieldPointerType;
    typedef typename VectorFieldType::ConstPointer            VectorFieldConstPointerType;
    typedef typename VectorFieldType::IndexType               VectortFieldIndexType;

    /**
     * Type relative to the field geometry
     */
    typedef typename VectorFieldType::PointType               OriginType;
    typedef typename VectorFieldType::SpacingType             SpacingType;
    typedef typename VectorFieldType::DirectionType           DirectionType;
    typedef typename VectorFieldType::RegionType              RegionType;

    /**
     * Interpolation type
     */
    typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<VectorFieldType, ScalarType>
                                                              InterpolateFunctionType;
    typedef typename InterpolateFunctionType::Pointer         InterpolateFunctionPointerType;

    /**
     * Object that computes the exponential of a staitionary velocity field
     */
    typedef itk::StationaryVelocityFieldExponential<VectorFieldType, VectorFieldType>
                                                              SVFExponentialType;

    /**
     * Dimension of the domain space
     */
    itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
    itkStaticConstMacro(ParametersDimension, unsigned int, NDimensions);

    /**
     * Generic constructors.
     */
    itkNewMacro( Self );
    itkTypeMacro( StationaryVelocityFieldTransform, Transform );

    /**
     * Method not implemented for stationary velocity fields.
     * @retun parameters
     */
    virtual const ParametersType&       GetParameters(void) const
    {
        itkExceptionMacro("This type of transform does not handle any parameters.");
    }

    /**
     * Method not implemented for stationary velocity fields.
     * @param param parameters
     */
    virtual void                        SetFixedParameters( const ParametersType & param )
    {
        itkExceptionMacro("This type of transform does not handle any parameters.");
    }

    /**
     * Method not implemented for stationary velocity fields.
     * @param param parameters
     */
    virtual void                        SetParameters( const ParametersType & param )
    {
        itkExceptionMacro("This type of transform does not handle any parameters.");
    }


    /**
     * Method not implemented for stationary velocity fields.
     * @param p value
     */
    virtual void                        SetParametersByValue( const ParametersType & p )
    {
        itkExceptionMacro("This type of transform does not handle any parameters.");
    }

    /**
     * Sets the transform to identity.
     */
    virtual void                        SetIdentity(void);

    /**
     * Gets the vector field containing the transformation parameters.
     * @return  vector  field
     */
    virtual const VectorFieldType *     GetParametersAsVectorField(void) const;

    /**
     * Sets the vector field containing the transformation parameters.
     * @param field vector field
     */
    virtual void                        SetParametersAsVectorField(const VectorFieldType * field);

    /**
     * Gets the inverse of the transformation.
     * The proposed implementation is not perfectly correct in a theoretical point of view.
     * This will be fixed in the near future.
     * @param  inverse  transformation to invert
     * @return true if the inversion was a success, false otherwise
     */
    virtual bool                        GetInverse(Self* inverse) const;

    /**
     * Gets the inverse of thes transformation.
     * The proposed implementation is not perfectly correct in a theoretical point of view.
     * This will be fixed in the near future.
     * @return inverted transformation
     */
    virtual InverseTransformBasePointer GetInverseTransform(void) const;

    /**
     * Transforms a point. This method should be used to transform one or a few number of points.
     * This method should not be used to resample an entire image using for instance a
     * itk::ResampleImageFilter. Instead, first estimate the corresponding displacement field using
     * GetDisplacementFieldAsVectorField, and second resample the image using this displacement field.
     * @param  point  point
     * @return transformed point
     */
    virtual OutputPointType             TransformPoint(const InputPointType  & point) const;

    /**
     * Transforms a vector. See TransformPoint method for more details.
     * @param  vector  vector
     * @return transformed vector
     */
    virtual OutputVectorType            TransformVector(const InputVectorType & vector) const;

    /**
     * Transforms a vnl_vector. See TransformPoint method for more details.
     * @param  vector  vector
     * @return transformed vector
     */
    virtual OutputVnlVectorType         TransformVector(const InputVnlVectorType & vector) const;

    /**
     * Method not implemented for stationary velocity fields.
     * @param  vector  vector
     * @return transformed vector
     */
    virtual OutputCovariantVectorType   TransformCovariantVector(const InputCovariantVectorType & vector) const
    {
        itkExceptionMacro("Cannot transform covariant vector!");
    }

    /**
     * Computes and returns the exponential of the stationary velocity field (i.e. displacement
     * field). The output field is actually not a transformation object but a simple itk::Image
     * (where displacements are stored into itk::Vector objects, hence "vector field") since there
     * is no displacement field transform object in ITK so far.
     * @param  scheme  numerical scheme used to estimate the exponential of the velocity field
     * @return vector field containing the displacements
     */
    VectorFieldPointerType              GetDisplacementFieldAsVectorField(typename SVFExponentialType::NumericalScheme scheme = SVFExponentialType::FORWARD_EULER) const;

    /**
     * Gets the origin of the field.
     * @retun origin
     */
    OriginType                          GetOrigin(void) const;

    /**
     * Gets the spacing of the field.
     * @retun spacing
     */
    SpacingType                         GetSpacing(void) const;

    /**
     * Gets the direction of the field.
     * @retun direction
     */
    DirectionType                       GetDirection(void) const;

    /**
     * Gets the region object that defines the size and starting index for the largest possible
     * region this image could represent.
     * @retun region
     */
    RegionType                          GetLargestPossibleRegion(void) const;

    /**
     * Method not implemented for stationary velocity fields.
     * @param point point
     * @retun Jacobian matrix
     */
    virtual const JacobianType &        GetJacobian(const InputPointType  & point) const
    {
        itkExceptionMacro("This type of transform does not handle Jacobian!");
    }





protected:

    /**
     * Prints contents.
     */
    void PrintSelf(std::ostream &os, Indent indent) const;

    /**
     * Default constructor.
     */
    StationaryVelocityFieldTransform(void);

    /**
     * Destructor.
     */
    virtual ~StationaryVelocityFieldTransform(){};



private:

    /**
     * Parameters = vector field
     */
    VectorFieldConstPointerType         m_VectorField;

    /**
     * Interpolation function
     */
    InterpolateFunctionPointerType      m_InterpolateFunction;

};


} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStationaryVelocityFieldTransform.txx"
#endif

#endif
