/*=========================================================================

Program:   vtkINRIA3D
Module:    $Id: itkDisplacementFieldTransform.h 1 2008-01-22 19:01:33Z ntoussaint $
Language:  C++
Author:    $Author: ntoussaint vgarcia $
Date:      $Date: 2008-01-22 20:01:33 +0100 (Tue, 22 Jan 2008) $
Version:   $Revision: 1 $

Copyright (c) 2007 INRIA - Asclepios Project. All rights reserved.
See Copyright.txt for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef _itkDisplacementFieldTransform_h_
#define _itkDisplacementFieldTransform_h_

#include <itkTransform.h>
#include <itkProcessObject.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkImage.h>
#include "itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h"


namespace itk
{


/**
 *
 * This class defines the non-linear displacement field transformation.
 *
 * There are two templates for this class:
 *
 *   TScalarType  Type of the transformation parameters. Can be "float" or "double".
 *
 *   NDimensions  Dimension of the transformation space.
 *
 * The parameters of a displacement field is a vector field. Since a vector field cannot be simply
 * represented by an simple linear array, we store this vector field into a VectorFieldType which is
 * simply an itk::Image where the type of voxels is itk::Vector. The methods
 * SetParametersAsVectorField and GetParametersAsVectorField respectively allow to set and get the
 * vector field.
 *
 * This transform is only able to transform points and vectors. Transforming a covariant vector has
 * no meaning in the sense that we need to transform the vector with the jacobian of the
 * transformation with respect to the parameters, which cannot be computed.
 *
 * The user can compute the Jacobian of the transformation with respect to the coordinates
 * at a specific location with the method GetJacobianWithRespectToCoordinates(), and its determinant
 * (scalar) with the method GetJacobianDeterminantWithRespectToCoordinates(). However, the
 * determinant of the transformation with respect its parameters can not be computed as there is no
 * parametric handling.
 *
 * There is no I/O support as there is no parameters, the user might save the
 * vector field instead.
 *
 * @brief  Displacement field transformation
 *
 * @class  DisplacementFieldTransform (itk)
 *
 * @author Nicolas Toussaint, Vincent Garcia INRIA
 */
template <class TScalarType=float, unsigned int NDimensions=3>
class ITK_EXPORT DisplacementFieldTransform : public Transform<TScalarType, NDimensions, NDimensions>
{


public:

    typedef DisplacementFieldTransform                        Self;
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
     * Dimension of the domain space
     */
    itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
    itkStaticConstMacro(ParametersDimension, unsigned int, NDimensions);

    /**
     * Generic constructors.
     */
    itkNewMacro( Self );
    itkTypeMacro( DisplacementFieldTransform, Transform );

    /**
     * Gets the transformation parameters.
     */
    virtual const ParametersType&       GetParameters(void) const
    {
        itkExceptionMacro("This type of transform does not handle any parameters.");
    }

    /**
     * Sets the transformation parameters.
     */
    virtual void                        SetParameters( const ParametersType & )
    {
        itkExceptionMacro("This type of transform does not handle any parameters.");
    }

    /**
     * Sets the transformation parameters.
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
     * Gets the parameters as vector field.
     */
    virtual const VectorFieldType *     GetParametersAsVectorField(void) const;

    /**
     * Sets the parameters as vector field.
     */
    virtual void                        SetParametersAsVectorField(const VectorFieldType * field);

    /**
     * Gets an inverse of this transformation.
     */
    virtual bool                        GetInverse(Self* inverse) const;

    /**
     * Gets an inverse of this transformation.
     */
    virtual InverseTransformBasePointer GetInverseTransform(void) const;

    /**
     * Transforms a point.
     */
    virtual OutputPointType             TransformPoint(const InputPointType  & point) const;

    /**
     * Transforms a vector.
     */
    virtual OutputVectorType            TransformVector(const InputVectorType & vector) const;

    /**
     * Transforms a vnl_vector.
     */
    virtual OutputVnlVectorType         TransformVector(const InputVnlVectorType & vector) const;

    /**
     * Transform a CovariantVector.
     */
    virtual OutputCovariantVectorType   TransformCovariantVector(const InputCovariantVectorType &) const
    {
        itkExceptionMacro("Cannot transform covariant vector!");
    }

    /**
     * Gets the origin of the field.
     */
    OriginType                          GetOrigin(void);

    /**
     * Gets the spacing of the field.
     */
    SpacingType                         GetSpacing(void);

    /**
     * Gets the direction of the field.
     */
    DirectionType                       GetDirection(void);

    /**
     * Gets the region object that defines the size and starting index for the largest possible
     * region this image could represent.
     */
    RegionType                          GetLargestPossibleRegion(void);

    /**
     * This method computes the Jacobian matrix of the transformation
     * at a given input point. The rank of the Jacobian will also indicate
     * if the transform is invertible at this point.
     */
    virtual const JacobianType &        GetJacobian(const InputPointType  & point) const
    {
        itkExceptionMacro("This type of transform does not handle Jacobian!");
    }

    /**
     * Gets the jacobian of the transformation with respect to the coordinates at a specific point.
     */
    virtual vnl_matrix_fixed<double,NDimensions,NDimensions>
                                        GetSpatialJacobian(const InputPointType  & point) const;
    /**
     * Gets the jacobian determinant of transformation with respect to the coordinates at a specific point.
     */
    virtual ScalarType                  GetSpatialJacobianDeterminant(const InputPointType  & point) const;



protected:

    /**
     * Prints contents.
     */
    void PrintSelf(std::ostream &os, Indent indent) const;

    /**
     * Default constructor.
     */
    DisplacementFieldTransform(void);

    /**
     * Destructor.
     */
    virtual ~DisplacementFieldTransform(){};



private:

    /**
     * Parameters = vector field
     */
    VectorFieldConstPointerType     m_VectorField;

    /**
     * Interpolation function
     */
    InterpolateFunctionPointerType  m_InterpolateFunction;

    /**
     * Weights used to scale partial derivatives during Jacobian computation
     */
    double                          m_DerivativeWeights[NDimensions];

};


} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDisplacementFieldTransform.txx"
#endif

#endif
