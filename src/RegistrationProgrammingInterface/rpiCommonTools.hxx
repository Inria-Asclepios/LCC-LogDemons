#ifndef _RPI_COMMON_TOOLS_HXX_
#define _RPI_COMMON_TOOLS_HXX_


#include <iostream>
#include <sstream>

#include <itkImageBase.h>
#include <itkImageIOBase.h>
#include <itkDisplacementFieldTransform.h>
#include <itkStationaryVelocityFieldTransform.h>
#include <itkAffineTransform.h>
#include <itkEuler3DTransform.h>

#include <rpiRegistrationMethod.hxx>


// Namespace RPI : Registration Programming Interface
namespace rpi
{


/**
 * Reads image information.
 * @param   fileName  name of the image to read
 * @return  image information
 */
inline itk::ImageIOBase::Pointer readImageInformation( std::string fileName );


/**
 * Reads image.
 * @param   fileName  name of the image to read
 * @return  image
 */
template<class TImage>
typename TImage::Pointer readImage( std::string fileName );


/**
 * Registers useful ITK transformations with the itk::TransformFactory object.
 */
void registerUsefulITKTransformations(void);


/**
 * Reads an ITK Euler3D transformation from an input file.
 * @param   filename  input file name
 * @return  Euler3D transformation
 */
template<class TTransformScalarType>
typename itk::Euler3DTransform<TTransformScalarType>::Pointer
readEuler3DTransformation( std::string fileName );


/**
 * Reads an ITK affine 3D transformation from an input file.
 * @param   filename  input file name
 * @return  Affine transformation
 */
template<class TTransformScalarType>
typename itk::AffineTransform<TTransformScalarType, 3>::Pointer
readAffineTransformation( std::string fileName );


/**
 * Reads an ITK linear 3D transformation from an input file.
 * @param   filename  input file name
 * @return  transformation
 */
template<class TTransformScalarType>
typename itk::Transform<TTransformScalarType, 3, 3>::Pointer
readLinearTransformation( std::string fileName );


/**
 * Reads an ITK displacement field 3D from an input image.
 * @param   filename  input image name
 * @return  displacement field
 */
template<class TTransformScalarType>
typename itk::DisplacementFieldTransform<TTransformScalarType, 3>::Pointer
readDisplacementField( std::string fileName );


/**
 * Reads an ITK stationary velocity field 3D from an input image.
 * @param   filename  input image name
 * @return  stationary velocity field
 */
template<class TTransformScalarType>
typename itk::StationaryVelocityFieldTransform<TTransformScalarType, 3>::Pointer
readStationaryVelocityField( std::string fileName );


/**
 * Converts a linear 3D transformation into a displacement field 3D transformation.
 * The geometry of the displacement field is taken from the input image.
 * @param   image            image
 * @param   linearTransform  linear transformation
 * @return  displacement field transformation
 */
template<class TLinearScalarType, class TFieldScalarType, class TImage>
typename itk::DisplacementFieldTransform<TFieldScalarType, TImage::ImageDimension>::Pointer
linearToDisplacementFieldTransformation(
        TImage * image,
        itk::Transform<TLinearScalarType, TImage::ImageDimension, TImage::ImageDimension> * transform );


/**
 * Converts a linear 3D transformation into a stationary velocity field 3D transformation.
 * The geometry of the displacement field is taken from the input image.
 * @param   image            image
 * @param   linearTransform  linear transformation
 * @return  stationary velocity field transformation
 */
template<class TLinearScalarType, class TFieldScalarType, class TImage>
typename itk::StationaryVelocityFieldTransform<TFieldScalarType, TImage::ImageDimension>::Pointer
linearToStationaryVelocityFieldTransformation(
        TImage * image,
        itk::Transform<TLinearScalarType, TImage::ImageDimension, TImage::ImageDimension> * transform );


/**
 * Gets and writes the output transformation into an output file.
 * @param  transformation  registration object
 * @param  fileName        name of the output file
 */
template<class TTransformScalarType, int TDimension>
void writeLinearTransformation(
        itk::Transform<TTransformScalarType, TDimension, TDimension> * transform,
        std::string fileName );


/**
 * Writes a displacement field into an output file.
 * @param  field     displacement field
 * @param  fileName  name of the output file
 */
template<class TTransformScalarType, int TDimension>
void writeDisplacementFieldTransformation(
        itk::Transform<TTransformScalarType, TDimension, TDimension> * field,
        std::string fileName );


/**
 * Writes a stationary velocity field into an output file.
 * @param  field     stationary velocity field
 * @param  fileName  name of the output file
 */
template<class TTransformScalarType, int TDimension>
void writeStationaryVelocityFieldTransformation(
        itk::Transform<TTransformScalarType, TDimension, TDimension> * field,
        std::string fileName );


/**
 * Image interpolator used for image resampling.
 */
enum ImageInterpolatorType{
    INTERPOLATOR_NEAREST_NEIGHBOR,  /** nearest neighbor */
    INTERPOLATOR_LINEAR,            /** linear           */
    INTERPOLATOR_BSLPINE,           /** b-spline         */
    INTERPOLATOR_SINUS_CARDINAL     /** sinus cardinal   */
};

/**
 * Return the input image interpolation type as string.
 * @param  type image interpolation type
 * @return string describing the image interpolation type
 */
inline std::string getImageInterpolatorTypeAsString(ImageInterpolatorType interpolator);


/**
 * Resamples an image given a transformation and a geometry and saves the resampled image into an output file.
 * @param  transform     transformation
 * @param  image         image to resample
 * @param  origin        origin of the ouptut image
 * @param  spacing       voxel size of the ouptut image
 * @param  size          size of the ouptut image
 * @param  direction     orientation of the ouptut image
 * @param  fileName      name of the output file
 * @param  interpolator  type of the image interpolation
 */
template<class TImage, class TTransformScalarType>
void resampleAndWriteImage(
    itk::Transform<TTransformScalarType, TImage::ImageDimension, TImage::ImageDimension> * transform,
    TImage * image,
    const typename TImage::PointType     & origin,
    const typename TImage::SpacingType   & spacing,
    const typename TImage::SizeType      & size,
    const typename TImage::DirectionType & direction,
    std::string fileName,
    ImageInterpolatorType interpolator = INTERPOLATOR_LINEAR
    );


/**
 * Resamples the moving image in the geometry of the fixed image using a given a transformation.
 * Then saves the resampled image into an output file.
 * @param  transform    transformation
 * @param  fixedImage   fixed image
 * @param  movingImage  moving image
 * @param  fileName     name of the output file
 * @param  interpolator type of the image interpolation
 */
template<class TFixedImage, class TMovingImage, class TTransformScalarType>
void resampleAndWriteImage(
    itk::Transform<TTransformScalarType, TFixedImage::ImageDimension, TFixedImage::ImageDimension> * transform,
    TFixedImage  * fixedImage,
    TMovingImage * movingImage,
    std::string fileName,
    ImageInterpolatorType interpolator = INTERPOLATOR_LINEAR
    );


/**
 * Resamples a tensor image given a transformation and a geometry and saves the resampled tensor image into an output file.
 * @param  transform     transformation
 * @param  tensor        tensor image to resample
 * @param  origin        origin of the ouptut tensor image
 * @param  spacing       voxel size of the ouptut tensor image
 * @param  size          size of the ouptut tensor image
 * @param  direction     orientation of the ouptut tensor image
 * @param  fileName      name of the output file
 */
template<class TTensor, class TTransformScalarType>
void resampleAndWriteTensor(
    itk::Transform<TTransformScalarType, TTensor::ImageDimension, TTensor::ImageDimension> * transform,
    TTensor * tensor,
    const typename TTensor::PointType     & origin,
    const typename TTensor::SpacingType   & spacing,
    const typename TTensor::SizeType      & size,
    const typename TTensor::DirectionType & direction,
    std::string fileName
    );


/**
 * Resamples the moving tensor image in the geometry of the fixed tensor image using a given a transformation.
 * Then saves the resampled tensor image into an output file.
 * @param  transform     transformation
 * @param  fixedTensor   fixed tensor image
 * @param  movingTensor  moving tensor image
 * @param  fileName      name of the output file
 */
template<class TFixedTensor, class TMovingTensor, class TTransformScalarType>
void resampleAndWriteTensor(
    itk::Transform<TTransformScalarType, TFixedTensor::ImageDimension, TFixedTensor::ImageDimension> * transform,
    TFixedTensor  * fixedTensor,
    TMovingTensor * movingTensor,
    std::string fileName
    );


/**
 * Gets the geometry of an image from its header. The geometry is writen into the input references
 * (origin, spacing, size, direction). The function has one template NDimension representing the
 * dimension of the geometry-related inputs.
 * @param  fileName   input image
 * @param  origin     image origin
 * @param  spacing    voxel size
 * @param  size       image size
 * @param  direction  image direction
 */
template <unsigned int NDimension>
void getGeometryFromImageHeader(
    std::string fileName,
    typename itk::ImageBase<NDimension>::PointType     & origin,
    typename itk::ImageBase<NDimension>::SpacingType   & spacing,
    typename itk::ImageBase<NDimension>::SizeType      & size,
    typename itk::ImageBase<NDimension>::DirectionType & direction
    );


/**
 * Gets the geometry from an input image. The geometry is writen into the input references
 * (origin, spacing, size, direction). The function has one template TImageType. The type of the
 * geometry-related inputs is deduced from this template.
 * @param  image      input image
 * @param  origin     image origin
 * @param  spacing    voxel size
 * @param  size       image size
 * @param  direction  image direction
 */
template <class TImageType>
void getGeometryFromImage(
    const TImageType * image,
    typename TImageType::PointType     & origin,
    typename TImageType::SpacingType   & spacing,
    typename TImageType::SizeType      & size,
    typename TImageType::DirectionType & direction
    );


/**
 * Parses a string and generate a vector of unsigned int.
 * @param   str  string to parse
 * @return  vector of unsigned int
 */
template<typename T>
std::vector<T> StringToVector( const std::string & str );


/**
 * Displays a std::vector as a string.
 * @param   vector  the std::vector to be displayed
 * @return  std::vector as a string
 */
template<typename T>
std::string VectorToString( std::vector<T> vector );


/**
 * Displays a boolean as a string.
 * @param   value  boolean value
 * @return  boolean as a string
 */
inline std::string BooleanToString( bool value );


} // End of namespace


#include "rpiCommonTools.cxx"


#endif // _RPI_COMMON_TOOLS_HXX_
