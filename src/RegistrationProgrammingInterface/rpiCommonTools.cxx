#ifndef _RPI_COMMON_TOOLS_CXX_
#define _RPI_COMMON_TOOLS_CXX_

#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <cerrno>
#include <climits>

#include <itkImage.h>
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkTransform.h>
#include <itkTransformFactory.h>
#include <itkMatrixOffsetTransformBase.h>
#include <itkAffineTransform.h>
#include <itkEuler2DTransform.h>
#include <itkEuler3DTransform.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkTransformToDisplacementFieldSource.h>
#include <itkTransformToVelocityFieldSource.h>
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <itkConstantBoundaryCondition.h>

#include "rpiCommonTools.hxx"



// Namespace RPI : Registration Programming Interface
namespace rpi
{


//--------------------------------------------------------------------------------------------------
//
// Local functions
//
//--------------------------------------------------------------------------------------------------


/**
 * Reads a transform from a text file and returns a TransformBase::Pointer object.
 * @param  fileName path to the input transformation file
 * @return transformation read as TransformBase::Pointer object
 */
inline
itk::TransformBase::Pointer
readTransformBase( std::string fileName )
{

    // Type definition
    typedef  itk::TransformFileReader                TransformReaderType;
    typedef  TransformReaderType::TransformListType  TransformListType;

    // Define transformation reader
    TransformReaderType::Pointer reader = TransformReaderType::New();
    reader->SetFileName( fileName );

    // Update the reader
    try
    {
        reader->Update();
    }
    catch( itk::ExceptionObject& )
    {
        throw std::runtime_error( "Could not read the input transformation." );
    }

    // Get the list of transformations
    TransformListType * list = reader->GetTransformList();
    if ( list->empty() )
        throw std::runtime_error( "The input file does not contain any transformation." );
    else if ( list->size()>1 )
        throw std::runtime_error( "The input file contains more than one transformation." );

    // Return the pointer on the transformation
    return list->front();

}



/**
 * Casts the scalar type of an input Euler3D transformation.
 * @param  input input transformation
 * @return transformation casted
 */
template<class TInputScalarType, class TOutputScalarType>
typename itk::Euler3DTransform<TOutputScalarType>::Pointer
castEuler3DTransform(itk::Euler3DTransform<TInputScalarType> * input)
{
    // Typedef for input transformation
    typedef  itk::Euler3DTransform<TInputScalarType>       InputTransformType;
    typedef  typename InputTransformType::ParametersType   InputParameterType;
    typedef  typename InputTransformType::InputPointType   InputCenterType;

    // Typedef for output transformation
    typedef  itk::Euler3DTransform<TOutputScalarType>      OutputTransformType;
    typedef  typename OutputTransformType::ParametersType  OutputParameterType;
    typedef  typename OutputTransformType::InputPointType  OutputCenterType;

    // Get input parameters and center
    InputParameterType iP = input->GetParameters();
    InputCenterType    iC = input->GetCenter();

    // Create output parameters and center
    OutputParameterType oP(InputTransformType::ParametersDimension);
    OutputCenterType    oC;
    for (int i=0; i<OutputTransformType::ParametersDimension; i++)
        oP[i] = iP[i];
    for (int i=0; i<OutputTransformType::OutputSpaceDimension; i++)
        oC[i] = iC[i];

    // Create the output transformation
    typename OutputTransformType::Pointer output = OutputTransformType::New();
    output->Register();
    output->SetParameters( oP );
    output->SetCenter(     oC );

    return output;
}



/**
 * Casts the scalar type of an input affine 3D transformation.
 * @param  input input transformation
 * @return transformation casted
 */
template<class TInputScalarType, class TOutputScalarType>
typename itk::AffineTransform<TOutputScalarType,3>::Pointer
castAffineTransform(itk::AffineTransform<TInputScalarType,3> * input)
{
    // Typedef for input transformation
    typedef  itk::AffineTransform<TInputScalarType,3>      InputTransformType;
    typedef  typename InputTransformType::ParametersType   InputParameterType;
    typedef  typename InputTransformType::InputPointType   InputCenterType;

    // Typedef for output transformation
    typedef  itk::AffineTransform<TOutputScalarType,3>     OutputTransformType;
    typedef  typename OutputTransformType::ParametersType  OutputParameterType;
    typedef  typename OutputTransformType::InputPointType  OutputCenterType;

    // Get input parameters and center
    InputParameterType iP = input->GetParameters();
    InputCenterType    iC = input->GetCenter();

    // Create output parameters and center
    OutputParameterType oP(InputTransformType::ParametersDimension);
    OutputCenterType    oC;
    for (int i=0; i<OutputTransformType::ParametersDimension; i++)
        oP[i] = iP[i];
    for (int i=0; i<OutputTransformType::OutputSpaceDimension; i++)
        oC[i] = iC[i];

    // Create the output transformation
    typename OutputTransformType::Pointer output = OutputTransformType::New();
    output->SetParameters( oP );
    output->SetCenter(     oC );

    return output;
}



/**
 * Casts the scalar type of an input matrix offset 3D transformation.
 * @param  input input transformation
 * @return transformation casted
 */
template<class TInputScalarType, class TOutputScalarType>
typename itk::MatrixOffsetTransformBase<TOutputScalarType,3,3>::Pointer
castMatrixOffsetTransform(itk::MatrixOffsetTransformBase<TInputScalarType,3,3> * input)
{
    // Typedef for input transformation
    typedef  itk::MatrixOffsetTransformBase<TInputScalarType,3,3>   InputTransformType;
    typedef  typename InputTransformType::ParametersType            InputParameterType;
    typedef  typename InputTransformType::InputPointType            InputCenterType;

    // Typedef for output transformation
    typedef  itk::MatrixOffsetTransformBase<TOutputScalarType,3,3>  OutputTransformType;
    typedef  typename OutputTransformType::ParametersType           OutputParameterType;
    typedef  typename OutputTransformType::InputPointType           OutputCenterType;

    // Get input parameters and center
    InputParameterType iP = input->GetParameters();
    InputCenterType    iC = input->GetCenter();

    // Create output parameters and center
    OutputParameterType oP(InputTransformType::ParametersDimension);
    OutputCenterType    oC;
    for (int i=0; i<OutputTransformType::ParametersDimension; i++)
        oP[i] = iP[i];
    for (int i=0; i<OutputTransformType::OutputSpaceDimension; i++)
        oC[i] = iC[i];

    // Create the output transformation
    typename OutputTransformType::Pointer output = OutputTransformType::New();
    output->SetParameters( oP );
    output->SetCenter(     oC );

    return output;
}



/**
 * Convert an Euler3D transformation into an affine 3D transformation
 * @param  euler euler transformation
 * @return affine transformation
 */
template<class TEuler3DScalarType, class TAffineScalarType>
typename itk::AffineTransform<TAffineScalarType,3>::Pointer
convertEuler3DToAffine(itk::Euler3DTransform<TEuler3DScalarType> * euler)
{
    // Typedef for euler3D transformation
    typedef  itk::Euler3DTransform<TEuler3DScalarType>  Euler3DType;
    typedef  typename Euler3DType::MatrixType           Euler3DMatrixType;
    typedef  typename Euler3DType::OutputVectorType     Euler3DTranslationType;
    typedef  typename Euler3DType::InputPointType       Euler3DCenterType;

    // Typedef for affine transformation
    typedef  itk::AffineTransform<TAffineScalarType,3>  AffineType;
    typedef  typename AffineType::MatrixType            AffineMatrixType;
    typedef  typename AffineType::OutputVectorType      AffineTranslationType;
    typedef  typename AffineType::InputPointType        AffineCenterType;

    // Get input matrix, translation, and center
    Euler3DMatrixType      iM = euler->GetMatrix();
    Euler3DTranslationType iT = euler->GetTranslation();
    Euler3DCenterType      iC = euler->GetCenter();

    // Create output matrix, translation, and center
    AffineMatrixType       oM;
    AffineTranslationType  oT;
    AffineCenterType       oC;
    for (int i=0; i<AffineType::OutputSpaceDimension; i++){
        oT[i] = iT[i];
        oC[i] = iC[i];
        for (int j=0; j<AffineType::OutputSpaceDimension; j++)
            oM(i,j) = iM(i,j);
    }

    // Create the output transformation
    typename AffineType::Pointer affine = AffineType::New();
    affine->SetMatrix(      oM );
    affine->SetTranslation( oT );
    affine->SetCenter(      oC );

    return affine;
}



/**
 * Convert a MatroixOffset transformation into an affine 3D transformation
 * @param  matrixOffset matrixOffset transformation
 * @return affine transformation
 */
template<class TMatrixOffsetScalarType, class TAffineScalarType>
typename itk::AffineTransform<TAffineScalarType,3>::Pointer
convertMatrixOffsetToAffine(itk::MatrixOffsetTransformBase<TMatrixOffsetScalarType,3,3> * input)
{
    // Typedef for euler3D transformation
    typedef  itk::MatrixOffsetTransformBase<TMatrixOffsetScalarType>  InputTransformType;
    typedef  typename InputTransformType::ParametersType              InputParameterType;
    typedef  typename InputTransformType::InputPointType              InputCenterType;

    // Typedef for output transformation
    typedef  itk::AffineTransform<TAffineScalarType,3>                OutputTransformType;
    typedef  typename OutputTransformType::ParametersType             OutputParameterType;
    typedef  typename OutputTransformType::InputPointType             OutputCenterType;

    // Get input parameters and center
    InputParameterType iP = input->GetParameters();
    InputCenterType    iC = input->GetCenter();

    // Create output parameters and center
    OutputParameterType oP(InputTransformType::ParametersDimension);
    OutputCenterType    oC;
    for (int i=0; i<OutputTransformType::ParametersDimension; i++)
        oP[i] = iP[i];
    for (int i=0; i<OutputTransformType::OutputSpaceDimension; i++)
        oC[i] = iC[i];

    // Create the output transformation
    typename OutputTransformType::Pointer output = OutputTransformType::New();
    output->SetParameters( oP );
    output->SetCenter(     oC );

    return output;
}



//--------------------------------------------------------------------------------------------------
//
// Functions of rpi::CommonTools
//
//--------------------------------------------------------------------------------------------------

inline
itk::ImageIOBase::Pointer
readImageInformation( std::string fileName )
{

    // Define image IO
    itk::ImageIOBase::Pointer imageIO  = itk::ImageIOFactory::CreateImageIO( fileName.c_str(), itk::ImageIOFactory::ReadMode );

    // Test if image exists
    if ( !imageIO )
        throw std::runtime_error( "Could not read image " + fileName + "." );

    // Read image information
    try
    {
        imageIO->SetFileName( fileName );
        imageIO->ReadImageInformation();
    }
    catch( itk::ExceptionObject& )
    {
        throw std::runtime_error( "Could not read image information." );
    }
    return imageIO;
}



template<typename TImage>
typename TImage::Pointer readImage( std::string fileName )
{
    typedef itk::ImageFileReader<TImage>  ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName( fileName );
    try
    {
        reader->Update();
    }
    catch( itk::ExceptionObject& )
    {
        throw std::runtime_error( "Could not read the input image." );
    }
    return reader->GetOutput();
}



void registerUsefulITKTransformations(void)
{
    // Type definition
    typedef  itk::AffineTransform< double, 3 >            AffineDouble;
    typedef  itk::AffineTransform< float, 3 >             AffineFloat;
    typedef  itk::MatrixOffsetTransformBase< double, 3 >  MatrixOffsetDouble;
    typedef  itk::MatrixOffsetTransformBase< float, 3 >   MatrixOffsetFloat;
    typedef  itk::Euler3DTransform<float>                 Euler3DFloat;
    typedef  itk::Euler3DTransform<double>                Euler3DDouble;

    // Register useful transformations
    itk::TransformFactory< AffineDouble >::RegisterTransform ();
    itk::TransformFactory< AffineFloat >::RegisterTransform ();
    itk::TransformFactory< MatrixOffsetDouble >::RegisterTransform ();
    itk::TransformFactory< MatrixOffsetFloat >::RegisterTransform ();
    itk::TransformFactory< Euler3DDouble >::RegisterTransform ();
    itk::TransformFactory< Euler3DFloat >::RegisterTransform ();
}



template<class TTransformScalarType>
typename itk::Euler3DTransform<TTransformScalarType>::Pointer
readEuler3DTransformation( std::string fileName )
{

    // Type definition
    typedef  itk::TransformBase             TransformBaseType;
    typedef  itk::Euler3DTransform<float>   Euler3DFloat;
    typedef  itk::Euler3DTransform<double>  Euler3DDouble;

    // Read the TransformBase object
    TransformBaseType * transform = readTransformBase(fileName).GetPointer();

    // Try to cast as Euler3DFloat
    Euler3DFloat *euler3DFloat = dynamic_cast<Euler3DFloat*>(transform);
    if(euler3DFloat != 0)
        return castEuler3DTransform<float,TTransformScalarType>(euler3DFloat);

    // Try to cast as Euler3DDouble
    Euler3DDouble *euler3DDouble = dynamic_cast<Euler3DDouble*>(transform);
    if(euler3DDouble != 0)
        return castEuler3DTransform<double,TTransformScalarType>(euler3DDouble);

    // Otherwise, throw an exception
    throw std::runtime_error( "Transformation type not supported." );
    return 0;
}



template<class TTransformScalarType>
typename itk::AffineTransform<TTransformScalarType, 3>::Pointer
readAffineTransformation( std::string fileName )
{

    // Type definition
    typedef  itk::TransformBase                             TransformBaseType;
    typedef  itk::AffineTransform<float, 3>                 AffineFloat;
    typedef  itk::AffineTransform<double, 3>                AffineDouble;
    typedef  itk::Euler3DTransform<float>                   Euler3DFloat;
    typedef  itk::Euler3DTransform<double>                  Euler3DDouble;
    typedef  itk::MatrixOffsetTransformBase<float,  3 ,3 >  MatrixOffsetFloat;
    typedef  itk::MatrixOffsetTransformBase<double, 3 ,3 >  MatrixOffsetDouble;

    // Read the TransformBase object
    TransformBaseType * transform = readTransformBase(fileName).GetPointer();

    // Try to cast as AffineFloat
    AffineFloat *affineFloat = dynamic_cast<AffineFloat*>(transform);
    if(affineFloat != 0)
        return castAffineTransform<float,TTransformScalarType>(affineFloat);

    // Try to cast as AffineDouble
    AffineDouble *affineDouble = dynamic_cast<AffineDouble*>(transform);
    if(affineDouble != 0)
        return castAffineTransform<double,TTransformScalarType>(affineDouble);

    // Try to cast as Euler3DFloat
    Euler3DFloat *euler3DFloat = dynamic_cast<Euler3DFloat*>(transform);
    if(euler3DFloat != 0)
        return convertEuler3DToAffine<float,TTransformScalarType>(euler3DFloat);

    // Try to cast as Euler3DDouble
    Euler3DDouble *euler3DDouble = dynamic_cast<Euler3DDouble*>(transform);
    if(euler3DDouble != 0)
        return convertEuler3DToAffine<double,TTransformScalarType>(euler3DDouble);

    // Try to cast as MatrixOffsetFloat
    MatrixOffsetFloat *matrixOffsetFloat = dynamic_cast<MatrixOffsetFloat*>(transform);
    if(matrixOffsetFloat != 0)
        return convertMatrixOffsetToAffine<float,TTransformScalarType>(matrixOffsetFloat);

    // Try to cast as MatrixOffsetDouble
    MatrixOffsetDouble *matrixOffsetDouble = dynamic_cast<MatrixOffsetDouble*>(transform);
    if(matrixOffsetDouble != 0)
        return convertMatrixOffsetToAffine<double,TTransformScalarType>(matrixOffsetDouble);

    // Otherwise, throw an exception
    throw std::runtime_error( "Transformation type not supported." );
    return 0;
}



template<class TTransformScalarType>
typename itk::Transform<TTransformScalarType, 3, 3>::Pointer
readLinearTransformation( std::string fileName )
{

    // Type definition
    typedef  itk::TransformBase                            TransformBaseType;
    typedef  itk::Euler3DTransform<float>                  Euler3DFloat;
    typedef  itk::Euler3DTransform<double>                 Euler3DDouble;
    typedef  itk::AffineTransform<float, 3>                AffineFloat;
    typedef  itk::AffineTransform<double, 3>               AffineDouble;
    typedef  itk::MatrixOffsetTransformBase<float, 3, 3>   MatrixOffsetFloat;
    typedef  itk::MatrixOffsetTransformBase<double, 3, 3>  MatrixOffsetDouble;

    // Read the TransformBase object
    TransformBaseType * transform = readTransformBase(fileName).GetPointer();

    // Note : In the following, we return the pointer (using GetPointer) instead of the itk::Pointer
    // to allow the automatic cast of the returned transform into a itk::Transform object.

    // Try to cast as Euler3DFloat
    Euler3DFloat *euler3DFloat = dynamic_cast<Euler3DFloat*>(transform);
    if(euler3DFloat != 0)
        return castEuler3DTransform<float,TTransformScalarType>(euler3DFloat).GetPointer();

    // Try to cast as Euler3DDouble
    Euler3DDouble *euler3DDouble = dynamic_cast<Euler3DDouble*>(transform);
    if(euler3DDouble != 0)
        return castEuler3DTransform<double,TTransformScalarType>(euler3DDouble).GetPointer();

    // Try to cast as AffineFloat
    AffineFloat *affineFloat = dynamic_cast<AffineFloat*>(transform);
    if(affineFloat != 0)
        return castAffineTransform<float,TTransformScalarType>(affineFloat).GetPointer();

    // Try to cast as AffineDouble
    AffineDouble *affineDouble = dynamic_cast<AffineDouble*>(transform);
    if(affineDouble != 0)
        return castAffineTransform<double,TTransformScalarType>(affineDouble).GetPointer();

    // Try to cast as MatrixOffsetFloat
    MatrixOffsetFloat *matrixOffsetFloat = dynamic_cast<MatrixOffsetFloat*>(transform);
    if(matrixOffsetFloat != 0)
        return castMatrixOffsetTransform<float,TTransformScalarType>(matrixOffsetFloat).GetPointer();

    // Try to cast as MatrixOffsetDouble
    MatrixOffsetDouble *matrixOffsetDouble = dynamic_cast<MatrixOffsetDouble*>(transform);
    if(matrixOffsetDouble != 0)
        return castMatrixOffsetTransform<double,TTransformScalarType>(matrixOffsetDouble).GetPointer();

    // Otherwise, throw an exception
    throw std::runtime_error( "Transformation type not supported." );
    return 0;
}



template<class TTransformScalarType>
typename itk::DisplacementFieldTransform<TTransformScalarType, 3>::Pointer
readDisplacementField( std::string fileName )
{

    typedef itk::DisplacementFieldTransform<TTransformScalarType, 3>
            FieldTransformType;

    typedef typename FieldTransformType::DisplacementFieldType
            DisplacementFieldType;

    typedef itk::ImageFileReader<DisplacementFieldType>
            ReaderType;

    // Read the displacement field (image)
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( fileName );
    try
    {
        reader->Update();
    }
    catch( itk::ExceptionObject& err )
    {
        throw std::runtime_error( "Could not read the input transformation." );
    }
    typename DisplacementFieldType::Pointer field = reader->GetOutput();

    // Create the displacement field (transformation)
    typename FieldTransformType::Pointer transform = FieldTransformType::New();
    transform->SetDisplacementField( static_cast<typename DisplacementFieldType::Pointer>( field.GetPointer() ) );
    return transform;
}



template<class TTransformScalarType>
typename itk::StationaryVelocityFieldTransform<TTransformScalarType, 3>::Pointer
readStationaryVelocityField( std::string fileName )
{

    typedef itk::StationaryVelocityFieldTransform<TTransformScalarType, 3>
            FieldTransformType;

    typedef typename FieldTransformType::DisplacementFieldType
            DisplacementFieldType;

    typedef itk::ImageFileReader<DisplacementFieldType>
            ReaderType;

    // Read the velocity field (image)
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( fileName );
    try
    {
        reader->Update();
    }
    catch( itk::ExceptionObject& err )
    {
        throw std::runtime_error( "Could not read the input transformation." );
    }
    typename DisplacementFieldType::Pointer field = reader->GetOutput();

    // Create the velocity field (transformation)
    typename FieldTransformType::Pointer transform = FieldTransformType::New();
    transform->SetParametersAsVectorField( static_cast<typename DisplacementFieldType::ConstPointer>( field.GetPointer() ) );
    return transform;
}



template<class TLinearScalarType, class TFieldScalarType, class TImage>
typename itk::DisplacementFieldTransform<TFieldScalarType, TImage::ImageDimension>::Pointer
linearToDisplacementFieldTransformation(
        TImage * image,
        itk::Transform<TLinearScalarType, TImage::ImageDimension, TImage::ImageDimension> * transform )
{
    typedef  itk::DisplacementFieldTransform<TFieldScalarType, TImage::ImageDimension>
            FieldTransformType;

    typedef  typename FieldTransformType::DisplacementFieldType
            DisplacementFieldType;

    typedef  itk::TransformToDisplacementFieldSource<DisplacementFieldType, TLinearScalarType>
            GeneratorType;

    // Create a field generator
    typename GeneratorType::Pointer fieldGenerator = GeneratorType::New();
    fieldGenerator->SetTransform(                 transform );
    fieldGenerator->SetOutputParametersFromImage( image );

    // Update the field generator
    try
    {
        fieldGenerator->Update();
    }
    catch( itk::ExceptionObject& err )
    {
        throw std::runtime_error( "Could not generate a displacement field from a linear transformation." );
    }

    // Gets the displacement field (image)
    typename DisplacementFieldType::Pointer container = fieldGenerator->GetOutput();
    container->DisconnectPipeline();

    // Create the displacement field (transformation)
    typename FieldTransformType::Pointer field = FieldTransformType::New();
    field->SetDisplacementField( static_cast<typename DisplacementFieldType::Pointer>( container.GetPointer() ) );
    return field;
}



template<class TLinearScalarType, class TFieldScalarType, class TImage>
typename itk::StationaryVelocityFieldTransform<TFieldScalarType, TImage::ImageDimension>::Pointer
linearToStationaryVelocityFieldTransformation(
        TImage * image,
        itk::Transform<TLinearScalarType, TImage::ImageDimension, TImage::ImageDimension> * transform )
{
    typedef  itk::StationaryVelocityFieldTransform<TFieldScalarType, TImage::ImageDimension>
            FieldTransformType;

    typedef  typename FieldTransformType::DisplacementFieldType
            DisplacementFieldType;

    typedef  itk::TransformToVelocityFieldSource<DisplacementFieldType, TLinearScalarType>
            GeneratorType;

    // Create a field generator
    typename GeneratorType::Pointer fieldGenerator = GeneratorType::New();
    fieldGenerator->SetTransform(                 transform );
    fieldGenerator->SetOutputParametersFromImage( image );

    // Update the field generator
    try
    {
        fieldGenerator->Update();
    }
    catch( itk::ExceptionObject& err )
    {
        throw std::runtime_error( "Could not generate a velocity field from a linear transformation." );
    }

    // Gets the displacement field (image)
    typename DisplacementFieldType::Pointer container = fieldGenerator->GetOutput();
    container->DisconnectPipeline();

    // Create the displacement field (transformation)
    typename FieldTransformType::Pointer field = FieldTransformType::New();
    field->SetParametersAsVectorField( static_cast<typename DisplacementFieldType::ConstPointer>( container.GetPointer() ) );
    return field;
}



template<class TTransformScalarType, int TDimension>
void writeLinearTransformation(
        itk::Transform<TTransformScalarType, TDimension, TDimension> * transform,
        std::string fileName )
{
    typedef itk::TransformFileWriter TrasformWriterType;
    typename TrasformWriterType::Pointer transformWriter = TrasformWriterType::New();
    transformWriter->SetFileName( fileName );
    transformWriter->SetInput( transform );
    transformWriter->Update();
}



template<class TTransformScalarType, int TDimension>
void writeDisplacementFieldTransformation(
        itk::Transform<TTransformScalarType, TDimension, TDimension> * field,
        std::string fileName )
{
    // Type definition
    typedef itk::DisplacementFieldTransform<TTransformScalarType, TDimension>  FieldTransformType;
    typedef typename FieldTransformType::DisplacementFieldType                       DisplacementFieldType;
    typedef itk::ImageFileWriter<DisplacementFieldType>                              FieldWriterType;

    // Cast the transformation into a displacement field
    FieldTransformType * transform = dynamic_cast<FieldTransformType *>(field);

    // Write the output field
    typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
    fieldWriter->SetInput(    transform->GetDisplacementField() );
    fieldWriter->SetFileName( fileName );
    fieldWriter->Update();
}



template<class TTransformScalarType, int TDimension>
void writeStationaryVelocityFieldTransformation(
        itk::Transform<TTransformScalarType, TDimension, TDimension> * field,
        std::string fileName )
{
    // Type definition
    typedef itk::StationaryVelocityFieldTransform<TTransformScalarType, TDimension>  FieldTransformType;
    typedef typename FieldTransformType::VectorFieldType                             VectorFieldType;
    typedef itk::ImageFileWriter<VectorFieldType>                              FieldWriterType;

    // Cast the transformation into a displacement field
    FieldTransformType * transform = dynamic_cast<FieldTransformType *>(field);

    // Write the output field
    typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
    fieldWriter->SetInput(    transform->GetParametersAsVectorField() );
    fieldWriter->SetFileName( fileName );
    fieldWriter->Update();
}



inline std::string getImageInterpolatorTypeAsString(ImageInterpolatorType interpolator)
{
    switch(interpolator) {

    case INTERPOLATOR_NEAREST_NEIGHBOR:
        return "nearest neighbor";

    case INTERPOLATOR_LINEAR:
        return "linear";

    case INTERPOLATOR_BSLPINE:
        return "b-spline";

    case INTERPOLATOR_SINUS_CARDINAL:
        return "sinus cardinal";

    default:
        throw std::runtime_error("Image interpolator not supported.");
    }
}



template<class TImage, class TTransformScalarType>
void resampleAndWriteImage(
    itk::Transform<TTransformScalarType, TImage::ImageDimension, TImage::ImageDimension> * transform,
    TImage * image,
    const typename TImage::PointType     & origin,
    const typename TImage::SpacingType   & spacing,
    const typename TImage::SizeType      & size,
    const typename TImage::DirectionType & direction,
    std::string fileName,
    ImageInterpolatorType interpolator )
{

    // Create and initialize the resample image filter
    typedef itk::ResampleImageFilter<TImage, TImage, TTransformScalarType> ResampleFilterType;
    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput(             image );
    resampler->SetTransform(         transform );
    resampler->SetSize(              size );
    resampler->SetOutputOrigin(      origin );
    resampler->SetOutputSpacing(     spacing );
    resampler->SetOutputDirection(   direction );
    resampler->SetDefaultPixelValue( 0 );

    // Set the image interpolator
    switch(interpolator) {

    case INTERPOLATOR_NEAREST_NEIGHBOR:
        resampler->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<TImage, TTransformScalarType>::New());
        break;

    case INTERPOLATOR_LINEAR:
        // Nothing to do ; linear interpolator by default
        break;

    case INTERPOLATOR_BSLPINE:
        resampler->SetInterpolator(itk::BSplineInterpolateImageFunction<TImage, TTransformScalarType>::New());
        break;

    case INTERPOLATOR_SINUS_CARDINAL:
        resampler->SetInterpolator(itk::WindowedSincInterpolateImageFunction<
                                   TImage,
                                   TImage::ImageDimension,
                                   itk::Function::HammingWindowFunction<TImage::ImageDimension>,
                                   itk::ConstantBoundaryCondition<TImage>,
                                   TTransformScalarType
                                   >::New());
        break;

    default:
        throw std::runtime_error("Image interpolator not supported.");
        break;
    }

    // Resample and write the output image
    typedef itk::ImageFileWriter<TImage> ImageWriterType;
    typename ImageWriterType::Pointer imageWriter = ImageWriterType::New();
    imageWriter->SetFileName( fileName );
    imageWriter->SetInput( resampler->GetOutput() );
    imageWriter->Update();
}



template<class TFixedImage, class TMovingImage, class TTransformScalarType>
void
resampleAndWriteImage(
    itk::Transform<TTransformScalarType, TFixedImage::ImageDimension, TFixedImage::ImageDimension> * transform,
    TFixedImage  * fixedImage,
    TMovingImage * movingImage,
    std::string fileName,
    ImageInterpolatorType interpolator )
{
    resampleAndWriteImage<TMovingImage,TTransformScalarType>(
                transform,
                movingImage,
                fixedImage->GetOrigin(),
                fixedImage->GetSpacing(),
                fixedImage->GetLargestPossibleRegion().GetSize(),
                fixedImage->GetDirection(),
                fileName,
                interpolator );
}



template<class TTensor, class TTransformScalarType>
void resampleAndWriteTensor(
    itk::Transform<TTransformScalarType, TTensor::ImageDimension, TTensor::ImageDimension> * transform,
    TTensor * tensor,
    const typename TTensor::PointType     & origin,
    const typename TTensor::SpacingType   & spacing,
    const typename TTensor::SizeType      & size,
    const typename TTensor::DirectionType & direction,
    std::string fileName)
{
    // Type definition
    typedef  itk::DisplacementFieldTransform<TTransformScalarType, TTensor::ImageDimension>
            DFTransformType;
    typedef  itk::StationaryVelocityFieldTransform<TTransformScalarType, TTensor::ImageDimension>
            SVFTransformType;

    // If the transformation is linear, generate a displacement field and use recursion
    if (transform->IsLinear())
    {
        typename DFTransformType::Pointer df = linearToDisplacementFieldTransformation<TTransformScalarType, TTransformScalarType, TTensor>(tensor,transform);
        resampleAndWriteTensor<TTensor,TTransformScalarType>(df, tensor, origin, spacing, size, direction, fileName);
        return;
    }

    // If the transformation is a stationary velocity field, generate a displacement field and use recursion
    SVFTransformType* svf = dynamic_cast<SVFTransformType*>(transform);
    if (svf!=0)
    {
        typename DFTransformType::Pointer df = DFTransformType::New();
        df->SetParametersAsVectorField(svf->GetDisplacementFieldAsVectorField());
        resampleAndWriteTensor<TTensor,TTransformScalarType>(df, tensor, origin, spacing, size, direction, fileName);
        return;
    }

    // If the transformation is a displacement field, warp the displacement field
    DFTransformType* df = dynamic_cast<DFTransformType*>(transform);
    if (df!=0)
    {
        std::cout << "Tensor resampling not implemented yet!" << std::endl;
        // WarpTensorImageFilter
        return;
    }
}



template<class TFixedTensor, class TMovingTensor, class TTransformScalarType>
void resampleAndWriteTensor(
    itk::Transform<TTransformScalarType, TFixedTensor::ImageDimension, TFixedTensor::ImageDimension> * transform,
    TFixedTensor  * fixedTensor,
    TMovingTensor * movingTensor,
    std::string fileName)
{
    resampleAndWriteTensor<TMovingTensor,TTransformScalarType>(
                transform,
                movingTensor,
                fixedTensor->GetOrigin(),
                fixedTensor->GetSpacing(),
                fixedTensor->GetLargestPossibleRegion().GetSize(),
                fixedTensor->GetDirection(),
                fileName);
}



template <unsigned int NDimension>
void
getGeometryFromImageHeader(
    std::string fileName,
    typename itk::ImageBase<NDimension>::PointType     & origin,
    typename itk::ImageBase<NDimension>::SpacingType   & spacing,
    typename itk::ImageBase<NDimension>::SizeType      & size,
    typename itk::ImageBase<NDimension>::DirectionType & direction
    )
{
    // Read image header
    typename itk::ImageIOBase::Pointer image = readImageInformation( fileName );
    if (image->GetNumberOfDimensions()!=NDimension)
        throw std::runtime_error("Input image dimension does not fit geometry dimension.");

    // Fill geometry from image header
    for (unsigned int i=0; i<NDimension; i++)
    {
        origin[i]  = image->GetOrigin(i);
        spacing[i] = image->GetSpacing(i);
        size[i]    = image->GetDimensions(i);
        for (unsigned int j=0; j<NDimension; j++)
            direction(i,j) = image->GetDirection(j)[i];
    }
}



template <class TImageType>
void
getGeometryFromImage(
    const TImageType * image,
    typename TImageType::PointType     & origin,
    typename TImageType::SpacingType   & spacing,
    typename TImageType::SizeType      & size,
    typename TImageType::DirectionType & direction
    )
{

    // Input geometry
    typename TImageType::PointType     input_origin    = image->GetOrigin();
    typename TImageType::SpacingType   input_spacing   = image->GetSpacing();
    typename TImageType::SizeType      input_size      = image->GetLargestPossibleRegion().GetSize();
    typename TImageType::DirectionType input_direction = image->GetDirection();

    // Output geometry
    unsigned int dim = TImageType::ImageDimension;
    for (unsigned int i=0; i<dim; i++)
    {
        origin[i]  = input_origin[i];
        spacing[i] = input_spacing[i];
        size[i]    = input_size[i];
        for (unsigned int j=0; j<dim; j++)
            direction(i,j) = input_direction(i,j);
    }
}



inline int atoi_check( const char * str )
{
    char *endptr;
    long val= strtol(str, &endptr, 0);

    // Check for various possible errors
    if ( (errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || val>=INT_MAX || val<=INT_MIN )
    {
        std::cout<<std::endl;
        std::cout<<"Cannot parse integer. Out of bound."<<std::endl;
        exit( EXIT_FAILURE );
    }

    if (endptr == str || *endptr!='\0')
    {
        std::cout<<std::endl;
        std::cout<<"Cannot parse integer. Contains non-digits or is empty."<<std::endl;
        exit( EXIT_FAILURE );
    }

    return val;
}



template<typename T>
std::vector<T> StringToVector( const std::string & str)
{

    std::vector<T> vect;

    std::string::size_type crosspos = str.find('x',0);

    if (crosspos == std::string::npos)
    {
        // only one uint
        vect.push_back( static_cast<T>( atoi_check(str.c_str()) ));
        return vect;
    }

    // first uint
    vect.push_back( static_cast<T>(atoi_check( (str.substr(0,crosspos)).c_str() ) ) );

    while(true)
    {
        std::string::size_type crossposfrom = crosspos;
        crosspos =  str.find('x',crossposfrom+1);

        if (crosspos == std::string::npos)
        {
            vect.push_back( static_cast<T>( atoi_check( (str.substr(crossposfrom+1,str.length()-crossposfrom-1)).c_str() ) ) );
            return vect;
        }

        vect.push_back( static_cast<T>( atoi_check( (str.substr(crossposfrom+1,crosspos-crossposfrom-1)).c_str() ) ) );
    }
}



template<typename T>
std::string VectorToString( std::vector<T> vector )
{
    unsigned int size = vector.size();
    std::ostringstream oss;
    oss << "[ ";
    if (size>0)
        for (unsigned int i=0; i<size-1; i++)
            oss << vector.at(i) << ", ";
    oss << vector[size-1] << " ]";
    return oss.str();
}



inline std::string BooleanToString( bool value )
{
    return (std::string)(value?"true":"false");
}



} // End of namespace


#endif // _RPI_COMMON_TOOLS_CXX_
