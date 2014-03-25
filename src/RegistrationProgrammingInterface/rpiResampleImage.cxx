#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <itkConstantBoundaryCondition.h>

#include <itkTransform.h>
#include <itkIdentityTransform.h>
#include <itkDisplacementFieldTransform.h>
#include <itkStationaryVelocityFieldTransform.h>

#include <rpiCommonTools.hxx>

#include <tclap/CmdLine.h>

#include <mipsInrimageImageIOFactory.h>


/**
 * @description  Resample an input scalar image given an input transformation.
 * @author       Vincent Garcia
 * @date         07/06/2011
 */


/**
 * Trasnformation enumeration.
 */
enum TransformationEnum {
    IDENTITY_TRANSFORMATION,
    LINEAR_TRANSFORMATION,
    DISPLACEMENT_FIELD_TRANSFORMATION,
    STATIONARY_VELOCITY_FIELD_TRANSFORMATION
};


/**
 * Structure containing parameters.
 */
struct Param
{
    std::string                inputPath;
    std::string                outputPath;
    std::string                geometryPath;
    rpi::ImageInterpolatorType interpType;
    std::string                transformPath;
    TransformationEnum         transformType;
    bool                       verbose;
};


/**
 * Parses the command line arguments and deduces the corresponding Param structure.
 * @param  argc   number of arguments
 * @param  argv   array containing the arguments
 * @param  param  structure of parameters
 */
void parseParameters(int argc, char** argv, struct Param & param)
{

    // Program description
    std::string description = "\b\b\bDESCRIPTION\n";
    description += "Resample an input scalar image given an input transformation. ";
    description += "Transformation supported are identity, linear transformations, displacement ";
    description += "field transformations, and stationary velocity field transformations. ";
    description += "The image specified using the \"--geometry\" option is used to set the resample ";
    description += "image geometry. If the \"--geometry\" option is not used, the geometry of the ";
    description += "resampled image will be either similar to the input image geometry if the input ";
    description += "transformation is a linear transformation, or will be similar to the geometry ";
    description += "of the input field is the input transformation is a displacement field or a ";
    description += "stationary velocity field. Several image interpolation methods are proposed.";
    description += "\nAuthors : Vincent Garcia";

    // Option description
    std::string dInput   = "Path to the input image.";
    std::string dOutput  = "Path to the output resampled image.";
    std::string dGeom    = "Path to the image containing to output image geometry.";
    std::string dIT      = "If activated, this options implies to use the identity transformation (i.e. no transformation).";
    std::string dLT      = "Path to the input linear transformation.";
    std::string dDFT     = "Path to the input displacement field transformation.";
    std::string dSVFT    = "Path to the input stationary velocity field transformation.";
    std::string dInterp  = "Interpolation mode : 0 = nearest neighbor, 1 = linear, 2 = b-spline, and 3 = sinus cardinal (default 1).";
    std::string dVerbose = "Verbose mode.";

    try {

        // Define command line parser
        TCLAP::CmdLine cmd( description, ' ', "1.0", true);

        // Set options
        TCLAP::SwitchArg              aVerbose(  "",  "verbose", dVerbose, cmd, true);
        TCLAP::ValueArg<unsigned int> aInterp( "", "interpolation", dInterp, false, 1, "uint", cmd );
        TCLAP::ValueArg<std::string>  aGeom( "g", "geometry", dGeom, false, "", "string", cmd );
        TCLAP::ValueArg<std::string>  aOutput( "o", "output", dOutput, true, "", "string", cmd );
        TCLAP::ValueArg<std::string>  aInput( "i", "input", dInput, true, "", "string", cmd );

        TCLAP::SwitchArg              aIT(  "",  "identity-transform", dIT, false);
        TCLAP::ValueArg<std::string>  aLT( "l", "linear-transform", dLT, false, "", "string" );
        TCLAP::ValueArg<std::string>  aDFT( "d", "displacement-field", dDFT, false, "", "string" );
        TCLAP::ValueArg<std::string>  aSVFT( "s", "velocity-field", dSVFT, false, "", "string" );
        std::vector< TCLAP::Arg* > xorlist;
        xorlist.push_back(&aIT);
        xorlist.push_back(&aLT);
        xorlist.push_back(&aDFT);
        xorlist.push_back(&aSVFT);
        cmd.xorAdd(xorlist);

        // Parse command line
        cmd.parse( argc, argv );

        // Set parameters
        param.inputPath    = aInput.getValue();
        param.outputPath   = aOutput.getValue();
        param.geometryPath = aGeom.getValue();
        param.verbose      = aVerbose.getValue();

        // Set transformation parameter
        if (aIT.getValue())
        {
            param.transformType = IDENTITY_TRANSFORMATION;
            param.transformPath = "";
        }
        if (aLT.getValue().compare("")!=0)
        {
            param.transformType = LINEAR_TRANSFORMATION;
            param.transformPath = aLT.getValue();
        }
        else if (aDFT.getValue().compare("")!=0)
        {
            param.transformType = DISPLACEMENT_FIELD_TRANSFORMATION;
            param.transformPath = aDFT.getValue();
        }
        else if (aSVFT.getValue().compare("")!=0)
        {
            param.transformType = STATIONARY_VELOCITY_FIELD_TRANSFORMATION;
            param.transformPath = aSVFT.getValue();
        }

        // Set interpolation parameter
        unsigned int interp = aInterp.getValue();
        if      (interp==0)
            param.interpType = rpi::INTERPOLATOR_NEAREST_NEIGHBOR;
        else if (interp==1)
            param.interpType = rpi::INTERPOLATOR_LINEAR;
        else if (interp==2)
            param.interpType = rpi::INTERPOLATOR_BSLPINE;
        else if (interp==3)
            param.interpType = rpi::INTERPOLATOR_SINUS_CARDINAL;
        else
            throw std::runtime_error("Interpolation mode not recognized.");
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
        throw std::runtime_error("Unable to parse the command line arguments.");
    }
}


/**
 * Prints parameters.
 * @param param parameters
 */
void printParameters(struct Param param)
{
    if (!param.verbose)
        return;

    std::cout << "Input image    : " << param.inputPath << std::endl;
    std::cout << "Output image   : " << param.outputPath << std::endl;

    if (param.geometryPath.compare(""))
        std::cout << "Geometry       : " << param.geometryPath << std::endl;

    std::cout << "Transformation : " << param.transformPath << std::endl;

    std::cout << "Interpolation  : ";
    if (param.interpType == rpi::INTERPOLATOR_NEAREST_NEIGHBOR)
        std::cout << "nearest neighbor" << std::endl;
    else if (param.interpType == rpi::INTERPOLATOR_LINEAR)
        std::cout << "linear" << std::endl;
    else if (param.interpType == rpi::INTERPOLATOR_BSLPINE)
        std::cout << "b-spline" << std::endl;
    else if (param.interpType == rpi::INTERPOLATOR_SINUS_CARDINAL)
        std::cout << "sinus cardinal" << std::endl;
}


/**
 * Reads the input image in the correct format, ressamples it, and save it.
 * @param  param  parameters
 */
template <class TPixelType, unsigned int NDimension>
void resample(Param param)
{

    // Type definition
    typedef itk::Image<TPixelType, NDimension>                        ImageType;
    typedef itk::ImageFileReader<ImageType>                           ReaderType;
    typedef itk::ImageFileWriter<ImageType>                           WriterType;
    typedef itk::ResampleImageFilter<ImageType, ImageType, double>    ResampleFilterType;
    typedef itk::IdentityTransform<double, NDimension>                IdentityTransformType;
    typedef itk::Transform<double, NDimension, NDimension>            LinearTransformType;
    typedef itk::DisplacementFieldTransform<double, NDimension>       DFTransformType;
    typedef itk::StationaryVelocityFieldTransform<double, NDimension> SVFTransformType;
    typedef typename DFTransformType::VectorFieldType                 VectorFieldType;

    try
    {

        // Output image geometry
        typename ImageType::PointType     origin;
        typename ImageType::SpacingType   spacing;
        typename ImageType::SizeType      size;
        typename ImageType::DirectionType direction;

        // Read image
        typename ImageType::Pointer input = rpi::readImage<ImageType>( param.inputPath );

        // Switch on transformation type
        switch (param.transformType)
        {

            case IDENTITY_TRANSFORMATION :
                {
                    // Read transformation
                    typename IdentityTransformType::Pointer identity = IdentityTransformType::New();

                    // Get geometry
                    if (param.geometryPath.compare("")!=0)
                        rpi::getGeometryFromImageHeader<NDimension>(param.geometryPath, origin, spacing, size, direction);
                    else
                        rpi::getGeometryFromImage<ImageType>(input, origin, spacing, size, direction);

                    // Resample and write image
                    rpi::resampleAndWriteImage<ImageType,double>(identity, input, origin, spacing, size, direction, param.outputPath, param.interpType );
                }
                break;

            case LINEAR_TRANSFORMATION :
                {
                    // Read transformation
                    typename LinearTransformType::Pointer linear = rpi::readLinearTransformation<double>(param.transformPath);

                    // Get geometry
                    if (param.geometryPath.compare("")!=0)
                        rpi::getGeometryFromImageHeader<NDimension>(param.geometryPath, origin, spacing, size, direction);
                    else
                        rpi::getGeometryFromImage<ImageType>(input, origin, spacing, size, direction);

                    // Resample and write image
                    rpi::resampleAndWriteImage<ImageType,double>(linear, input, origin, spacing, size, direction, param.outputPath, param.interpType );
                }
                break;


            case DISPLACEMENT_FIELD_TRANSFORMATION :
                {
                    // Read transformation
                    typename DFTransformType::Pointer df = rpi::readDisplacementField<double>(param.transformPath);

                    // Get geometry
                    if (param.geometryPath.compare("")!=0)
                        rpi::getGeometryFromImageHeader<NDimension>(param.geometryPath, origin, spacing, size, direction);
                    else
                        rpi::getGeometryFromImage<VectorFieldType>(df->GetParametersAsVectorField(), origin, spacing, size, direction);

                    // Resample and write image
                    rpi::resampleAndWriteImage<ImageType,double>(df, input, origin, spacing, size, direction, param.outputPath, param.interpType );
                }
                break;


            case STATIONARY_VELOCITY_FIELD_TRANSFORMATION :
                {
                    // Read transformation
                    typename SVFTransformType::Pointer svf = rpi::readStationaryVelocityField<double>(param.transformPath);

                    // Get geometry
                    if (param.geometryPath.compare("")!=0)
                        rpi::getGeometryFromImageHeader<NDimension>(param.geometryPath, origin, spacing, size, direction);
                    else
                        rpi::getGeometryFromImage<VectorFieldType>(svf->GetParametersAsVectorField(), origin, spacing, size, direction);

                    // Create a displacement field from the stationary velocity field
                    typename DFTransformType::Pointer df = DFTransformType::New();
                    df->SetParametersAsVectorField( svf->GetDisplacementFieldAsVectorField() );

                    // Resample and write image
                    rpi::resampleAndWriteImage<ImageType,double>(df, input, origin, spacing, size, direction, param.outputPath, param.interpType );
                }
                break;

            default :
            {
                throw std::runtime_error("Transformation type not supported.");
            }
        }


    }
    catch( itk::ExceptionObject& e )
    {
        std::cerr << e << std::endl;
        throw std::runtime_error("Image resampling fail.");
    }
}


/**
 * Main function.
 */
int main(int argc, char** argv)
{

    // Allow to read and write Inrimage
    itk::InrimageImageIOFactory::RegisterOneFactory();

    // Parameters
    struct Param param;

    try
    {
        // Parse command line options
        parseParameters( argc, argv, param );
        printParameters( param );

        // Read ImageIO
        itk::ImageIOBase::Pointer         imageIO        = rpi::readImageInformation(param.inputPath);
        unsigned int                      dim            = imageIO->GetNumberOfDimensions();
        itk::ImageIOBase::IOPixelType     pixel_type     = imageIO->GetPixelType();
        itk::ImageIOBase::IOComponentType component_type = imageIO->GetComponentType();

        // Switch on pixel type
        if ( pixel_type==itk::ImageIOBase::SCALAR )
        {
            // Switch on pixel type and dimension
            if ( dim==3 )
            {
                if      ( component_type == itk::ImageIOBase::UCHAR  )
                    resample<unsigned char,  3>(param);
                else if ( component_type == itk::ImageIOBase::CHAR   )
                    resample<char,           3>(param);
                else if ( component_type == itk::ImageIOBase::USHORT )
                    resample<unsigned short, 3>(param);
                else if ( component_type == itk::ImageIOBase::SHORT  )
                    resample<short,          3>(param);
                else if ( component_type == itk::ImageIOBase::UINT   )
                    resample<unsigned int,   3>(param);
                else if ( component_type == itk::ImageIOBase::INT    )
                    resample<int,            3>(param);
                else if ( component_type == itk::ImageIOBase::ULONG  )
                    resample<unsigned long,  3>(param);
                else if ( component_type == itk::ImageIOBase::LONG   )
                    resample<long,           3>(param);
                else if ( component_type == itk::ImageIOBase::FLOAT  )
                    resample<float,          3>(param);
                else if ( component_type == itk::ImageIOBase::DOUBLE )
                    resample<double,         3>(param);
                else
                    throw std::runtime_error( "Component type not supported supported." );
            }
            else
                throw std::runtime_error( "Dimension not supported yet." );
        }
        else
            throw std::runtime_error( "Pixel type not supported. Only scalar images are supported yet." );
    }
    catch( std::exception& e )
    {
        std::cerr << std::endl << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    };

    return EXIT_SUCCESS;
}

