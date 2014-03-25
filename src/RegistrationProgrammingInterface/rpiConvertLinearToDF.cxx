#include <iostream>
#include <cstdlib>

#include <itkImage.h>
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>

#include <tclap/CmdLine.h>

#include <mipsInrimageImageIOFactory.h>

#include "rpiCommonTools.hxx"



/**
 * Converts a linear transformation into a displacement field transformation.
 * The geometry of the displacement field is taken from the input image.
 * The initial linear transformation file must contain an ITK 3D transformation where elements
 * are coded as 'double'. The output displacement field is saved into a 3D vector image where each
 * element contains a 3D vector (of 'float') describing a displacement along the X-, Y-, and Z-axes.
 * @author Vincent Garcia
 * @date   2011/05/13
 */



/**
 * Structure containing the IO parameters.
 */
struct Param{
    std::string inputTransformPath;
    std::string inputImagePath;
    std::string outputTransformPath;
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
    description += "Generates a displacement field transformation from the input 3D linear ";
    description += "transformation. The displacement field is stored into a 3D image where each ";
    description += "element is a 3D vector containing the displacement along the X-, Y-, and Z-axes. ";
    description += "The geometry of the displacement field is taken from the input 3D image. ";
    description += "\nAuthor : Vincent Garcia";

    // Option description
    std::string dInputTransform  = "Path to the input linear transformation. ";
    dInputTransform             += "The dimension of the transformation must be 3. ";
    dInputTransform             += "The transformation must be an ITK transformation where element ";
    dInputTransform             += "type is 'double'.";
    std::string dInputImage      = "Path to the input 3D Image.";
    std::string dOutputTransform = "Path to the output displacement field transformation.";

    try {

        // Define the command line parser
        TCLAP::CmdLine cmd( description, ' ', "1.0", true);
        TCLAP::ValueArg<std::string> aOutputTransform( "o", "output-transform", dOutputTransform, true, "", "string", cmd );
        TCLAP::ValueArg<std::string> aInputImage(      "i", "input-image",      dInputImage,      true, "", "string", cmd );
        TCLAP::ValueArg<std::string> aInputTransform(  "t", "input-transform",  dInputTransform,  true, "", "string", cmd );

        // Parse the command line
        cmd.parse( argc, argv );

        // Set the parameters
        param.inputTransformPath  = aInputTransform.getValue();
        param.inputImagePath      = aInputImage.getValue();
        param.outputTransformPath = aOutputTransform.getValue();

    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
        throw std::runtime_error("Unable to parse the command line arguments.");
    }
}



/**
 * Main function.
 */
int main(int argc, char** argv)
{

    // Allows the executable to read and write Inrimage
    itk::InrimageImageIOFactory::RegisterOneFactory();

    // Type definition
    typedef  itk::ImageIOBase                                       ImageIOType;
    typedef  itk::Image<float,3>                                    ImageType;
    typedef  double                                                 LinearScalarType;
    typedef  itk::Transform< LinearScalarType, 3, 3 >               LinearTransformType;
    typedef  float                                                  FieldScalarType;
    typedef  itk::DisplacementFieldTransform< FieldScalarType, 3 >  FieldTransformType;

    try
    {
        // Parse parameters
        struct Param param;
        parseParameters(argc, argv, param);

        // Read linear transformation
        LinearTransformType::Pointer linear = rpi::readLinearTransformation<LinearScalarType>( param.inputTransformPath );

        // Read image information
        ImageIOType::Pointer imageIO;
        imageIO = rpi::readImageInformation( param.inputImagePath );
        if (  imageIO->GetNumberOfDimensions()!=3 )
            throw std::runtime_error("Only images of dimension 3 are supported yet.");

        // Read image
        ImageType::Pointer image = rpi::readImage<ImageType>( param.inputImagePath );

        // Generate displacement field from linear transformation
        FieldTransformType::Pointer field = rpi::linearToDisplacementFieldTransformation<LinearScalarType, FieldScalarType, ImageType>( image, linear );

        // Write displacement field
        rpi::writeDisplacementFieldTransformation<FieldScalarType,3>(field, param.outputTransformPath );
    }
    catch( std::exception& e )
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;

}
