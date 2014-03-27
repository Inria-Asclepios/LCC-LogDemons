#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include <tclap/CmdLine.h>
#include <rpiCommonTools.hxx>
#include "rpiLCClogDemons.hxx"



/**
 * LCC log-Demons executable.
 * @author Marco Lorenzi, Vincent Garcia, Florence Dru, and Tom Vercauteren
 * @date   2013/09/04
 */



/**
 * Structure containing the parameters.
 */
struct Param{
    std::string  fixedImagePath;
    std::string  movingImagePath;
    std::string  outputImagePath;
    std::string  outputTransformPath;
    std::string  MaskImagePath;
    std::string  trueField;
    std::string  outputDisplacementFieldPath;
    std::string  intialLinearTransformPath;
    std::string  intialFieldTransformPath;
    std::string  iterations;
    unsigned int updateRule;
    unsigned int gradientType;
    double       maximumUpdateStepLength;
    double       updateFieldStandardDeviation;
    double       stationaryVelocityFieldStandardDeviation;
    double       SimilarityCriteriaStandardDeviation;
    double       SigmaI;
    bool         useHistogramMatching;
    unsigned int RegularizationType;
    double        HarmonicWeight;
    double       BendingWeight;
    bool         verbose;
    unsigned int BCHExpansion;
    rpi::ImageInterpolatorType interpolatorType;
    bool         BoundaryCheck;

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
    description += "Local correlation coefficient for the diffeomorphic demons registration method. ";
    description += "\nAuthors : Marco Lorenzi, Vincent Garcia and Tom Vercauteren";

    // Option description
    std::string des_interpolatorType;
    des_interpolatorType  = "Type of image interpolator used for resampling the moving image at the ";
    des_interpolatorType += "very end of the program. Must be an integer value among 0 = nearest ";
    des_interpolatorType += "neighbor, 1 = linear, 2 = b-spline, and 3 = sinus cardinal (default 1).";
    
    std::string des_RegularizationType      = "Type of Regularization: 0 = Gaussian convolution (classical Demons), 1 = Harmonic + Bending Energy (default 1)";

    std::string des_HarmonicWeight          = "Weight (will be scaled by 1e-3) of the penalization of the Harmonic Energy (default= 0.0)";

    std::string des_BendingWeight           = "Weight (will be scaled by 1e-6) of the penalization of the Bending Energy (default= 1)";

    std::string des_BCHExpansion            = "Number of terms in the BCH expansion (default 2). ";

    std::string des_useHistogramMatching    = "Use histogram matching before processing? (default false, not used with LCC criterion). ";

    std::string des_verbose   		        = "Algorithm verbosity (default NO verbosity). ";

    std::string des_BoundaryCheck           = "Boundary checking for the LCC (to be used when registering segmented images with constant 0 background. Default: true. Set this option for disabling boundary checking).";

    std::string des_velFieldSigma           = "Standard deviation of the Gaussian smoothing of the stationary velocity field (world units). ";
    des_velFieldSigma                      += "Setting it below 0.1 means no smoothing will be performed (default 1.5).";

    std::string des_upFieldSigma            = "Standard deviation of the Gaussian smoothing of the update field (world units). ";
    des_upFieldSigma                       += "Setting it below 0.1 means no smoothing will be performed (default 0.0).";
    
    std::string des_SimCrSigma              = "Standard deviation of the Gaussian smoothing of the similarity criteria (world units, default 3 mm).  ";
    
    std::string des_SigmaI                  = "Trade-off between similarity and regularization (sigma_i ^2): 0 (sharper but unregular deformations) < sigma_i <= 1 (smoother deformations but weaker correspondencies). Default: 0.15    ";

    std::string des_gradientType            = "Type of gradient used for computing the demons force. ";
    des_gradientType                       += "0 is symmetrized, 1 is fixed image, 2 is warped moving image, 3 is mapped moving image (default 0, not used with LCC criterion which is symmetric by default).";

    std::string des_maxStepLength           = "Maximum length of an update vector (world units). ";
    des_maxStepLength                      += "Setting it to 0 implies no restrictions will be made on the step length.";

    std::string des_updateRule              = "Update rule:  0 : SSD similarity - exp(v)<-exp(v)oexp(u) (log-domain), ";
    des_updateRule                          += "1 : SSD similarity - exp(v)<-symmetrized(exp(v)oexp(u)) (symmetric log-domain), 2 : LCC similarity (default 2).";


    std::string des_iterations              = "Number of iterations per level of resolution (from coarse to fine levels). ";
    des_iterations                         += "Levels must be separated by \"x\" (default 30x20x10).";

    std::string des_initLinearTransform     = "Path to the initial linear transformation.";

    std::string des_initFieldTransform      = "Path to the initial stationary velocity field transformation.";

    std::string des_outputImage             = "Path to the output image (default output_image.nii.gz).";

    std::string des_outputDisplacementField = "Path of the output displacement field transformation (default output_displacement_field.mha).";

    std::string des_MaskImage               = "Path of the mask image.";

    std::string des_trueField		         = "Path of the true deformation field (default none).";

    std::string des_outputTransform         = "Path of the output stationary velocity field transformation (default output_stationary_velocity_field.mha).";

    std::string des_movingImage             = "Path to the moving image.";

    std::string des_fixedImage              = "Path to the fixed image.";


    try {

        // Define the command line parser
        TCLAP::CmdLine cmd( description, ' ', "1.0", true );

        // Set options
        TCLAP::ValueArg<unsigned int>  arg_interpolatorType( "", "interpolator-type", des_interpolatorType, false, 1, "int", cmd);
        TCLAP::ValueArg<unsigned int>  arg_BCHExpansion( "", "bch-expansion", des_BCHExpansion, false, 2, "uint", cmd );
        TCLAP::SwitchArg               arg_useHistogramMatching( "",  "use-histogram-matching", des_useHistogramMatching, cmd, false);
        TCLAP::SwitchArg               arg_verbose( "V",  "verbosity", des_verbose, cmd, false);
        TCLAP::ValueArg<double>        arg_velFieldSigma( "d", "velocity-field-sigma", des_velFieldSigma, false, 1.5, "double", cmd );
        TCLAP::ValueArg<double>        arg_upFieldSigma( "u", "update-field-sigma", des_upFieldSigma, false, 0.0, "double", cmd );
        TCLAP::ValueArg<double>        arg_SimCrSigma( "C", "Sim-Cr-sigma", des_SimCrSigma, false, 3.0, "double", cmd );
        TCLAP::ValueArg<unsigned int>  arg_gradientType( "g", "gradient-type", des_gradientType, false, 0, "uint", cmd );
        TCLAP::ValueArg<double>        arg_maxStepLength( "l", "max-step-length", des_maxStepLength, false, 2.0, "double", cmd );
        TCLAP::ValueArg<unsigned int>  arg_updateRule( "r", "update-rule", des_updateRule, false, 1, "uint", cmd );
        TCLAP::ValueArg<std::string>   arg_iterations( "a", "iterations", des_iterations, false, "30x20x10", "uintxuintx...xuint", cmd );
        TCLAP::ValueArg<std::string>   arg_initLinearTransform( "", "initial-linear-transform", des_initLinearTransform, false, "", "string", cmd );
        TCLAP::ValueArg<std::string>   arg_initFieldTransform( "", "initial-transform", des_initFieldTransform,  false, "", "string", cmd );
        TCLAP::ValueArg<std::string>   arg_trueField( "T", "true-field", des_trueField,  false, "", "string", cmd );
        TCLAP::ValueArg<std::string>   arg_outputImage( "i", "output-image", des_outputImage, false, "output_image.nii.gz", "string", cmd );
        TCLAP::ValueArg<std::string>   arg_outputDisplacementField( "", "output-displacement-field", des_outputDisplacementField, false, "output_displacement_field.mha", "string", cmd );
        TCLAP::ValueArg<std::string>   arg_outputTransform( "t", "output-transform", des_outputTransform, false, "output_stationary_velocity_field.mha", "string", cmd );
        TCLAP::ValueArg<std::string>   arg_movingImage( "m", "moving-image", des_movingImage, true, "", "string", cmd );
        TCLAP::ValueArg<std::string>   arg_fixedImage( "f", "fixed-image", des_fixedImage, true, "", "string", cmd );
        TCLAP::ValueArg<std::string>   arg_MaskImage( "M", "mask-image", des_MaskImage, false, "", "string", cmd );
        TCLAP::ValueArg<double>        arg_SigmaI( "S", "sigma-I", des_SigmaI, false, 0.15, "double", cmd );
        TCLAP::ValueArg<unsigned int>  arg_RegularizationType( "R",  "regularization", des_RegularizationType, false, 0,"uint", cmd);
        TCLAP::ValueArg<double>        arg_HarmonicWeight( "x", "Harmonic-weight", des_HarmonicWeight, false, 0.0, "double", cmd );
        TCLAP::ValueArg<double>        arg_BendingWeight( "b", "bending-weight", des_BendingWeight, false, 1.0, "double", cmd );
        TCLAP::SwitchArg               arg_BoundaryCheck( "B", "boundary-check", des_BoundaryCheck, cmd, true);
        // Parse the command line
        cmd.parse( argc, argv );

        // Set the parameters
        param.fixedImagePath                           = arg_fixedImage.getValue();
        param.movingImagePath                          = arg_movingImage.getValue();
        param.MaskImagePath                            = arg_MaskImage.getValue();
        param.outputTransformPath                      = arg_outputTransform.getValue();
        param.trueField		                           = arg_trueField.getValue();
        param.outputDisplacementFieldPath              = arg_outputDisplacementField.getValue();
        param.outputImagePath                          = arg_outputImage.getValue();
        param.intialLinearTransformPath                = arg_initLinearTransform.getValue();
        param.intialFieldTransformPath                 = arg_initFieldTransform.getValue();
        param.iterations                               = arg_iterations.getValue();
        param.updateRule                               = arg_updateRule.getValue();
        param.maximumUpdateStepLength                  = arg_maxStepLength.getValue();
        param.gradientType                             = arg_gradientType.getValue();
        param.updateFieldStandardDeviation             = arg_upFieldSigma.getValue();
        param.SimilarityCriteriaStandardDeviation      = arg_SimCrSigma.getValue();
        param.SigmaI 				                   = arg_SigmaI.getValue();
        param.stationaryVelocityFieldStandardDeviation = arg_velFieldSigma.getValue();
        param.useHistogramMatching                     = arg_useHistogramMatching.getValue();
        param.verbose 	 		                       = arg_verbose.getValue();
        param.BCHExpansion                             = arg_BCHExpansion.getValue();
        param.HarmonicWeight                           = arg_HarmonicWeight.getValue();
        param.BendingWeight                            = arg_BendingWeight.getValue();
        param.BoundaryCheck                            = arg_BoundaryCheck.getValue();

	// Set the interpolator type
        unsigned int interpolator_type = arg_interpolatorType.getValue();
        if      ( interpolator_type==0 )
            param.interpolatorType = rpi::INTERPOLATOR_NEAREST_NEIGHBOR;
        else if ( interpolator_type==1 )
            param.interpolatorType = rpi::INTERPOLATOR_LINEAR;
       else if  ( interpolator_type==2 )
            param.interpolatorType = rpi::INTERPOLATOR_BSLPINE;
        else if ( interpolator_type==3 )
            param.interpolatorType = rpi::INTERPOLATOR_SINUS_CARDINAL;
        else
            throw std::runtime_error("Image interpolator not supported.");

    // Check the regularization type
       unsigned int regularization_type = arg_RegularizationType.getValue();
       if  ( regularization_type==0 )
            param.RegularizationType = 0;
        else if ( interpolator_type==1 )
            param.RegularizationType = 1;
        else
            throw std::runtime_error("Regularization type not supported.");
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
        throw std::runtime_error("Unable to parse the command line arguments.");
    }
}



/**
  * Prints parameters.
  * @param  fixedImagePath               path to the fixed image
  * @param  movingImagePath              path to the moving image
  * @param  outputImagePath              path to the output file containing the resampled image
  * @param  outputTransformPath          path to the output file containing the transformation
  * @param  outputDisplacementFieldPath  path to the output file containing the displacement field
  * @param  initialLinearTransformPath   path to the output file containing the transformation
  * @param  initialFieldTransformPath    path to the output file containing the transformation
  * @param  interpolatorType             interpolator type
  * @param  registration                 registration object
  * @param  regularization               regularization type
  */
template< class TFixedImage, class TMovingImage, class TTransformScalarType >
void PrintParameters( std::string fixedImagePath,
                      std::string movingImagePath,
                      std::string outputImagePath,
                      std::string MaskImagePath,
                      std::string outputTransformPath,
                      std::string outputDisplacementFieldPath,
                      std::string initialLinearTransformPath,
                      std::string initialFieldTransformPath,
                      unsigned int RegularizationType,
                      rpi::ImageInterpolatorType interpolatorType,
                      rpi::LCClogDemons<TFixedImage, TMovingImage, TTransformScalarType> * registration )
{
    // Print I/O parameters
    std::cout << std::endl;
    std::cout << "I/O PARAMETERS"                             << std::endl;
    std::cout << "  Fixed image path                             : " << fixedImagePath              << std::endl;
    std::cout << "  Moving image path                            : " << movingImagePath             << std::endl;
    std::cout << "  Output image path                            : " << outputImagePath             << std::endl;
    std::cout << "  Mask image path                            : " << MaskImagePath             << std::endl;
    std::cout << "  Output transformation path                   : " << outputTransformPath         << std::endl;
    std::cout << "  Output displacement field path               : " << outputDisplacementFieldPath << std::endl;
    if ( initialLinearTransformPath.compare("")!=0 )
        std::cout << "  Initial linear transform                     : " << initialLinearTransformPath << std::endl;
    if ( initialFieldTransformPath.compare("")!=0 )
        std::cout << "  Initial displacement field transform         : " << initialFieldTransformPath  << std::endl;
    std::cout << std::endl;

    // Print method parameters
    std::cout << "METHOD PARAMETERS"                          << std::endl;
    std::cout << "  Iterations                                   : " << rpi::VectorToString<unsigned int>( registration->GetNumberOfIterations() )     << std::endl;


    if(registration->GetUpdateRule() ==2)
      {
       std::cout << "  Similarity metric:                           : " << "LCC"                                                  << std::endl;
       std::cout << "  Similarity Criterion standard deviation      : " << registration->GetSimilarityCriteriaStandardDeviation()	    << " (world unit)" << std::endl;
       std::cout << "  Trade-off parameter                          : " << registration->GetSigmaI()	    << std::endl;
       std::cout << "  Boundary Checking                            : " << rpi::BooleanToString(registration->GetBoundaryCheck())	                     << std::endl;
      }
    else
      {
        std::cout << "  Similarity metric:                           : " << "SSD"                                                  << std::endl;
        std::cout << "  Maximum step length                          : " << registration->GetMaximumUpdateStepLength()                  << " (voxel unit)" << std::endl;
        std::cout << "  Gradient type                                : " << registration->GetGradientType()                                                << std::endl;
        std::cout << "  Use histogram matching?                      : " << rpi::BooleanToString( registration->GetUseHistogramMatching() )                << std::endl;
      }
    if (RegularizationType==1)
      {
        //THIS OPTION REQUIRES FFTW WHICH IS CURRENTLY NOT INSTALLED 
        RegularizationType=0;
        std::cout << "WARNING: This option is currently not available (requires fftw)" <<std::endl;
       /* std::cout << "Regularizer: Harmonic + Bending energy           "<<std::endl;
        std::cout << "  Weight of harmonic energy term               :" << registration->GetHarmonicWeight() << "" << std::endl;
        std::cout << "  Weight of bending energy term                :" << registration->GetBendingWeight()	    << "" << std::endl;
        */
      }
    if (RegularizationType==0)
      {
       std::cout << "Regularizer: Gaussian Convolution "<<std::endl;
       std::cout << "  Stationary velocity field standard deviation : " << registration->GetStationaryVelocityFieldStandardDeviation() << " (world unit)" << std::endl;
       std::cout << "  Update field standard deviation              : " << registration->GetUpdateFieldStandardDeviation()             << " (world unit)" << std::endl;
       
     }

    std::cout << "  Terms in BCH expansion                       : " << registration->GetNumberOfTermsBCHExpansion()                                   << std::endl << std::endl;
    std::cout << "  Interpolator type                     : " << rpi::getImageInterpolatorTypeAsString(interpolatorType)                    << std::endl;

}



/**
  * Starts the image registration.
  * @param   param  parameters needed for the image registration process
  * @return  EXIT_SUCCESS if the registration succeded, EXIT_FAILURE otherwise
  */
template< class TFixedImage, class TMovingImage >
int StartMainProgram(struct Param param)
{

    typedef double
            TransformScalarType;

    typedef rpi::LCClogDemons< TFixedImage, TMovingImage, TransformScalarType >
            RegistrationMethod;

    typedef itk::Transform<double, 3, 3>
            LinearTransformType;

    typedef rpi::DisplacementFieldTransform< TransformScalarType, 3 >
            FieldTransformType;


    // Creation of the registration object
    RegistrationMethod * registration = new RegistrationMethod();

    try
    {

        // Read input images
        typename TFixedImage::Pointer  fixedImage  = rpi::readImage< TFixedImage >(  param.fixedImagePath );
        typename TMovingImage::Pointer movingImage = rpi::readImage< TMovingImage >( param.movingImagePath );

        if (param.MaskImagePath.compare("")!=0)
          {
           typename TMovingImage::Pointer MaskImage = rpi::readImage< TMovingImage >( param.MaskImagePath );
           registration->UseMask(true);
           registration->SetMaskImage(MaskImage );
          }
        else
          registration->UseMask(false);

        // Set parameters
        registration->SetFixedImage(                               fixedImage );
        registration->SetMovingImage(                              movingImage );


        registration->SetNumberOfIterations(                       rpi::StringToVector<unsigned int>( param.iterations ) );
        registration->SetMaximumUpdateStepLength(                  param.maximumUpdateStepLength );
        registration->SetSimilarityCriteriaStandardDeviation(      param.SimilarityCriteriaStandardDeviation );
        registration->SetSigmaI(                                   param.SigmaI );
        registration->SetRegularizationType(                       param.RegularizationType);
        registration->SetBoundaryCheck(                            param.BoundaryCheck );

        switch (param.RegularizationType)
          {
           case 0:
            registration->SetUpdateFieldStandardDeviation(             param.updateFieldStandardDeviation );
            registration->SetStationaryVelocityFieldStandardDeviation( param.stationaryVelocityFieldStandardDeviation ); break;
           case 1:
            registration->SetHarmonicWeight(                           param.HarmonicWeight );
            registration->SetBendingWeight(                            param.BendingWeight ); break;
           }

        registration->SetVerbosity(                     	       param.verbose);
        registration->SetNumberOfTermsBCHExpansion(                param.BCHExpansion );



	 if ( param.trueField.compare("")!=0)
	 	{
         typename rpi::DisplacementFieldTransform<TransformScalarType,3>::Pointer TrueFieldImage = rpi::readDisplacementField<TransformScalarType>( param.trueField );
         registration->SetTrueField( TrueFieldImage);
		}

        // Set update rule
        switch( param.updateRule )
        {
        case 0:
            registration->SetUpdateRule(                               RegistrationMethod::UPDATE_LOG_DOMAIN );
            registration->SetUseHistogramMatching(                     param.useHistogramMatching );
            // Set gradient type
            switch( param.gradientType )
            {
            case 0:
                registration->SetGradientType( RegistrationMethod::GRADIENT_SYMMETRIZED );         break;
            case 1:
                registration->SetGradientType( RegistrationMethod::GRADIENT_FIXED_IMAGE );         break;
            case 2:
                registration->SetGradientType( RegistrationMethod::GRADIENT_WARPED_MOVING_IMAGE ); break;
            case 3:
                registration->SetGradientType( RegistrationMethod::GRADIENT_MAPPED_MOVING_IMAGE ); break;
            default:
                throw std::runtime_error( "Gradient type must fit in the range [0,3]." );
            } break;
        case 1:
            registration->SetUpdateRule(                               RegistrationMethod::UPDATE_SYMMETRIC_LOG_DOMAIN );
            registration->SetUseHistogramMatching(                     param.useHistogramMatching );
            // Set gradient type
            switch( param.gradientType )
            {
            case 0:
                registration->SetGradientType( RegistrationMethod::GRADIENT_SYMMETRIZED );         break;
            case 1:
                registration->SetGradientType( RegistrationMethod::GRADIENT_FIXED_IMAGE );         break;
            case 2:
                registration->SetGradientType( RegistrationMethod::GRADIENT_WARPED_MOVING_IMAGE ); break;
            case 3:
                registration->SetGradientType( RegistrationMethod::GRADIENT_MAPPED_MOVING_IMAGE ); break;
            default:
                throw std::runtime_error( "Gradient type must fit in the range [0,3]." );
            } break;
        case 2:
            registration->SetUpdateRule(   RegistrationMethod::UPDATE_SYMMETRIC_LOCAL_LOG_DOMAIN ); break;
        default:
            throw std::runtime_error( "Update rule must fit in the range [0,2]." );
        }





        // Initialize transformation
        if ( param.intialFieldTransformPath.compare("")!=0  &&  param.intialLinearTransformPath.compare("")!=0 )
        {
            throw std::runtime_error( "Cannot initialize with a stationary velocity field and a linear transformation." );
        }
        else if ( param.intialFieldTransformPath.compare("")!=0 )
        {
            typename FieldTransformType::Pointer field = rpi::readDisplacementField<TransformScalarType>( param.intialFieldTransformPath );
            registration->SetInitialTransformation( field );
        }
        else if ( param.intialLinearTransformPath.compare("")!=0 )
        {
            typename LinearTransformType::Pointer linear = rpi::readLinearTransformation<double>( param.intialLinearTransformPath );
            typename FieldTransformType::Pointer  field  = rpi::linearToDisplacementFieldTransformation<double, TransformScalarType, TFixedImage>( fixedImage, linear );
            registration->SetInitialTransformation( field );
        }


        // Print parameters
        PrintParameters< TFixedImage, TMovingImage, TransformScalarType >(
                param.fixedImagePath,
                param.movingImagePath,
                param.outputImagePath,
                param.MaskImagePath,
                param.outputTransformPath,
                param.outputDisplacementFieldPath,
                param.intialLinearTransformPath,
                param.intialFieldTransformPath,
                param.RegularizationType,
                param.interpolatorType,
                registration );


        // Display
        std::cout << "STARTING MAIN PROGRAM" << std::endl;


        // Start registration process
        std::cout << "  Registering images                    : " << std::flush;
        registration->StartRegistration();
        std::cout << "OK" << std::endl;


        // Write stationary velocity field
        std::cout << "  Writing stationary velocity field     : " << std::flush;
        rpi::writeStationaryVelocityFieldTransformation<TransformScalarType, TFixedImage::ImageDimension>(
                registration->GetTransformation(),
                param.outputTransformPath );
        std::cout << "OK" << std::endl;


        // Write displacement field
        std::cout << "  Writing displacement field            : " << std::flush;
        rpi::writeDisplacementFieldTransformation<TransformScalarType, TFixedImage::ImageDimension>(
                registration->GetDisplacementFieldTransformation(),
                param.outputDisplacementFieldPath );
        std::cout << "OK" << std::endl;


        // Write the output image
        std::cout << "  Writing image                        : " << std::flush;
        rpi::resampleAndWriteImage<TFixedImage, TMovingImage, TransformScalarType>(
		    fixedImage,
                    movingImage,
                    registration->GetDisplacementFieldTransformation(),
                    param.outputImagePath,
                    param.interpolatorType	);
        std::cout << "OK" << std::endl << std::endl;
    }
    catch( std::exception& e )
    {
        std::cerr << e.what() << std::endl;
        delete registration;
        return EXIT_FAILURE;
    };


    delete registration;
    return EXIT_SUCCESS;
}



/**
 * Main function.
 */
int main(int argc, char** argv)
{
   // itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
    // Parse parameters
    struct Param param;
    parseParameters( argc, argv, param);


    // Read image information
    itk::ImageIOBase::Pointer fixed_imageIO;
    itk::ImageIOBase::Pointer moving_imageIO;
    try
    {
        fixed_imageIO  = rpi::readImageInformation( param.fixedImagePath );
        moving_imageIO = rpi::readImageInformation( param.movingImagePath );
    }
    catch( std::exception& e )
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }


    // Only 3D images are supported yet
    if (  fixed_imageIO->GetNumberOfDimensions()!=3  &&  moving_imageIO->GetNumberOfDimensions()!=3  )
    {
        std::cerr << "Error: Only images of dimension 3 are supported yet." << std::endl;
        return EXIT_FAILURE;
    }


    // Only scalar images are supported yet
    if (  fixed_imageIO->GetPixelType() != itk::ImageIOBase::SCALAR  ||  moving_imageIO->GetPixelType() != itk::ImageIOBase::SCALAR  )
    {
        std::cerr << "Error: Only scalar images are supported yet." << std::endl;
        return EXIT_FAILURE;
    }


    // Start main program
    return StartMainProgram< itk::Image<double,3> , itk::Image<double,3> >( param );
}
