#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkDerivativeImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include <vnl/vnl_math.h>
#include <itkImageRegionIterator.h>
#include <itkDivideByConstantImageFilter.h>
#include <itkWarpVectorImageFilter.h>
#include <itkWarpImageFilter.h>
#include <itkWarpVectorImageFilter.h>
#include "string.h"
#include <tclap/CmdLine.h>
#include <itkExpImageFilter.h>
#include <itkLogImageFilter.h>

/*
 * The program implements the iterative computation of the logJacobian scalar map of a deformation field 
 * parametrized by a stationary velocity field (SVF). Given the one paramter subgroups property 
 * exp(2v)=exp(v) o exp(v), the jacobian computation is iteratively computed through the following scheme:
 *
 * Step1. Find n such that v0=inputSVF/n "small" (scaling step) 
 * Step2. Approximate log[det(Jac(v0))]=Div(v0)
 * Step3. Iterate n times log[det(Jac(vi))]=log[det(Jac(vi-1))|vi-1]+log[det(Jac(vi-1))] 
 *
 * for further details see

 * LCC-Demons: A robust and accurate symmetric diffeomorphic registration algorithm.
 * Lorenzi M, Ayache N, Frisoni GB, Pennec X; Neuroimage. 2013 Nov 1;81:470-83.
 *
 * The program requires the input SVF image path and the output logJacobian map path. 
 * It is possible to provide a Mask for the computation of scaling step (useful to control boundaries),
 * and to provide a scaling factor for the input SVF.
 */


/**
 * Structure containing the input parameters.
 */
struct Param{
    std::string  SVFImage;
    std::string  Mask;
    std::string  OutputImage;
    float ScalingFactor;
    bool  NumericalScheme;
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
    description += "Iterative Log Jacobian Map Computation";
    description += "\nAuthor : Marco Lorenzi";

    try {

        // Define the command line parser
        TCLAP::CmdLine cmd( description, ' ', "1.0", true);

        TCLAP::ValueArg<std::string>  arg_SVFImage( "i", "input-svf", "Path to the input stationary velocity field", true, "", "string", cmd );
        TCLAP::ValueArg<std::string>  arg_OutputImage( "o", "output-svf", "Path of the output LogJacobian map (default LogJacobian.mha).", false, "LogJacobian.mha", "string", cmd );
        TCLAP::ValueArg<std::string>  arg_Mask( "m", "mask", "Path to the mask (default whole image)", false, "null", "string", cmd );
        TCLAP::ValueArg<double>  arg_ScalingFactor( "s", "scaling-factor", "Scaling factor for the input velocity field (default 1)", false, 1.0, "double", cmd );
	TCLAP::ValueArg<double>  arg_NumericalScheme( "z", "numerical-scheme", "Numerical scheme for the exponential: 0 Scaling and squarings (default), 1 Forward Euler", false, 0.0, "bool", cmd );


        // Parse the command line
        cmd.parse( argc, argv );

        // Set the parameters
        param.SVFImage                     = arg_SVFImage.getValue();
        param.OutputImage                  = arg_OutputImage.getValue();
        param.Mask                         = arg_Mask.getValue();
        param.ScalingFactor                = arg_ScalingFactor.getValue();
	param.NumericalScheme              = arg_NumericalScheme.getValue();
        }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
        throw std::runtime_error("Unable to parse the command line arguments.");
    }
}


/**
  * Prints parameters.
  * @param  SVFImage                        path to the input SVF
  * @param  OutputImage      		    path to the output log jacobian map
  * @param  Mask                 	    path to the mask
  * @param  Scaling Factor                  ScalingFactor
  * @param  Numerical Scheme                Numerical Scheme
  */
void PrintParameters( std::string SVFImage,
                      std::string OutputImage,
                      std::string Mask,
                      double      ScalingFactor,
		      bool NumericalScheme
                    )
{

    // Print I/O parameters
    std::cout << std::endl;
    std::cout << "I/O PARAMETERS"                  << std::endl;
    std::cout << "  Input SVF image path           : " << SVFImage      << std::endl;
    std::cout << "  Output image path              : " << OutputImage     << std::endl;
    std::cout << "  Mask image path                : " << Mask     << std::endl;    std::cout << "  Scaling Factor                 : " << ScalingFactor << std::endl;
    std::cout << "  Numerical Scheme               : " << NumericalScheme << std::endl;
    std::cout << std::endl;
}


int main( int argc, char ** argv )
{
  
  //Parsing initial parameters
  struct Param param;
  parseParameters( argc, argv, param);

  //Definition of the type:  
  typedef itk::Vector<float,3>  VectorPixelType;
  typedef itk::Image<VectorPixelType,3> VectorImageType;	
  typedef itk::Image<float,3> ImageType;

  //Definition of Reader and Writer Classes:
  typedef itk::ImageFileReader< ImageType >  ScalarReaderType;
  typedef itk::ImageFileReader< VectorImageType >  VectorReaderType;

  ScalarReaderType::Pointer readerMask = ScalarReaderType::New();
  VectorReaderType::Pointer reader1  = VectorReaderType::New();

  float mult=param.ScalingFactor;
  try
  {
   reader1->SetFileName(param.SVFImage);
   reader1->Update();
  }
  catch 
   (TCLAP::ArgException &e)
  {
   std::cerr << "Error: " << e.what() << std::endl;
   return EXIT_FAILURE;
  }

  typedef itk::ImageRegionIterator<VectorImageType> IteratorType;  
  
  VectorImageType::Pointer Out1=VectorImageType::New();
  Out1->SetRegions(reader1->GetOutput()->GetLargestPossibleRegion() );
  Out1->SetSpacing( reader1->GetOutput()->GetSpacing() );
  Out1->SetOrigin( reader1->GetOutput()->GetOrigin() );
  Out1->SetDirection( reader1->GetOutput()->GetDirection() );
  Out1->Allocate();
  
  IteratorType out1(Out1,reader1->GetOutput()->GetLargestPossibleRegion());
 
  typedef itk::Index<3> IndexType;
  itk::Vector<float,3> Vect;  
  
  for (out1.GoToBegin();!out1.IsAtEnd();++out1 ) 
    {
      Vect.Fill(0);
      const IndexType index=out1.GetIndex();
      Vect= mult*(reader1->GetOutput()->GetPixel( index ));
      out1.Set( Vect);
    }

  double minpixelspacing = Out1->GetSpacing()[0];
  for (unsigned int i = 0; i<3; ++i)
  {
   if ( Out1->GetSpacing()[i] < minpixelspacing )
     minpixelspacing = Out1->GetSpacing()[i];
  }

   
  IteratorType InputIt (Out1, Out1->GetRequestedRegion());

  float norm2,maxnorm2=0;
   
  if (param.Mask!="null")
   try
    {
     readerMask->SetFileName(param.Mask);
     readerMask->Update();	
    }
   catch(TCLAP::ArgException &e)
    {
     std::cerr << "Error: " << e.what() << std::endl;
     return EXIT_FAILURE;
    }

  /* 
   *   Evaluate the maximum norm in the region of interest
   */
  
  for( InputIt.GoToBegin(); !InputIt.IsAtEnd(); ++InputIt )
  {
   if ((param.Mask!="null"))
   { if (readerMask->GetOutput()->GetPixel(InputIt.GetIndex())>0)
     {
      norm2 = InputIt.Get().GetSquaredNorm();
      if (maxnorm2<norm2) 
       maxnorm2=norm2;
     }
   }  
   else 
   {
     norm2 = InputIt.Get().GetSquaredNorm();
     if (maxnorm2<norm2) 
       maxnorm2=norm2;
    }
   }

  maxnorm2 /= (minpixelspacing * minpixelspacing);


  /**
    *    Evaluate the number numiter of iterations required 
   **/

 float divider;
 unsigned int numiter;

  if (param.NumericalScheme==1)
    {
  /**
   * Forward Euler: maxnorm(v)/numiter<0.5
   **/
      numiter = floor(2*sqrt(maxnorm2)+0.5);
      divider=numiter;
    }
  else
   {
  /**
   *  Scaling and Squaring: maxnorm(v)/2^numiter<0.5
   **/
    numiter = static_cast<unsigned int> (2.0 + 0.5 * std::log(maxnorm2)/vnl_math::ln2+0.5);
    divider = static_cast<float>(1<<numiter);
   }

  
  /**
    *     Scale the input SVF 
   **/
  typedef itk::DivideByConstantImageFilter<VectorImageType,float,VectorImageType> DividerType;
  DividerType::Pointer Divider=DividerType::New();
  Divider->SetInput(Out1);
  Divider->SetConstant(divider );
  Divider->Update();


  /**
    *     Initialize the iterative computation 
   **/


  typedef itk::AddImageFilter<VectorImageType,VectorImageType,VectorImageType> AdderType;
  AdderType::Pointer Adder=AdderType::New();  
  
  VectorImageType::Pointer Vect1=VectorImageType::New();
  Vect1=Divider->GetOutput();

  VectorImageType::Pointer warpedIm=VectorImageType::New();
   
  typedef itk::AddImageFilter<ImageType,ImageType,ImageType> AddImgFilterType;
  AddImgFilterType::Pointer Add=AddImgFilterType::New();

  /**
    *     Initial step: compute Div(v0)~logJac(v0) 
   **/

  
  typedef itk::VectorIndexSelectionCastImageFilter< VectorImageType,ImageType> FilterType;
  FilterType::Pointer componentExtractor0 = FilterType::New();
  componentExtractor0->SetIndex( 0 );
  FilterType::Pointer componentExtractor1 = FilterType::New();
  componentExtractor1->SetIndex( 1);
  FilterType::Pointer componentExtractor2 = FilterType::New();
  componentExtractor2->SetIndex( 2 );

  componentExtractor0->SetInput( Vect1 );
  componentExtractor1->SetInput( Vect1 );
  componentExtractor2->SetInput( Vect1 );

  // Compute partial derivatives

  typedef itk::DerivativeImageFilter<
             ImageType, ImageType > DerivativeFilterType;

  DerivativeFilterType::Pointer derivative0=DerivativeFilterType::New();
  DerivativeFilterType::Pointer derivative1=DerivativeFilterType::New();
  DerivativeFilterType::Pointer derivative2=DerivativeFilterType::New();

  derivative0->SetOrder(1);
  derivative0->SetDirection(0);

  derivative1->SetOrder(1);
  derivative1->SetDirection(1);

  derivative2->SetOrder(1);
  derivative2->SetDirection(2);

  derivative0->SetInput(componentExtractor0->GetOutput());
  derivative1->SetInput(componentExtractor1->GetOutput());
  derivative2->SetInput(componentExtractor2->GetOutput());

  derivative0->Update();
  derivative1->Update();
  derivative2->Update();
  
  AddImgFilterType::Pointer addFilter = AddImgFilterType::New();
  addFilter->SetInput1( derivative0->GetOutput() );
  addFilter->SetInput2( derivative1->GetOutput() );
  addFilter->Update();

  AddImgFilterType::Pointer addFilter1 = AddImgFilterType::New();
  addFilter1->SetInput1( derivative2->GetOutput() );
  addFilter1->SetInput2(addFilter->GetOutput());
  addFilter1->Update();

  ImageType::Pointer Image=addFilter1->GetOutput();

    typedef itk::AddImageFilter<VectorImageType,VectorImageType,VectorImageType> AdderType;
  AdderType::Pointer VectorAdder=AdderType::New();  


typedef itk::WarpVectorImageFilter<VectorImageType,VectorImageType,VectorImageType> WarperType;

  WarperType::Pointer VectorWarper=WarperType::New();
  VectorWarper->SetOutputOrigin(Out1->GetOrigin());
  VectorWarper->SetOutputSpacing(Out1->GetSpacing());
  VectorWarper->SetOutputDirection(Out1->GetDirection());
//  VectorWarper->SetEdgePaddingValue(1); 

  typedef itk::ExpImageFilter<ImageType,ImageType> ExpImageFilterType;
  ExpImageFilterType::Pointer ExpImage=ExpImageFilterType::New();

  typedef itk::LogImageFilter<ImageType,ImageType> LogImageFilterType;
  LogImageFilterType::Pointer LogImage=LogImageFilterType::New();


  ExpImage->SetInput(Image);
  ExpImage->Update();

  typedef itk::WarpImageFilter<ImageType,ImageType,VectorImageType> WarpImgType;
  WarpImgType::Pointer WarpImg=WarpImgType::New();

  WarpImg->SetOutputOrigin(Vect1->GetOrigin());
  WarpImg->SetOutputSpacing(Vect1->GetSpacing());
  WarpImg->SetOutputDirection(Vect1->GetDirection());
  WarpImg->SetEdgePaddingValue(1); 

   VectorImageType::Pointer WarpedIm=VectorImageType::New();
  
   if (param.NumericalScheme==1)
     {  
      /**
        * Iterative step Forward Euler: compute numiter times logJac(exp(vi))=logJac(exp(vi-1))|exp(v0)+logJac(exp(vi-1)) 
       **/
      for( int i=0; i<numiter; i++ )
        {
         WarpImg->SetInput(Image);
         VectorWarper->SetInput(Vect1);
         
         if (i==0)
	   {
	    WarpImg->SetDisplacementField(Vect1);
	   }
	  else
           {
	    WarpImg->SetDisplacementField(VectorAdder->GetOutput());
	   }

	  WarpImg->SetOutputOrigin(Vect1->GetOrigin());
          WarpImg->SetOutputSpacing(Vect1->GetSpacing());
          WarpImg->SetOutputDirection(Vect1->GetDirection());
          WarpImg->UpdateLargestPossibleRegion();
          
          Add->SetInput1(WarpImg->GetOutput());
 
	  if (i==0)
	    {
	     Add->SetInput2(Image);
	     VectorWarper->SetDisplacementField(Vect1);
	    }
          else
	    {
	     Add->SetInput2(Add->GetOutput());
	     VectorWarper->SetDisplacementField(VectorAdder->GetOutput());
            }
	  
	  Add->Update();
          VectorWarper->GetOutput()->SetRequestedRegion(Vect1->GetRequestedRegion());
          VectorWarper->Update();
          WarpedIm = VectorWarper->GetOutput();
          WarpedIm->DisconnectPipeline();
           
          VectorAdder->SetInput1(Vect1);
	  VectorAdder->SetInput2(WarpedIm);
	  VectorAdder->GetOutput()->SetRequestedRegion(Vect1->GetRequestedRegion() );
	  VectorAdder->Update(); 

	}
      }
   else
      {
  /**
    *     Iterative step Scaling and Squaring: compute numiter times logJac(exp(vi))=log[det(Jac(exp(vi-1)))|exp(vi-1)]+log[det(Jac(exp(vi-1)))] 
   **/
      for( int i=0; i<numiter; i++ )
	{
         WarpImg->SetInput(ExpImage->GetOutput());
         WarpImg->SetDisplacementField(Vect1);
         WarpImg->Update();

         LogImage->SetInput(WarpImg->GetOutput());
         LogImage->Update();

         Add->SetInput1(LogImage->GetOutput());
         Add->SetInput2(Image);
         Add->Update();

         VectorWarper->SetInput(Vect1); 
         VectorWarper->SetDisplacementField(Vect1); 
         VectorWarper->GetOutput()->SetRequestedRegion(Vect1->GetRequestedRegion());
         VectorWarper->Update();
   
         WarpedIm= VectorWarper->GetOutput(); 
         WarpedIm->DisconnectPipeline();

         VectorAdder->SetInput1(Vect1);
         VectorAdder->SetInput2(WarpedIm);
         VectorAdder->GetOutput()->SetRequestedRegion(Vect1->GetRequestedRegion());
         VectorAdder->Update();
     

         Vect1=VectorAdder->GetOutput();
         Vect1->DisconnectPipeline();
         Image=Add->GetOutput();
         Image->DisconnectPipeline();
         ExpImage->SetInput(Add->GetOutput());
         ExpImage->Update();
        }
      }

  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer WriterImg=ImageWriterType::New();

  WriterImg->SetInput(Add->GetOutput());
  WriterImg->SetFileName(param.OutputImage);
  WriterImg->Update();


  return EXIT_SUCCESS;
}
