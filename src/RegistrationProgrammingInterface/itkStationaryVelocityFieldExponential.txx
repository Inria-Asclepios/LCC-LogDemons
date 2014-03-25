#ifndef _itkStationaryVelocityFieldExponential_txx_
#define _itkStationaryVelocityFieldExponential_txx_

#include "itkStationaryVelocityFieldExponential.h"
#include "itkImageRegionConstIterator.h"
#include "itkDerivativeImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"


/**
 * Given vector field V, the function Divergence(V) compute the divergence scalar map d/dx Vx + d/dy Vy + d/dz Vz. 
 * This is an external function which is templated on the input vector Image type and the otput scalar image type
 */

template<class TVectorType,class TScalarType>
itk::SmartPointer<TScalarType> Divergence(itk::SmartPointer<TVectorType> VectorField)
  {

  /** Extract the three component of the vector field  **/
  typedef typename itk::VectorIndexSelectionCastImageFilter< TVectorType,TScalarType> FilterType;
  typename FilterType::Pointer componentExtractor0 = FilterType::New();
  componentExtractor0->SetIndex( 0 );
  typename FilterType::Pointer componentExtractor1 = FilterType::New();
  componentExtractor1->SetIndex( 1);
  typename FilterType::Pointer componentExtractor2 = FilterType::New();
  componentExtractor2->SetIndex( 2 );

  componentExtractor0->SetInput( VectorField );
  componentExtractor1->SetInput( VectorField );
  componentExtractor2->SetInput( VectorField );

  /** Compute the first order derivative for each component along the correspondent direction **/ 
  typedef itk::DerivativeImageFilter<TScalarType, TScalarType > DerivativeFilterType;

  typename DerivativeFilterType::Pointer derivative0=DerivativeFilterType::New();
  typename DerivativeFilterType::Pointer derivative1=DerivativeFilterType::New();
  typename DerivativeFilterType::Pointer derivative2=DerivativeFilterType::New();

  derivative0->SetOrder(1);
  derivative0->SetDirection(0);

  derivative1->SetOrder(1);
  derivative1->SetDirection(1);

  derivative2->SetOrder(1);
  derivative2->SetDirection(2);

  derivative0->SetInput(componentExtractor0->GetOutput());
  derivative1->SetInput(componentExtractor1->GetOutput());
  derivative2->SetInput(componentExtractor2->GetOutput());

  /** Compute the divergence  **/
  typedef typename itk::AddImageFilter<TScalarType,TScalarType,TScalarType> AddImgFilterType;

  typename  AddImgFilterType::Pointer addFilter = AddImgFilterType::New();
  addFilter->SetInput1( derivative0->GetOutput() );
  addFilter->SetInput2( derivative1->GetOutput() );
  addFilter->Update();
  

  typename AddImgFilterType::Pointer addFilter1 = AddImgFilterType::New();
  addFilter1->SetInput1( derivative2->GetOutput() );
  addFilter1->SetInput2(addFilter->GetOutput());
  addFilter1->Update();


  return(addFilter1->GetOutput());

};


/**
 * 
 *  Definitions of the methods for itkVelocityFieldExponential
 *
**/

namespace itk
{
/**
 * Constructor
 */
template <class TInputImage, class TOutputImage> StationaryVelocityFieldExponential<TInputImage, TOutputImage>::StationaryVelocityFieldExponential()
{
  m_Divider 	    = 		DivideByConstantType::New();
  m_Multiplier 	    = 		MultiplyByConstantType::New();

  m_VectorWarper    = 		VectorWarperType::New();
  m_ScalarWarper    = 		ScalarWarperType::New();

  m_IterativeScheme =	   	SCALING_AND_SQUARING;  /** Default: Scaling and squaring **/
  m_MultiplicativeFactor = 	1.0;  /** Default: No scaling factor **/	
  m_LogJacobianDeterminantComputation=	0;  /** Default: No Log-Jacobian determinant computation**/

  m_LogJacobianDetImage=0;

  m_ScalarAdder = ScalarAdderType::New();
  m_VectorAdder = VectorAdderType::New();

}

/**
 * Print out a description of self
 */
template <class TInputImage, class TOutputImage>
void
StationaryVelocityFieldExponential<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "IterativeScheme: "
     << m_IterativeScheme << std::endl;
  os << indent << "LogJacobianComputation:   "
     << (m_LogJacobianDeterminantComputation?"On":"Off") << std::endl;

  return;
}


/**
 * GenerateData
 */
template <class TInputImage, class TOutputImage>
void
StationaryVelocityFieldExponential<TInputImage,TOutputImage>
::GenerateData()
{
  itkDebugMacro(<<"Actually executing");
 
  InputImageConstPointer inputPtr;

  /*Scaling the input velocity field*/
  if (m_MultiplicativeFactor !=1.0)
    {
    m_Multiplier->SetInput(this->GetInput());
    m_Multiplier->SetConstant(m_MultiplicativeFactor); 
    inputPtr=m_Multiplier->GetOutput();
    }

  else inputPtr=this->GetInput();

  unsigned int numiter = 0;

  //* TODO: Add different scaling methods based on the iterative scheme **/
    
  // Compute a good number of iterations based on the rationale
  // that the initial first order approximation,
  // exp(Phi/2^N) = Phi/2^N,
  // needs to be diffeomorphic. For this we simply impose to have
  // max(norm(Phi)/2^N) < 0.5*pixelspacing

  InputPixelRealValueType maxnorm2 = 0.0;

  double minpixelspacing = inputPtr->GetSpacing()[0];
  for (unsigned int i = 1; i<itkGetStaticConstMacro(ImageDimension); ++i)
    {
    if ( inputPtr->GetSpacing()[i] < minpixelspacing )
      {
      minpixelspacing = inputPtr->GetSpacing()[i];
      }
    }

  typedef ImageRegionConstIterator<InputImageType> InputConstIterator;
  InputConstIterator InputIt = InputConstIterator(
  inputPtr, inputPtr->GetRequestedRegion());

  for( InputIt.GoToBegin(); !InputIt.IsAtEnd(); ++InputIt )
    {
    InputPixelRealValueType norm2 = InputIt.Get().GetSquaredNorm();
    if (norm2>maxnorm2) maxnorm2=norm2;
    }

  // Divide the norm by the minimum pixel spacing
  maxnorm2 /= vnl_math_sqr(minpixelspacing);

  InputPixelRealValueType numiterfloat = 2.0 +  0.5 * vcl_log(maxnorm2)/vnl_math::ln2;


  // take the ceil and threshold
  numiter = static_cast<unsigned int>(numiterfloat + 1.0);
    
  std::cout<<"numiter "<<numiter<<" -- numiterfloat "<<numiterfloat<<" -- divider "<<static_cast<InputPixelRealValueType>(1<<numiter) <<std::endl;

  // Get the first order approximation (division by 2^numiter)
  m_Divider->SetInput(inputPtr);
  m_Divider->GraftOutput( this->GetOutput() );
  m_Divider->SetConstant( static_cast<InputPixelRealValueType>(1<<numiter) );

  m_Divider->Update();

  // Region passing stuff
  this->GraftOutput( m_Divider->GetOutput() );
  this->GetOutput()->Modified();

  if (m_LogJacobianDeterminantComputation)
    {
    m_LogJacobianDetImage=Divergence<InputImageType,ScalarImageType>(m_Divider->GetOutput());
    m_ScalarWarper->SetOutputOrigin(this->GetOutput()->GetOrigin());
    m_ScalarWarper->SetOutputSpacing(this->GetOutput()->GetSpacing());
    m_ScalarWarper->SetOutputDirection(this->GetOutput()->GetDirection());
    }

  /** Initialize the warper for the iterative composition of the vector field **/
  m_VectorWarper->SetOutputOrigin(inputPtr->GetOrigin());
  m_VectorWarper->SetOutputSpacing(inputPtr->GetSpacing());
  m_VectorWarper->SetOutputDirection(inputPtr->GetDirection());

  /** Initialize the update for the exponential computation **/
  InputImagePointer ActualVectorField=m_Divider->GetOutput();

  /** First step of the exponential computation  **/
  this->GraftOutput( m_Divider->GetOutput() );
  this->GetOutput()->Modified();



  /**  Define the number of iteration based on the iterative scheme: scaling and squaring: N, small steps: 2^N **/
  if (m_IterativeScheme==FORWARD_EULER)
    numiter=static_cast<InputPixelRealValueType>(1<<numiter)-1;

  for( unsigned int i=0; i<numiter; i++ )
    {
    /** Compute the Log-Jacobian determinant inside whitin the exponential computation **/
    if (m_LogJacobianDeterminantComputation)
      {
      if (i==0)  /**  The initial approximation for the Jacobian is the divergence **/
        m_ScalarWarper->SetInput(m_LogJacobianDetImage);
      else
        m_ScalarWarper->SetInput(m_ScalarAdder->GetOutput());	

      m_ScalarWarper->SetDeformationField(ActualVectorField);
      m_ScalarWarper->UpdateLargestPossibleRegion();
      
      m_ScalarAdder->SetInput1(m_ScalarWarper->GetOutput());

      if (i==0||m_IterativeScheme==FORWARD_EULER)
        m_ScalarAdder->SetInput2(m_LogJacobianDetImage);
      else 
        m_ScalarAdder->SetInput2(m_ScalarAdder->GetOutput());

      m_ScalarAdder->Update();
      }

    /** Compute the exponential by iterative composition **/
    
    m_VectorWarper->SetInput(this->GetOutput());
    
    m_VectorWarper->SetDeformationField(ActualVectorField);
    
    m_VectorWarper->GetOutput()->SetRequestedRegion(this->GetOutput()->GetRequestedRegion() );

    m_VectorWarper->Update();

    OutputImagePointer warpedIm = m_VectorWarper->GetOutput();
    warpedIm->DisconnectPipeline();

    m_VectorAdder->SetInput1(ActualVectorField);

    m_VectorAdder->SetInput2(warpedIm);
    m_VectorAdder->GetOutput()->SetRequestedRegion(
    this->GetOutput()->GetRequestedRegion() );

    m_VectorAdder->Update();

    // Region passing stuff
    this->GraftOutput( m_VectorAdder->GetOutput() );

    // Make a call to modified. This seems only necessary for
    // a non-inplace adder but it doesn't hurt anyhow.
    this->GetOutput()->Modified();
	
    if (m_IterativeScheme==SCALING_AND_SQUARING)	
      ActualVectorField=this->GetOutput();
		
    }
   
  if (m_LogJacobianDeterminantComputation) 
    m_LogJacobianDetImage=m_ScalarAdder->GetOutput();

  }


} // end namespace itk

#endif



