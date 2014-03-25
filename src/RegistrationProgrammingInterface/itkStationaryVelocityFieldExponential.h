#ifndef _itkStationaryVelocityFieldExponential_h_
#define _itkStationaryVelocityFieldExponential_h_

#include "itkImageToImageFilter.h"
#include "itkDivideByConstantImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkOppositeImageFilter.h"
#include "itkWarpVectorImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h"
#include "itkAddImageFilter.h"
#include <itkMultiplyByConstantImageFilter.h>

namespace itk
{

 /**
  * @description StationaryVelocityFieldExponential (itk)
  * The Filter computes the Lie group exponential Exp(V) of the input vector field V.
  * The computation is performed by using the one-parameter subgroup property: Exp(V) = Exp(V/2) o Exp(V/2)
  * 
  * The computation is implemented with two iterative schemes ( SetIterativeScheme() method ): 
  * 
  * Scaling and squaring (faster, less accurate)
  * - find N such that Exp(V)[0]~ V[0]= V/2^N is "small"
  * - iteratively compute Exp(V)[i] = Exp(V)[i-1] o Exp(V)[i-i] (N steps)
  *
  *
  * Forward Euler (slower, more accuracy)
  * - find N such that Exp(V)[0]~ V[0]= V/2^N is "small"
  * - iteratively compute Exp(V)[i] = Exp(V)[i-1] o Exp(V)[0] (2^N steps)
  *
  * 
  * The filter implements also the computation of the log-Jacobian determinant log|Exp(V)| of the deformation field Exp(V).  
  * From Exp(V) = Exp(V/2) o Exp(V/2) we have  log|Exp(V)|= log|Exp(V/2)|o Exp(V/2) + log|Exp(V/2)|.
  * The iterative rule can be applied to both the "Scaling and squaring" and "Small steps update" schemes within the exponential computation cycle. 
  * A first approximation for log(Exp(V))[0] is Div(V/2^N).
  *
  * The filter is templated over the input and output vector fields
  */

template <class TInputImage, class TOutputImage>
class ITK_EXPORT StationaryVelocityFieldExponential : public ImageToImageFilter<TInputImage, TOutputImage>
{

public:

  typedef  StationaryVelocityFieldExponential        Self;
  typedef  ImageToImageFilter<TInputImage, TOutputImage>     Superclass;
  typedef  SmartPointer<Self>                                Pointer;
  typedef  SmartPointer<const Self>                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  itkTypeMacro (StationaryVelocityFieldExponential, ImageToImageFilter);

  typedef TInputImage 					InputImageType;
  typedef typename InputImageType::Pointer		InputImagePointer;
  typedef typename InputImageType::ConstPointer		InputImageConstPointer;
  typedef typename InputImageType::PixelType		InputPixelType;
  typedef typename InputPixelType::RealValueType	InputPixelRealValueType;
    

  typedef TOutputImage					OutputImageType;
  typedef typename OutputImageType::Pointer		OutputImagePointer;
  typedef typename OutputImageType::PixelType		OutputPixelType;

   
   
  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,TInputImage::ImageDimension );
  itkStaticConstMacro(PixelDimension, unsigned int,InputPixelType::Dimension );
  itkStaticConstMacro(OutputPixelDimension, unsigned int,OutputPixelType::Dimension );
    
  /** Define the scalar image type from the vector type **/
  typedef typename itk::Image<InputPixelRealValueType,ImageDimension>		ScalarImageType;
  typedef typename ScalarImageType::Pointer					ScalarImagePointer;

  /** Enable the computation of the Log Jacobian Scalar Image	0:No (default), 1:Yes      **/
  itkSetMacro(LogJacobianDeterminantComputation, bool);
  itkGetConstMacro(LogJacobianDeterminantComputation, bool);
  itkBooleanMacro(LogJacobianDeterminantComputation);

  /** Set multiplicative factor for the input Velocity **/
  itkSetMacro(MultiplicativeFactor, InputPixelRealValueType);
  itkGetConstMacro(MultiplicativeFactor, InputPixelRealValueType);


  /** Get the Log Jacobian Determinant Scalar Map **/

  ScalarImagePointer GetLogJacobianDeterminant(void)
  {
    if (this->m_LogJacobianDeterminantComputation==0)
    {itkExceptionMacro("Log Jacobian Computation was not enabled");}
    else 
      return m_LogJacobianDetImage;
  }

  /** Iterative scheme for the exponential computation **/
  enum NumericalScheme 
  {
    SCALING_AND_SQUARING,  /** Scaling and squaring  **/
    FORWARD_EULER	   /** Forward Euler  **/
  };


  /** Set the numerical scheme for the exponential computation **/
  void SetIterativeScheme(NumericalScheme scheme)
  {
    m_IterativeScheme = scheme;
  }

protected:
  StationaryVelocityFieldExponential();
  virtual ~StationaryVelocityFieldExponential(){};

  void PrintSelf(std::ostream& os,Indent indent) const;
  void GenerateData();

      
  typedef DivideByConstantImageFilter<InputImageType, InputPixelRealValueType, InputImageType>		DivideByConstantType;
  typedef MultiplyByConstantImageFilter<InputImageType, InputPixelRealValueType, InputImageType>	MultiplyByConstantType;

  typedef WarpVectorImageFilter< InputImageType, InputImageType, InputImageType>	                VectorWarperType;
  typedef VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<InputImageType, double>	FieldInterpolatorType;

  typedef WarpImageFilter< ScalarImageType, ScalarImageType, InputImageType>				ScalarWarperType;
  typedef AddImageFilter <InputImageType, InputImageType, InputImageType>				VectorAdderType;
  typedef AddImageFilter <ScalarImageType, ScalarImageType, ScalarImageType>				ScalarAdderType;

  typedef typename MultiplyByConstantType::Pointer							MultiplyByConstantPointer;
    
  typedef typename DivideByConstantType::Pointer							DivideByConstantPointer;
  typedef typename VectorWarperType::Pointer								VectorWarperPointer;
  typedef typename ScalarWarperType::Pointer								ScalarWarperPointer;
  typedef typename FieldInterpolatorType::OutputType							FieldInterpolatorOutputType;


  typedef typename VectorAdderType::Pointer								VectorAdderPointer;	
  typedef typename ScalarAdderType::Pointer								ScalardAdderPointer;

private:

  StationaryVelocityFieldExponential(const Self&);
  void operator=(const Self&);

  NumericalScheme				m_IterativeScheme;
  InputPixelRealValueType			m_MultiplicativeFactor;
  bool						m_LogJacobianDeterminantComputation;

  ScalarImagePointer    			m_LogJacobianDetImage;
  
  MultiplyByConstantPointer			m_Multiplier;
  DivideByConstantPointer			m_Divider;
  VectorWarperPointer				m_VectorWarper;
  ScalarWarperPointer				m_ScalarWarper;
  VectorAdderPointer				m_VectorAdder;
  ScalardAdderPointer				m_ScalarAdder;

};


} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStationaryVelocityFieldExponential.txx"
#endif


#endif
