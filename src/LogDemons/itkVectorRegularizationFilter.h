#ifndef __itkVectorRegularizationFilter_h
#define __itkVectorRegularizationFilter_h
 
#include "itkImageToImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include <itkComposeImageFilter.h>

#include "itkFFTWForwardFFTImageFilter.h"
#include "itkFFTWInverseFFTImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "fftw3.h"


namespace itk
{
/** \class VectorRegularizationFilter
 **/
template< class TRealVectorImage>
class VectorRegularizationFilter:public ImageToImageFilter< TRealVectorImage, TRealVectorImage >
{
public:
  /** Standard class typedefs. */
  typedef VectorRegularizationFilter             Self;
  typedef ImageToImageFilter< TRealVectorImage, TRealVectorImage > Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef typename NeighborhoodIterator<TRealVectorImage>::RadiusType RadiusType;

  typedef TRealVectorImage                                     RealVectorImageType;
  typedef typename TRealVectorImage::PixelType::RealValueType  RealPixelType;
  typedef typename TRealVectorImage::Pointer                   RealVectorImagePointer;
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TRealVectorImage::ImageDimension);

  typedef typename itk::Image<RealPixelType,ImageDimension>    RealImageType;
  typedef typename RealImageType::Pointer                      RealImagePointer;

  typedef typename itk::ComposeImageFilter<RealImageType,RealVectorImageType> ComposeVectorType;
  typedef typename ComposeVectorType::Pointer        ComposeVectorPointer;

  typedef typename itk::VectorIndexSelectionCastImageFilter< RealVectorImageType,RealImageType> ExtractorFilterType;
  typedef typename ExtractorFilterType::Pointer        ExtractorFilterPointer;




  typedef typename  itk::FFTWForwardFFTImageFilter<RealImageType >  FFTFilterType;
  typedef typename  FFTFilterType::Pointer            FFTFilterPointer;

  
  typedef typename  FFTFilterType::OutputImageType    ComplexImageType;
  typedef typename  ComplexImageType::Pointer         ComplexImagePointer;
  typedef typename  ComplexImageType::PixelType       ComplexPixelType;

  typedef typename itk::Vector<ComplexPixelType, ImageDimension> ComplexVectorType;
  typedef typename itk::Image<ComplexVectorType,ImageDimension>  ComplexVectorImageType;
  typedef typename ComplexVectorImageType::Pointer               ComplexVectorImagePointer;
  
  typedef itk::FFTWInverseFFTImageFilter< ComplexImageType, RealImageType >  FFTInvFilterType;
  typedef typename  FFTInvFilterType::Pointer            FFTInvFilterPointer;

  typedef typename itk::ComposeImageFilter<ComplexImageType,ComplexVectorImageType> ComposeComplexVectorFilterType;
  typedef typename ComposeComplexVectorFilterType::Pointer        ComposeComplexVectorPointer;

  typedef typename itk::VectorIndexSelectionCastImageFilter< ComplexVectorImageType,ComplexImageType> ExtractorComplexFilterType;
  typedef typename ExtractorComplexFilterType::Pointer        ExtractorComplexFilterPointer;
  
  typedef typename itk::FFTShiftImageFilter< ComplexImageType, ComplexImageType > ShiftFilterType;
  typedef typename ShiftFilterType::Pointer           ShiftFilterPointer;


  /** Method for creation through the object factory. */
  itkNewMacro(Self);
 
  /** Run-time type information (and related methods). */
  itkTypeMacro(VectorRegularizationFilter, ImageToImageFilter);
 
  void SetFactor(float factor[]);


protected:
  VectorRegularizationFilter();
  ~VectorRegularizationFilter(){}
 
  //virtual void BeforeThreadedGenerateData();
 
  virtual void GenerateData();

  
 
private:
  VectorRegularizationFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
 
  FFTFilterPointer                   m_FFTFilter;
  FFTInvFilterPointer                m_FFTInvFilter;
  ShiftFilterPointer                 m_ShiftFilter;

  ComposeComplexVectorPointer        m_ComplexComposer;
  ExtractorComplexFilterPointer      m_ComplexExtractor;

  ComposeVectorPointer               m_RealComposer;
  ExtractorFilterPointer             m_RealExtractor;

  ComplexImagePointer                m_ComplexImageArray[ImageDimension];

  RealImagePointer                   m_RealImageArray[ImageDimension];

  float 			     m_factor[2];
  
  int   			     m_dimensions[ImageDimension];

  ComplexVectorImagePointer          m_ComplexSolution;
  
};
} //namespace ITK
 
 
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorRegularizationFilter.txx"
#endif
 
#endif // __itkVectorRegularizationFilter_h
