#ifndef __itkBendingEnergyRegularizationFilter_h
#define __itkBendingEnergyRegularizationFilter_h
 
#include "itkImageToImageFilter.h"
#include "itkNeighborhoodIterator.h"

namespace itk
{
/** \class BendingEnergyRegularizationFilter
 **/
template< class TImage>
class BendingEnergyRegularizationFilter:public ImageToImageFilter< TImage, TImage >
{
public:
  /** Standard class typedefs. */
  typedef BendingEnergyRegularizationFilter             Self;
  typedef ImageToImageFilter< TImage, TImage > Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef typename NeighborhoodIterator<TImage>::RadiusType RadiusType;

  typedef typename TImage::PixelType::RealValueType ImagePixelType;
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TImage::ImageDimension);
  
  typedef TImage                                ImageType;
  typedef typename ImageType::Pointer           TImagePointer;
  
  typedef typename TImage::PixelType                  ScalarPixelType;
  typedef typename itk::Image<ScalarPixelType,ImageDimension> ScalarImageType;
  typedef typename ScalarImageType::Pointer  ScalarImagePointer;

 
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
 
  /** Run-time type information (and related methods). */
  itkTypeMacro(BendingEnergyRegularizationFilter, ImageToImageFilter);
 
  void SetVelocity(TImagePointer V);
  
  void SetCurrentUpdate(TImagePointer U); 

  void SetSqLaplacian(TImagePointer Lp );

  void SetRadius(RadiusType R);

  void SetFactor(ImagePixelType factor);


protected:
  BendingEnergyRegularizationFilter();
  ~BendingEnergyRegularizationFilter(){}
 
  virtual void BeforeThreadedGenerateData();
 
  /** Does the real work. */
  virtual void ThreadedGenerateData(const typename Superclass::OutputImageRegionType& outputRegionForThread, int threadId);


 
private:
  BendingEnergyRegularizationFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
  
  TImagePointer         m_Velocity;
  TImagePointer         m_CurrentUpdate;
  TImagePointer         m_SqLaplacian;
  RadiusType            m_Radius;
  ImagePixelType        m_Factor;

};
} //namespace ITK
 
 
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBendingEnergyRegularizationFilter.txx"
#endif
 
#endif // __itkBendingEnergyRegularizationFilter_h
