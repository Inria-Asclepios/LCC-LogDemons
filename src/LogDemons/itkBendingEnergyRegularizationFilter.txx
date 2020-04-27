#ifndef __itkBendingEnergyRegularizationFilter_txx
#define __itkBendingEnergyRegularizationFilter_txx
 
#include "itkBendingEnergyRegularizationFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"

namespace itk
{
template<class TImage>
BendingEnergyRegularizationFilter<TImage>::BendingEnergyRegularizationFilter()
{
  m_Velocity=ImageType::New();
  m_CurrentUpdate=ImageType::New();
  m_SqLaplacian=ImageType::New();
  
  for (int i=0;i<ImageDimension;++i)
     m_Radius[i]=5;
  
  m_Factor=120.0;
}
 
 
template<class TImage>
void
BendingEnergyRegularizationFilter<TImage>::SetVelocity(TImagePointer V)
{  m_Velocity=V;
    }

template<class TImage>
void
BendingEnergyRegularizationFilter<TImage>::SetCurrentUpdate(TImagePointer U)
{  m_CurrentUpdate=U;
    }

template<class TImage>
void
BendingEnergyRegularizationFilter<TImage>::SetSqLaplacian(TImagePointer Sq)
{  m_SqLaplacian=Sq;
    }

template<class TImage>
void
BendingEnergyRegularizationFilter<TImage>::SetRadius(RadiusType R)
{  m_Radius=R;
    }

 template<class TImage>
void
BendingEnergyRegularizationFilter<TImage>::SetFactor(ImagePixelType factor)
{
 m_Factor=factor;
}

template<class TImage>
void BendingEnergyRegularizationFilter<TImage>::BeforeThreadedGenerateData()
{
    }



template<class TImage>
void BendingEnergyRegularizationFilter<TImage>
::DynamicThreadedGenerateData(const typename Superclass::OutputImageRegionType& outputRegionForThread)
{
  typename TImage::ConstPointer input = this->GetInput();
  typename TImage::Pointer output = this->GetOutput();
 
  typedef itk::NeighborhoodAlgorithm
    ::ImageBoundaryFacesCalculator< TImage > FaceCalculatorType;  
 
  FaceCalculatorType faceCalculator;
  typename FaceCalculatorType::FaceListType faceList;
  faceList = faceCalculator(output, outputRegionForThread,
                            m_Radius);
    
  typename FaceCalculatorType::FaceListType::iterator fit;

  typename TImage::IndexType indexV;
  typename TImage::PixelType Vector;

  for ( fit=faceList.begin(); fit != faceList.end(); ++fit)
    {
     itk::ImageRegionIterator<TImage> out(output, *fit);
     itk::ConstNeighborhoodIterator<TImage> it(m_Radius, input, *fit);

   
     for( it.GoToBegin(), out.GoToBegin(); ! it.IsAtEnd(); ++it, ++out )
      {
       indexV=it.GetIndex(); 	
       float norm=m_SqLaplacian->GetPixel(indexV).GetSquaredNorm();
       for(int i=0;i<ImageDimension;++i)
        Vector[i]=0.5* (m_Factor*(m_Velocity->GetPixel(indexV)[i]-m_CurrentUpdate->GetPixel(indexV)[i])-m_SqLaplacian->GetPixel(indexV)[i])/(norm+m_Factor/4)  ;
        out.Set(Vector);
 
      }
   }
}
 
}// end namespace
 
#endif //__itkBendingEnergyRegularizationFilter_txx
