#ifndef __itkVectorRegularizationFilter_txx
#define __itkVectorRegularizationFilter_txx
 
#include "itkVectorRegularizationFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNeighborhoodIterator.h"
#include"itkNeighborhoodInnerProduct.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkImageRegionIterator.h"

#include "itkImageFileWriter.h"

namespace itk
{
template<class TRealVectorImage>
VectorRegularizationFilter<TRealVectorImage>::VectorRegularizationFilter()
{
  m_ShiftFilter  = ShiftFilterType::New();
  m_FFTFilter    = FFTFilterType::New();
  m_FFTInvFilter = FFTInvFilterType::New();
  
  m_ComplexComposer  = ComposeComplexVectorFilterType::New();
 
  m_ComplexExtractor = ExtractorComplexFilterType::New();
  m_ComplexComposer  = ComposeComplexVectorFilterType::New();

  m_RealExtractor = ExtractorFilterType::New();
  m_RealComposer  = ComposeVectorType::New();
 
  for (int i=0;i<ImageDimension;++i)
    {
     m_ComplexImageArray[i] = ComplexImageType::New();
     m_RealImageArray[i] = RealImageType::New();
      }

  m_factor[0]=1e-3;
  m_factor[1]=1e-3;
  m_ComplexSolution = ComplexVectorImageType::New();

}
 
 
template<class TRealVectorImage>
void VectorRegularizationFilter<TRealVectorImage>::SetFactor(float factor[])
  {  
    for (int i=0;i<2;++i)
      m_factor[i]=factor[i];
  }
 
template<class TRealVectorImage>
void VectorRegularizationFilter<TRealVectorImage>
::GenerateData()
{
  for (int i=0;i<ImageDimension;++i)
    {
     m_RealExtractor->SetInput(this->GetInput());
     m_RealExtractor->SetIndex(i);
     m_RealExtractor->Update();
     m_FFTFilter->SetInput(m_RealExtractor->GetOutput());

     if (this->GetOutput()->GetLargestPossibleRegion().GetSize()[0]%2==1)
      {
    //   m_FFTInvFilter->SetActualXDimensionIsOdd(true);
      }
     m_FFTFilter->Update();
     m_ShiftFilter->SetInput(m_FFTFilter->GetOutput());
     m_ShiftFilter->SetInverse(false);
     m_ShiftFilter->Update();
     m_ComplexImageArray[i]=m_ShiftFilter->GetOutput();
     m_ComplexImageArray[i]->DisconnectPipeline();
     m_ComplexComposer->SetInput(i,m_ComplexImageArray[i]);
    }

  m_ComplexComposer->Update();
  m_ComplexSolution=m_ComplexComposer->GetOutput();  
 
  for (int i=0;i<ImageDimension;++i)
          m_dimensions[i]=m_ComplexSolution->GetLargestPossibleRegion().GetSize()[i];

  typename itk::ImageRegionIterator<ComplexVectorImageType> out(m_ComplexSolution,m_ComplexSolution->GetLargestPossibleRegion());

  typename itk::Index<ImageDimension> Index;
  typename ComplexImageType::PixelType  Vector;

  float SqPosition;
  float scaling;

    for( out.GoToBegin(); ! out.IsAtEnd(); ++out )
      {
        SqPosition=0;
        Index= out.GetIndex();     
        for (int i=0;i<ImageDimension;++i)
          SqPosition+=(Index[i]-m_dimensions[i]/2)*(Index[i]-m_dimensions[i]/2);

        scaling=1+m_factor[0]*SqPosition + m_factor[1]*SqPosition*SqPosition;     

        out.Set(out.Get()*(1/scaling));
 
      }


  
  for (int i=0;i<ImageDimension;++i)
    {
     m_ComplexExtractor->SetInput(m_ComplexSolution);
     m_ComplexExtractor->SetIndex(i);
     m_ComplexExtractor->Update();
     m_ShiftFilter->SetInput(m_ComplexExtractor->GetOutput());
     m_ShiftFilter->SetInverse(true);
     m_ShiftFilter->Update();
     m_FFTInvFilter->SetInput(m_ShiftFilter->GetOutput());

     if (!m_ShiftFilter->GetOutput()->GetLargestPossibleRegion().GetSize()[0]%2)
      {
      // m_FFTInvFilter->SetActualXDimensionIsOdd(true);
      }

     m_FFTInvFilter->Update();
     m_RealImageArray[i]=m_FFTInvFilter->GetOutput();
     m_RealImageArray[i]->DisconnectPipeline();
     m_RealComposer->SetInput(i,m_RealImageArray[i]);
    }   

 m_RealComposer->Update();
 this->GraftOutput(m_RealComposer->GetOutput());
}
 
}// end namespace
 
#endif //__itkVectorRegularizationFilter_txx
