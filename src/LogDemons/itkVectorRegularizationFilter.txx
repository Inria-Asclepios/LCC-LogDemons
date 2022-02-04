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

      try
      {
          m_RealExtractor->Update();
      }
      catch (itk::ExceptionObject & err)
      {
          std::cerr << "ExceptionObject caught !" << std::endl;
          std::cerr << err << std::endl;
      }
     m_FFTFilter->SetInput(m_RealExtractor->GetOutput());

      try
      {
          m_FFTFilter->Update();
      }
      catch (itk::ExceptionObject & err)
      {
          std::cerr << "ExceptionObject caught !" << std::endl;
          std::cerr << err << std::endl;
      }
      m_ShiftFilter->SetInput(m_FFTFilter->GetOutput());
      m_ShiftFilter->SetInverse(false);

      try
      {
          m_ShiftFilter->Update();
      }
      catch (itk::ExceptionObject & err)
      {
          std::cerr << "ExceptionObject caught !" << std::endl;
          std::cerr << err << std::endl;
      }
      m_ComplexImageArray[i]=m_ShiftFilter->GetOutput();
      m_ComplexImageArray[i]->DisconnectPipeline();
      m_ComplexComposer->SetInput(i,m_ComplexImageArray[i]);
    }

  try
  {
      m_ComplexComposer->Update();
  }
  catch (itk::ExceptionObject & err)
  {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
  }
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

      try
      {
          m_ComplexExtractor->Update();
      }
      catch (itk::ExceptionObject & err)
      {
          std::cerr << "ExceptionObject caught !" << std::endl;
          std::cerr << err << std::endl;
      }
      m_ShiftFilter->SetInput(m_ComplexExtractor->GetOutput());
      m_ShiftFilter->SetInverse(true);

      try
      {
          m_ShiftFilter->Update();
      }
      catch (itk::ExceptionObject & err)
      {
          std::cerr << "ExceptionObject caught !" << std::endl;
          std::cerr << err << std::endl;
      }
      m_FFTInvFilter->SetInput(m_ShiftFilter->GetOutput());

      try
      {
          m_FFTInvFilter->Update();
      }
      catch (itk::ExceptionObject & err)
      {
          std::cerr << "ExceptionObject caught !" << std::endl;
          std::cerr << err << std::endl;
      }
     m_RealImageArray[i]=m_FFTInvFilter->GetOutput();
     m_RealImageArray[i]->DisconnectPipeline();
     m_RealComposer->SetInput(i,m_RealImageArray[i]);
    }   

  try
  {
      m_RealComposer->Update();
  }
  catch (itk::ExceptionObject & err)
  {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
  }
 this->GraftOutput(m_RealComposer->GetOutput());
}
 
}// end namespace
 
#endif //__itkVectorRegularizationFilter_txx
