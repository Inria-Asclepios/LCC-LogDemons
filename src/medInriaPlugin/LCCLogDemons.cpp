/*=========================================================================

 medInria

 Copyright (c) INRIA 2013. All rights reserved.
 See LICENSE.txt for details.
 
  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

=========================================================================*/

#include "LCCLogDemons.h"

#include <dtkCore/dtkAbstractData.h>
#include <dtkCore/dtkAbstractDataFactory.h>
#include <dtkCore/dtkAbstractProcessFactory.h>

// /////////////////////////////////////////////////////////////////
//
// /////////////////////////////////////////////////////////////////

#include "itkImageRegistrationMethod.h"

#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"


#include "time.h"

#include <rpiLCClogDemons.hxx>
#include <rpiCommonTools.hxx>

// /////////////////////////////////////////////////////////////////
// LCCLogDemonsPrivate
// /////////////////////////////////////////////////////////////////
typedef itk::Image< float, 3 >  RegImageType;
typedef double TransformScalarType;

class LCCLogDemonsPrivate
{
public:
    LCCLogDemons * proc;
    template <class PixelType>
    int update();
    template <typename PixelType>
    bool writeTransform(const QString& file);
    
    rpi::LCClogDemons< RegImageType, RegImageType, float > * registrationMethod;
    rpi::LCClogDemons< RegImageType, RegImageType, float >::UpdateRule updateRule;
    rpi::LCClogDemons< RegImageType, RegImageType, float >::GradientType gradientType;
    std::vector<unsigned int> iterations;
    double maxStepLength, similarityCriteriaSigma, sigmaI, updateFieldSigma, velocityFieldSigma;
    unsigned int BCHExpansion;
    bool verbose, boundaryCheck, useHistogramMatching;
};

// /////////////////////////////////////////////////////////////////
// LCCLogDemons
// /////////////////////////////////////////////////////////////////

LCCLogDemons::LCCLogDemons() : itkProcessRegistration(), d(new LCCLogDemonsPrivate)
{
    d->proc = this;
    d->registrationMethod = NULL;
}

LCCLogDemons::~LCCLogDemons()
{
    d->proc = NULL;

    
    if (d->registrationMethod)
        delete d->registrationMethod;

    d->registrationMethod = NULL;
    
    delete d;
    d = 0;
}

bool LCCLogDemons::registered()
{
    return dtkAbstractProcessFactory::instance()->registerProcessType("LCCLogDemons",
                                                                      createLCCLogDemons);
}

void LCCLogDemons::setUpdateRule(int rule)
{
    d->updateRule = (rpi::LCClogDemons< RegImageType, RegImageType, float >::UpdateRule)rule;
}

void LCCLogDemons::setVerbosity(bool verbose)
{
    d->verbose = verbose;
}

void LCCLogDemons::setNumberOfIterations(std::vector<unsigned int> iterations)
{
    d->iterations = iterations;
}

void LCCLogDemons::setMaximumUpdateStepLength(double value)
{
    d->maxStepLength = value;
}

void LCCLogDemons::setGradientType(int value)
{
    d->gradientType = (rpi::LCClogDemons< RegImageType, RegImageType, float >::GradientType)value;;
}

void LCCLogDemons::setUseHistogramMatching(bool flag)
{
    d->useHistogramMatching = flag;
}

void LCCLogDemons::setBoundaryCheck(bool flag)
{
    d->boundaryCheck = flag;
}

void LCCLogDemons::setSigmaI(double sigmaI)
{
    d->sigmaI = sigmaI;
}

void LCCLogDemons::setSimilarityCriteriaSigma(double std)
{
    d->similarityCriteriaSigma = std;
}

void LCCLogDemons::setUpdateFieldSigma(double sigma)
{
    d->updateFieldSigma = sigma;
}

void LCCLogDemons::setVelocityFieldSigma(double sigma)
{
    d->velocityFieldSigma = sigma;
}

void LCCLogDemons::setNumberOfTermsBCHExpansion(unsigned int number)
{
    d->BCHExpansion = number;
}

QString LCCLogDemons::description() const
{
    return "LCCLogDemons";
}



// /////////////////////////////////////////////////////////////////
// Templated Version of update
// /////////////////////////////////////////////////////////////////


template <typename PixelType>
int LCCLogDemonsPrivate::update()
{
    registrationMethod = new rpi::LCClogDemons<RegImageType,RegImageType,float> ();
    proc->fixedImage().GetPointer()->SetSpacing(proc->movingImages()[0].GetPointer()->GetSpacing());
    typedef itk::Image<unsigned short, 3> OldImageType;
    typename OldImageType::Pointer bloup = (OldImageType*)proc->fixedImage().GetPointer();
    typename OldImageType::Pointer bloup2 = (OldImageType*)proc->movingImages()[0].GetPointer();
    typedef itk::CastImageFilter <OldImageType, RegImageType> CastFilter;
    typename CastFilter::Pointer castFilter1 = CastFilter::New();
    typename CastFilter::Pointer castFilter2 = CastFilter::New();

    castFilter1->SetInput(bloup);
    castFilter2->SetInput(bloup2);
    typename RegImageType::Pointer input1 = castFilter1->GetOutput();
    typename RegImageType::Pointer input2 = castFilter2->GetOutput();

    registrationMethod->SetFixedImage((const RegImageType*) input1);
    registrationMethod->SetMovingImage((const RegImageType*) input2);
    
    registrationMethod->SetUpdateRule(updateRule);
    registrationMethod->SetVerbosity(verbose);

    //Log
    registrationMethod->SetMaximumUpdateStepLength(maxStepLength);
    registrationMethod->SetGradientType(gradientType);
    registrationMethod->SetUseHistogramMatching(useHistogramMatching);

    registrationMethod->SetNumberOfIterations(iterations);
    registrationMethod->SetBoundaryCheck(boundaryCheck);
    registrationMethod->SetSimilarityCriteriaStandardDeviation(similarityCriteriaSigma);
    registrationMethod->SetSigmaI(sigmaI);
    registrationMethod->SetUpdateFieldStandardDeviation(updateFieldSigma);
    registrationMethod->SetStationaryVelocityFieldStandardDeviation(velocityFieldSigma);
    registrationMethod->SetNumberOfTermsBCHExpansion(BCHExpansion);
    
    // Run the registration
    time_t t1 = clock();
    try {
        registrationMethod->StartRegistration();
    }
    catch( std::exception & err )
    {
        qDebug() << "ExceptionObject caught ! (startRegistration)" << err.what();
        return 1;
    }
    
    time_t t2 = clock();
    
    qDebug() << "Elasped time: " << (double)(t2-t1)/(double)CLOCKS_PER_SEC;
    
    typedef itk::ResampleImageFilter< RegImageType,RegImageType, float>    ResampleFilterType;
    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetTransform(registrationMethod->GetTransformation());
    resampler->SetInput((const RegImageType*)proc->movingImages()[0].GetPointer());
    resampler->SetSize( proc->fixedImage()->GetLargestPossibleRegion().GetSize() );
    resampler->SetOutputOrigin( proc->fixedImage()->GetOrigin() );
    resampler->SetOutputSpacing( proc->fixedImage()->GetSpacing() );
    resampler->SetOutputDirection( proc->fixedImage()->GetDirection() );
    resampler->SetDefaultPixelValue( 0 );

    
    try {
        resampler->Update();
    }
    catch (itk::ExceptionObject &e) {
        qDebug() << e.GetDescription();
        return 1;
    }
    
    itk::ImageBase<3>::Pointer result = resampler->GetOutput();
    result->DisconnectPipeline();
    
    if (proc->output())
        proc->output()->setData (result);
    return 0;
}

int LCCLogDemons::update(itkProcessRegistration::ImageType imgType)
{
    if(fixedImage().IsNull() || movingImages()[0].IsNull())
        return 1;

    return d->update<float>();
}


template <typename PixelType>
bool LCCLogDemonsPrivate::writeTransform(const QString& file)
{

    
    if (registrationMethod)
    {
        try{
            rpi::writeDisplacementFieldTransformation<float, 3>(registrationMethod->GetTransformation(),
                                                                file.toStdString());
        }
        catch (std::exception)
        {
            return false;
        }
        return true;
    }
    else
    {
        return false;
    }
    
}

bool LCCLogDemons::writeTransform(const QString& file)
{
    if(d->registrationMethod == NULL)
        return 1;
    
    return d->writeTransform<float>(file);
}

itk::Transform<double,3,3>::Pointer LCCLogDemons::getTransform()
{
    return NULL;
}

QString LCCLogDemons::getTitleAndParameters()
{
    return QString();
}

// /////////////////////////////////////////////////////////////////
// Type instanciation
// /////////////////////////////////////////////////////////////////

dtkAbstractProcess *createLCCLogDemons()
{
    return new LCCLogDemons;
}
