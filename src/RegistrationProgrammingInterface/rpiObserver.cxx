#ifndef _RPI_OBSERVER_CXX_
#define _RPI_OBSERVER_CXX_


#include "rpiObserver.hxx"


// Namespace RPI : Registration Programming Interface
namespace rpi
{


template <class TFixedImage, class TMovingImage, class TTransformScalarType>
Observer<TFixedImage,TMovingImage,TTransformScalarType>::
Observer(void)
{
    this->m_registrationMethod = NULL;
}


template <class TFixedImage, class TMovingImage, class TTransformScalarType>
Observer<TFixedImage,TMovingImage,TTransformScalarType>::
~Observer(void)
{
    // Do nothing
}


template <class TFixedImage, class TMovingImage, class TTransformScalarType>
RegistrationMethod<TFixedImage, TMovingImage, TTransformScalarType> *
Observer<TFixedImage,TMovingImage,TTransformScalarType>::
GetRegistrationMethod(void)
{
    return this->m_registrationMethod;
}


template <class TFixedImage, class TMovingImage, class TTransformScalarType>
void
Observer<TFixedImage,TMovingImage,TTransformScalarType>::
SetRegistrationMethod(RegistrationMethod<TFixedImage, TMovingImage, TTransformScalarType> * method)
{
    this->m_registrationMethod = method;
}


} // End of namespace

#endif // _RPI_OBSERVER_CXX_
