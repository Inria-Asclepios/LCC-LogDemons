#ifndef _RPI_OBSERVER_HXX_
#define _RPI_OBSERVER_HXX_


// Namespace RPI : Registration Programming Interface
namespace rpi
{


// Forward declaration to avoid circular references.
template < class TFixedImage, class TMovingImage, class TTransformScalarType >
class RegistrationMethod;




/**
 * Observer for RegistrationMethod.
 * Because it contains an instance of the RegistrationMethod class, the Observer class is templated.
 * Because an instance of the class Observer is supposed to be attached to an instance of the class
 * Registration method,
 * Template must fit the template of the registration method to whom the observe is attached.
 *
 *   TFixedImage           Type of the fixed image.
 *
 *   TMovingImage          Type of the moving image.
 *
 *   TTransformScalarType  Type of the transformation parameters.
 *
 * @author Vincent Garcia
 * @date   2011/06/24
 */
template < class TFixedImage, class TMovingImage, class TTransformScalarType=double >
class Observer
{


public:

    /**
     * Class constructor.
     */
    Observer();

    /**
     * Class destructor.
     */
    virtual ~Observer(void);

    /**
     * Updates the observer.
     */
    virtual void Update(void) = 0;

    /**
     * Gets registration method.
     * @param registration method
     */
    RegistrationMethod<TFixedImage, TMovingImage, TTransformScalarType> * GetRegistrationMethod(void);

    /**
     * Sets registration method.
     * @return registration method
     */
     void SetRegistrationMethod(RegistrationMethod<TFixedImage, TMovingImage, TTransformScalarType> * method);


protected:

    /**
     * Registration method.
     */
    RegistrationMethod<TFixedImage, TMovingImage, TTransformScalarType> * m_registrationMethod;

};

} // End of namespace


/** Add the source code file (template) */
#include "rpiObserver.cxx"

#endif // _RPI_OBSERVER_HXX_
