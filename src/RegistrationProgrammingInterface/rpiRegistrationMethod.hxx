#ifndef _RPI_REGISTRATION_METHOD_HXX_
#define _RPI_REGISTRATION_METHOD_HXX_


#include <itkImage.h>
#include <itkTransform.h>

#include "rpiObserver.hxx"


// Namespace RPI : Registration Programming Interface
namespace rpi
{


/**
 * Abstract class implementing a generic image registration method.
 *   TFixedImage           Type of the fixed image. Should be "itk::Image".
 *   TMovingImage          Type of the moving image. Should be "itk::Image".
 *   TTransformScalarType  Type of the scalar value of the transformation. Should be "float" or "double".
 *
 * @author Vincent Garcia
 * @date   2010/04/21
 */
template < class TFixedImage, class TMovingImage, class TTransformScalarType=double >
class RegistrationMethod {


public:

    typedef itk::Image< typename TFixedImage::PixelType, TFixedImage::ImageDimension >
            FixedImageType;

    typedef typename FixedImageType::ConstPointer
            FixedImageConstPointerType;

    typedef itk::Image< typename TMovingImage::PixelType, TMovingImage::ImageDimension >
            MovingImageType;

    typedef typename MovingImageType::ConstPointer
            MovingImageConstPointerType;

    typedef itk::Transform< TTransformScalarType, TFixedImage::ImageDimension, TMovingImage::ImageDimension >
            TransformType;

    typedef typename TransformType::Pointer
            TransformPointerType;

    typedef Observer<TFixedImage, TMovingImage, TTransformScalarType>
            ObserverType;

    /**
     * Registration status.
     */
    enum RegistrationStatus {
        REGISTRATION_STATUS_PROCESSING, /** Processing registration */
        REGISTRATION_STATUS_STOP        /** Stop: no activity       */
    };


protected:

    /**
     * Fixed image
     */
    FixedImageConstPointerType   m_fixedImage;

    /**
     * Moving image
     */
    MovingImageConstPointerType  m_movingImage;

    /**
     * Transformation estimated by the registration method
     */
    TransformPointerType         m_transform;

    /**
     * Observer
     */
    std::vector<ObserverType*>   m_observers;

    /**
     * Registration status
     */
    RegistrationStatus           m_registrationStatus;


public :

    /**
     * Class constructor
     */
    RegistrationMethod(void);

    /**
     * Class destructor
     */
    virtual ~RegistrationMethod(void);

    /**
     * Gets the fixed image.
     * @return  fixed image
     */
    FixedImageConstPointerType   GetFixedImage(void) const;

    /**
     * Sets the fixed image.
     * @param  image  fixed image
     */
    void                         SetFixedImage(const FixedImageType * image);

    /**
     * Gets the moving image.
     * @return  moving image
     */
    MovingImageConstPointerType  GetMovingImage(void) const;

    /**
     * Sets the moving image.
     * @param  image  moving image
     */
    void                         SetMovingImage(const MovingImageType * image);

    /**
     * Gets the transformation estimated by the registration method.
     * If the computeRegistration function has not been called yet, the return transformation
     * is equal to the initial transformation if set, otherwise, the transformation is the identity.
     * @return  estimated transformation
     */
    TransformPointerType         GetTransformation(void) const;

    /**
     * Attachs an observer to the registration method
     * @param observer observer
     */
    void                         AttachObserver(ObserverType * observer);

    /**
     * Notify observer that the registration status changed.
     */
    void                         Notify(void);

    /**
     * Gets the registration status.
     * @return registration status
     */
    RegistrationStatus           GetRegistrationStatus(void);

    /**
     * Sets the registration status.
     * @param status registration status
     */
    void                         SetRegistrationStatus(RegistrationStatus status);

    /**
     * Registers the moving image on the fixed image.
     */
    virtual void                 StartRegistration(void) = 0;

};


} // End of namespace


/** Add the source code file (template) */
#include "rpiRegistrationMethod.cxx"


#endif // _RPI_REGISTRATION_METHOD_HXX_
