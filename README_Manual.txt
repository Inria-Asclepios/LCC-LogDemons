----------Overview------------

LCClogDemons provides an implementation of the LCC log-Demons algorithm for the diffeomorphic registration parameterized by stationary velocity fields presented in the research article 

LCC-Demons: A robust and accurate symmetric diffeomorphic registration algorithm. Lorenzi M, Ayache N, Frisoni GB, Pennec X. Neuroimage. 2013 Nov 1;81:470-83. 

The software inherits from the Symmetric Log-Demons algorithm implemented by Tom Vercauteren and Florence Dru (Asclepios INRIA and Mauna Kea Technologies - http://www.insight-journal.org/browse/publication/644). 
LCClogDemons implements the symmetric Local Correlation Coefficient (LCC) as a similarity measure, and thus it is unbiased with respect to local linear intensity bias of the images. Moreover, in addition to the standard Gaussian convolution for the regularization of the stationary velocity field, the current implementation enables Harmonic and/or Curvature penalization.  LCClogDemons is suited for both inter and intra-subject registration, and compares positively with respect to state-of-art methods.

For further informations and feedbacks please contact the following address:
 marco.lorenzi@inria.fr

----------Installation------------

LCClogDemons requires ITK4 (current version 4.4.1) compiled by linking against fftw libraries (current version 3.3.3).
This means that ITK must be compiled with ccmake by enabling the flags ITK_USE_SYSTEM_FFTW, ITK_USE_FFTWD, and ITK_USE_FFTWF.

To ease the process we provide here a bash installer (setup.sh) which should automatically download and locally install the required dependencies. On the terminal just type ./setup.sh

----------Usage------------------

LCClogDemons enables diffeormorphic registration of images based on different
similarity criteria and regularization methods.

***Similarity metric: LCC 

LCC is enabled by setting the option -r 2
LCC kernel size (in mm) is defined by the option -C <kernel_size>. Typical
values range from 3 to 6, but may vary depending on the application.

By using the LCC you may want to set the following parameter:
-S <similarity_over_regularization_trade-off>

***Similarity metric: SSD

SSD is enabled by setting the option -r 1 (SSD-based symmetric log-domain -suggested), or -r 0 (SSD-based forward log-domain).

By using the SSD you may want to set the following parameters:
-g <gradient_type>
-l <max_step_length>
--use-histogram-matching 

***Regularization scheme: Harmonic + Curvature energy 

Harmonic + Curvature regularization is enabled by the option -R 1. The weight
for the penalization of the harmonic energy is set by the option -x
<harmonic_weight>, while the weight for the penalization of the curvature
energy is set by the option -b <curvature_weight>

***Regularization scheme: Gaussian convolution

Gaussian convolution is the classical log-Demons regularization scheme. It is
equivalent to the penalization of the harmonic energy associated to the velocity field. Gaussian convolution is enabled by the option -R 0. 
The kernel size for the smoothing of the velocity field is set by the option
-d <velocity_regularization>, while the smoothing of the update field (fluid-like regularization) is set by the option -u <update_regularization>.   
  
***Other options.
LCC-Demons enables the following main options
-V verbosity
-a number of iterations for the multiresolution scheme

------------Examples------------
Inter-subject registration of brain images.

cd LCClogDemons/example

***Typical parameters for the registration of brain images based on LCC and Curvature regularization are

../build/bin/rpiLCClogDemons -f Image1.nii.gz -m Image2.nii.gz -r 2 -R 1 -C 3 -a 30x20x10 -x 0 -b 0.2 -S 0.15 -V

***Typical parameters for the registration of brain images based on LCC and Gaussian regularization are

../build/bin/rpiLCClogDemons -f Image1.nii.gz -m Image2.nii.gz -r 2 -R 0 -C 3 -a 30x20x10 -u 0.5 -d 1.0 -S 0.15

------------log-Jacobian determinant computation------------

The present suite provides a tool for the consistent and robust computation of
the log-Jacobian determinant of a given stationary velocity field provided for
instance by the LCC-Demons.
The syntax is the following

./SVFLogJacobian -i <input_stationary_velocity_field_path> -o <output_logJacobian_determinant_path>

Other available options are

-s <scaling_factor> : optional scaling factor for the stationary velocity
field
-m <mask_image> : restricts into the mask the calculation of the step-size for the iterative computation 












 



 
