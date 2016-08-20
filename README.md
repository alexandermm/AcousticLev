# AcousticLev
C++ header only library to calculate ultrasonic phase array/s settings to levitate particles using automatic differentiation.

The library is an implementation of the algorithm found in this [paper] (http://www.nature.com/ncomms/2015/151027/ncomms9661/full/ncomms9661.html). 
All that is needed to use the library is to put the folder acousticHologram in the same folder with the C++ file that uses it. An example of the library being used is particleLev.cpp. One can use the python file showData.py to see the results of the calculations. The python file also writes a .vtk file that can be viewed in [ParaView] (http://www.paraview.org), that shows the magnitude of the pressure, potential and force fields due to the phased array/s.

The library uses the following dependencies:
* [Eigen] (http://eigen.tuxfamily.org) (for matrix operations)
* [NLopt] (http://ab-initio.mit.edu/nlopt/) (for optimization)
* [MayaVi] (http://code.enthought.com/projects/mayavi/) (in the python file for visualization)

The library uses forward mode automatic differentiation to quickly compute the pressure field at any point due to the phased array/s without any approximation errors. It uses modified [FADBAD++] (http://www.fadbad.com/fadbad.html) files that include Bessel functions to model the sound transducer's pressure field.



## How to use the library (following the code in particleLev.cpp)
* The first step is to declare all the geometric, physical and optimization constants used for the calculation. The variable sf is used to scale the physical values so that they are all closer in magnitude. opt_scale is used to scale the optimization function used by NLopt so that the solution converges faster. The w Eigen matrix is used to weight different directions 
(x,y or z) so that the shape of the sound field around the levitation point is adecuate (see paper). When using one phased array, the w should be weighted towards the direction the array is pointing to.

* The second step is to declare rectangular phased array objects. When the arrays are declared, their center is at the origin with the each of the sound transducers in the array pointing in the same direction as the z-axis. Each declared array can be translated or rotated. After moving all the needed arrays to their desired locations, one adds them together (using the overloaded + sign) to form a single fixed array. This new array now contains all the transducer locations and normal vectors which are used to calculate the transducer settings (see paper). The class declaration is in the file transducerArray.hpp.

* The third step is to declare the desired levitation points, and the Mj Eigen matrix. The Mj matrix stores all the sound pressure derivatives (up to third order if needed) with respect to a fixed levitation point, for each of the transducers in the final single phased array. This values are used thorughout the calculations needed to optimize the array and to compute field values (pressure, potential or force) at a given point (see paper). The values are exact and they are calculated using forward mode automatic differentiation using an analytic model of the sound field for each transducer and FADBAD++. The functions doing this are in the file transducerDerivatives.cpp. 

* The final step is to call NLopt and the optimization function and its gradient in optFunction.hpp to optimize the phases of each of the transducers. This function is given as U_aa (eq. 10) in the paper. This function should have a 2 in front of the second term and was derived by using (mag(p))^2 = p.p instead of mag(p) = p.p (eq. 11) as shown in the paper.

* For post-processing, one can calculate the imaginary pressure, the potential and force due to the pressure field produced by all of the transducers to check that all the forces are pointing towards the desired levitation point. This data can be accessed by the showData.py script to produce a .vtk file. The functions that can do this are in the file soundField.hpp.

* Also one can compute particle paths that a particle would take from a given point to a levitation point. This includes air resistance. This is done using a [4th order Rugge-Kutta scheme] (http://lpsa.swarthmore.edu/NumInt/NumIntFourth.html) with a given fixed time step. The functions doing this are found in particlePath.hpp.



## Details about how the automatic differentiation is done
The [automatic differentiation] (https://en.wikipedia.org/wiki/Automatic_differentiation) is done in forward mode three times to get the third order derivatives needed for the optimization function, using the forward mode file of FADBAD++ (fadiff.h together with the main header file fadbad.h). In theory, one should perform the differentiation in reverse mode the first time since the number of dependent variables (1) is lower than the number of independadent variables (3). That is the behabiour for a large number of independent variables. However, when one programs both alternatives, using the reverse mode differentiation once first and then forward differentiation once (for calculating forces for example) is 35% slower than using forward differentiation twice. This is due to the more complex implementation needed for reverse differentiation, and the fact that the ratio between dependent and independent variables is low. By declaring the number of variables being differentitated each time forward mode is used (using the stack-based forward method implementation of FADBAD++), the run time was reduced a further 33% (again for calculating forces). The real and complex components of the pressure field are differentiated separately so one can use FADBAD++, since the transducer sound field (eq. 7) can be separated into real and complex terms. 



## Details on how the functions in the paper are derived
The U_aa function can be derived from the Gorkov potential (eq. 3), by noting that (mag(p))^2 = p.p, and by appling the product rule on eq. 12. This yields for the derivative with respect to variable c:

(a.b)_c = (a.b_c) + (a_c.b)

Substituting this equation twice in (eq. 3) yields the U_aa function.

