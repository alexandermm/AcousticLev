# AcousticLev
C++ header only library to calculate ultrasonic phase array/s settings to levitate particles using automatic differentiation.

The library is an implementation of the algorithm found in this [paper] (http://www.nature.com/ncomms/2015/151027/ncomms9661/full/ncomms9661.html) 
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

* 


