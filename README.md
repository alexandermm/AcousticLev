# AcousticLev
C++ header only library to calculate ultrasonic phase array setting to levitate particles using automatic differentiation.

The library is an implementation of the algorithm found in this [paper] (http://www.nature.com/ncomms/2015/151027/ncomms9661/full/ncomms9661.html) 
All that is needed to use the library is to put the folder acousticHologram in the same folder with the C++ file that uses it. An example of the library being used is particleLev.cpp. One can use the python file showData.py to see the results of the calculations.

The library uses the following dependencies:
* Eigen (for matrix operations)
* NLopt (for optimization)
* MayaVi (in the python file for visualization)

The library uses forward mode automatic differentiation to quickly compute the pressure field at any point due to the phased array/s without any approximation errors. It uses modified FADBAD++ files that include Bessel functions to model the sound transducer's pressure field.
