
# AcousticLev
AcousticLev is a C++ header only library to calculate ultrasonic phase array/s settings to levitate particles using automatic differentiation (AD).

The library is an implementation of the algorithm found in this [paper] (http://www.nature.com/ncomms/2015/151027/ncomms9661/full/ncomms9661.html). 
All that is needed to use the library is to put the folder acousticHologram in the same folder with the C++ file that uses it. An example of the library being used is particleLev.cpp. One can use the python file showData.py to see the results of the calculations. The python file also writes a .vtk file that can be viewed in [ParaView] (http://www.paraview.org), that shows the magnitude of the pressure, potential and force fields due to the phased array/s.

The library uses the following dependencies:
* [Eigen] (http://eigen.tuxfamily.org) (for matrix operations)
* [NLopt] (http://ab-initio.mit.edu/nlopt/) (for optimization)
* [MayaVi] (http://code.enthought.com/projects/mayavi/) (in the python file for visualization)

The library uses forward mode automatic differentiation to quickly compute the pressure field at any point due to the phased array/s without any approximation errors. It uses modified [FADBAD++] (http://www.fadbad.com/fadbad.html) files that include Bessel functions to model the sound transducer's pressure field.

More information on the libabry can be found in the manual (AcousticLev.pdf).
