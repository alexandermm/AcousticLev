

style TYPE="text/css">
    code.has-jax {font: inherit; font-size: 100%; background: inherit; border: inherit;}
</style>
<script type="text/x-mathjax-config">
    MathJax.Hub.Config({
        tex2jax: {
            inlineMath: [['$','$'], ['\\(','\\)']],
            skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'] // removed 'code' entry
        }
    });
    MathJax.Hub.Queue(function() {
        var all = MathJax.Hub.getAllJax(), i;
        for(i = 0; i < all.length; i += 1) {
            all[i].SourceElement().parentNode.className += ' has-jax';
        }
    });
</script>
<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>



# AcousticLev
C++ header only library to calculate $x^2$ ultrasonic phase array/s settings to levitate particles using automatic differentiation.

The library is an implementation of the algorithm found in this [paper] (http://www.nature.com/ncomms/2015/151027/ncomms9661/full/ncomms9661.html). 
All that is needed to use the library is to put the folder acousticHologram in the same folder with the C++ file that uses it. An example of the library being used is particleLev.cpp. One can use the python file showData.py to see the results of the calculations. The python file also writes a .vtk file that can be viewed in [ParaView] (http://www.paraview.org), that shows the magnitude of the pressure, potential and force fields due to the phased array/s.

The library uses the following dependencies:
* [Eigen] (http://eigen.tuxfamily.org) (for matrix operations)
* [NLopt] (http://ab-initio.mit.edu/nlopt/) (for optimization)
* [MayaVi] (http://code.enthought.com/projects/mayavi/) (in the python file for visualization)

The library uses forward mode automatic differentiation to quickly compute the pressure field at any point due to the phased array/s without any approximation errors. It uses modified [FADBAD++] (http://www.fadbad.com/fadbad.html) files that include Bessel functions to model the sound transducer's pressure field.



## How to use the library (following the code in particleLev.cpp)
* The first step is to declare all the geometric, physical and optimization constants used for the calculation. The variable sf is used to scale the physical values so that they are all closer in magnitude. opt_scale is used to scale the optimization function used by NLopt so that the solution converges faster. The w Eigen matrix is used to weight different directions 
(x,y or z) so that the shape of the sound field around the levitation point is adequate (see paper). When using one phased array, the w should be weighted towards the direction the array is pointing to.

* The second step is to declare rectangular phased array objects. When the arrays are declared, their center is at the origin with the each of the sound transducers in the array pointing in the same direction as the z-axis. Each declared array can be translated or rotated. After moving all the needed arrays to their desired locations, one adds them together (using the overloaded + sign) to form a single fixed array. This new array now contains all the transducer locations and normal vectors which are used to calculate the transducer settings (see paper). The class declaration is in the file transducerArray.hpp.

* The third step is to declare the desired levitation points, and the Mj Eigen matrix. The Mj matrix stores all the sound pressure derivatives (up to third order if needed) with respect to a fixed levitation point, for each of the transducers in the final single phased array. This values are used throughout the calculations needed to optimize the array and to compute field values (pressure, potential or force) at a given point (see paper). The values are exact and they are calculated using forward mode automatic differentiation using an analytic model of the sound field for each transducer and FADBAD++. The functions doing this are in the file transducerDerivatives.cpp. 

* The final step is to call NLopt and the optimization function and its gradient in optFunction.hpp to optimize the phases of each of the transducers. This function is given as U_aa (eq. 10) in the paper. This function should have a 2 in front of the second term and was derived by using (mag(p))^2 = p.p instead of mag(p) = p.p (eq. 11) as shown in the paper.

* For post-processing, one can calculate the imaginary pressure, the potential and force due to the pressure field produced by all of the transducers to check that all the forces are pointing towards the desired levitation point. This data can be accessed by the showData.py script to produce a .vtk file. The functions that can do this are in the file soundField.hpp.

* Also one can compute particle paths that a particle would take from a given point to a levitation point. This includes air resistance. This is done using a [4th order Rugge-Kutta scheme] (http://lpsa.swarthmore.edu/NumInt/NumIntFourth.html) with a given fixed time step. The functions doing this are found in particlePath.hpp.


The following sections go into more detail on how the library is designed and what equations are used. One does not need to read them in order to use the library. 



## Details about how the automatic differentiation is done
The [automatic differentiation] (https://en.wikipedia.org/wiki/Automatic_differentiation) is done in forward mode three times to get the third order derivatives needed for the optimization function, using the forward mode file of FADBAD++ (fadiff.h together with the main header file fadbad.h). **The advantage of this method compared to finite differences is that the result is exact up to machine precision, and there is no error due to the step size or the frequency of the sound from the transducer being modelled.** In theory, one should perform the differentiation in reverse mode the first time since the number of dependent variables (1) is lower than the number of independent variables (3). That is the behaviour for a large number of independent variables. However, when one programs both alternatives, using the reverse mode differentiation once first and then forward differentiation once for calculating forces, the run time is 35% slower than using forward differentiation twice. This is due to the more complex implementation needed for reverse differentiation, and the fact that the ratio between dependent and independent variables is low. **This illustrates the fact that even if the time complexity in big O notation is lower for one of two algorithms, one should still run both algorithms if the number of variables is low, because this type of time complexity excludes constants and lower order terms.** 

By declaring the number of variables being differentiated each time forward mode is used (using the stack-based forward method implementation of FADBAD++), the run time was reduced a further 33% (again for calculating forces). The real and complex components of the pressure field are differentiated separately so one can easily use FADBAD++, since the transducer sound field (eq. 7) can be separated into real and complex terms with type double. 



## Details on how the equations in the paper are derived
The U_aa function can be derived from the Gorkov potential (eq. 3), by noting that (mag(p))^2 = p.p, and by applying the product rule on eq. 12. This yields for the derivative with respect to variable c:

(a.b)_c = (a.b_c) + (a_c.b)

Substituting this equation twice in (eq. 3) yields the U_aa function (with an extra 2 in front of the second term).

To find eq. 13 used for the gradient calculation which is then used by NLopt, one again uses the product rule on eq. 12. Note that the sound field equation for each transducer can be separated into a real and an imaginary component. The terms are in the form k * cos(phasej) and k * sin(phasej). For example, since (cos (x))_x = -sin(x), When one does the derivative of Real(P_g)_phasej one gets a term equal to the negative imaginary component, -Imag(P^j_g), hence the negative terms in eq. 13. No complex calculus needs to be used.

The paper uses acoustic the radiation force equation originally derived in the paper: [Acoustofluidics 7: the acoustic radiation force on small particles] (http://web-files.ait.dtu.dk/bruus/TMF/publications/pub2011/Bruus_Acoustofluidics_Tutorial_07_Lab_Chip_12_1014_2012.pdf). 



## Library architecture
The library was designed to be accurate but also fast, hence the decision to use automatic differentiation. The library was also designed to be easy to read. This was done by naming the data types in a similar way as the variables in the paper, and using the inline function id() found in common.hpp.

The Eigen matrix Mj for example, is based on the variable M^j in the paper, and it contains all the transducer pressure derivative terms (without the exp(i * phase^j) factor), where the rows are for each transducer and each column is for each derivative. The function id() is used in the code and with Mj in particular as follows: when Mj is filled with the output of the functions using FADBAD++ (in transducerDerivatives.hpp) the function id(a,b,c) corresponds to the location to store the first derivative (a) with respect to variable x,y,z or none, the second derivative (b) with respect to the same variables and so on. So for example, if one wants to get or write the value of the pressure for transducer j (M^j in the paper), one uses Mj(j, id(0,0,0)). If one wants M^j_yx, one uses Mj(j, id(2,1,0)) and if one wants M^j_zzz, one uses Mj(n, id(3,3,3)).

When one needs to calculate only the pressure, the potential and force at a given point, one only needs Mj values up to the second derivative. Therefore there are two functions to fill Mj, one for the optimization function used by NLopt that fills Mj up to third derivatives using forward mode AD three times (getThirdDerivsMjs()), and one for calculating pressure, potential and force that fills Mj up to second derivatices using forward mode AD two times (getSecondDerivsMjs()). Both functions are found in transducerDerivatives.hpp.

For calculating pressure, potential and force, there is a common initial loop that sums all the pressures from each transducer, however one only needs to sum up to second order derivatives for calculating the force (see force equation), up to first order derivatives for the potential (see eq. 3 in paper) and only the pressure terms for the pressure. This loop can be taken out of the three corresponding functions (inline pressureDerivatives() in transducerDerivatives.hpp) and when summing over pressure terms, their first, or their second derivatives, the summation will be over id() = 0, id() = 0 to 3 and id() = 0 to 15, since id(a,b,c) is equal to a + b * MAX_ID + c * MAX_ID * MAX_ID with MAX_ID = 4. This can be done because the right derivatives are stored contiguously in the columns of the Mj matrix.
