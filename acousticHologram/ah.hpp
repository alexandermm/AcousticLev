#ifndef AH_HPP
#define AH_HPP


//Type definitions
//typedef long double ld
typedef Matrix<int, Dynamic, 1      > ivec;

typedef Matrix<double, Dynamic, Dynamic> dmat;
typedef Matrix<double, Dynamic, 1      > dvec;

typedef Matrix<complex<double>, Dynamic, Dynamic> cmat;
typedef Matrix<complex<double>, Dynamic, 1      > cvec;


//Constants
double PI=acos(-1.0);
complex<double> i(0.0,1.0);
//long double pi=acosl(-1.0L)
const int MAX_ID  = 4;
const int NUM_IDS = 3;


#include "common.hpp"
#include "transducerArray.hpp"
#include "transducerDerivatives.hpp"
#include "soundField.hpp"
#include "optFunction.hpp"
#include "ahFileIO.hpp"
#include "fieldSampling.hpp"
#include "particlePath.hpp"


#endif
