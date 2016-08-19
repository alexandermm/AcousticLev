#ifndef TRANSDUCERDERIVATIVES_HPP
#define TRANSDUCERDERIVATIVES_HPP

#include "fadiff.h"
using namespace fadbad;



template <typename DerivedA, typename DerivedB>
F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS> 
realFunc(const F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS>& x, const F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS>& y, 
			     const F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS>& z, 
	                     const double &scaledP_0, const double &k, const double &r, 
	                     const DenseBase<DerivedA>& transP, const DenseBase<DerivedB>& transN)
{
	F< F< F<double,NUM_IDS> ,NUM_IDS > ,NUM_IDS > dx = x-transP(0), dy = y-transP(1), dz = z-transP(2);
	F< F< F<double,NUM_IDS> ,NUM_IDS > ,NUM_IDS > dj = sqrt( dx*dx+dy*dy+dz*dz );
	F< F< F<double,NUM_IDS> ,NUM_IDS > ,NUM_IDS > sin_angle = sqrt( (dy*transN(2)-dz*transN(1))*(dy*transN(2)-dz*transN(1)) + 
			    (dz*transN(0)-dx*transN(2))*(dz*transN(0)-dx*transN(2)) + (dx*transN(1)-dy*transN(0))*(dx*transN(1)-dy*transN(0)))/dj;

        return scaledP_0*j0( k*r*sin_angle )*cos(k*dj)/dj;
}

template <typename DerivedA, typename DerivedB>
F< F<double,NUM_IDS> ,NUM_IDS> 
realDFunc(F< F<double,NUM_IDS> ,NUM_IDS>& o_dfdx, F< F<double,NUM_IDS> ,NUM_IDS>& o_dfdy, F< F<double,NUM_IDS> ,NUM_IDS>& o_dfdz,
			 const F< F<double,NUM_IDS> ,NUM_IDS>&i_x,const F< F<double,NUM_IDS> ,NUM_IDS>&i_y,const F< F<double,NUM_IDS> ,NUM_IDS>&i_z, 
			 const double &scaledP_0, const double &k, const double &r,
			 const DenseBase<DerivedA>& transP, const DenseBase<DerivedB>& transN)
{
	F< F< F<double,NUM_IDS> ,NUM_IDS > ,NUM_IDS > x(i_x),y(i_y),z(i_z);   // Initialize arguments
        
	x.diff(0); y.diff(1); z.diff(2);
	
	F< F< F<double,NUM_IDS> ,NUM_IDS > ,NUM_IDS > f(realFunc(x,y,z, scaledP_0,k,r, transP, transN));    // Evaluate function and record DAG
	
	o_dfdx=f.d(0); o_dfdy=f.d(1); o_dfdz=f.d(2);
	return f.x();
}

template <typename DerivedA, typename DerivedB>
F<double,NUM_IDS> realDDFunc(F<double,NUM_IDS>& o_dfdxdx, F<double,NUM_IDS>& o_dfdxdy, F<double,NUM_IDS>& o_dfdxdz,
		     F<double,NUM_IDS>& o_dfdydx, F<double,NUM_IDS>& o_dfdydy, F<double,NUM_IDS>& o_dfdydz,
		     F<double,NUM_IDS>& o_dfdzdx, F<double,NUM_IDS>& o_dfdzdy, F<double,NUM_IDS>& o_dfdzdz,
		     F<double,NUM_IDS>& o_dfdx, F<double,NUM_IDS>& o_dfdy, F<double,NUM_IDS>& o_dfdz,
		     const F<double,NUM_IDS>& i_x, const F<double,NUM_IDS>& i_y, const F<double,NUM_IDS>& i_z, 
		     const double &scaledP_0, const double &k, const double &r, 
		     const DenseBase<DerivedA>& transP, const DenseBase<DerivedB>& transN)
{
	F< F<double,NUM_IDS> ,NUM_IDS > x(i_x),y(i_y),z(i_z), dfdx,dfdy,dfdz;// Initialize arguments
	
	x.diff(0); y.diff(1); z.diff(2);

	F< F<double,NUM_IDS> ,NUM_IDS > f(realDFunc(dfdx,dfdy,dfdz, x,y,z, scaledP_0,k,r, transP, transN));// Evaluate function and derivatives
	
	o_dfdxdx=dfdx.d(0); o_dfdxdy=dfdx.d(1); o_dfdxdz=dfdx.d(2);	
	o_dfdydx=dfdy.d(0); o_dfdydy=dfdy.d(1); o_dfdydz=dfdy.d(2);
	o_dfdzdx=dfdz.d(0); o_dfdzdy=dfdz.d(1); o_dfdzdz=dfdz.d(2);
	o_dfdx=dfdx.x(); o_dfdy=dfdy.x(); o_dfdz=dfdz.x();	
	return f.x();                   
}

///////////////////////////////////////////////////////////////////////////////////////////////////////


template <typename DerivedA, typename DerivedB>
F< F<double,NUM_IDS> ,NUM_IDS> realFunc(const F< F<double,NUM_IDS> ,NUM_IDS>& x, const F< F<double,NUM_IDS> ,NUM_IDS>& y, 
			     const F< F<double,NUM_IDS> ,NUM_IDS>& z, 
	                     const double &scaledP_0, const double &k, const double &r, 
	                     const DenseBase<DerivedA>& transP, const DenseBase<DerivedB>& transN)
{
	F< F<double,NUM_IDS> ,NUM_IDS> dx = x-transP(0), dy = y-transP(1), dz = z-transP(2);
	F< F<double,NUM_IDS> ,NUM_IDS> dj = sqrt( dx*dx+dy*dy+dz*dz );
	F< F<double,NUM_IDS> ,NUM_IDS> sin_angle = sqrt( (dy*transN(2)-dz*transN(1))*(dy*transN(2)-dz*transN(1)) + 
			    (dz*transN(0)-dx*transN(2))*(dz*transN(0)-dx*transN(2)) + (dx*transN(1)-dy*transN(0))*(dx*transN(1)-dy*transN(0)))/dj;

        return scaledP_0*j0( k*r*sin_angle )*cos(k*dj)/dj;
}
template <typename DerivedA, typename DerivedB>
F<double,NUM_IDS> realDFunc(F<double,NUM_IDS>& o_dfdx, F<double,NUM_IDS>& o_dfdy, F<double,NUM_IDS>& o_dfdz,
			 const F<double,NUM_IDS>& i_x, const F<double,NUM_IDS>& i_y, const F<double,NUM_IDS>& i_z, 
			 const double &scaledP_0, const double &k, const double &r,
			 const DenseBase<DerivedA>& transP, const DenseBase<DerivedB>& transN)
{
	F< F<double,NUM_IDS> ,NUM_IDS> x(i_x),y(i_y),z(i_z);   // Initialize arguments
        
	x.diff(0); y.diff(1); z.diff(2);
	
	F< F<double,NUM_IDS> ,NUM_IDS> f(realFunc(x,y,z, scaledP_0,k,r, transP, transN));    // Evaluate function and record DAG
	
	o_dfdx=f.d(0); o_dfdy=f.d(1); o_dfdz=f.d(2);
	return f.x();
}
template <typename DerivedA, typename DerivedB, typename DerivedC>
void realDDFunc(const double& i_x, const double& i_y, const double& i_z, 
		     const int& n, const double &scaledP_0, const double &k, const double &r, 
		     DenseBase<DerivedA>& Mj, const DenseBase<DerivedB>& transP, const DenseBase<DerivedC>& transN)
{
	F<double,NUM_IDS> x(i_x),y(i_y),z(i_z), dfdx,dfdy,dfdz;// Initialize arguments
	
	x.diff(0); y.diff(1); z.diff(2);

	F<double,NUM_IDS> f(realDFunc(dfdx,dfdy,dfdz, x,y,z, scaledP_0,k,r, transP, transN));// Evaluate function and derivatives
	
	Mj(n, id(0,0,0)).real(f.x());
		
	//First derivatives
	Mj(n, id(1,0,0)).real(dfdx.x());    Mj(n, id(2,0,0)).real(dfdy.x());    Mj(n, id(3,0,0)).real(dfdz.x());

	//Second derivatives
	Mj(n, id(1,1,0)).real(dfdx.d(0));   Mj(n, id(2,1,0)).real(dfdy.d(0));   Mj(n, id(3,1,0)).real(dfdz.d(0));
	Mj(n, id(1,2,0)).real(dfdx.d(1));   Mj(n, id(2,2,0)).real(dfdy.d(1));   Mj(n, id(3,2,0)).real(dfdz.d(1));
	Mj(n, id(1,3,0)).real(dfdx.d(2));   Mj(n, id(2,3,0)).real(dfdy.d(2));   Mj(n, id(3,3,0)).real(dfdz.d(2));      
}


template <typename DerivedA, typename DerivedB>
F< F<double,NUM_IDS> ,NUM_IDS> 
imagFunc(const F< F<double,NUM_IDS> ,NUM_IDS>& x, const F< F<double,NUM_IDS> ,NUM_IDS>& y, const F< F<double,NUM_IDS> ,NUM_IDS>& z, 
	                     const double &scaledP_0, const double &k, const double &r, 
	                     const DenseBase<DerivedA>& transP, const DenseBase<DerivedB>& transN)
{
	F< F<double,NUM_IDS> ,NUM_IDS> dx = x-transP(0), dy = y-transP(1), dz = z-transP(2);
	F< F<double,NUM_IDS> ,NUM_IDS> dj = sqrt( dx*dx+dy*dy+dz*dz );
	F< F<double,NUM_IDS> ,NUM_IDS> sin_angle = sqrt( (dy*transN(2)-dz*transN(1))*(dy*transN(2)-dz*transN(1)) + 
			    (dz*transN(0)-dx*transN(2))*(dz*transN(0)-dx*transN(2)) + (dx*transN(1)-dy*transN(0))*(dx*transN(1)-dy*transN(0)))/dj;

        return scaledP_0*j0( k*r*sin_angle )*sin(k*dj)/dj;
}
template <typename DerivedA, typename DerivedB>
F<double,NUM_IDS> imagDFunc(F<double,NUM_IDS>& o_dfdx, F<double,NUM_IDS>& o_dfdy, F<double,NUM_IDS>& o_dfdz,
			 const F<double,NUM_IDS>& i_x, const F<double,NUM_IDS>& i_y, const F<double,NUM_IDS>& i_z, 
			 const double &scaledP_0, const double &k, const double &r,
			 const DenseBase<DerivedA>& transP, const DenseBase<DerivedB>& transN)
{
	F< F<double,NUM_IDS> ,NUM_IDS> x(i_x),y(i_y),z(i_z);   // Initialize arguments
        
	x.diff(0); y.diff(1); z.diff(2);
	
	F< F<double,NUM_IDS> ,NUM_IDS> f(imagFunc(x,y,z, scaledP_0,k,r, transP, transN));    // Evaluate function and record DAG
	
	o_dfdx=f.d(0); o_dfdy=f.d(1); o_dfdz=f.d(2);
	return f.x();
}
template <typename DerivedA, typename DerivedB, typename DerivedC>
void imagDDFunc(const double& i_x, const double& i_y, const double& i_z, 
		     const int& n, const double &scaledP_0, const double &k, const double &r, 
		     DenseBase<DerivedA>& Mj, const DenseBase<DerivedB>& transP, const DenseBase<DerivedC>& transN)
{
	F<double,NUM_IDS> x(i_x),y(i_y),z(i_z), dfdx,dfdy,dfdz;// Initialize arguments
	
	x.diff(0); y.diff(1); z.diff(2);

	F<double,NUM_IDS> f(imagDFunc(dfdx,dfdy,dfdz, x,y,z, scaledP_0,k,r, transP, transN));// Evaluate function and derivatives
	
	Mj(n, id(0,0,0)).imag(f.x());
		
	//First derivatives
	Mj(n, id(1,0,0)).imag(dfdx.x());    Mj(n, id(2,0,0)).imag(dfdy.x());    Mj(n, id(NUM_IDS,0,0)).imag(dfdz.x());

	//Second derivatives
	Mj(n, id(1,1,0)).imag(dfdx.d(0));   Mj(n, id(2,1,0)).imag(dfdy.d(0));   Mj(n, id(3,1,0)).imag(dfdz.d(0));
	Mj(n, id(1,2,0)).imag(dfdx.d(1));   Mj(n, id(2,2,0)).imag(dfdy.d(1));   Mj(n, id(3,2,0)).imag(dfdz.d(1));
	Mj(n, id(1,3,0)).imag(dfdx.d(2));   Mj(n, id(2,3,0)).imag(dfdy.d(2));   Mj(n, id(3,3,0)).imag(dfdz.d(2));                 
}




















///////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename DerivedA, typename DerivedB, typename DerivedC>
void realDDDFunc(const double& i_x, const double& i_y, const double& i_z, 
		 const int& n, const double& scaledP_0, const double& k, const double& r, 
		 DenseBase<DerivedA>& Mj, const DenseBase<DerivedB>& transP, const DenseBase<DerivedC>& transN)
{
	F<double,NUM_IDS> x(i_x),y(i_y),z(i_z), dfdx,dfdy,dfdz, dfdxdx,dfdydx,dfdzdx, dfdxdy,dfdydy,dfdzdy, dfdxdz,dfdydz,dfdzdz;// Initialize arguments
	
	x.diff(0); y.diff(1); z.diff(2);                    // Second order wrt. y

	F<double,NUM_IDS> 
	f(realDDFunc(dfdxdx,dfdxdy,dfdxdz, dfdydx,dfdydy,dfdydz, dfdzdx,dfdzdy,dfdzdz, dfdx,dfdy,dfdz, x,y,z, scaledP_0,k,r, transP, transN));
	
	
	Mj(n, id(0,0,0)).real(f.x());
		
	//First derivatives
	Mj(n, id(1,0,0)).real(dfdx.x());    Mj(n, id(2,0,0)).real(dfdy.x());    Mj(n, id(3,0,0)).real(dfdz.x());

	//Second derivatives
	Mj(n, id(1,1,0)).real(dfdx.d(0));   Mj(n, id(2,1,0)).real(dfdy.d(0));   Mj(n, id(3,1,0)).real(dfdz.d(0));
	Mj(n, id(1,2,0)).real(dfdx.d(1));   Mj(n, id(2,2,0)).real(dfdy.d(1));   Mj(n, id(3,2,0)).real(dfdz.d(1));
	Mj(n, id(1,3,0)).real(dfdx.d(2));   Mj(n, id(2,3,0)).real(dfdy.d(2));   Mj(n, id(3,3,0)).real(dfdz.d(2));

	//Third derivatives
	Mj(n, id(1,1,1)).real(dfdxdx.d(0)); Mj(n, id(2,1,1)).real(dfdydx.d(0)); Mj(n, id(3,1,1)).real(dfdzdx.d(0)); 
	Mj(n, id(1,2,1)).real(dfdxdy.d(0)); Mj(n, id(2,2,1)).real(dfdydy.d(0)); Mj(n, id(3,2,1)).real(dfdzdy.d(0));
	Mj(n, id(1,3,1)).real(dfdxdz.d(0)); Mj(n, id(2,3,1)).real(dfdydz.d(0)); Mj(n, id(3,3,1)).real(dfdzdz.d(0));

	Mj(n, id(1,1,2)).real(dfdxdx.d(1)); Mj(n, id(2,1,2)).real(dfdydx.d(1)); Mj(n, id(3,1,2)).real(dfdzdx.d(1)); 
	Mj(n, id(1,2,2)).real(dfdxdy.d(1)); Mj(n, id(2,2,2)).real(dfdydy.d(1)); Mj(n, id(3,2,2)).real(dfdzdy.d(1));
	Mj(n, id(1,3,2)).real(dfdxdz.d(1)); Mj(n, id(2,3,2)).real(dfdydz.d(1)); Mj(n, id(3,3,2)).real(dfdzdz.d(1));

	Mj(n, id(1,1,3)).real(dfdxdx.d(2)); Mj(n, id(2,1,3)).real(dfdydx.d(2)); Mj(n, id(3,1,3)).real(dfdzdx.d(2)); 
	Mj(n, id(1,2,3)).real(dfdxdy.d(2)); Mj(n, id(2,2,3)).real(dfdydy.d(2)); Mj(n, id(3,2,3)).real(dfdzdy.d(2));
	Mj(n, id(1,3,3)).real(dfdxdz.d(2)); Mj(n, id(2,3,3)).real(dfdydz.d(2)); Mj(n, id(3,3,3)).real(dfdzdz.d(2));	
}




template <typename DerivedA, typename DerivedB>
F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS> 
imagFunc(const F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS>& x, const F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS>& y, 
			     const F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS>& z, 
			     const double &scaledP_0, const double &k, const double &r, 
			     const DenseBase<DerivedA>& transP, const DenseBase<DerivedB>& transN)
{
	F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS> dx = x-transP(0), dy = y-transP(1), dz = z-transP(2);
	F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS> dj = sqrt( dx*dx+dy*dy+dz*dz );
	F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS> sin_angle = sqrt( (dy*transN(2)-dz*transN(1))*(dy*transN(2)-dz*transN(1)) + 
			    (dz*transN(0)-dx*transN(2))*(dz*transN(0)-dx*transN(2)) + (dx*transN(1)-dy*transN(0))*(dx*transN(1)-dy*transN(0)))/dj;

        return scaledP_0*j0( k*r*sin_angle )*sin(k*dj)/dj;
}

template <typename DerivedA, typename DerivedB>
F< F<double,NUM_IDS> ,NUM_IDS> 
imagDFunc(F< F<double,NUM_IDS> ,NUM_IDS>& o_dfdx, F< F<double,NUM_IDS> ,NUM_IDS>& o_dfdy, F< F<double,NUM_IDS> ,NUM_IDS>& o_dfdz,
			 const F< F<double,NUM_IDS> ,NUM_IDS>& i_x, const F< F<double,NUM_IDS> ,NUM_IDS>& i_y, 
			 const F< F<double,NUM_IDS> ,NUM_IDS>& i_z, 
			 const double &scaledP_0, const double &k, const double &r,
			 const DenseBase<DerivedA>& transP, const DenseBase<DerivedB>& transN)
{
	F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS> x(i_x),y(i_y),z(i_z);   // Initialize arguments
        
	x.diff(0); y.diff(1); z.diff(2);
	
	F< F< F<double,NUM_IDS> ,NUM_IDS> ,NUM_IDS> f(imagFunc(x,y,z, scaledP_0,k,r, transP, transN));    // Evaluate function and record DAG
	
	o_dfdx=f.d(0); o_dfdy=f.d(1); o_dfdz=f.d(2);
	return f.x();
}

template <typename DerivedA, typename DerivedB>
F<double,NUM_IDS> imagDDFunc(F<double,NUM_IDS>& o_dfdxdx, F<double,NUM_IDS>& o_dfdxdy, F<double,NUM_IDS>& o_dfdxdz,
		     F<double,NUM_IDS>& o_dfdydx, F<double,NUM_IDS>& o_dfdydy, F<double,NUM_IDS>& o_dfdydz,
		     F<double,NUM_IDS>& o_dfdzdx, F<double,NUM_IDS>& o_dfdzdy, F<double,NUM_IDS>& o_dfdzdz,
		     F<double,NUM_IDS>& o_dfdx, F<double,NUM_IDS>& o_dfdy, F<double,NUM_IDS>& o_dfdz,
		     const F<double,NUM_IDS>& i_x, const F<double,NUM_IDS>& i_y, const F<double,NUM_IDS>& i_z, 
		     const double &scaledP_0, const double &k, const double &r, 
		     const DenseBase<DerivedA>& transP, const DenseBase<DerivedB>& transN)
{
	F< F<double,NUM_IDS> ,NUM_IDS> x(i_x),y(i_y),z(i_z), dfdx,dfdy,dfdz;// Initialize arguments
	
	x.diff(0); y.diff(1); z.diff(2);

	F< F<double,NUM_IDS> ,NUM_IDS> f(imagDFunc(dfdx,dfdy,dfdz, x,y,z, scaledP_0,k,r, transP, transN));// Evaluate function and derivatives
	
	o_dfdxdx=dfdx.d(0); o_dfdxdy=dfdx.d(1); o_dfdxdz=dfdx.d(2);	
	o_dfdydx=dfdy.d(0); o_dfdydy=dfdy.d(1); o_dfdydz=dfdy.d(2);
	o_dfdzdx=dfdz.d(0); o_dfdzdy=dfdz.d(1); o_dfdzdz=dfdz.d(2);
	o_dfdx=dfdx.x(); o_dfdy=dfdy.x(); o_dfdz=dfdz.x();	
	return f.x();                   
}

template <typename DerivedA, typename DerivedB, typename DerivedC>
void imagDDDFunc(const double& i_x, const double& i_y, const double& i_z, 
		 const int& n, const double& scaledP_0, const double& k, const double& r, 
		 DenseBase<DerivedA>& Mj, const DenseBase<DerivedB>& transP, const DenseBase<DerivedC>& transN)
{
	F<double,NUM_IDS> x(i_x),y(i_y),z(i_z), dfdx,dfdy,dfdz, dfdxdx,dfdydx,dfdzdx, dfdxdy,dfdydy,dfdzdy, dfdxdz,dfdydz,dfdzdz;
	
	x.diff(0); y.diff(1); z.diff(2);                    // Second order wrt. y

	F<double,NUM_IDS> 
	f(imagDDFunc(dfdxdx,dfdxdy,dfdxdz, dfdydx,dfdydy,dfdydz, dfdzdx,dfdzdy,dfdzdz, dfdx,dfdy,dfdz, x,y,z, scaledP_0,k,r, transP, transN));
	

	Mj(n, id(0,0,0)).imag(f.x());
		
	//First derivatives
	Mj(n, id(1,0,0)).imag(dfdx.x());    Mj(n, id(2,0,0)).imag(dfdy.x());    Mj(n, id(3,0,0)).imag(dfdz.x());

	//Second derivatives
	Mj(n, id(1,1,0)).imag(dfdx.d(0));   Mj(n, id(2,1,0)).imag(dfdy.d(0));   Mj(n, id(3,1,0)).imag(dfdz.d(0));
	Mj(n, id(1,2,0)).imag(dfdx.d(1));   Mj(n, id(2,2,0)).imag(dfdy.d(1));   Mj(n, id(3,2,0)).imag(dfdz.d(1));
	Mj(n, id(1,3,0)).imag(dfdx.d(2));   Mj(n, id(2,3,0)).imag(dfdy.d(2));   Mj(n, id(3,3,0)).imag(dfdz.d(2));

	//Third derivatives
	Mj(n, id(1,1,1)).imag(dfdxdx.d(0)); Mj(n, id(2,1,1)).imag(dfdydx.d(0)); Mj(n, id(3,1,1)).imag(dfdzdx.d(0)); 
	Mj(n, id(1,2,1)).imag(dfdxdy.d(0)); Mj(n, id(2,2,1)).imag(dfdydy.d(0)); Mj(n, id(3,2,1)).imag(dfdzdy.d(0));
	Mj(n, id(1,3,1)).imag(dfdxdz.d(0)); Mj(n, id(2,3,1)).imag(dfdydz.d(0)); Mj(n, id(3,3,1)).imag(dfdzdz.d(0));

	Mj(n, id(1,1,2)).imag(dfdxdx.d(1)); Mj(n, id(2,1,2)).imag(dfdydx.d(1)); Mj(n, id(3,1,2)).imag(dfdzdx.d(1)); 
	Mj(n, id(1,2,2)).imag(dfdxdy.d(1)); Mj(n, id(2,2,2)).imag(dfdydy.d(1)); Mj(n, id(3,2,2)).imag(dfdzdy.d(1));
	Mj(n, id(1,3,2)).imag(dfdxdz.d(1)); Mj(n, id(2,3,2)).imag(dfdydz.d(1)); Mj(n, id(3,3,2)).imag(dfdzdz.d(1));

	Mj(n, id(1,1,3)).imag(dfdxdx.d(2)); Mj(n, id(2,1,3)).imag(dfdydx.d(2)); Mj(n, id(3,1,3)).imag(dfdzdx.d(2)); 
	Mj(n, id(1,2,3)).imag(dfdxdy.d(2)); Mj(n, id(2,2,3)).imag(dfdydy.d(2)); Mj(n, id(3,2,3)).imag(dfdzdy.d(2));
	Mj(n, id(1,3,3)).imag(dfdxdz.d(2)); Mj(n, id(2,3,3)).imag(dfdydz.d(2)); Mj(n, id(3,3,3)).imag(dfdzdz.d(2));	
}


template <typename DerivedA>
void getThirdDerivMjs(const double& i_x, const double& i_y, const double& i_z, const double& P_0, const double& k, const double& r, 
		 DenseBase<DerivedA>& Mj, RectangularTransducerArray& tArray)
{
	dvec transP(3); 
	dvec transN(3);
	dmat transPs = tArray.transPoints();
	dmat transNs = tArray.transNormals();
	int N_trans = tArray.nTrans();

	for (int n = 0; n < N_trans; n++) {
    
		transP = transPs.row(n);
		transN = transNs.row(n);

		// Get real components
		realDDDFunc(i_x,i_y,i_z, n,P_0,k,r, Mj,transP,transN);
		//Get imaginary components
		imagDDDFunc(i_x,i_y,i_z, n,P_0,k,r, Mj,transP,transN);
  	}
}


template <typename DerivedA>
void getSecondDerivMjs(const double& i_x, const double& i_y, const double& i_z, const double& P_0, const double& k, const double& r, 
		 DenseBase<DerivedA>& Mj, RectangularTransducerArray& tArray)
{
	dvec transP(3); 
	dvec transN(3);
	dmat transPs = tArray.transPoints();
	dmat transNs = tArray.transNormals();
	int N_trans = tArray.nTrans();

	for (int n = 0; n < N_trans; n++) {
    
		transP = transPs.row(n);
		transN = transNs.row(n);

		// Get real components
		realDDFunc(i_x,i_y,i_z, n,P_0,k,r, Mj,transP,transN);
		//Get imaginary components
		imagDDFunc(i_x,i_y,i_z, n,P_0,k,r, Mj,transP,transN);
  	}
}

#endif
