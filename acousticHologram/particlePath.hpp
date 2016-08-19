#ifndef PARTICLEPATH_HPP
#define PARTICLEPATH_HPP


inline double getMaxRe(dmat& path, const double ReFactor)
{
	dmat Vs = path.block(0,3, path.rows(),3);
	
	dmat Vx = Vs.col(0);
	dmat Vy = Vs.col(1);
	dmat Vz = Vs.col(2);

	dmat VMags = Vx.cwiseProduct(Vx) + Vy.cwiseProduct(Vy) + Vz.cwiseProduct(Vz);		
	VMags = VMags.cwiseSqrt();

	return (VMags.maxCoeff())*ReFactor;
}



// Right hand side of the system of first order ODEs to solve
inline dmat RHS( const dmat& X , dvec& totalPhases, dmat& K, 
		const double P_0, const double k, const double r,		
		dmat& transPs, dmat& transNs, const int N_trans, 
		const double ReFactor, const double dragFactor, const double mass_p)
{
  	// Right hand side of the system of first order ODEs to solve, in this case:
  	// dx,y,z/dt = v;
  	// dx,y,z/dt(v) = F/m;
  	dmat output(1, NUM_IDS*2);

	//Positions
	dmat pt = X.block(0,0, 1,3);

	//Velocities
	dmat V(1,3);
	V = X.block(0,3, 1,3);
  	output.block(0,0, 1,3) = V;
  	

	//Acceleration
	//Get Mjs for each node using transducer info
	cmat Mj(N_trans, powInt(MAX_ID, NUM_IDS));
	dvec transP(3); 
	dvec transN(3);

	for (int n = 0; n < N_trans; n++) 
	{
		transP = transPs.row(n);
		transN = transNs.row(n);

		// Get real components (remember: with (real/imag)DDFunc only 1st and 2nd derivatives of Mj get filled in) 
		realDDFunc(pt(0),pt(1),pt(2), n,P_0,k,r, Mj,transP,transN);
		//Get imaginary components
		imagDDFunc(pt(0),pt(1),pt(2), n,P_0,k,r, Mj,transP,transN);
  	}	
	//Get acoustic force	 
	dmat acousticForce = particleForce(totalPhases, K,Mj);


	//Get drag force
	double vMag = V.norm();
	double Re_p = ReFactor*vMag;
	double C_D;

	if (Re_p <= 1.0e-8)                   C_D = 0.0;
	if ((Re_p > 1.0e-8) && (Re_p <= 1.0)) C_D = 24.0/Re_p;
	if ((Re_p > 1.0) && (Re_p <= 1.0e3))  C_D = 24.0/Re_p + 4.0/pow(Re_p, 0.33); 
	if (Re_p > 1.0e3)                     C_D = 0.4; 

	dmat dragForce(1,3);	
	dragForce = -C_D*dragFactor*vMag*V;


	//Total acceleration
	output.block(0,3, 1,3) = (acousticForce + dragForce)/mass_p;

  	return output;
}



//Function to get particle path using a 4th order Rugge-Kutta scheme (RK4)
inline dmat getParticlePathUsingRK4(const double dt,const double tMax, const dmat& particleP0, 
				dvec& totalPhases, dmat& K, 
				const double P_0, const double k, const double r,
				RectangularTransducerArray& tArray, 
				const double ReFactor, const double dragFactor, const double mass_p)
{
  	//initialize values
	dmat transPs = tArray.transPoints();
	dmat transNs = tArray.transNormals();
	int N_trans = tArray.nTrans();

  	double t   = 0.0;
  	int nSteps = round( ( tMax - t ) / dt );
	int NDIM   = NUM_IDS*2;

	// State variable X is [position, velocity]
  	dmat X(nSteps+1, NDIM); X.fill(0.0);
	X.block(0,0, 1,3) = particleP0;


	dmat k1(1,NDIM);
	dmat k2(1,NDIM);
	dmat k3(1,NDIM);
	dmat k4(1,NDIM);
	
  	for (int stepNumber = 0; stepNumber < nSteps; stepNumber++)
  	{
    		k1 = RHS( X.row(stepNumber)                 , totalPhases,K,P_0,k,r,transPs,transNs,N_trans, ReFactor,dragFactor,mass_p);
    		k2 = RHS( X.row(stepNumber) + dt / 2.0 * k1 , totalPhases,K,P_0,k,r,transPs,transNs,N_trans, ReFactor,dragFactor,mass_p);
    		k3 = RHS( X.row(stepNumber) + dt / 2.0 * k2 , totalPhases,K,P_0,k,r,transPs,transNs,N_trans, ReFactor,dragFactor,mass_p);
    		k4 = RHS( X.row(stepNumber) + dt * k3       , totalPhases,K,P_0,k,r,transPs,transNs,N_trans, ReFactor,dragFactor,mass_p);

    		X.row(stepNumber+1) = X.row(stepNumber) + dt / 6.0 * ( k1 + 2.0 * k2 + 2.0 * k3 + k4 );
  	} 

	return X;
}

#endif




