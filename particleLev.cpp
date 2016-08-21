#include <stdio.h>
#include <iostream>
#include <complex>
#include <ctime>

#include <Eigen/Dense>
#include <nlopt.hpp>

//Namespaces
using namespace std;
using namespace Eigen;

#include "./acousticHologram/ah.hpp"
//Can use #include <acousticHologram/ah> if folder found in the default include path




int main()
{
	//To get exactly same results as in original paper
	//All in SI units except: m -> mm
	double sf = 1.0e3;                                //Scale factor to change from m to mm in all units used

	//K_1 and K_2 constants
	double diam_p = (1.0e-3)*sf;                      //Particle diameter
	double V      = PI*diam_p*diam_p*diam_p/6.0;      //Particle volume
	double rho_p  = 29.36/(sf*sf*sf);                 //Particle density
	double c_p    = 1950.0*sf;                        //Speed of sound for expanded polystyrine (styrofoam)
	double rho_0  = 1.21/(sf*sf*sf);                  //Density of air
	double mu_0   = 1.8e-5/(sf*sf);			  //Air viscosity at 20 C
	double c_0    = 340.29*sf;                        //Speed of sound at sea level
	double freq   = 4.0e4;
	double omega  = 2.0*PI*freq;                      //Angular frequency of transducers

	dmat K(1,2);
	K << 0.25*V*(1.0/c_0/c_0/rho_0 - 1.0/c_p/c_p/rho_p), 0.75*V*(rho_0-rho_p)/omega/omega/rho_0/(rho_0+2.0*rho_p);

	//Optimization weights
	//w_p, w_x, w_y, w_z
	dmat w(1,4);
	w << 1.0, 1.0, 1.0, 1.0e2; //Bottle trap

	//Wave number
	double k = omega/c_0;

	//Wavelength
	double wl = c_0/freq;
	printf ("\nWavelength of transducer array sound waves in air: %3.2e m\n\n", wl/sf);


	//Transducer power constant (in Pa using mm)
	double P_0 = 1.0e4/sf;
	
	//Constant used to get good magnitude of derivatives for nlopt optimization function
	const double opt_scale = 5e-4;

	//Transducer piston radius
	double r = (3e-3)*sf;


	//Get transducer locations and normal vectors
	double dx_trans = (7e-3)*sf;
	double dy_trans = (7e-3)*sf;
	int N_transx = 20;
	int N_transy = 20;
	
	RectangularTransducerArray tArray(N_transx,N_transy, dx_trans,dy_trans);
	tArray.translate((N_transx-1)*dx_trans/2.0, (N_transy-1)*dy_trans/2.0, 0.0);


	//Levitation points
	int N_levPoints = 2, someN = 13;
	dmat levPoints(N_levPoints, 3);
	levPoints.row(0) << dx_trans*(N_transx-1)/2.0, dy_trans*(N_transy-1)/2.0, 0.20*sf;
	levPoints.row(1) << dx_trans*someN/2.0, dy_trans*someN/2.0, 0.18*sf;


	//Mj matrix (that will use id function)
	cmat Mj(tArray.nTrans(), powInt(MAX_ID, NUM_IDS));

	//Matrices used to pass levitation points
	dvec pLev(3);	

	dvec phases(tArray.nTrans());
	vector<double> phasesv; phasesv.resize(phases.size());

	dvec totalPhases(tArray.nTrans()); totalPhases.fill(0.0); 
	clock_t begin, end;
	double Mj_elapsed_secs = 0.0, opt_elapsed_secs = 0.0;

	//For every levitation point:
	for (int m = 0; m < N_levPoints; m++) {
		//Get point
		pLev = levPoints.row(m);		
		
		//Getting Mj complex values
		begin = clock();
		getThirdDerivMjs(pLev(0),pLev(1),pLev(2), P_0,k,r, Mj,tArray);
		end = clock();
  		Mj_elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;
	
		//Optimize phases for levitation point using l-bfgs	
		nlopt::opt opt(nlopt::LD_LBFGS, tArray.nTrans());

		//Set up objective function
		my_func_data func_data = {tArray.nTrans(),opt_scale, w,K,Mj};
		opt.set_min_objective(&myfunc, &func_data);

		//Set optimization parameters
		opt.set_xtol_rel(1e-10);

		//Initial guess (initialize in Eigen vector and change to c++ vector to go into nlopt function)
		phases.fill(0.0);
		dvec::Map(&phasesv[0], phases.size()) = phases;

		//Optimize and store minimum value
		double minf;

		begin = clock();
		nlopt::result result = opt.optimize(phasesv, minf);
		end = clock();
  		opt_elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;

		printf("Minimum value for point %d is %8.10f\n", m+1,minf);

		//Change result back to an Eigen vector
		phases = dvec::Map(phasesv.data(), phasesv.size());

		//Add to total phase looping back to -pi if total higher than pi
		for (int n = 0; n < tArray.nTrans(); n++) {
    
			totalPhases(n) += phases(n);
			if (totalPhases(n) > PI) totalPhases(n) -= 2.0*PI;
  		}

	}
	//Report time
	cout << endl;
	printf("Time to calculate Mj values: %4.3f secs\n", Mj_elapsed_secs);
	printf("Time to optimize points:     %4.3f secs\n\n", opt_elapsed_secs);
	//cout << Mj.block(0,0, 10,3) << endl << endl;
	
	//Write information
	writeTransducerAndPhaseData("levData.txt", sf,k,P_0,r, K,levPoints,tArray,totalPhases);

	
	//Compute pressure field in some volume surrounding levitation points
	dmat pCenter(1,3);
	pCenter << (6.64e-2)*sf, (6.64e-2)*sf, (19.01e-2)*sf;

	double dx    = (4e-3)*sf, dy  = (4e-3)*sf, dz  = (4e-3)*sf;
	int nDx      = 3;
	int numSide  = nDx*2+1;
	int numNodes = powInt(numSide, 3);

	//Make grid points
	dmat pts = get3Dgrid(pCenter, dx,dy,dz, nDx);
		

	//Calculate pressure magnitude, Gorkov potential and force experienced by particles for each node
	dvec pt(3);
	double elapsed_secs;
	dvec pressureMags(numNodes);	
	dvec potential(numNodes);
	dmat forces(numNodes, 3);
	
	begin = clock();
	for (int m = 0; m < numNodes; m++) {
		pt = pts.row(m);

		//Get Mjs for each node using transducer info
		getSecondDerivMjs(pt(0),pt(1),pt(2), P_0,k,r, Mj,tArray);

		pressureMags(m) = abs(imagPressure(totalPhases, Mj));
		potential(m)    = gorkovPotential(totalPhases, K,Mj);
		forces.row(m)   = particleForce(totalPhases, K,Mj);
	}
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("Time to calculate P mag, Gorkov pot and force at nodes: %4.3f secs\n", elapsed_secs);

	//Write pressure, potential and force field data
	writeFieldData("fieldData.txt", pts,pressureMags,potential,forces); 
	


	//Get particle paths for two particles at specified locations
	int numParticles = 2;

	dmat particleP0s(numParticles,3);
	particleP0s.row(0) << (5.51e-2)*sf, (7.01e-2)*sf, (20.51e-2)*sf;
	particleP0s.row(1) << (4.01e-2)*sf, (6.01e-2)*sf, (18.41e-2)*sf;


	//Keep particles in cube
	dmat pCloudCenter(1,3); pCloudCenter = pCenter;
	dmat halfSide(1,3); halfSide << 0.032*sf, 0.032*sf, 0.032*sf;
	
	double xMin = pCloudCenter(0) - halfSide(0); double xMax = pCloudCenter(0) + halfSide(0);
	double yMin = pCloudCenter(1) - halfSide(1); double yMax = pCloudCenter(1) + halfSide(1);
	double zMin = pCloudCenter(2) - halfSide(2); double zMax = pCloudCenter(2) + halfSide(2);

	dmat particlePFs(numParticles,3);

	double tMax     = 1.0;
  	double dt       = 0.01;

	double ReFactor   = rho_0*diam_p/mu_0;
	double dragFactor = 0.5*rho_0*PI*0.25*diam_p*diam_p;
	double mass_p     = V*rho_p;
	double maxRe;

	int nSteps = round( (tMax) / dt );
	int NDIM   = NUM_IDS*2;
  	
	dmat particlePath(nSteps+1, NDIM);
	dmat particlePaths(3*numParticles, nSteps+1);
	
	begin = clock();
	for (int m = 0; m < numParticles; m++)
	{
		particlePath = 
		getParticlePathUsingRK4(dt,tMax,particleP0s.row(m), totalPhases,K,P_0,k,r,tArray, ReFactor,dragFactor,mass_p);

		particlePFs.row(m) = particlePath.block(nSteps,0, 1,3);
		particlePaths.block(3*m,0, 3,nSteps+1) = (particlePath.transpose()).block(0,0, 3,nSteps+1);

		if (m == 0) maxRe = getMaxRe(particlePath, ReFactor);
	}

	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("Time to calculate particle paths:                      %4.3f secs\n\n", elapsed_secs);	

	printf("Max Reynolds number of particle 1: %4.1f \n", maxRe);
	printf("Final location of particle 1:      (%3.2e, %3.2e, %3.2e) m\n\n",
		particlePFs(0,0)/sf,particlePFs(0,1)/sf,particlePFs(0,2)/sf);

	//Write particle path data
	writeParticlePathData("pathData.txt", dt,particlePaths, pCloudCenter,halfSide);
	
	//Write final particle position file
	writeLastParticlePositions("positionData.txt", diam_p,sf, pCloudCenter,halfSide, levPoints,particlePFs);
	

	return 0;
}
