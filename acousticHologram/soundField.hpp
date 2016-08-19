#ifndef SOUNDFIELD_HPP
#define SOUNDFIELD_HPP


//Complex dot product
inline double compDot(complex<double> a, complex<double> b)
{
	return a.real()*b.real() + a.imag()*b.imag();
}

//Gradient of complex dot product
inline double gradCompDot(complex<double> pf, complex<double> pfj,  complex<double> pg, complex<double> pgj)
{
	return pf.imag()*pgj.real() + pfj.real()*pg.imag() - pf.real()*pgj.imag() - pfj.imag()*pg.real();
}


//Imaginary pressure for each node
template <typename DerivedA, typename DerivedB, typename DerivedC>
inline complex<double> imagPressure(DenseBase<DerivedA>& phases, DenseBase<DerivedB>& K, DenseBase<DerivedC>& Mj)
{
	int N_trans = phases.size();	
	complex<double> pressureComp = 0.0 + 0.0i;
	
	for (int n = 0; n < N_trans; n++) {
		pressureComp += exp(i*phases(n)) * Mj(n, id(0,0,0));
	}
	 
	return pressureComp;
}

//Presure derivatives for Gorkov potential and Force calculations
template <typename DerivedA, typename DerivedB, typename DerivedC>
inline void pressureDerivatives(DenseBase<DerivedA>& p, DenseBase<DerivedB>& phases, DenseBase<DerivedC>& Mj, int N_trans,int pSize)
{
	p.fill(0.0 + 0.0i);

	for (int x = 0; x < N_trans; x++) {
        	for (int n = 0; n < pSize; n++) p(n) += exp(i*phases(x))*Mj(x, n);
	}
}

//Gorkov potential for each node
template <typename DerivedA, typename DerivedB, typename DerivedC>
inline double gorkovPotential(DenseBase<DerivedA>& phases, DenseBase<DerivedB>& K, DenseBase<DerivedC>& Mj)
{
	int N_trans = phases.size();	
	int pSize = powInt(MAX_ID, NUM_IDS);
	cvec p(pSize);

	pressureDerivatives(p, phases, Mj, N_trans, pSize);

	//Calculate potential
	double P   = compDot(p(id(0,0,0)), p(id(0,0,0)));
	double P_x = compDot(p(id(1,0,0)), p(id(1,0,0)));
	double P_y = compDot(p(id(2,0,0)), p(id(2,0,0)));
	double P_z = compDot(p(id(3,0,0)), p(id(3,0,0)));
	 
	return K(0)*P - K(1)*(P_x + P_y + P_z);
}



//Forces experienced by particle at each node
template <typename DerivedA, typename DerivedB, typename DerivedC>
inline dmat particleForce(DenseBase<DerivedA>& phases, DenseBase<DerivedB>& K, DenseBase<DerivedC>& Mj)
{
	int a, n;	
	int pSize = powInt(MAX_ID, NUM_IDS);
	cvec p(pSize);	
	dmat forceComponents(1,3);
	int N_trans = phases.size();

	pressureDerivatives(p, phases, Mj, N_trans, pSize);

	double sumDerivs;
	for (a = 1; a < MAX_ID; a++) {
		sumDerivs = 0.0;

		for (n = 1; n < MAX_ID; n++) sumDerivs += compDot(p(id(n,a,0)), p(id(n,0,0)));

		//Changed sign of expression since: F = -grad(U)
		forceComponents(a-1) = -2.0*(K(0)*compDot(p(id(a,0,0)), p(id(0,0,0))) - K(1)*sumDerivs);
	}
	 
	return forceComponents;
}

#endif
