#ifndef OPTFUNCTION_HPP
#define OPTFUNCTION_HPP


typedef struct {
	int N_trans;
	
	double opt_scale;

	dmat w;
	dmat K;
	cmat Mj;
} my_func_data;



//Objective function and gradient
double myfunc(const std::vector<double> &phasesv, std::vector<double> &gradv, void *func_data)
{
	//Get function input data
	my_func_data *d = (my_func_data*) func_data;
    	int N_trans = d->N_trans;
	
	double opt_scale = d->opt_scale;

	dvec w; w = dvec::Map(d->w.data(),d->w.size());
	dvec K; K = dvec::Map(d->K.data(),d->K.size());
	
	cmat Mj = (d->Mj)*opt_scale;		


	//Declarations
	double f;
	int x,l,m,n,a; 
	int pSize = powInt(MAX_ID, NUM_IDS);

	cvec p(pSize); p.fill(0.0 + 0.0i);
	cmat pj(N_trans, pSize); pj.fill(0.0 + 0.0i);
	
	dvec Uderiv(3);
	dvec df(N_trans);
	dvec phases(N_trans);
	phases = dvec::Map(&phasesv[0], phasesv.size());

	
	//Objective function
	for (x = 0; x < N_trans; x++) {
        	for (n = 0; n < pSize; n++) {
			pj(x, n)  = exp(i*phases(x))*Mj(x, n);
			p(n)     += pj(x, n);
		}
	}

	double sumDerivs;	

	for (a = 1; a < MAX_ID; a++) {
		sumDerivs = 0.0;

		for (n = 1; n < MAX_ID; n++) sumDerivs += compDot(p(id(n,a,0)), p(id(n,a,0))) + compDot(p(id(a,0,0)), p(id(n,a,a)));

		Uderiv(a-1) = 2.0*(K(0)*(compDot(p(id(a,0,0)), p(id(a,0,0))) + compDot(p(id(0,0,0)), p(id(a,a,0)))) - K(1)*sumDerivs);
	}
	
	f  = w(0)*compDot(p(id(0,0,0)), p(id(0,0,0))) - Uderiv.dot(w.tail(3));
	

	//Gradient
	for (x = 0; x < N_trans; x++) {
		for (a = 1; a < MAX_ID; a++) {
			sumDerivs = 0.0;

			for (n = 1; n < MAX_ID; n++) {         
				sumDerivs += gradCompDot(p(id(n,a,0)),pj(x,id(n,a,0)),p(id(n,a,0)),pj(x,id(n,a,0))) + 
					     gradCompDot(p(id(a,0,0)),pj(x,id(a,0,0)),p(id(n,a,a)),pj(x,id(n,a,a)));
			}

			Uderiv(a-1) = 2.0*(K(0)*(gradCompDot(p(id(a,0,0)),pj(x,id(a,0,0)), p(id(a,0,0)),pj(x,id(a,0,0))) + 
					gradCompDot(p(id(0,0,0)),pj(x,id(0,0,0)), p(id(a,a,0)),pj(x,id(a,a,0)))) - K(1)*sumDerivs);
		}

		df(x) = w(0)*gradCompDot(p(id(0,0,0)),pj(x,id(0,0,0)), p(id(0,0,0)),pj(x,id(0,0,0))) - Uderiv.dot(w.tail(3));
	}

	//Map Eigen vector to c++ vector 
	dvec::Map(&gradv[0], df.size()) = df;

	return f;
}

#endif
