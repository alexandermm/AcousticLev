#ifndef TRANSDUCERARRAY_HPP
#define TRANSDUCERARRAY_HPP

#include <Eigen/Geometry>


class RectangularTransducerArray 
{
	private:
		//Variables
		int N_trans;
		ivec sizeSubTransducerArrays;
	
		dmat transPs; //Eigen matrix representing 3d locations of transducers in a column
		dmat transNs; //Eigen matrix representing their corresponding normals
	
	public:
		//Constructor
		RectangularTransducerArray() {};
		RectangularTransducerArray(int, int, double, double);

		void rotate(double, double, double);		
		void translate(double, double, double);
		int nTrans() {return N_trans;}	
		dmat transPoints()   {return transPs;}
		dmat transNormals()  {return transNs;}
		ivec sizeSubArrays() {return sizeSubTransducerArrays;}

		//Class overloaded operators
		RectangularTransducerArray operator+(const RectangularTransducerArray&);
};

//Constructor
RectangularTransducerArray::RectangularTransducerArray(const int N_transx,const int N_transy, const double dx_trans,const double dy_trans)
{
	N_trans  = N_transx*N_transy;
	ivec intHolder(1); intHolder << N_trans;
	sizeSubTransducerArrays = intHolder;

	//Transducer array locations
	dvec xHolder(1,1); xHolder.setLinSpaced(N_transx, -dx_trans*(N_transx-1)/2.0,dx_trans*(N_transx-1)/2.0);
	dvec yHolder(1,1); yHolder.setLinSpaced(N_transy, -dy_trans*(N_transy-1)/2.0,dy_trans*(N_transy-1)/2.0);

	dvec x = xHolder.replicate(N_transy,1);

	dmat yHolder2 = yHolder.replicate(1, N_transx);
	yHolder2.transposeInPlace();
	dvec y(Map<VectorXd>(yHolder2.data(), N_trans));

	dvec z(N_trans); z.fill(0.0);

	transPs.resize(N_trans, 3);
	transPs.col(0) = x; transPs.col(1) = y; transPs.col(2) = z;

	//Transducer unit normal vectors
	dmat vec(1,3); vec << 0.0, 0.0, 1.0;
	transNs.resize(N_trans, 3);
	transNs = vec.replicate(N_trans,1);
}

//Rotate transducer array about origin using Euler angles in radians
void RectangularTransducerArray::rotate(double xAngle, double yAngle, double zAngle)
{
	//Make rotation matrix
	Matrix3d rot;
	rot = AngleAxisd(xAngle, Vector3d::UnitX()) * AngleAxisd(yAngle, Vector3d::UnitY()) * AngleAxisd(zAngle, Vector3d::UnitZ());	
	
	//Rotate transducer locations
	transPs = (rot*transPs.transpose()).transpose();

	//Rotate normals
	transNs = (rot*transNs.transpose()).transpose();
}

//Translate transducer array (just need to move points)
void RectangularTransducerArray::translate(double dx, double dy, double dz)
{
	//Translate transducer locations
	dvec oneVec = dvec::Ones(N_trans);

	transPs.col(0) += oneVec*dx; 
	transPs.col(1) += oneVec*dy; 
	transPs.col(2) += oneVec*dz;
}

//Class overloaded operators
RectangularTransducerArray RectangularTransducerArray::operator+(const RectangularTransducerArray& right)
{
	//Increase total number of transducer	
	this->N_trans += right.N_trans;

	//Write number of transducers in each sub array 
	this->sizeSubTransducerArrays.conservativeResize(this->sizeSubTransducerArrays.size()+1);	
	this->sizeSubTransducerArrays(this->sizeSubTransducerArrays.size()-1) = right.N_trans;

	//Append points and normals to first transducer array object
	int oldNumRows = this->transPs.rows();
	this->transPs.conservativeResize(oldNumRows+right.transPs.rows(), this->transPs.cols());
	this->transPs.block(oldNumRows,0, right.transPs.rows(),this->transPs.cols()) = right.transPs;

	this->transNs.conservativeResize(oldNumRows+right.transNs.rows(), this->transNs.cols());
	this->transNs.block(oldNumRows,0, right.transNs.rows(),this->transNs.cols()) = right.transNs;

	return *this;
}

#endif

