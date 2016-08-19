#ifndef AHFILEIO_HPP
#define AHFILEIO_HPP

#include "EigenFileIO.hpp"

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline void writeTransducerAndPhaseData(const char* filename, const double sf, const double k, const double P_0, const double r,
			const DenseBase<DerivedA>& K, const DenseBase<DerivedB>& levPoints, 
			RectangularTransducerArray& tArray, const DenseBase<DerivedC>& totalPhases)
{
	dmat transPs = tArray.transPoints();
	dmat transNs = tArray.transNormals();	
	int N_levPoints = levPoints.rows();
	int N_trans     = transPs.rows();	

	//Store all information in a single Eigen matrix
	dmat toWrite(7 + N_levPoints+N_trans*3,3); toWrite.fill(0.0);

	toWrite(0,0) = N_levPoints;
	toWrite(1,0) = N_trans;
	toWrite(2,0) = sf;
	toWrite(3,0) = k;
	toWrite(4,0) = P_0;
	toWrite(5,0) = r;
	toWrite(6,0) = K(0);    toWrite(6,1) = K(1);
	
	toWrite.block(7,                      0, N_levPoints,3) = levPoints;
	toWrite.block(7+N_levPoints,          0, N_trans    ,3) = transPs;
	toWrite.block(7+N_levPoints+N_trans,  0, N_trans    ,3) = transNs;

	toWrite.block(7+N_levPoints+N_trans*2,0, N_trans    ,1) = totalPhases;

	//Write information
	write_matrix(filename, toWrite);
}




template <typename DerivedA, typename DerivedB, typename DerivedC, typename DerivedD>
inline void writeFieldData(const char* filename,  
			const DenseBase<DerivedA>& pts, const DenseBase<DerivedB>& pressureMags, 
			const DenseBase<DerivedC>& potential, const DenseBase<DerivedD>& forces)
{
	int numNodes = pts.rows();
	dmat toWritePF(numNodes*4,3); toWritePF.fill(0.0);

	toWritePF.block(0,         0, numNodes,3) = pts;
	toWritePF.block(numNodes  ,0, numNodes,1) = pressureMags;
	toWritePF.block(numNodes*2,0, numNodes,1) = potential;
	toWritePF.block(numNodes*3,0, numNodes,3) = forces;
	
	//Write information
	write_matrix(filename, toWritePF);
}



template <typename DerivedA, typename DerivedB, typename DerivedC>
inline void writeParticlePathData(const char* filename, const double dt, const DenseBase<DerivedA>& paths, const DenseBase<DerivedB>& pCloudCenter,
					const DenseBase<DerivedC>& halfSide)
{
	int numRows = paths.rows();
	int numCols = paths.cols();

	dmat toWrite(numRows+3, numCols); toWrite.fill(0.0);

	toWrite(0,0) = dt;
	toWrite.block(1,0, 1,3)             = pCloudCenter;
	toWrite.block(2,0, 1,3)             = halfSide;
	toWrite.block(3,0, numRows,numCols) = paths;
	
	
	//Write information
	write_matrix(filename, toWrite);
}




template <typename DerivedA, typename DerivedB, typename DerivedC, typename DerivedD>
inline void writeLastParticlePositions(const char* filename, const double diam_p,const double sf, 
DenseBase<DerivedA>& pCloudCenter, DenseBase<DerivedB>& halfSide, DenseBase<DerivedC>& levPoints, DenseBase<DerivedD>& particlePFs)
{
	int numPaths      = particlePFs.rows();
	int numLevPoints  = levPoints.rows();


	dmat toWrite(3+numLevPoints+numPaths, 3); toWrite.fill(0.0);

	toWrite(0,0) = diam_p; toWrite(0,1) = sf; toWrite(0,2) = numLevPoints;
	
	toWrite.block(1,0, 1,3)             = pCloudCenter;
	toWrite.block(2,0, 1,3)             = halfSide;
	
	toWrite.block(3,0,              numLevPoints, 3) = levPoints;
	toWrite.block(3+numLevPoints,0, numPaths,3)      = particlePFs;
	
	
	//Write information
	write_matrix(filename, toWrite);
}

#endif
