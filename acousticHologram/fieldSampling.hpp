#ifndef FIELDSAMPLING_HPP
#define FIELDSAMPLING_HPP


inline dmat get3Dgrid(dmat &pCenter, const double dx,const double dy,const double dz, const int nDx)
{
	int numSide  = nDx*2+1;
	int numNodes = powInt(numSide, 3);

	dvec xLine(1,1); xLine.setLinSpaced(numSide, pCenter(0)-nDx*dx, pCenter(0)+nDx*dx);
	dvec yLine(1,1); yLine.setLinSpaced(numSide, pCenter(1)-nDx*dy, pCenter(1)+nDx*dy);
	dvec zLine(1,1); zLine.setLinSpaced(numSide, pCenter(2)-nDx*dz, pCenter(2)+nDx*dz);

	dvec xi = xLine.replicate(numSide, 1);
	dmat yLine2 = yLine.replicate(1, numSide);
	yLine2.transposeInPlace();
	dvec yi(Map<VectorXd>(yLine2.data(), numSide*numSide));

	dmat pts(numNodes, 3);
	pts.col(0) = xi.replicate(numSide, 1);
	pts.col(1) = yi.replicate(numSide, 1);
	
	dmat zLine2 = zLine.replicate(1, numSide*numSide);
	zLine2.transposeInPlace();
	dvec zi(Map<VectorXd>(zLine2.data(), numNodes));
	pts.col(2) = zi;

	return pts;
}


#endif
