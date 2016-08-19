#ifndef COMMON_HPP
#define COMMON_HPP


//Power function for integers
inline int powInt(int x, unsigned n)
{
        int y = 1;
        for (unsigned i = 0; i < n; i++)
            y *= x;
        return y;
}
//Index function for Mj matrix
inline int id(int x, int y, int z)
{
	return x + y*MAX_ID + z*MAX_ID*MAX_ID;
}


//Uniformly distributed double type random number generator
inline double doubleUniformRand(double fMin, double fMax)
{
    double f = ((double)rand() / (double)(RAND_MAX));
    return fMin + f * (fMax - fMin);
}



//Extra Eigen operations

//Remove rows or columns (using matrix row/column numbering starting at 0)
inline void removeRow(dmat& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

inline void removeColumn(dmat& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

#endif
