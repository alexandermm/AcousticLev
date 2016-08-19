#ifndef EIGENFILEIO_HPP
#define EIGENFILEIO_HPP

#include <fstream>
#include <string>
#include <iomanip>

template <typename Derived>
inline void write_matrix(const char* filename, const MatrixBase<Derived>& matrix)
{
	int numRows = matrix.rows();
	int numCols = matrix.cols();
	
	int r,c;

	ofstream myfile(filename);
	
	if (myfile.is_open())
	{
		myfile << scientific << showpoint << setprecision(16);		

  		for (r=0; r < numRows; r++) {
			for (c=0; c < numCols; c++) {
				myfile << " " << matrix(r,c);
			}		
			myfile << "\n";
		}    
	}
	else cout << "Unable to open file for writing Eigen matrix";

    	myfile.close();
}



inline dmat read_matrix(const char* filename)
{
	int i=0, maxj=0, j;

	dmat result(1,1);
	double value;
	
	string line;
  	ifstream myfile(filename);

	if (myfile.is_open())
  	{
    		while (getline(myfile, line))
    		{
			istringstream iss(line);

			j=0;
			while (iss >> value)
   			{					      						
				result.conservativeResize(i+1, maxj+1);				
				result(i,j) = value;
				
				if(i == 0) maxj++;
				j++;				
   			}

			if(i == 0) maxj--;		
			i++;
    		}

    		myfile.close();
  	}
  	else cout << "Unable to open file for reading Eigen matrix";
	
	return result;
}

#endif
