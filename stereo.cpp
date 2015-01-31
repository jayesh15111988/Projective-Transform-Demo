#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#define PI 3.14
#define E 2.718281828

using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//



typedef SDoublePlane SDoubleMatrix;

//Function for Homographic transformation
SDoublePlane projective_transformation(const SDoublePlane &input, const SDoubleMatrix &homography)
{
double interpolated;
int rows = input.rows(), cols = input.cols(),i,j,k;
SDoublePlane output(rows,cols);

for(i=0;i<rows;i++)
{
	for(j=0;j<cols;j++)
	{
		double x1[3] = {j,i,1};
		double x[3] = {0,0,0};
		x[0] = homography[0][0]*x1[0] + homography[0][1]*x1[1] + homography[0][2]*x1[2];
		x[1] = homography[1][0]*x1[0] + homography[1][1]*x1[1] + homography[1][2]*x1[2];
		x[2] = homography[2][0]*x1[0] + homography[2][1]*x1[1] + homography[2][2]*x1[2];
		double temp = x[0];
		x[0] = x[1]/x[2]; //rows temp1 Y coordinate
		x[1] = temp/x[2]; //cols temp2 X coordinate
		int t1=x[0];
int t2=x[1];
		if(x[0] <= 0 || x[1] <= 0 || x[0] >= rows || x[1] >= cols)
		{
			output[i][j] = 0;
			continue;
		}
		else
		{
		int r1=floor(x[0]),r2=r1+1; 
		if(r2>=rows) {
			r2=rows-1;
		}
		
		int c1=floor(x[1]), c2=c1+1; 
		if(c2>=cols) {
			c2=cols-1;
		}
	
		//We are interpolating resulting image just to avoid any jagged edges or steps in the output. Application of interpolation is visibly strong to 
		//create a smoother image as an output

		interpolated = (input[r1][c1])*((x[1]-c1)*(r2-x[0]))/((r2-r1)*(c2-c1)) + (input[r2][c1])*(((x[0]-r1)*(c2-x[1]))/((r2-r1)*(c2-c1))) +
		(input[r1][c1])*(((c2-x[1])*(r2-x[0]))/((r2-r1)*(c2-c1))) + (input[r2][c2])*((x[0]-r1)*(x[1]-c1))/((r2-r1)*(c2-c1));
		output[i][j] =interpolated;
       }
	}
}

return output;
}




int main(int argc, char *argv[])
{
  if(argc < 2)
    {
      cerr << "usage: " << argv[0] << " image_file1 " << endl;
      return 1;
    }

  string input_filename1 = argv[1];
  string output_filename = "Projective_Transform.png";

//Third parameter is a output file name
	if(argc > 2) {
  		output_filename = argv[2];
	}


// apply a projection matrix
  // Given projection matrix
  SDoubleMatrix proj_matrix(3,3);
  proj_matrix[0][0] = 0.907; proj_matrix[0][1] = 0.258; proj_matrix[0][2] = -182;
  proj_matrix[1][0] = -0.153; proj_matrix[1][1] = 1.44; proj_matrix[1][2] = 58;
  proj_matrix[2][0] = -0.000306; proj_matrix[2][1] =  0.000731; proj_matrix[2][2] = 1;
  
//References
//http://www.bluebit.gr/matrix-calculator/calculate.aspx to get inverse of given matrix

  //inverse matrix
  SDoublePlane inv_proj_matrix(3,3);
  inv_proj_matrix[0][0]=1.124668581; inv_proj_matrix[0][1]=-0.314676504; inv_proj_matrix[0][2]=222.940924688;
  inv_proj_matrix[1][0]=0.108839051; inv_proj_matrix[1][1]=0.685058665; inv_proj_matrix[1][2]=-19.924695338;
  inv_proj_matrix[2][0]=0.000264587; inv_proj_matrix[2][1]=-0.00059706; inv_proj_matrix[2][2]=1.082784875;
  
  SDoublePlane image1 = SImageIO::read_png_file(input_filename1.c_str());

  SDoublePlane projected = projective_transformation(image1, inv_proj_matrix);
  SImageIO::write_png_file(output_filename.c_str(), projected, projected, projected);
	return 0;
}
