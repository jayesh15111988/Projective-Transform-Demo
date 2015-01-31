#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#define PI 3.14
#define E 2.718281828

//prototype of function to find minimum of d value

int minvalue(int,int);
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


// Rotate an image about its center by the specified angle, using a nearest-neighbor approach
//  (i.e., no interpolation...)
//
//Function to calculate simple disparity map. This function calculates the disparity map for each and every pixel in image
//After calculating disparity it sets the each pixel equal to its disparity value and gives result

SDoublePlane disparity_map(const SDoublePlane &input1, const SDoublePlane &input2)
{
int temp=0;
int sum=0;
int min=30000;
SDoublePlane dup(input1.rows(),input1.cols());
SDoublePlane result(input1.rows(), input1.cols());

int count=0;

for(int i=0;i<input1.rows();i++)
{
    min=30000;
    for(int j=0;j<input1.cols();j++)
    {
		min=30000;
        sum=0;
		for(int d=0;d<50;d++)
		{       
			sum=0;
			for(int k=i-1;k<i+2;k++)
			{
				for(int l=j-1;l<j+2;l++)
				{
					if(k>=0 && k<input1.rows() && l+d>=0 && l+d<input1.cols())
					{
						sum+=pow((input1[k][l]-input2[k][l+d]),2);
					}
				}
			}
			//Finding value of d corresponding to minimum "sum"
			if(sum<min)
			{
				min=sum;
				dup[i][j]=abs(input1[i][j]-input2[i][j+d]);
				result[i][j]=d;
			}
		}
	}
}

return result;
}

//Given the two intergers, function returns the minimum of two

int minvalue(int a,int b)
{
	if(a>b)
	{
		return b;
	}
	else
	{
		return a;
	}
}

//Evaluating stereo using Hidden Markov Model 
int fun(int d1,int d2)
{
if(pow((d1-d2),2)==0)
{
return 2000;
}
else
{
return 0;
}

}
SDoublePlane hmm_stereo(const SDoublePlane &input1, const SDoublePlane &input2)
{
int minimum;
int sum=0;
int min30=30000;

SDoublePlane result(input1.rows(), input1.cols());

int count=0;

//Calculating disparity value for each pixel in image

for(int i=0;i<input1.rows();i++)
{
    for(int j=0;j<input1.cols();j++)
    {
		min30=30000;
		for(int d=0;d<50;d++)
		{       
			sum=0;
			for(int k=i-1;k<i+2;k++)
			{
				for(int l=j-1;l<j+2;l++)
				{
					if(k>=0 && k<input1.rows() && l+d>=0 && l+d<input1.cols())
					{
						sum+=pow((input1[k][l]-input2[k][l+d]),2);
					}
				}
			}
			if(sum<min30)
			{
			min30=sum;
			result[i][j]=d;
			}
		}
	}
}


SDoublePlane messages1(input1.rows(),input1.cols());
SDoublePlane messages2(input2.rows(),input2.cols());

SDoublePlane store(input1.rows(),input1.cols());
SDoublePlane store1(input2.rows(),input2.cols());

SDoublePlane output(input1.rows(),input1.cols());

SDoublePlane saved1(input1.rows(),input1.cols());
SDoublePlane saved2(input1.rows(),input1.cols());

SDoublePlane sum1(input1.rows(),input1.cols());
SDoublePlane sum2(input2.rows(),input2.cols());

//calculating sum of neighbours for each pixel in both images

for(int i=0;i<input1.rows();i++)
{
	for(int j=0;j<input1.cols();j++)
	{
		sum1[i][j]=0;
		for(int k=i-1;k<i+2;k++)
		{
			for(int l=j-1;l<j+2;l++)
			{
				if(k>=0 && l>=0 && k<input1.rows() && l<input1.cols())
				{
					sum1[i][j]+=input1[k][l];
				}
			}
		}
	}
}


for(int i=0;i<input2.rows();i++)
{
	for(int j=0;j<input2.cols();j++)
	{
		sum2[i][j]=0;
		for(int k=i-1;k<i+2;k++)
		{
			for(int l=j-1;l<j+2;l++)
			{
				if(k>=0 && l>=0 && k<input2.rows() && l<input2.cols())
				{
					sum2[i][j]+=input2[k][l];
				}
			}
		}
	}
}

int start,term;
int min;

//Sending message from pixel(i,j) to the pixel (i,j+1) at the right

for(int i=0;i<input1.rows();i++)
{
	start=i;
	for(int k=0;k<input1.cols()-1;k++)
	{
		min=30000;
		for(int d1=0;d1<50;d1++)
		{
			if(k==0)
			{
				term=0;
			}
			else
			{
				term=store[start][k-1];
			}
			for(int d2=0;d2<50;d2++)
			{
				messages1[start][k]=pow((sum1[start][k]-sum2[start][k+d2]),2)+abs(d1-d2)+term;
				if(messages1[start][k]<min)
				{
					saved1[start][k]=d1;
					min=messages1[start][k];
					store[start][k]=min;
				}
			}
		}	      
		saved1[start][k]=saved1[start][k];
	}
}

//Sending message to the pixel at left 

for(int i=0;i<input1.rows();i++)
{
	start=i;
	for(int k=input1.cols()-1;k>=1;k--)
	{
		min=30000;
		for(int d1=0;d1<50;d1++)
		{
			if(k==input1.cols()-1)
			{
				term=0;
			}
			else
			{
				term=store1[start][k+1];
			}
			
			for(int d2=0;d2<50;d2++)
			{
				messages2[start][k]=pow((sum1[start][k]-sum2[start][k+d1]),2)+abs(d1-d2)+term;
				if(messages2[start][k]<min)
				{
					saved2[start][k]=d2;
					min=messages2[start][k];
					store1[start][k]=min;
				}
			}
		}
		saved2[start][k]=saved2[start][k];
	}
}

//Calculating minimum value of d among the three values viz, disparity value of pixel, message sent by previous pixel and message sent by next pixel

for(int i=0;i<input1.rows();i++)
{
	for(int j=0;j<input1.cols();j++)
	{
		if(j+1<input1.cols()&& j-1>=0)
		{
			output[i][j]=minvalue(saved1[i][j],minvalue(saved2[i][j],result[i][j]));
		}
	}
}

return output;
}

//Stereo matching using Markov Random fields for better result

SDoublePlane mrf_stereo(const SDoublePlane &input1, const SDoublePlane &input2)
{

int sum=0;
int min30=30000;
SDoublePlane dup(input1.rows(),input1.cols());
SDoublePlane result(input1.rows(), input1.cols());

int count=0;

//Calculating disparity value for each pixel

for(int i=0;i<input1.rows();i++)
{
    for(int j=0;j<input1.cols();j++)
    {
		min30=30000;
		sum=0;
		for(int d=0;d<50;d++)
		{       
			sum=0;
			for(int k=i-1;k<i+2;k++)
			{
				for(int l=j-1;l<j+2;l++)
				{
					if(k>=0 && k<input1.rows() && l+d>=0 && l+d<input1.cols())
					{
						sum+=pow((input1[k][l]-input2[k][l+d]),2);
					}
				}
			}
			if(sum<min30)
			{
			min30=sum;
			dup[i][j]=abs(input1[i][j]-input2[i][j+d]);
			result[i][j]=d;
			}
		}
	}
}


SDoublePlane sum1(input1.rows(),input1.cols());
SDoublePlane sum2(input2.rows(),input2.cols());

SDoublePlane messages1(input1.rows(),input1.cols());
SDoublePlane messages2(input2.rows(),input2.cols());

SDoublePlane store(input1.rows(),input1.cols());
SDoublePlane store1(input2.rows(),input2.cols());
SDoublePlane store2(input2.rows(),input2.cols());
SDoublePlane store3(input2.rows(),input2.cols());

SDoublePlane output(input1.rows(),input1.cols());

SDoublePlane saved1(input1.rows(),input1.cols());
SDoublePlane saved2(input1.rows(),input1.cols());
SDoublePlane saved3(input1.rows(),input1.cols());
SDoublePlane saved4(input1.rows(),input1.cols());



for(int i=0;i<input1.rows();i++)
{
	for(int j=0;j<input1.cols();j++)
	{
		sum1[i][j]=0;
		for(int k=i-1;k<i+2;k++)
		{
			for(int l=j-1;l<j+2;l++)
			{
				if(k>=0 && l>=0 && k<input1.rows() && l<input1.cols())
				{
					sum1[i][j]+=input1[k][l];
				}
			}
		}
	}
}


for(int i=0;i<input2.rows();i++)
{
	for(int j=0;j<input2.cols();j++)
	{
		sum2[i][j]=0;
		for(int k=i-1;k<i+2;k++)
		{
			for(int l=j-1;l<j+2;l++)
			{
				if(k>=0 && l>=0 && k<input2.rows() && l<input2.cols())
				{
					sum2[i][j]+=input2[k][l];
				}
			}
		}
	}
}

int start,term;
int min;

//Sending message to the right neighbour

for(int i=0;i<input1.rows();i++)
{
	start=i;
	for(int k=0;k<input1.cols()-1;k++)
	{
		min=30000;
		for(int d1=0;d1<50;d1++)
		{
			if(k==0)
			{
				term=0;
			}
			else
			{
				term=store[start][k-1];
			}
			for(int d2=0;d2<50;d2++)
			{
if(k+d2<input1.cols())
{
				messages1[start][k]=pow((sum1[start][k]-sum2[start][k+d2]),2)+pow((d1-d2),2)+term;
				if(messages1[start][k]<min)
				{
					saved1[start][k]=d1;
					min=messages1[start][k];
					store[start][k]=min;
				}
			}}
		}
	}
}

//Sending message to the neighbour at bottom

for(int i=0;i<input1.cols();i++)
{
	start=i;
	for(int k=0;k<input1.rows()-1;k++)
	{
		min=30000;
		for(int d1=0;d1<50;d1++)
		{
			if(k==0)
			{
				term=0;
			}
			else
			{
				term=store2[k-1][start];
			}
			for(int d2=0;d2<50;d2++)
			{
			if(k+d2<input1.rows())
			{
				messages1[k][start]=pow((sum1[k][start]-sum2[k+d2][start]),2)+pow((d1-d2),2)+term;
				if(messages1[k][start]<min)
				{
					saved3[k][start]=d1;
					min=messages1[k][start];
					store2[k][start]=min;
				}
			}
			}
		}
		saved3[k][start]=saved3[k][start];
	}
}

//Sending message to the neighbour to the left

for(int i=0;i<input1.rows();i++)
{
	start=i;
	for(int k=input1.cols()-1;k>=1;k--)
	{
		min=30000;
		for(int d1=0;d1<50;d1++)
		{
			if(k==input1.cols()-1)
			{
				term=0;
			}
			else
			{
				term=store1[start][k+1];
			}
			for(int d2=0;d2<50;d2++)
			{
if(k+d2<input1.cols())
{
				messages2[start][k]=pow((sum1[start][k]-sum2[start][k+d2]),2)+pow((d1-d2),2)+term;
				if(messages2[start][k]<min)
				{
					saved2[start][k]=d2;
					min=messages2[start][k];
					store1[start][k]=min;
				}
			}}
		}
		saved2[start][k]=saved2[start][k];
	}
}

//Sending message to the neighbour at the top

for(int i=0;i<input1.cols();i++)
{
	start=i;
	for(int k=input1.rows()-1;k>=1;k--)
	{
		min=30000;
		for(int d1=0;d1<50;d1++)
		{
			if(k==input1.rows()-1)
			{
				term=0;
			}
			else
			{   
				term=store3[k+1][start];
			}
			for(int d2=0;d2<50;d2++)
			{
			if(k+d2<input1.rows())
			{
				messages2[k][start]=pow((sum1[k][start]-sum2[k+d2][start]),2)+pow((d1-d2),2)+term;
				if(messages2[k][start]<min)
				{
					saved4[k][start]=d2;
					min=messages2[k][start];
					store3[k][start]=min;
				}
			}			
			}	
		}
		saved4[k][start]=saved4[k][start];
	}
}


int minimum;
for(int i=0;i<input1.rows();i++)
{
	for(int j=0;j<input1.cols();j++)
	{
		if(i-1>=0 && i+1<input1.rows())
		{
			minimum=minvalue(minvalue(saved1[i][j],saved3[i][j]),minvalue(saved4[i][j],result[i][j]));
			output[i][j]=minimum;
		}
	}
}

return output;
}
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
		if(r2>=rows) r2=rows-1;
		int c1=floor(x[1]), c2=c1+1; 
		if(c2>=cols) c2=cols-1;
	
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
  if(argc != 4 && argc != 3)
    {
      cerr << "usage: " << argv[0] << " image_file1 image_file2 [gt_file]" << endl;
      return 1;
    }

  string input_filename1 = argv[1], input_filename2 = argv[2];
  string gt_filename;
  if(argc == 4)
	gt_filename = argv[3];

  // read in images and gt
  SDoublePlane image1 = SImageIO::read_png_file(input_filename1.c_str());
  SDoublePlane image2 = SImageIO::read_png_file(input_filename2.c_str());
  SDoublePlane gt;
  
  if(gt_filename != "")
  {
    gt = SImageIO::read_png_file(gt_filename.c_str());
    // gt maps are scaled by a factor of 3, undo this...
  /* for(int i=0; i<gt.rows(); i++)
     for(int j=0; j<gt.cols(); j++)
      gt[i][j] = gt[i][j] / 3.0;*/
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
  
  SDoublePlane projected = projective_transformation(image1, inv_proj_matrix);
  SImageIO::write_png_file("Projective_Transform.png", projected, projected, projected);


// do stereo using simple technique
 SDoubleMatrix disp = disparity_map(image1, image2);
  SImageIO::write_png_file("disp_simple.png", disp, disp, disp);

  // do stereo using hmm
  SDoubleMatrix disp2 = hmm_stereo(image1, image2);
  SImageIO::write_png_file("disp_hmm.png", disp2, disp2, disp2);
  
 //  do stereo using mrf
  SDoubleMatrix disp3 = mrf_stereo(image1, image2);
  SImageIO::write_png_file("disp_mrf.png", disp3, disp3, disp3);

  // Measure error with respect to ground truth, if we have it...
 if(gt_filename != "")
    {
	double err=0;
    	for(int i=0; i<gt.rows(); i++)
	for(int j=0; j<gt.cols(); j++)
	err += sqrt((disp[i][j] - gt[i][j])*(disp[i][j] - gt[i][j]));

    cout << "Simple stereo technique mean error = " << err/gt.rows()/gt.cols() << endl;

    err=0;
    for(int i=0; i<gt.rows(); i++)
	for(int j=0; j<gt.cols(); j++)
	err += sqrt((disp2[i][j] - gt[i][j])*(disp2[i][j] - gt[i][j]));

    cout<<"HMM stereo technique mean error = "<< err/gt.rows()/gt.cols() << endl;


    err=0;
	for(int i=0; i<gt.rows(); i++)
	for(int j=0; j<gt.cols(); j++)
	err += sqrt((disp3[i][j] - gt[i][j])*(disp3[i][j] - gt[i][j]));

    cout << "MRF stereo technique mean error = " << err/gt.rows()/gt.cols() << endl;

    }
	return 0;
}
