/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #1:
 * svmmath.cpp
 *		a routine for intersecting >2 lines (for vanishing point
 *		computation);
 *		routines for computing the homography for the reference
 *		plane and arbitrary polygons
 **************************************************************/

#pragma warning(disable : 4996)

#include "svmmath.h"
#include "jacob.h"
#include "vec.h"
#include <cstring>
#include <cstdio>
#include <assert.h>
#include <iostream>

#include "Eigen/Core"
#include "MinEig.h"

using namespace Eigen;
using namespace std;

//
// TODO 1: BestFitIntersect()
//		Given lines, the list of 3 or more lines to be intersected,
//		find the best fit intersection point.
//		See http://www-2.cs.cmu.edu/~ph/869/www/notes/vanishing.txt.
//
SVMPoint BestFitIntersect(const std::list<SVMLine> &lines, int imgWidth, int imgHeight)
{
  // check
  if (lines.size() < 2)
    {
      fprintf(stderr, "Not enough lines to compute the best fit.");
      abort();
    }

  SVMPoint bestfit;
  list<SVMLine>::const_iterator iter;

  // To accumulate stuff
  typedef Matrix<double, Dynamic, 3, RowMajor> Matrix3;

  int numLines = (int) lines.size();
  Matrix3 A = Matrix3::Zero(numLines, 3);
  
  // Transformation for numerical stability

  // Note: iterate through the lines list as follows:
  //		for (iter = lines.begin(); iter != lines.end(); iter++) {
  //			...iter is the pointer to the current line...
  //		}
  // Note: Function to find eigenvector with smallest eigenvalue is MinEig(A, eval, evec)
  //
  /******** BEGIN TODO ********/

  
printf("TODO: svmmath.cpp:61\n"); 
//fl_message("TODO: svmmath.cpp:61\n");

  double eigenvec[3];
  double eval;
int w_set_by_user=(imgWidth+imgHeight)/2;
vector<Vec3<double>> temp_line;
int i=0;
for (iter = lines.begin(); iter != lines.end(); iter++) {

	Vec3<double> a1;
    Vec3<double> a2;
	

	a1[0]=(iter->pnt1->u)-imgWidth/2;
	a1[1]=(iter->pnt1->v)-imgHeight/2;
	a1[2]=w_set_by_user;

	a2[0]=(iter->pnt2->u)-imgWidth/2;
	a2[1]=(iter->pnt2->v)-imgHeight/2;
	a2[2]=w_set_by_user;
	
	//a1[0]/=a1[2];
	//a1[1]/=a1[2];
 //   a1[2]/=a1[2];
	//
	//a2[0]/=a2[2];
	//a2[1]/=a2[2];
 //   a2[2]/=a2[2];
	
	Vec3<double> myline=cross(a1,a2);
	temp_line.push_back(myline);
	A(i,0)=myline[0];
	A(i,1)=myline[1];
	A(i,2)=myline[2];
	i++;
}
	
	MinEig(A, eval, eigenvec);
	bestfit.u=(eigenvec[0]/eigenvec[2])*w_set_by_user+imgWidth/2;
	bestfit.v=(eigenvec[1]/eigenvec[2])*w_set_by_user+imgHeight/2;
	bestfit.w=1;	

/******** END TODO ********/
 return bestfit;
}


//
// TODO 2: ConvertToPlaneCoordinate()
//		Given a plane defined by points, converts their coordinates into
//		a plane coordinate of your choise.
//              See the pdf titled "Homography from Polygon in R^3 to Image Plane",
//              whose link can be found from the project page.
//
//      The final divisors you apply to the u and v coordinates should be saved uScale and vScale
//
void ConvertToPlaneCoordinate(const vector<SVMPoint>& points, vector<Vec3d>& basisPts, double &uScale, double &vScale)
{
  int numPoints = points.size();
  //int numBasis = basisPts.size();
  /******** BEGIN TODO ********/
printf("numPoints=%d\n",numPoints);  

Vec4<double> p;
Vec4<double> q;
Vec4<double> r;
Vec4<double> ex;
Vec4<double> ey;
Vec4<double> s;
Vec4<double> t;

p[0]=points[0].X;
p[1]=points[0].Y;
p[2]=points[0].Z;
p[3]=points[0].W;

r[0]=points[1].X;
r[1]=points[1].Y;
r[2]=points[1].Z;
r[3]=points[1].W;

ex=p-r;
ex.normalize();

double exy;
double best_exy=100.0;
Vec4<double> temp_ey;
Vec4<double> temp_q;

for(int k=2;k<numPoints;k++){
temp_q[0]=points[k].X;
temp_q[1]=points[k].Y;
temp_q[2]=points[k].Z;
temp_q[3]=points[k].W;
temp_ey=temp_q-r;
temp_ey.normalize();
//exy=ex*temp_ey;
exy=ex[0]*temp_ey[0]+ex[1]*temp_ey[1]+ex[2]*temp_ey[2]+ex[3]*temp_ey[3];
if(exy<best_exy){
ey=temp_ey;
best_exy=exy;
q=temp_q;
}
}

Vec4<double> temp_qr;
double qrlength;
temp_qr=q-r;
qrlength=temp_qr[0]*ex[0]+temp_qr[1]*ex[1]+temp_qr[2]*ex[2]+temp_qr[3]*ex[3];
s[0]=qrlength*ex[0];
s[1]=qrlength*ex[1];
s[2]=qrlength*ex[2];
s[3]=qrlength*ex[3];
//s=((q-r)*ex)*ex;
ey=(q-r)-s;


double u_max=0;
double v_max=0;
double u_min=1e10;
double v_min=1e10;
double u;
double v;

for(int k=0;k<numPoints;k++){
SVMPoint temp2=points[k];
Vec4<double> tansferpoint;

tansferpoint[0]=temp2.X;
tansferpoint[1]=temp2.Y;
tansferpoint[2]=temp2.Z;
tansferpoint[3]=temp2.W;
tansferpoint=tansferpoint-r;
//if((tansferpoint-r)*ex>u_max){
//	u_max=(tansferpoint-r)*ex;
//}
//if((tansferpoint-r)*ex<u_min){
//	u_min=(tansferpoint-r)*ex;
//}
//if((tansferpoint-r)*ey>v_max){
//	v_max=(tansferpoint-r)*ey;
//}
//if((tansferpoint-r)*ey<v_min){
//	v_min=(tansferpoint-r)*ey;
//}
u=tansferpoint[0]*ex[0]+tansferpoint[1]*ex[1]+tansferpoint[2]*ex[2]+tansferpoint[3]*ex[3];
v=tansferpoint[0]*ey[0]+tansferpoint[1]*ey[1]+tansferpoint[2]*ey[2]+tansferpoint[3]*ey[3];
//u=tansferpoint*ex;
//v=tansferpoint*ey;
basisPts.push_back(Vec3d(u,v,1));
}

//uScale=u_max-u_min;
//vScale=v_max-v_min;

//for(int k=0;k<numPoints;k++){
//SVMPoint temp2=points[k];
//Vec4<double> tansferpoint;
//tansferpoint[0]=temp2.X;
//tansferpoint[1]=temp2.Y;
//tansferpoint[2]=temp2.Z;
//double ui=(tansferpoint-r)*ex;
//double vi=(tansferpoint-r)*ey;
//basisPts.push_back(Vec3d(ui,vi,1));
//}
	/******** END TODO ********/
}



//
// TODO 3: ComputeHomography()
//		Computes the homography H from the plane specified by "points" to the image plane,
//		and its inverse Hinv.
//		If the plane is the reference plane (isRefPlane == true), don't convert the
//		coordinate system to the plane. Only do this for polygon patches where
//		texture mapping is necessary.
//		Coordinate system conversion is to be implemented in a separate routine
//		ConvertToPlaneCoordinate.
//		For more detailed explaination, see the pdf titled
//              "Homography from Polygon in R^3 to Image Plane", whose link can be found from
//              the project page.
//
void ComputeHomography(CTransform3x3 &H, CTransform3x3 &Hinv, const vector<SVMPoint> &points, vector<Vec3d> &basisPts, bool isRefPlane)
{
  int i;
  int numPoints = (int) points.size();
  assert( numPoints >= 4 );

  basisPts.clear();
  if (isRefPlane) // reference plane
    {
      for (i=0; i < numPoints; i++)
        {
          Vec3d tmp = Vec3d(points[i].X, points[i].Y, points[i].W); // was Z, not W
          basisPts.push_back(tmp);
        }
    }
  else // arbitrary polygon
    {
      double uScale, vScale; // unused in this function
      ConvertToPlaneCoordinate(points, basisPts, uScale, vScale);
    }

  // A: 2n x 9 matrix where n is the number of points on the plane
  //    as discussed in lecture
  int numRows = 2 * numPoints;
  const int numCols = 9;

  typedef Matrix<double, Dynamic, 9, RowMajor> MatrixType;
  MatrixType A = MatrixType::Zero(numRows, numCols);

  /******** BEGIN TODO ********/
  
printf("TODO: svmmath.cpp:187\n"); 

for(int i=0;i<numPoints;i++){
	Vec3d points_prime=basisPts[i];
	A(2*i,0)=points_prime[0];
	A(2*i,1)=points_prime[1];
	A(2*i,2)=1;

	A(2*i+1,3)=points_prime[0];
	A(2*i+1,4)=points_prime[1];
	A(2*i+1,5)=1;

	A(2*i,6)=-points[i].u*points_prime[0];
	A(2*i,7)=-points[i].u*points_prime[1];
	A(2*i,8)=-points[i].u;
	A(2*i+1,6)=-points[i].v*points_prime[0];
	A(2*i+1,7)=-points[i].v*points_prime[1];
	A(2*i+1,8)=-points[i].v;
}

 double eval, h[9];
 MinEig(A, eval, h);

 H[0][0] = h[0];
 H[0][1] = h[1];
 H[0][2] = h[2];

 H[1][0] = h[3];
 H[1][1] = h[4];
 H[1][2] = h[5];

 H[2][0] = h[6];
 H[2][1] = h[7];
 H[2][2] = h[8];

 /******** END TODO ********/

 // compute inverse of H
 if (H.Determinant() == 0)
   fl_alert("Computed homography matrix is uninvertible \n");
 else
 Hinv = H.Inverse();

 int ii;
 printf("\nH=[\n");
 for (ii=0; ii<3; ii++)
   printf("%e\t%e\t%e;\n", H[ii][0]/H[2][2], H[ii][1]/H[2][2], H[ii][2]/H[2][2]);
 printf("]\nHinv=[\n");

 for (ii=0; ii<3; ii++)
   printf("%e\t%e\t%e;\n", Hinv[ii][0]/Hinv[2][2], Hinv[ii][1]/Hinv[2][2], Hinv[ii][2]/Hinv[2][2]);

 printf("]\n\n");
}
