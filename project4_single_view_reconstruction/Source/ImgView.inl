/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #2:
 * ImgView.inl (included from ImgView.cpp)
 *		contains routines for computing the 3D position of points
 ***************************************************************/

//
// TODO 4: sameXY()
//		Computes the 3D position of newPoint using knownPoint
//		that has the same X and Y coordinate, i.e. is directly
//		below or above newPoint.
//		See lecture slide on measuring heights.
//
// HINT1: make sure to dehomogenize points when necessary
// HINT2: there is a degeneracy that you should look out for involving points already in line with the reference
// HINT3: make sure to get the sign of the result right, i.e. whether it is above or below ground
void ImgView::sameXY()
{
  if (pntSelStack.size() < 2)
    {
        printf("Not enough points on the stack.");
		//fl_alert("Not enough points on the stack.");
      return;
    }
  
  SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
  SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];
  
  if( !knownPoint.known() )
    {
        printf("Can't compute relative values for unknown point.");
		//fl_alert("Can't compute relative values for unknown point.");
      return;
    }
  
  if( refPointOffPlane == NULL )
    {
        printf("Need to specify the reference height first.");
		//fl_alert("Need to specify the reference height first.");
      return;
    }
  
  /******** BEGIN TODO ********/
  
  // See the lecture note on measuring heights
  // using a known point directly below the new point.
  
  // printf("sameXY() to be implemented!\n");
  
  
  
printf("TODO: ImgView.inl:49\n"); 

    newPoint.X = knownPoint.X;
	newPoint.Y = knownPoint.Y;
	newPoint.W = 1.0;

	Mat3d homographyMatrix = Mat3d(H[0][0], H[0][1], H[0][2], H[1][0], H[1][1], H[1][2], H[2][0], H[2][1], H[2][2]);

	double height_mag;
	Vec3d horizon;
	
	Vec3d b, r, t;
	Vec3d t_0, b_0; 
	Vec3d v, t_1, v_x, v_y, v_z;

	Vec3d unconverted_b, unconverted_b_0;
	Vec3d l1, l2, l3, l4;
	double dp_parallel_check;

	r = Vec3d(refPointOffPlane->u/refPointOffPlane->w, refPointOffPlane->v/refPointOffPlane->w, 1.0);  
	t_0 = Vec3d(newPoint.u/newPoint.w, newPoint.v/newPoint.w, 1.0);

	v_x = Vec3d(xVanish.u, xVanish.v, xVanish.w);
	v_y = Vec3d(yVanish.u, yVanish.v, yVanish.w);
	v_z = Vec3d(zVanish.u, zVanish.v, zVanish.w);

	unconverted_b_0 = Vec3d(knownPoint.X, knownPoint.Y, 1.0);
	unconverted_b = Vec3d(refPointOffPlane->X, refPointOffPlane->Y, 1.0);

	b_0 = homographyMatrix*unconverted_b_0;
	b = homographyMatrix*unconverted_b;

	v = cross(cross(b, b_0), cross(v_x, v_y));
	t = cross(cross(v, t_0), cross(r, b));

	if (v[0] == 0 && v[1] == 0 && v[2] == 0){
		t = t_0 + b - b_0;
	}


	t[0] = t[0] / t[2]; 
	r[0] = r[0] / r[2]; 
	v_z[0] = v_z[0] / v_z[2]; 
	b[0] = b[0] / b[2]; 

	t[1] = t[1] / t[2]; 
	r[1] = r[1] / r[2]; 
	v_z[1] = v_z[1] / v_z[2]; 
	b[1] = b[1] / b[2]; 

	t[2] = 1.0;
	r[2] = 1.0;
	v_z[2] = 1.0;
	b[2] = 1.0;

	l1 = t - b;
	l2 = v_z - r;
	l3 = r - b;
	l4 = v_z - t;

	height_mag = ( l1.length()*l2.length() )/( l3.length()*l4.length() )*referenceHeight; 

	dp_parallel_check = (t_0 - b_0)[0] * (v_z - b)[0] + (t_0 - b_0)[1] * (v_z - b)[1] + (t_0 - b_0)[2] * (v_z - b)[2]; 
	if (dp_parallel_check >= 0){
		newPoint.Z = -height_mag;
	}
	else{
		newPoint.Z = height_mag;
	}
	/******** END TODO ********/
 
 newPoint.known(true);
 
 printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );
 
 redraw();
}



//
// TODO 5: sameZPlane()
//		Compute the 3D position of newPoint using knownPoint
//		that lies on the same plane and whose 3D position is known.
//		See the man on the box lecture slide.
//		If newPoint is on the reference plane (Z==0), use homography (this->H, or simply H) directly.
//
// HINT: For this function, you will only need to use the three vanishing points and the reference homography 
//       (in addition to the known 3D location of knownPoint, and the 2D location of newPoint)
void ImgView::sameZPlane()
{
  if (pntSelStack.size() < 2)
    {
      fl_alert("Not enough points on the stack.");
      return;
    }
  
  SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
  SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];
  
  if( !knownPoint.known() )
    {
      fl_alert("Can't compute relative values for unknown point.");
      return;
    }
  
  /******** BEGIN TODO ********/
  
printf("TODO: ImgView.inl:142\n"); 

newPoint.Z = knownPoint.Z;
newPoint.W = 1.0;

Mat3d homographyMatrix = Mat3d(H[0][0], H[0][1], H[0][2], H[1][0], H[1][1], H[1][2], H[2][0], H[2][1], H[2][2]);
Mat3d homographyInverseMatrix = Mat3d(Hinv[0][0], Hinv[0][1], Hinv[0][2], Hinv[1][0], Hinv[1][1], Hinv[1][2], Hinv[2][0], Hinv[2][1], Hinv[2][2]);

Vec3d t1,m0,b0;
Vec3d b1, v, vx, vy, vz;

t1 = Vec3d(knownPoint.u, knownPoint.v, knownPoint.w);
m0 = Vec3d(newPoint.u, newPoint.v, newPoint.w);
b1 = homographyMatrix * Vec3d(knownPoint.X, knownPoint.Y, 1.0);

vx = Vec3d(xVanish.u, xVanish.v, xVanish.w);
vy = Vec3d(yVanish.u, yVanish.v, yVanish.w);
vz = Vec3d(zVanish.u, zVanish.v, zVanish.w);

v = cross(cross(m0, t1), cross(vx, vy));

if (knownPoint.Z == 0){
	b0 = Vec3d(newPoint.u, newPoint.v, newPoint.w);
}
else{
	if (v[0] == 0 && v[1] == 0 && v[2] == 0){
		b0 = cross(b1, cross(m0, t1));
	}
	else{
		b0 = cross(cross(b1, v), cross(m0, vz));
	}
}

b0[0] = b0[0] / b0[2];
b0[1] = b0[1] / b0[2];
b0[2] = 1.0;

m0 = homographyInverseMatrix*b0;
m0[0] = m0[0]/m0[2];
m0[1] = m0[1]/m0[2];

newPoint.X = m0[0];
newPoint.Y = m0[1];

	/******** END TODO ********/
 
 newPoint.known(true);
 
 printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );
 
 redraw();
}

