#include<stdio.h>
#include<stdbool.h>
#include<math.h>


//COLLIDE_SATURNE: Checks if the asteroid collides with Saturne.
int collide_saturne(double time_step, double x, double y, double Vx, double Vy,double ast_r, double x0, double y0, double r){

  /*
  *Let's check if the segment (distance traveled by the ast.) is long enought so that it CAN (maybe) cross the circle.
  *To do so we are going to compare the length of the segment with the distance between
  *the EDGE of the circle and the "beginig point" of the segment
  */

  double segLen = time_step*sqrt(pow(Vx,2)+pow(Vy,2));
  double circleDist = sqrt(pow(x0-x,2)+pow(y0-y,2))-r-ast_r;

  if (segLen >= circleDist){

    /*
    *Let's consider the equation of the right going trough the center of the circle & perpendicular to the first segment.
    It's length is equal to dist (see below, proof on the internet)
    If this distance is smaller than the radius, it means that they collided.
    */

    //Check if the asteroid is going straight down/up (just in case)
    if (fabs(Vx)<1e-200){

      printf("Cannot divide by 0!\n");
      return -1;

    }

    //Determine the equation of the right, associated to the first segment (the ax+by+c=0 right)
    double a = Vy/Vx;
    double b = -1;
    double c = y;

    //Dist is the distance between the center of the circle and this right.
    double dist = fabs(a*x0 + b*y0 + c)/sqrt(pow(a,2)+pow(b,2));

    //If dist. is smaller than the sum of the radiuses, they collided!
    if (dist <= r + ast_r) return 1;
    else return 0;

  }

  else return 0;

}


/*
*To check if two asteroids are going to collide (between towe points), we are goign to "project" the steroid
*in ttime (projected_x = x + Vx*projected_t), and we are going to check if they end up at the same position,
*at the same time. With some reflexion + internet, we concluded that we needed to solve a polynomial equation (2 degree):
* if there is 0/1 root -> there is no collision, if there are 2, the smallest one is the time at which they collide.
*To solve this equation we are not going to use any num. method as we already know how to (easilly) find the root.
*
*To calculate the resulting velocities of the collision between two asteroids, we are going to use
*the conservation of momentum. This will use a system of 4 equation: we will solve it using the Gauss-Jordan method.
*/

//POLY_ROOT: Returns -1 if there is no/one root, returns the lowest of the two roots otherwise.
double poly_root(double a, double b, double c){

  double delta = pow(b,2)-a*c;

  if (delta <= 0) return -1;
  else{

    double dr = sqrt(delta);
    double t = fmin((-b-dr)/a,(-b+dr)/a);

    if(t<0) return -1;
    else return t;

  }

}

//ROW_MULTIPLY: Multiplies a row by another one with a coef. alpha.
void row_multiply(int N, double M[N][N], int R1, int R2, double alpha){

  if (R1 < N && R2 < N && R1 >= 0 && R2 >= 0){

    for (int i=0; i<N; i++){

      M[R1][i] += alpha * M[R2][i];

    }

  }

  else {

    return;

  }

}

//ROW_SWAP: Swaps two row in a Matrice
void row_swap(int N, double M[N][N], int R1, int R2){

  //VARIABLES
  double aux[N];

  if (R1 < N && R2 < N && R1 >= 0 && R2 >= 0){

    for (int i=0; i<N; i++){

      aux[i] = M[R1][i];
      M[R1][i] = M[R2][i];
      M[R2][i] = aux[i];

    }

  }

  else {

    return;

  }

}

//SOLVE: Solves a system of N equations/N variables (for this I just did the same thing as in the TP4)
void solve(int N, double A[N][N], double b[N]){

  //VARIABLE
  double c, sum=0;
  int tmp=0, tmpRow=0;

  //SIMPLIFICATION
  for (int k=0; k<N; k++){

	  if (A[k][k]==0){

	    //Looking for the greatest A[i][k]
	    for (int i=k; i<N;i++){

        if(A[i][k]>tmp){

          tmp=A[i][k];
          tmpRow = i;

        }

	    }

    //At this point tmpRow is the row at which A[i][k] is the greatest > we swap both row
    row_swap(N, A, tmpRow, k);

	  }

	  else{

      for (int i=k+1; i<N; i++){

        c = A[i][k]/A[k][k];
        //Changing the row n°i
        row_multiply(N, A, i, k, -c);
        //Chaging the b vector as well
        b[i]-=c*b[k];

      }

    }

  }

  //SOLUTION
  for (int i=N-1; i>=0; i--){

    //Sum o f row i
    sum = 0;
    for (int j=i+1; j<N; j++){

      sum += A[i][j]*b[j];

    }

    //Next, we solve
    b[i]=b[i]-sum;
    b[i]=b[i]/A[i][i];

  }

}

//COLLIDE_AST: Check if an asteroid collides with another. If yes, it calculates the new position/velocities
double collide_ast(double time_step, double A1[4], double r1, double m1,double A2[4], double r2, double m2){

  //Variables
  //Restitution coef.
  double e = 1;
  //Position/Velocity at the current time
  double DVx = A1[2] - A2[2];
  double DVy = A1[3] - A2[3];
  double Dx = A1[0] - A2[0];
  double Dy = A1[1] - A2[1];

  //First let's check if they collide (if yes we create a new point)
  double a = pow(DVx,2)+pow(DVy,2);
  double b = Dx*DVx+Dy*DVy;
  double dmin = r1 +r2;
  double c = pow(Dx,2)+pow(Dy,2)-pow(dmin,2);
  double collision_time = poly_root(a, b, c);

  //If there is no collision in the "current" time step, no need to continue
  if (collision_time < 0 || collision_time > time_step) return -1;

  //If we are still here, they collide at t = curr_t + collision_time, this means we can deduce a new position
  A1[0] += A1[2]*collision_time;
  A1[1] += A1[3]*collision_time;
  A2[0] += A2[2]*collision_time;
  A2[1] += A2[3]*collision_time;

  //Tengential and perbendiculara directions
  double nx = A2[0]-A1[0];
  double ny = A2[1]-A1[1];
  double tx = ny;
  double ty = -nx;
  //The equation system
  double A[4][4] = {{m1,0,m2,0},{0,m1,0,m2},{tx,ty,0,0},{nx,ny,-nx,-ny}};
  double s[4] = {m2*A2[2]+m1*A1[2], m2*A2[3]+m1*A1[3], A1[2]*tx+A1[3]*ty, -e*((A1[2]-A2[2])*nx+(A1[3]-A2[3])*ny)};

  //Solving the system
  solve(4,A,s);

  //At this point s sould contain the new speeds of the asteroids, let's update them
  A1[2] = s[0];
  A1[3] = s[1];
  A2[2] = s[2];
  A2[3] = s[3];

  //Return the (delta) time before the collision
  return collision_time;

}
