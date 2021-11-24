#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL2_gfxPrimitives.h>
#include "Collision.h"



//Acceleration (Saturne) : https://astronomia.fr/seminaires/annee1314/mecaniqueCeleste.php
//Collision between asteroid: https://www.f-legrand.fr/scidoc/docimg/sciphys/meca/collidisques/collidisques.html

//SDL (Display)
#define MAX_FPS 60
#define SCREEN_WIDTH 1000
#define SCREEN_HEIGHT 1000

//VARIABLES
double G = 6.67430e-11;
double M = 5.683e26;
double R = 5.8232e7;
unsigned long int Ttotal = 300e3;
unsigned long int nb_step = 1e3;
int NB_AST = 300;
double MAX_DIST = 1000e7;
bool DO_COLLISION = true;
bool DO_RESPAWN = false;

//ASTEROID DISTRIB VARIABLES
double MIN_AST_RADIUS = 5e3;
double MAX_AST_RADIUS = 10e4;

double MIN_AST_DIST = 90e6;
double MAX_AST_DIST = 130e6;

double AST_SHOWER_SPEED = 18995;
double AST_SHOWER_Y_SPEED_VARIATION = 1000;
double MIN_AST_SHOWER_Y = 101.5e6;
double MAX_AST_SHOWER_Y = 102e6;
double AST_SHOWER_X = 0;

double AST_SPEED_VARIATION = 0;


//AST: Structure representing an asteroid
typedef struct ast_s{

  double mass;
  double radius;
  double x, y;
  double Vx, Vy;
  double local_time;
  double time_before_next_collision;
  bool isEnabled;
  bool isMoon;

} ast_t;


//INIT_AST: Intialise a new ast with a given Radius, mass, x, y, Vx, Vy
ast_t init_ast(double m, double r, double x, double y, double Vx, double Vy){

  ast_t ast;
  ast.mass = m;
  ast.radius = r;
  ast.x = x;
  ast.Vx = Vx;
  ast.y = y;
  ast.Vy = Vy;
  ast.local_time = 0;
  ast.time_before_next_collision = -1;
  ast.isEnabled = true;
  ast.isMoon = false;

  return ast;

}

//RAND_INIT
ast_t rand_init_ast(int min_radius, int max_radius, int min_dist, int max_dist, double speed_variation){


  //Variables
  double mass, radius, x, y , Vx, Vy;
  double dist, speed;
  double mass_coef = 900;

  //Pick a radius
  radius = min_radius+rand()%(max_radius-min_radius);

  //Pick a mass
  mass = mass_coef*(4/3)*M_PI*pow(radius,3);

  //Pick a dist.
  dist = min_dist+rand()%(max_dist-min_dist);
  /*Now thanks to the circle equation, we deduce a x/y couple
  We pick a random x between [-dist; dist], thus deducing the corresponding y*/
  x = -dist+rand()%(2*(int)dist);
  /*Now if we want the y to have the same prob. of being >0 or <0, we pick a rand between [0;1]
  and choose a sign depending on the result*/
  if(rand()%2==1){
    y=sqrt(pow(dist,2)-pow(x,2));
  }
  else{
    y=-sqrt(pow(dist,2)-pow(x,2));
  }

  /*
  * Now we pick a velocity , but as we want it to be the tangent velocity, we are going to apply it
  * along the right perpendicular to the one between the centers of the ast. & Sarturn
  */

  speed = sqrt(G*M/dist);
  //Let's randomize the speed a bit
  //Randomize the speed a little bit
  if(2*speed_variation*speed > 1) speed = (speed-speed*speed_variation)+ rand() % (int)(2*speed_variation*speed);
  Vx = (-y/sqrt(pow(x,2)+pow(y,2)))*speed;
  Vy = (x/sqrt(pow(x,2)+pow(y,2)))*speed;

  //Now we initialize this asteroid.
  return init_ast(mass, radius, x, y, Vx, Vy);

}

//INIT_AST_SHOWER
ast_t init_ast_shower(int min_radius, int max_radius, int min_y, int max_y, double speed, double speed_variation, double ySpeed_variation){


  //Variables
  double mass, radius, x, y , Vx, Vy;
  double mass_coef = 900;

  //Pick a radius
  radius = min_radius+rand()%(max_radius-min_radius);

  //Pick a mass
  mass = mass_coef*(4/3)*M_PI*pow(radius,3);

  x = AST_SHOWER_X;
  y = min_y+rand()%(max_y - min_y);

  //Let's randomize the speed a bit
  if(2*speed_variation*speed > 1) speed = (speed-speed*speed_variation)+ rand() % (int)(2*speed_variation*speed);

  Vx = -speed;

  if (2*ySpeed_variation > 1) Vy = rand() % (int)(ySpeed_variation);
  if (rand()%2==1) Vy = -Vy;


  //Now we initialize this asteroid.
  return init_ast(mass, radius, x, y, Vx, Vy);

}

//INIT_DISTRIB: Intializes all the asteroids arround saturne with a uniform distribution
/*
 * n: the number of asteroids
 * ast: the array containing the asteroids
 * min_radius & max_radius: min and max radius of an asteroid
 * min_mass & max_mass: the min and max mass of an asteroids
 * min_speed & max_speed: the min and max speed of an asteroid
 * min_dist & max_dist: the min and max distance to the center of saturn for an asteroid (width of the ring)
 */
void init_distrib(int n , ast_t ast[n], int min_radius, int max_radius, int min_dist, int max_dist, double speed_variation){

  //Variables
  srand(time(NULL));
  double x,y,Vx,Vy, speed;

/***********************  SATELLITES IMPLEMENTATION   *******************************/

  // we decide to spawn the satellite on the x-axis, so y=0
  y=0;


  /***************** Titan, the biggest satellite of Saturn  *************************/
  x=2e8;
  speed = sqrt(G*M/x);
  Vx = (-y/sqrt(pow(x,2)+pow(y,2)))*speed;
  Vy = (x/sqrt(pow(x,2)+pow(y,2)))*speed;
  ast[0] = init_ast(6e25, 3.575e6, x, y, Vx, Vy);
  ast[0].isMoon = true;

  /***************  Rhéa   *******************************************/

  x=-2.5e8;
  speed = -sqrt(G*M/abs(x));
  Vx = (-y/sqrt(pow(x,2)+pow(y,2)))*speed;
  Vy = (x/sqrt(pow(x,2)+pow(y,2)))*speed;
  //ast[1] = init_ast(3.07e20, 7.643e5, x, y, Vx, Vy);
  ast[1] = init_ast(6e25, 7.643e6, x, y, Vx, Vy);
  ast[1].isMoon = true;

  /************  Japet  ***************************************/

  x=3.56e8;
  speed = sqrt(G*M/x);
  Vx = (-y/sqrt(pow(x,2)+pow(y,2)))*speed;
  Vy = (x/sqrt(pow(x,2)+pow(y,2)))*speed;
  //ast[2] = init_ast(8.06e20,7.356e5, x, y, Vx, Vy);
  ast[2].isMoon = true;


  /***************  Dioné   *******************************************/

  x=3.77e8;
  speed = sqrt(G*M/x);
  Vx = (-y/sqrt(pow(x,2)+pow(y,2)))*speed;
  Vy = (x/sqrt(pow(x,2)+pow(y,2)))*speed;
  //ast[3] = init_ast(9.55e19, 5.617e5, x, y, Vx, Vy);
  ast[3].isMoon = true;

  /***************  Théthys   *******************************************/

  x=2.95e8;
  speed = sqrt(G*M/x);
  Vx = (-y/sqrt(pow(x,2)+pow(y,2)))*speed;
  Vy = (x/sqrt(pow(x,2)+pow(y,2)))*speed;
  //ast[4] = init_ast(6.18e20, 5.33e5, x, y, Vx, Vy);
  ast[4].isMoon = true;


  /*************************************************************************************/

  //For every asteroid
  for(int i=5; i<n; i++){

    ast[i] = rand_init_ast(min_radius, max_radius, min_dist, max_dist, speed_variation);
    //ast[i] = init_ast_shower(min_radius, max_radius, MIN_AST_SHOWER_Y, MAX_AST_SHOWER_Y, AST_SHOWER_SPEED, speed_variation, AST_SHOWER_Y_SPEED_VARIATION);

  }

}

//ALL_DISABLED: returns true if all the asteroids are disabled
bool all_disabled(ast_t ast[], int size){

  for (int i=1; i<size; i++){

    if (ast[i].isEnabled==true) return false;

  }

  return true;

}


/*
* We found the expression of x'' thanks to the Second Newton law
* The m masses of the asteroid are cancelling each other
* x[0]: x and x[1]: x'
*/

//FX: returns the derivative of the x vector
void fx(double y[2], double x[2], double xp[2], ast_t ast[], int nb_ast, int curr_ast){

  //The distance tween two asteroids
  double dist;

  xp[0] = x[1];
  xp[1] = -G*M*x[0]/pow(sqrt(pow(x[0],2)+pow(y[0],2)),3);

  //For each asteroid, we add its corresponding force to the acceleration
  for (int i=0; i<nb_ast; i++){

    dist = sqrt(pow(x[0]-ast[i].x,2)+pow(y[0]-ast[i].y,2));

    //We don't want the asteroid to calculate the force it exerces on itself (or if they are too far appart )
    if (ast[i].isEnabled==false || curr_ast == i /*|| (dist >10e10 && ast[i].isMoon == false)*/) continue;

    //Otherwise we add the force (from the other asteroids) to the acceleration
    xp[1]+= -G*(ast[i].mass)*(x[0]-ast[i].x)/pow(dist,3);

  }

}

//FY: returns the derivative of the Y vector
void fy(double y[2], double x[2], double yp[2], ast_t ast[], int nb_ast, int curr_ast){

  //The distance between two asteroids
  double dist;

  yp[0] = y[1];
  yp[1] = -G*M*y[0]/pow(sqrt(pow(x[0],2)+pow(y[0],2)),3);

  //For each asteroid
  for (int i=0; i<nb_ast; i++) {

    dist = sqrt(pow(x[0]-ast[i].x,2)+pow(y[0]-ast[i].y,2));

    //It is not influenced by itself or asteroids beyond 5000 km
    if(ast[i].isEnabled==false || curr_ast == i /*|| (dist >10e10 && ast[i].isMoon == false)*/) continue;

    //Add the force to the acceleration
    yp[1]+= -G*(ast[i].mass)*(y[0]-ast[i].y)/pow(dist,3);


  }

}

//RK4: Give an estimation of the new x/y/Vx/Vy for a given asteroid and time step
ast_t RK4(ast_t ast, int j, ast_t ast_list[], int n, double time_step){

  //Variables
  double xTmp[2];
  double yTmp[2];
  double x[2], y[2];
  double kx1[2], kx2[2], kx3[2], kx4[2], ky1[2], ky2[2], ky3[2], ky4[2];


  //We take its x,y,Vx,Vy
  y[0]=ast.y;
  y[1]=ast.Vy;
  x[0]=ast.x;
  x[1]=ast.Vx;

  //Temporary variables (for the different k variables)
  xTmp[0]=x[0];
  xTmp[1]=x[1];
  yTmp[0]=y[0];
  yTmp[1]=y[1];

  //Estimation of KX1
  fx(yTmp, xTmp,kx1, ast_list, n, j);

  //Estimation of KY1
  fy(yTmp, xTmp, ky1, ast_list, n, j);

  xTmp[0]=x[0]+kx1[0]*time_step/2;
  xTmp[1]=x[1]+kx1[1]*time_step/2;
  yTmp[0]=y[0]+ky1[0]*time_step/2;
  yTmp[1]=y[1]+ky1[1]*time_step/2;

  //Estimation of KX2
  fx(yTmp, xTmp, kx2, ast_list, n, j);

  //Estimation of KY2
  fy(yTmp, xTmp, ky2, ast_list, n, j);

  xTmp[0]=x[0]+kx2[0]*time_step/2;
  xTmp[1]=x[1]+kx2[1]*time_step/2;
  yTmp[0]=y[0]+ky2[0]*time_step/2;
  yTmp[1]=y[1]+ky2[1]*time_step/2;

  //Estimation of KX3
  fx(yTmp, xTmp, kx3, ast_list, n, j);

  //Estimation of KY3
  fy(yTmp, xTmp, ky3, ast_list, n, j);

  xTmp[0]=x[0]+kx3[0]*time_step;
  xTmp[1]=x[1]+kx3[1]*time_step;
  yTmp[0]=y[0]+ky3[0]*time_step;
  yTmp[1]=y[1]+ky3[1]*time_step;

  //Estimation of KX4
  fx(yTmp, xTmp, kx4, ast_list, n, j);

  //Estimation of KY4
  fy(yTmp, xTmp, ky4, ast_list, n, j);

  //Estimation of X[i]
  ast.x = x[0]+time_step*(kx1[0]+kx2[0]*2+kx3[0]*2+kx4[0])/6;
  ast.Vx = x[1]+time_step*(kx1[1]+kx2[1]*2+kx3[1]*2+kx4[1])/6;

  //Estimation of Y[i]
  ast.y = y[0]+time_step*(ky1[0]+ky2[0]*2+ky3[0]*2+ky4[0])/6;
  ast.Vy = y[1]+time_step*(ky1[1]+ky2[1]*2+ky3[1]*2+ky4[1])/6;

  return ast;

}

//ERROR_OCCURED: Prints the error and exit the program when called (SDL)
void error_occured(char* mess){

    SDL_Log("ERROR: %s: %s\n", mess, SDL_GetError());
    SDL_Quit();
    exit(EXIT_FAILURE);

}

//RENDER_DISPLAY: Dispaly Saturn + every asteroid at a given frame (SDL)
void render_display(SDL_Renderer *rend,double scale, double xShift, double yShift, int n, ast_t ast[n]){

  //Clear the previous frame
  if (SDL_SetRenderDrawColor(rend, 0,0,0, SDL_ALPHA_OPAQUE) != 0) error_occured("Draw Color setup");
  if(SDL_RenderClear(rend)!=0) error_occured("Clearing previous frame");

  /*- - - - - - - - - - - - - DRAWING - - - - - - - - - - - - - - - */

  //Draw Saturne (Yellow)
  if (filledCircleRGBA(rend, SCREEN_WIDTH/2 + xShift, SCREEN_HEIGHT/2 + yShift, R*scale, 255, 250, 0, SDL_ALPHA_OPAQUE)!=0) error_occured ("Drawing Saturne");

  //Draw all the asteroids
  for (int i=0; i<n; i++){

    //We only display it if the asteroid is enabled
    if (ast[i].isEnabled == true){

      if (ast[i].isMoon) {

        if (filledCircleRGBA(rend, ast[i].x*scale+SCREEN_WIDTH/2+xShift, -ast[i].y*scale+SCREEN_HEIGHT/2+yShift, ast[i].radius*scale, 255, 50, (i*50)%250,SDL_ALPHA_OPAQUE)!=0) error_occured ("Drawing Circle");

      }

      else if (filledCircleRGBA(rend, ast[i].x*scale+SCREEN_WIDTH/2+xShift, -ast[i].y*scale+SCREEN_HEIGHT/2+yShift, ast[i].radius*scale, 255, 255, 255,SDL_ALPHA_OPAQUE)!=0) error_occured ("Drawing Circle");

    }

  }

  /*- - - - - - - - - - - - END DRAWING - - - - - - - - - - - - - - - */

  //Render the frame
  SDL_RenderPresent(rend);

}


//TRAJECTORY: Estimation of the solution of the differential equations taking collision & the different interactions into account
void trajectory(long unsigned int nb_step, ast_t ast[], int nb_ast){

  //Variables
  FILE * file = fopen("Hist.txt", "w+");
  double h=(double)Ttotal/nb_step;
  int nb_bin = 50, R_min=R, R_max= R+200e6;
  int hist[nb_bin];
  double dx = (R_max-R_min)/nb_bin;
  //(SDL)
  SDL_Window *window = NULL;
  SDL_Renderer *renderer = NULL;
  bool program_lauched = true;
  unsigned int frame_start = 0, frame_time = 0;
  double scale = 0.000002;
  int camera_xShift=0, camera_yShift=0;

  //Initialisation
  //(SDL)
  if(SDL_Init(SDL_INIT_VIDEO)!=0) error_occured("SDL initialization");
  if(SDL_CreateWindowAndRenderer(SCREEN_HEIGHT, SCREEN_WIDTH, 0, &window, &renderer)!=0) error_occured("Window/Renderer initialization");

  /*- - - - - - - - - - - - - - - - RUNNING - - - - - - - - - - - - - - - - - - - -*/

  //If all the asteroids are Disabled (aka they crashed on Saturne), no need to continue.
  while (all_disabled(ast, nb_ast)==false && program_lauched == true){

    //Variables (SDL)
    SDL_Event event;
    frame_start = SDL_GetTicks();

    //Check for events (SDL)
    while (SDL_PollEvent(&event)){

      //Gestion of the events
      switch (event.type){

        //If the user presses on the keyboard (any key)
        case SDL_KEYDOWN : {

          switch (event.key.keysym.sym)
          {

            //Zoom in
            case SDLK_a :
              if (scale <= 1e-3){

                scale += scale/10;
                camera_xShift += scale*camera_xShift;

              }
              continue;

            //Zoom out
            case SDLK_z :
              if(scale >= 1e-8){

                scale -= scale/10;
                camera_xShift -= scale*camera_xShift;

              }
              continue;

            //Move up
            case SDLK_UP :
              camera_yShift+=10e6*scale;
              continue;

            //Move down
            case SDLK_DOWN :
              camera_yShift-=10e6*scale;
              continue;

            //Move left
            case SDLK_LEFT :
              camera_xShift+=10e6*scale;
              continue;

            //Move left
            case SDLK_RIGHT :
              camera_xShift-=10e6*scale;
              continue;

            default:
              continue;

            }

        }

        //If the user Quits the program
        case SDL_QUIT : {

          program_lauched = false;
          break;

        }

        //Default (Any other event)
        default: break;

      }

    }

    //Let's reset the histogram for this time nb_step
    for (int i=0; i<nb_bin; i++){
      hist[i]=0;
    }

    //For every asteroid
    for (int j=0; j<nb_ast; j++){

      //We add one time step to everybody (after all the calculation at a given t)
      ast[j].local_time += h;

      //If the asteroid has reached the max time :we disable it
      if (ast[j].local_time > Ttotal) ast[j].isEnabled = false;

      //If the asteroid has colided with a planet or if it already reached the max time :we ignore it.
      if (ast[j].isEnabled==false) continue;

    }


    /*Let's check if the asteroid will collide with another one (if yes update information & plot a new point).
    dt represent the time before a given collision happens*/
    double dt;
    //These arrays represent an asteroid ([x, y, Vx, Vy]): This allows us to update the pos/speed IF they collide
    double A1[4];
    double A2[4];

    //For every asteroid, we check if there is a collsion
    for (int j=0; j<nb_ast;j++){

      for (int k=0; k<nb_ast;k++){

        //An asteroid never collides with itself (and if it already collided, we don't bother checking again)
        if (j==k || ast[j].time_before_next_collision >= 0 || ast[k].time_before_next_collision >= 0) continue;

        //Otherwise, we check if the asteroid collide with someone else
        else{

          A1[0] = ast[j].x;
          A1[1] = ast[j].y;
          A1[2] = ast[j].Vx;
          A1[3] = ast[j].Vy;

          A2[0] = ast[k].x;
          A2[1] = ast[k].y;
          A2[2] = ast[k].Vx;
          A2[3] = ast[k].Vy;

          //Check for collision
          dt = collide_ast(h, A1, ast[j].radius, ast[j].mass, A2, ast[k].radius, ast[k].mass);
          if (DO_COLLISION == false) dt=-1;
          //If there is a collsion, we update its x/y/Vx/Vy
          if (dt>=0){

            //If the asteroid j is a moon, we disable the ast. k
            if (ast[j].isMoon==true){

              ast[k].isEnabled=false;
              continue;

            }

            ast[j].time_before_next_collision = dt;
            ast[j].x = A1[0];
            ast[j].y = A1[1];
            ast[j].Vx = A1[2];
            ast[j].Vy = A1[3];

            ast[k].time_before_next_collision = dt;
            ast[k].x = A2[0];
            ast[k].y = A2[1];
            ast[k].Vx = A2[2];
            ast[k].Vy = A2[3];

            //Change the local time of that asteroid to set it when the collision occured (go back in time)
            ast[j].local_time -= h-ast[j].time_before_next_collision;
            ast[k].local_time -= h-ast[k].time_before_next_collision;

          }

        }

      }

    }

    //Now that we checked for collisions, we estimates the new positions/Speed
    for (int i=0; i<nb_ast; i++){

      //If they are enabled
      if (ast[i].isEnabled == true){

        //If it collided, the estimation is made over h-time_before_next_collision and not h
        if(ast[i].time_before_next_collision>=0){

          //New estimation for x/y/Vx/Vy @ t=t+h
          ast[i] = RK4(ast[i], i, ast, nb_ast, h-ast[i].time_before_next_collision);

          //Re-align the local time (t=t+time_before_next_collision -> t=t+h)
          ast[i].local_time += h-ast[i].time_before_next_collision;

          //Reset the dt for the next time step
          ast[i].time_before_next_collision = -1;

        }

        else {

          //Estimates (RK4) the new x/y/Vx/Vy of each asteroids
          ast[i] = RK4(ast[i], i, ast, nb_ast, h);

        }

        //Now, last thing we need to do is to check is they collide with Saturm / Or if it goes too far in space
        if (collide_saturne(h, ast[i].x, ast[i].y, ast[i].Vx, ast[i].Vy, ast[i].radius, 0, 0, R)==1 || sqrt(pow(ast[i].x,2) + pow(ast[i].y,2))>= MAX_DIST ) {

          printf("Ast. n°%d despawned! @time n°%.2lf\n", i, ast[i].local_time);

          //If there is a respawn
          if (DO_RESPAWN == true){

            //Let's respawn the asteroid
            ast[i]=rand_init_ast(MIN_AST_RADIUS, MAX_AST_RADIUS, MIN_AST_DIST, MAX_AST_DIST, AST_SPEED_VARIATION);
            //ast[i]=init_ast_shower(MIN_AST_RADIUS, MAX_AST_RADIUS, MIN_AST_SHOWER_Y, MAX_AST_SHOWER_Y, AST_SHOWER_SPEED, AST_SPEED_VARIATION, AST_SHOWER_Y_SPEED_VARIATION);

          }

          //If the asteroïd chrashes into the planet: we disable that asteroid
          else {

            ast[i].isEnabled = false;
            continue;

          }

        }

      }

    }

    //Render this time step
    render_display(renderer, scale, camera_xShift, camera_yShift, nb_ast, ast);

    //At the end of the frame, we check if we need to delay to adjust the FrameRate (SDL)
    frame_time = SDL_GetTicks() - frame_start;
    if ((int)(1000/MAX_FPS)>frame_time) SDL_Delay((int)(1000/MAX_FPS)-frame_time);

    //Now, let's fill the histogram
    double dist;
    int indice;
    double sum = 0;
    int nb=0;
    double moyenne=0;
    int isIn = 0;

    for (int n=0; n<nb_ast; n++){

      dist = sqrt(pow(ast[n].x,2)+pow(ast[n].y,2));

      /*
      if (dist<=R_max && dist>= R_min && ast[n].isEnabled==true){

        indice = (int)((dist-R_min)/dx);
        hist[indice]++;

      }
      */
      if (ast[n].y>=-R && ast[n].y<=R && ast[n].x<0 && ast[n].isMoon){
        isIn = 10;
      }


      if (ast[n].y>=-R && ast[n].y<=R && ast[n].x<0){

        sum += dist;
        nb ++;

      }

    }

    if (nb !=0) moyenne = sum/nb;


    //Let's print the Histogram in the file
    for (int m=0; m<nb_bin; m++){

      //fprintf(file, "%d\t%d\t%d\n", (int)ast[0].local_time, (int)(R_min+dx*m), hist[m]);
      fprintf(file, "%d\t%lf\t%d\n", (int)ast[0].local_time, moyenne, isIn);

    }

    fprintf(file, "\n");

  }


  /*- - - - - - - - - - - - - - - - END - - - - - - - - - - - - - - - - - -*/

  //Close everything after we are done (SDL)
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();
  fclose(file);

}

//MAIN
int main(int argc, char const *argv[]){

  //Variables
  int nb_ast = NB_AST;
  ast_t ast[nb_ast];

  //We initialize the rings
  init_distrib(nb_ast, ast, MIN_AST_RADIUS, MAX_AST_RADIUS, MIN_AST_DIST, MAX_AST_DIST, AST_SPEED_VARIATION);

  //Calculate the orbits
  trajectory(nb_step, ast, nb_ast);

	return 0;

}
