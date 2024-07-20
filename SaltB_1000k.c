/*=====================================================================
  mps.c   

    Sample program of the MPS method

    Moving Particle Semi-implicit Method, Academic Press, 2018
    ISBN: 9780128127797

    SALT B, 1000K DENSITY 2045.6, VISCOSITY 4.13E-6
=======================================================================*/
#include <stdio.h>
#include <math.h>

#define DIM                  2
#define PARTICLE_DISTANCE    0.01
#define DT                   0.0001
#define OUTPUT_INTERVAL      20

/* for three-dimensional simulation */
/*
#define DIM                  3
#define PARTICLE_DISTANCE    0.075
#define DT                   0.003
#define OUTPUT_INTERVAL      2 
*/

#define ARRAY_SIZE           5000
#define FINISH_TIME          4
#define KINEMATIC_VISCOSITY  (4.13E-6)
#define FLUID_DENSITY        2045.6 
#define GRAVITY_X  0      
#define GRAVITY_Y  -9.8
#define GRAVITY_Z  0.0      
#define RADIUS_FOR_NUMBER_DENSITY  (2.1*PARTICLE_DISTANCE) 
#define RADIUS_FOR_GRADIENT        (2.1*PARTICLE_DISTANCE) 
#define RADIUS_FOR_LAPLACIAN       (3.1*PARTICLE_DISTANCE) 
#define COLLISION_DISTANCE         (0.5*PARTICLE_DISTANCE)
#define THRESHOLD_RATIO_OF_NUMBER_DENSITY  0.97   
#define COEFFICIENT_OF_RESTITUTION 0.2
#define COMPRESSIBILITY (0.45E-9)
#define EPS             (0.01 * PARTICLE_DISTANCE)     
#define ON              1
#define OFF             0
#define RELAXATION_COEFFICIENT_FOR_PRESSURE 0.2
#define GHOST  -1
#define FLUID   0
#define WALL    2
#define DUMMY_WALL  3
#define GHOST_OR_DUMMY  -1
#define SURFACE_PARTICLE 1    
#define INNER_PARTICLE   0      
#define DIRICHLET_BOUNDARY_IS_NOT_CONNECTED 0 
#define DIRICHLET_BOUNDARY_IS_CONNECTED     1 
#define DIRICHLET_BOUNDARY_IS_CHECKED       2 
#define PI 3.14159265358979323846
void initializeParticlePositionAndVelocity_for2dim( void );
void initializeParticlePositionAndVelocity_for3dim( void );
void calculateConstantParameter( void );
void calculateNZeroAndLambda( void );
double weight( double distance, double re );
void mainLoopOfSimulation( void );
void calculateGravity( void );
void calculateViscosity( void );
void moveParticle( void );
void collision( void );
void calculatePressure( void );
void calculateParticleNumberDensity( void );
void setBoundaryCondition( void );
void setSourceTerm( void );
void setMatrix( void );
void exceptionalProcessingForBoundaryCondition( void );
void checkBoundaryCondition( void );
void increaseDiagonalTerm( void );
void solveSimultaniousEquationsByGaussEliminationMethod( void );
void removeNegativePressure( void );
void setMinimumPressure( void );
void calculatePressureGradient( void );
void moveParticleUsingPressureGradient( void );
void writeData_inProfFormat( void );
void writeData_inVtuFormat( void );
static double AccelerationX[ARRAY_SIZE];
static double AccelerationY[ARRAY_SIZE];
static double AccelerationZ[ARRAY_SIZE];
static int    ParticleType[ARRAY_SIZE];
static double PositionX[ARRAY_SIZE];
static double PositionY[ARRAY_SIZE];
static double PositionZ[ARRAY_SIZE];
static double VelocityX[ARRAY_SIZE];
static double VelocityY[ARRAY_SIZE];
static double VelocityZ[ARRAY_SIZE];
static double Pressure[ARRAY_SIZE];
static double ParticleNumberDensity[ARRAY_SIZE];
static int    BoundaryCondition[ARRAY_SIZE];
static double SourceTerm[ARRAY_SIZE];
static int    FlagForCheckingBoundaryCondition[ARRAY_SIZE];
static double CoefficientMatrix[ARRAY_SIZE * ARRAY_SIZE];
static double MinimumPressure[ARRAY_SIZE];
int    FileNumber;
double Time;  
int    NumberOfParticles;
double Re_forParticleNumberDensity,Re2_forParticleNumberDensity; 
double Re_forGradient,     Re2_forGradient; 
double Re_forLaplacian,    Re2_forLaplacian; 
double N0_forParticleNumberDensity;
double N0_forGradient;
double N0_forLaplacian;
double Lambda;
double collisionDistance,collisionDistance2;
double FluidDensity;
int nilai_y;

int main( void ) {

  printf("\n*** START MPS-SIMULATION ***\n");
  if( DIM == 2 ){
    initializeParticlePositionAndVelocity_for2dim();
  }else{
    initializeParticlePositionAndVelocity_for3dim();
  }
  calculateConstantParameter();
  mainLoopOfSimulation();
  printf("*** END ***\n\n");
  return 0;
  
}


void initializeParticlePositionAndVelocity_for2dim( void ){
  
  int iX, iY;
  int nX, nY;
  double x, y, z;
  int i = 0;
  int flagOfParticleGeneration;

  nX = (int)(1.0/PARTICLE_DISTANCE+120);  
  nY = (int)(0.6/PARTICLE_DISTANCE+120);
  for(iY= -95;iY<nY;iY++){
    for(iX= -95;iX<nX;iX++){
      x = PARTICLE_DISTANCE * (double)(iX);
      y = PARTICLE_DISTANCE * (double)(iY);
      z = 0.0; 
      flagOfParticleGeneration = OFF;

  /* dummy wall region, 2 lapis*/
      /*Active Core*/
      if( ((x>=-45*PARTICLE_DISTANCE+EPS))&&((x<=-11*PARTICLE_DISTANCE+EPS)) &&((y>=90*PARTICLE_DISTANCE+EPS))&&((y<=125*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=DUMMY_WALL;
	      flagOfParticleGeneration = ON;
      }
      /*Pipa vertikal active core*/
      if( ((x>=-35*PARTICLE_DISTANCE+EPS))&&((x<=-21*PARTICLE_DISTANCE+EPS)) &&((y>=56*PARTICLE_DISTANCE+EPS))&&((y<=94*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=DUMMY_WALL;
	      flagOfParticleGeneration = ON;
      }
      /*Pipa horizontal*/
      if( ((x>=-35*PARTICLE_DISTANCE+EPS))&&((x<=7*PARTICLE_DISTANCE+EPS)) &&((y>=54*PARTICLE_DISTANCE+EPS))&&((y<=68*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=DUMMY_WALL;
	      flagOfParticleGeneration = ON;
      }
      /*Pipa vertikal drain tank*/
      if( ((x>=-7*PARTICLE_DISTANCE+EPS))&&((x<=7*PARTICLE_DISTANCE+EPS)) &&((y>=28*PARTICLE_DISTANCE+EPS))&&((y<=67*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=DUMMY_WALL;
	      flagOfParticleGeneration = ON;
      }
      /*Drain Tank*/
      if( ((x>=-29*PARTICLE_DISTANCE+EPS))&&((x<=28*PARTICLE_DISTANCE+EPS)) &&((y>=-2*PARTICLE_DISTANCE+EPS))&&((y<=30*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=DUMMY_WALL;
	      flagOfParticleGeneration = ON;
      }



  /* wall region, 3 lapis */
      /*Active Core*/
      if( ((x>=-43*PARTICLE_DISTANCE+EPS))&&((x<=-13*PARTICLE_DISTANCE+EPS)) &&((y>=92*PARTICLE_DISTANCE+EPS))&&((y<=125*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=WALL;
	      flagOfParticleGeneration = ON;
      }
      /*Pipa vertikal active core*/
      if( ((x>=-33*PARTICLE_DISTANCE+EPS))&&((x<=-23*PARTICLE_DISTANCE+EPS)) &&((y>=56*PARTICLE_DISTANCE+EPS))&&((y<=94*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=WALL;
	      flagOfParticleGeneration = ON;
      }
      /*Pipa horizontal*/
      if( ((x>=-23*PARTICLE_DISTANCE+EPS))&&((x<=5*PARTICLE_DISTANCE+EPS)) &&((y>=56*PARTICLE_DISTANCE+EPS))&&((y<=66*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=WALL;
	      flagOfParticleGeneration = ON;
      }
      /*Pipa vertikal drain tank*/
      if( ((x>=-5*PARTICLE_DISTANCE+EPS))&&((x<=5*PARTICLE_DISTANCE+EPS)) &&((y>=28*PARTICLE_DISTANCE+EPS))&&((y<=56*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=WALL;
	      flagOfParticleGeneration = ON;
      }
      /*Drain Tank*/
      //Kotak atas
      if( ((x>=-27*PARTICLE_DISTANCE+EPS))&&((x<=26*PARTICLE_DISTANCE+EPS)) &&((y>=0*PARTICLE_DISTANCE+EPS))&&((y<=28*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=WALL;
	      flagOfParticleGeneration = ON;
      }
      //Silinder bawah
      if( ((x<=(1)*sqrt(y*70*PARTICLE_DISTANCE+EPS)))&&((x>=(-1)*sqrt(y*70*PARTICLE_DISTANCE+EPS))) &&((y>=0*PARTICLE_DISTANCE+EPS))&&((y<=10*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=WALL;
	      flagOfParticleGeneration = ON;
      }
      //Alas
      if( ((x>=-5*PARTICLE_DISTANCE+EPS))&&((x<=4*PARTICLE_DISTANCE+EPS)) &&((y>=0*PARTICLE_DISTANCE+EPS))&&((y<=3*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=WALL;
	      flagOfParticleGeneration = ON;
      }



  /* empty region */
    /*Active Core*/
    //atas
      if( ((x>=-40*PARTICLE_DISTANCE+EPS))&&((x<=-16*PARTICLE_DISTANCE+EPS)) &&((y>=112*PARTICLE_DISTANCE+EPS))&&((y<=125*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
    /*Pipa vertical active core*/
    //kosong tengah
      if( ((x>=-30*PARTICLE_DISTANCE+EPS))&&((x<=-26*PARTICLE_DISTANCE+EPS)) &&((y<=88*PARTICLE_DISTANCE+EPS))&&((y>=62*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
      //bawah lengkungan
      if( ((x>=-30*PARTICLE_DISTANCE+EPS))&&((x<=-23*PARTICLE_DISTANCE+EPS)) &&((y<=62*PARTICLE_DISTANCE+EPS))&&((y>=61*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
      if( ((x>=-29*PARTICLE_DISTANCE+EPS))&&((x<=-23*PARTICLE_DISTANCE+EPS)) &&((y<=61*PARTICLE_DISTANCE+EPS))&&((y>=60*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
      if( ((x>=-28*PARTICLE_DISTANCE+EPS))&&((x<=-23*PARTICLE_DISTANCE+EPS)) &&((y<=60*PARTICLE_DISTANCE+EPS))&&((y>=59*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
      //atas lengkungan
      if( ((x>=-26*PARTICLE_DISTANCE+EPS))&&((x<=-23*PARTICLE_DISTANCE+EPS)) &&((y<=63*PARTICLE_DISTANCE+EPS))&&((y>=62*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
      if( ((x>=-26*PARTICLE_DISTANCE+EPS))&&((x<=-24*PARTICLE_DISTANCE+EPS)) &&((y<=64*PARTICLE_DISTANCE+EPS))&&((y>=63*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
      if( ((x>=-26*PARTICLE_DISTANCE+EPS))&&((x<=-25*PARTICLE_DISTANCE+EPS)) &&((y<=65*PARTICLE_DISTANCE+EPS))&&((y>=64*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }

    /*Pipa horizontal*/
      if( ((x>=-23*PARTICLE_DISTANCE+EPS))&&((x<=0*PARTICLE_DISTANCE+EPS)) &&((y>=59*PARTICLE_DISTANCE+EPS))&&((y<=63*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
      //atas lengkungan
      if( ((x>=0*PARTICLE_DISTANCE+EPS))&&((x<=2*PARTICLE_DISTANCE+EPS)) &&((y>=58*PARTICLE_DISTANCE+EPS))&&((y<=61*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
      if( ((x>=0*PARTICLE_DISTANCE+EPS))&&((x<=1*PARTICLE_DISTANCE+EPS)) &&((y>=61*PARTICLE_DISTANCE+EPS))&&((y<=62*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
      //bawah lengkungan
      if( ((x>=-4*PARTICLE_DISTANCE+EPS))&&((x<=-2*PARTICLE_DISTANCE+EPS)) &&((y>=58*PARTICLE_DISTANCE+EPS))&&((y<=61*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
      if( ((x>=-3*PARTICLE_DISTANCE+EPS))&&((x<=-2*PARTICLE_DISTANCE+EPS)) &&((y>=57*PARTICLE_DISTANCE+EPS))&&((y<=58*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
    /*Pipa Vertikal drain tank*/
      if( ((x>=-2*PARTICLE_DISTANCE+EPS))&&((x<=2*PARTICLE_DISTANCE+EPS)) &&((y>=24*PARTICLE_DISTANCE+EPS))&&((y<=59*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
    /*Drain Tank*/
      //Kotak atas
      if( ((x>=-24*PARTICLE_DISTANCE+EPS))&&((x<=23*PARTICLE_DISTANCE+EPS)) &&((y>=10*PARTICLE_DISTANCE+EPS))&&((y<=25*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }
    //Silinder bawah
      if( ((x<=(1)*sqrt(y*45*PARTICLE_DISTANCE+EPS)))&&((x>=(-1)*sqrt(y*45*PARTICLE_DISTANCE+EPS))) &&((y>=3*PARTICLE_DISTANCE+EPS))&&((y<=10*PARTICLE_DISTANCE+EPS))){ 
	      flagOfParticleGeneration = OFF;
      }


  /* fluid region */
   /*Active Core*/
      if( ((x>=-40*PARTICLE_DISTANCE+EPS))&&((x<=-16*PARTICLE_DISTANCE+EPS)) &&((y>=101*PARTICLE_DISTANCE+EPS))&&((y<=112*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=FLUID;
	      flagOfParticleGeneration = ON;
      }
      if( ((x>=-39*PARTICLE_DISTANCE+EPS))&&((x<=-17*PARTICLE_DISTANCE+EPS)) &&((y>=100*PARTICLE_DISTANCE+EPS))&&((y<=101*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=FLUID;
	      flagOfParticleGeneration = ON;
      }
      if( ((x>=-38*PARTICLE_DISTANCE+EPS))&&((x<=-18*PARTICLE_DISTANCE+EPS)) &&((y>=99*PARTICLE_DISTANCE+EPS))&&((y<=100*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=FLUID;
	      flagOfParticleGeneration = ON;
      }
      if( ((x>=-37*PARTICLE_DISTANCE+EPS))&&((x<=-19*PARTICLE_DISTANCE+EPS)) &&((y>=98*PARTICLE_DISTANCE+EPS))&&((y<=99*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=FLUID;
	      flagOfParticleGeneration = ON;
      }
      if( ((x>=-36*PARTICLE_DISTANCE+EPS))&&((x<=-20*PARTICLE_DISTANCE+EPS)) &&((y>=97*PARTICLE_DISTANCE+EPS))&&((y<=98*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=FLUID;
	      flagOfParticleGeneration = ON;
      }
      if( ((x>=-35*PARTICLE_DISTANCE+EPS))&&((x<=-21*PARTICLE_DISTANCE+EPS)) &&((y>=96*PARTICLE_DISTANCE+EPS))&&((y<=97*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=FLUID;
	      flagOfParticleGeneration = ON;
      }
      if( ((x>=-34*PARTICLE_DISTANCE+EPS))&&((x<=-22*PARTICLE_DISTANCE+EPS)) &&((y>=95*PARTICLE_DISTANCE+EPS))&&((y<=96*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=FLUID;
	      flagOfParticleGeneration = ON;
      }
      if( ((x>=-30*PARTICLE_DISTANCE+EPS))&&((x<=-26*PARTICLE_DISTANCE+EPS)) &&((y>=88*PARTICLE_DISTANCE+EPS))&&((y<=95*PARTICLE_DISTANCE+EPS))){ 
	      ParticleType[i]=FLUID;
	      flagOfParticleGeneration = ON;
      }

      if( flagOfParticleGeneration == ON){
	      PositionX[i]=x; PositionY[i]=y; PositionZ[i]=z;	
	      i++;
      }
    }
  }
  NumberOfParticles = i;
  for(i=0;i<NumberOfParticles;i++) { VelocityX[i]=0.0; VelocityY[i]=0.0; VelocityZ[i]=0.0; }
}


void initializeParticlePositionAndVelocity_for3dim( void ){
  int iX, iY, iZ;
  int nX, nY, nZ;
  double x, y, z;
  int i = 0;
  int flagOfParticleGeneration;

  nX = (int)(1.0/PARTICLE_DISTANCE)+5;  
  nY = (int)(0.6/PARTICLE_DISTANCE)+5; 
  nZ = (int)(0.3/PARTICLE_DISTANCE)+5;
  for(iX= -4;iX<nX;iX++){
    for(iY= -4;iY<nY;iY++){
      for(iZ= -4;iZ<nZ;iZ++){
	x = PARTICLE_DISTANCE * iX;
	y = PARTICLE_DISTANCE * iY;
	z = PARTICLE_DISTANCE * iZ;
	flagOfParticleGeneration = OFF;

	/* dummy wall region */
	if( (((x>-4.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+4.0*PARTICLE_DISTANCE+EPS))&&( (y>0.0-4.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)))&&( (z>0.0-4.0*PARTICLE_DISTANCE+EPS)&&(z<=0.3+4.0*PARTICLE_DISTANCE+EPS ))){  
	  ParticleType[i]=DUMMY_WALL;
	  flagOfParticleGeneration = ON;
	}

	/* wall region */
	if( (((x>-2.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+2.0*PARTICLE_DISTANCE+EPS))&&( (y>0.0-2.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)))&&( (z>0.0-2.0*PARTICLE_DISTANCE+EPS)&&(z<=0.3+2.0*PARTICLE_DISTANCE+EPS ))){  
	  ParticleType[i]=WALL;
	  flagOfParticleGeneration = ON;
	}

	/* wall region */
	if( (((x>-4.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+4.0*PARTICLE_DISTANCE+EPS))&&( (y>0.6-2.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)))&&( (z>0.0-4.0*PARTICLE_DISTANCE+EPS)&&(z<=0.3+4.0*PARTICLE_DISTANCE+EPS ))){  
	  ParticleType[i]=WALL;
	  flagOfParticleGeneration = ON;
	}

	/* empty region */
	if( (((x>0.0+EPS)&&(x<=1.00+EPS))&&( y>0.0+EPS ))&&( (z>0.0+EPS )&&(z<=0.3+EPS ))){  
	  flagOfParticleGeneration = OFF;
	}

	/* fluid region */
	if( (((x>0.0+EPS)&&(x<=4.0+EPS))&&( (y>0.0+EPS)&&(y<0.5+EPS) ))&&( (z>0.0+EPS )&&(z<=0.3+EPS ))){
	  ParticleType[i]=FLUID;
	  flagOfParticleGeneration = ON;
	}

	if( flagOfParticleGeneration == ON){
	  PositionX[i]=x; 
	  PositionY[i]=y; 
	  PositionZ[i]=z;
	  i++;
	}
      }
    }
  }
  NumberOfParticles = i;
  for(i=0;i<NumberOfParticles;i++) { VelocityX[i]=0.0; VelocityY[i]=0.0; VelocityZ[i]=0.0; }
}


void calculateConstantParameter( void ){

  Re_forParticleNumberDensity  = RADIUS_FOR_NUMBER_DENSITY;  
  Re_forGradient       = RADIUS_FOR_GRADIENT;  
  Re_forLaplacian      = RADIUS_FOR_LAPLACIAN;  
  Re2_forParticleNumberDensity = Re_forParticleNumberDensity*Re_forParticleNumberDensity;
  Re2_forGradient      = Re_forGradient*Re_forGradient;
  Re2_forLaplacian     = Re_forLaplacian*Re_forLaplacian;
  calculateNZeroAndLambda();
  FluidDensity       = FLUID_DENSITY;
  collisionDistance  = COLLISION_DISTANCE; 
  collisionDistance2 = collisionDistance*collisionDistance;
  FileNumber=0;
  Time=0.0;
}


void calculateNZeroAndLambda( void ){
  int iX, iY, iZ;
  int iZ_start, iZ_end;
  double xj, yj, zj, distance, distance2;
  double xi, yi, zi;

  if( DIM == 2 ){
    iZ_start = 0; iZ_end = 1;
  }else{
    iZ_start = -4; iZ_end = 5;
  }

  N0_forParticleNumberDensity = 0.0;
  N0_forGradient      = 0.0;
  N0_forLaplacian     = 0.0;
  Lambda              = 0.0;
  xi = 0.0;  yi = 0.0;  zi = 0.0;

  for(iX= -95;iX<100;iX++){
    for(iY= -95;iY<100;iY++){
      for(iZ= iZ_start;iZ<iZ_end;iZ++){
	if( ((iX==0)&&(iY==0)) && (iZ==0) )continue;
	xj = PARTICLE_DISTANCE * (double)(iX);
	yj = PARTICLE_DISTANCE * (double)(iY);
	zj = PARTICLE_DISTANCE * (double)(iZ);
	distance2 = (xj-xi)*(xj-xi)+(yj-yi)*(yj-yi)+(zj-zi)*(zj-zi);
	distance = sqrt(distance2);
	N0_forParticleNumberDensity += weight(distance, Re_forParticleNumberDensity);
	N0_forGradient      += weight(distance, Re_forGradient);
	N0_forLaplacian     += weight(distance, Re_forLaplacian);
	Lambda              += distance2 * weight(distance, Re_forLaplacian);
      }
    }
  }
  Lambda = Lambda/N0_forLaplacian;
}


double weight( double distance, double re ){
  double weightIJ;

  if( distance >= re ){
    weightIJ = 0.0;
  }else{
    weightIJ = (re/distance) - 1.0;
  }
  return weightIJ;
}


void mainLoopOfSimulation( void ){
  int iTimeStep = 0;

  writeData_inVtuFormat();
  writeData_inProfFormat();

  while(1){
    calculateGravity();
    calculateViscosity();
    moveParticle();
    collision();
    calculatePressure();
    calculatePressureGradient();
    moveParticleUsingPressureGradient();
    iTimeStep++;
    Time += DT;
    if( (iTimeStep % OUTPUT_INTERVAL) == 0 ){
      printf("TimeStepNumber: %4d   Time: %lf(s)   NumberOfParticles: %d\n", iTimeStep, Time, NumberOfParticles);
      writeData_inVtuFormat();
      writeData_inProfFormat();
    }
    if( Time >= FINISH_TIME ){break;}
  }
}


void calculateGravity( void ){
  int i;

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      AccelerationX[i]=GRAVITY_X;
      AccelerationY[i]=GRAVITY_Y;
      AccelerationZ[i]=GRAVITY_Z;
    }else{
      AccelerationX[i]=0.0;
      AccelerationY[i]=0.0;
      AccelerationZ[i]=0.0;
    }
  }
}


void calculateViscosity( void ){
  int i,j;
  double viscosityTermX, viscosityTermY, viscosityTermZ;
  double distance, distance2;
  double w;
  double xij, yij, zij;
  double a;

  a = (KINEMATIC_VISCOSITY)*(2.0*DIM)/(N0_forLaplacian*Lambda);
  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] != FLUID) continue;
    viscosityTermX = 0.0;  viscosityTermY = 0.0;  viscosityTermZ = 0.0;

    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;
      xij = PositionX[j] - PositionX[i];
      yij = PositionY[j] - PositionY[i];
      zij = PositionZ[j] - PositionZ[i];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);
      if(distance<Re_forLaplacian){
	w =  weight(distance, Re_forLaplacian);
	viscosityTermX +=(VelocityX[j]-VelocityX[i])*w;
	viscosityTermY +=(VelocityY[j]-VelocityY[i])*w;
	viscosityTermZ +=(VelocityZ[j]-VelocityZ[i])*w;
      }
    }
    viscosityTermX = viscosityTermX * a;
    viscosityTermY = viscosityTermY * a;
    viscosityTermZ = viscosityTermZ * a;
    AccelerationX[i] += viscosityTermX;
    AccelerationY[i] += viscosityTermY;
    AccelerationZ[i] += viscosityTermZ;
  }
}


void moveParticle( void ){
  int i;

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      VelocityX[i] += AccelerationX[i]*DT; 
      VelocityY[i] += AccelerationY[i]*DT; 
      VelocityZ[i] += AccelerationZ[i]*DT;

      PositionX[i] += VelocityX[i]*DT; 
      PositionY[i] += VelocityY[i]*DT; 
      PositionZ[i] += VelocityZ[i]*DT;
    }
    AccelerationX[i]=0.0;
    AccelerationY[i]=0.0;
    AccelerationZ[i]=0.0;
  }
}


void collision( void ){
  int    i,j;
  double xij, yij, zij;
  double distance,distance2;
  double forceDT; /* forceDT is the impulse of collision between particles */
  double mi, mj;
  double velocity_ix, velocity_iy, velocity_iz;
  double e = COEFFICIENT_OF_RESTITUTION;
  static double VelocityAfterCollisionX[ARRAY_SIZE];
  static double VelocityAfterCollisionY[ARRAY_SIZE];
  static double VelocityAfterCollisionZ[ARRAY_SIZE];  

  for(i=0;i<NumberOfParticles;i++){ 
    VelocityAfterCollisionX[i] = VelocityX[i];
    VelocityAfterCollisionY[i] = VelocityY[i];
    VelocityAfterCollisionZ[i] = VelocityZ[i];     
  }
  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      mi = FluidDensity;
      velocity_ix = VelocityX[i];  
      velocity_iy = VelocityY[i];  
      velocity_iz = VelocityZ[i];
      for(j=0;j<NumberOfParticles;j++){
	if( (j==i) || (ParticleType[j]==GHOST) ) continue;
	xij = PositionX[j] - PositionX[i];
	yij = PositionY[j] - PositionY[i];
	zij = PositionZ[j] - PositionZ[i];
	distance2 = (xij*xij) + (yij*yij) + (zij*zij);
	if(distance2<collisionDistance2){
	  distance = sqrt(distance2);
	  forceDT = (velocity_ix-VelocityX[j])*(xij/distance)
	           +(velocity_iy-VelocityY[j])*(yij/distance)
	           +(velocity_iz-VelocityZ[j])*(zij/distance);
	  if(forceDT > 0.0){
	    mj = FluidDensity;
	    forceDT *= (1.0+e)*mi*mj/(mi+mj);
	    velocity_ix -= (forceDT/mi)*(xij/distance); 
	    velocity_iy -= (forceDT/mi)*(yij/distance); 
	    velocity_iz -= (forceDT/mi)*(zij/distance);
	    /*
	    if(j>i){ fprintf(stderr,"WARNING: Collision occured between %d and %d particles.\n",i,j); }
	    */
	  }
	}
      }
      VelocityAfterCollisionX[i] = velocity_ix; 
      VelocityAfterCollisionY[i] = velocity_iy; 
      VelocityAfterCollisionZ[i] = velocity_iz;
    }
  }
  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      PositionX[i] += (VelocityAfterCollisionX[i]-VelocityX[i])*DT; 
      PositionY[i] += (VelocityAfterCollisionY[i]-VelocityY[i])*DT; 
      PositionZ[i] += (VelocityAfterCollisionZ[i]-VelocityZ[i])*DT;
      VelocityX[i] = VelocityAfterCollisionX[i]; 
      VelocityY[i] = VelocityAfterCollisionY[i]; 
      VelocityZ[i] = VelocityAfterCollisionZ[i];
    }
  }
}


void calculatePressure( void ){
  calculateParticleNumberDensity();
  setBoundaryCondition();
  setSourceTerm();
  setMatrix();
  solveSimultaniousEquationsByGaussEliminationMethod();
  removeNegativePressure();
  setMinimumPressure();
}


void calculateParticleNumberDensity( void ){
  int    i,j;
  double xij, yij, zij;
  double distance, distance2;
  double w;

  for(i=0;i<NumberOfParticles;i++){
    ParticleNumberDensity[i] = 0.0;
    if(ParticleType[i] == GHOST) continue;
    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;
      xij = PositionX[j] - PositionX[i];
      yij = PositionY[j] - PositionY[i];
      zij = PositionZ[j] - PositionZ[i];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);
      w =  weight(distance, Re_forParticleNumberDensity);
      ParticleNumberDensity[i] += w;
    }
  }
}


void setBoundaryCondition( void ){
  int i;
  double n0 = N0_forParticleNumberDensity;
  double beta = THRESHOLD_RATIO_OF_NUMBER_DENSITY;

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i]==GHOST || ParticleType[i]== DUMMY_WALL ){
      BoundaryCondition[i]=GHOST_OR_DUMMY;
    }else if( ParticleNumberDensity[i] < beta * n0 ){
      BoundaryCondition[i]=SURFACE_PARTICLE;
    }else{
      BoundaryCondition[i]=INNER_PARTICLE;
    }
  }
}


void setSourceTerm( void ){
  int i;
  double n0    = N0_forParticleNumberDensity;
  double gamma = RELAXATION_COEFFICIENT_FOR_PRESSURE;

  for(i=0;i<NumberOfParticles;i++){
    SourceTerm[i]=0.0;
    if(ParticleType[i]==GHOST || ParticleType[i]== DUMMY_WALL ) continue;
    if(BoundaryCondition[i]==INNER_PARTICLE){
      SourceTerm[i] = gamma * (1.0/(DT*DT))*((ParticleNumberDensity[i]-n0)/n0);
    }else if(BoundaryCondition[i]==SURFACE_PARTICLE){
      SourceTerm[i]=0.0;
    }
  }
}


void setMatrix( void ){
  double xij, yij, zij;
  double distance, distance2;
  double coefficientIJ;
  double n0 = N0_forLaplacian;
  int    i,j;
  double a;
  int n = NumberOfParticles;

  for(i=0;i<NumberOfParticles;i++){
    for(j=0;j<NumberOfParticles;j++){
      CoefficientMatrix[i*n+j] = 0.0;
    }
  }

  a = 2.0*DIM/(n0*Lambda);
  for(i=0;i<NumberOfParticles;i++){
    if(BoundaryCondition[i] != INNER_PARTICLE) continue;
    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (BoundaryCondition[j]==GHOST_OR_DUMMY) ) continue;
      xij = PositionX[j] - PositionX[i];
      yij = PositionY[j] - PositionY[i];
      zij = PositionZ[j] - PositionZ[i];
      distance2 = (xij*xij)+(yij*yij)+(zij*zij);
      distance  = sqrt(distance2);
      if(distance>=Re_forLaplacian)continue;
      coefficientIJ = a * weight(distance, Re_forLaplacian)/FluidDensity;
      CoefficientMatrix[i*n+j]  = (-1.0)*coefficientIJ;
      CoefficientMatrix[i*n+i] += coefficientIJ;
    }
    CoefficientMatrix[i*n+i] += (COMPRESSIBILITY)/(DT*DT);
  }
  exceptionalProcessingForBoundaryCondition();
}


void exceptionalProcessingForBoundaryCondition( void ){
  /* If tere is no Dirichlet boundary condition on the fluid, 
     increase the diagonal terms of the matrix for an exception. This allows us to solve the matrix without Dirichlet boundary conditions. */
  checkBoundaryCondition();
  increaseDiagonalTerm();
}


void checkBoundaryCondition( void ){
  int i,j,count;
  double xij, yij, zij, distance2;

  for(i=0;i<NumberOfParticles;i++){
    if (BoundaryCondition[i]==GHOST_OR_DUMMY){
      FlagForCheckingBoundaryCondition[i]=GHOST_OR_DUMMY;
    }else if (BoundaryCondition[i]==SURFACE_PARTICLE){
      FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_CONNECTED;
    }else{
      FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_NOT_CONNECTED;
    }
  }

  do {
    count=0;
    for(i=0;i<NumberOfParticles;i++){
      if(FlagForCheckingBoundaryCondition[i]==DIRICHLET_BOUNDARY_IS_CONNECTED){
	for(j=0;j<NumberOfParticles;j++){
	  if( j==i ) continue;
	  if((ParticleType[j]==GHOST) || (ParticleType[j]== DUMMY_WALL)) continue;
	  if(FlagForCheckingBoundaryCondition[j]==DIRICHLET_BOUNDARY_IS_NOT_CONNECTED){
	    xij = PositionX[j] - PositionX[i];
	    yij = PositionY[j] - PositionY[i];
	    zij = PositionZ[j] - PositionZ[i];
	    distance2 = (xij*xij)+(yij*yij)+(zij*zij);
	    if(distance2>=Re2_forLaplacian)continue;
	    FlagForCheckingBoundaryCondition[j]=DIRICHLET_BOUNDARY_IS_CONNECTED;
	  }
	}
	FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_CHECKED;
	count++;
      }
    }
  } while (count!=0); /* This procedure is repeated until the all fluid or wall particles (which have Dirhchlet boundary condition in the particle group) are in the state of "DIRICHLET_BOUNDARY_IS_CHECKED".*/

  for(i=0;i<NumberOfParticles;i++){
    if(FlagForCheckingBoundaryCondition[i]==DIRICHLET_BOUNDARY_IS_NOT_CONNECTED){
      fprintf(stderr,"WARNING: There is no dirichlet boundary condition for %d-th particle.\n",i );
    }
  }
}


void increaseDiagonalTerm( void ){
  int i;
  int n = NumberOfParticles;

  for(i=0;i<n;i++) {
    if(FlagForCheckingBoundaryCondition[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED ){
      CoefficientMatrix[i*n+i] = 2.0 * CoefficientMatrix[i*n+i];
    }
  }
}


void solveSimultaniousEquationsByGaussEliminationMethod( void ){
  int    i,j,k;
  double c;
  double sumOfTerms;
  int    n = NumberOfParticles;

  for(i=0; i<n; i++){ 
    Pressure[i] = 0.0; 
  }
  for(i=0; i<n-1; i++){
    if ( BoundaryCondition[i] != INNER_PARTICLE ) continue;
    for(j=i+1; j<n; j++){
      if(BoundaryCondition[j]==GHOST_OR_DUMMY) continue;
      c = CoefficientMatrix[j*n+i]/CoefficientMatrix[i*n+i];
      for(k=i+1; k<n; k++){
	CoefficientMatrix[j*n+k] -= c * CoefficientMatrix[i*n+k];
      }
      SourceTerm[j] -= c*SourceTerm[i];
    }
  }
  for( i=n-1; i>=0; i--){
    if ( BoundaryCondition[i] != INNER_PARTICLE ) continue;
    sumOfTerms = 0.0;
    for( j=i+1; j<n; j++ ){
      if(BoundaryCondition[j]==GHOST_OR_DUMMY) continue;
      sumOfTerms += CoefficientMatrix[i*n+j] * Pressure[j];
    }
    Pressure[i] = (SourceTerm[i] - sumOfTerms)/CoefficientMatrix[i*n+i];
  }
}


void removeNegativePressure( void ){
  int i;

  for(i=0;i<NumberOfParticles;i++) {
    if(Pressure[i]<0.0)Pressure[i]=0.0;
  }
}


void setMinimumPressure( void ){
  double xij, yij, zij, distance2;
  int i,j;

  for(i=0;i<NumberOfParticles;i++) {
    if(ParticleType[i]==GHOST || ParticleType[i]==DUMMY_WALL)continue;
    MinimumPressure[i]=Pressure[i];
    for(j=0;j<NumberOfParticles;j++) {
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;
      if(ParticleType[j]==DUMMY_WALL) continue;
      xij = PositionX[j] - PositionX[i];
      yij = PositionY[j] - PositionY[i];
      zij = PositionZ[j] - PositionZ[i];
      distance2 = (xij*xij)+(yij*yij)+(zij*zij);
      if(distance2>=Re2_forGradient)continue;
      if( MinimumPressure[i] > Pressure[j] ){
	MinimumPressure[i] = Pressure[j];
      }
    }
  }
}


void calculatePressureGradient( void ){
  int    i,j;
  double gradient_x, gradient_y, gradient_z;
  double xij, yij, zij;
  double distance, distance2;
  double w,pij;
  double a;

  a =DIM/N0_forGradient;
  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] != FLUID) continue;
    gradient_x = 0.0;  gradient_y = 0.0;  gradient_z = 0.0;
    for(j=0;j<NumberOfParticles;j++){
      if( j==i ) continue;
      if( ParticleType[j]==GHOST ) continue;
      if( ParticleType[j]==DUMMY_WALL ) continue;
      xij = PositionX[j] - PositionX[i];
      yij = PositionY[j] - PositionY[i];
      zij = PositionZ[j] - PositionZ[i];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);
      if(distance<Re_forGradient){
	w =  weight(distance, Re_forGradient);
	pij = (Pressure[j] - MinimumPressure[i])/distance2;
	gradient_x += xij*pij*w;
	gradient_y += yij*pij*w;
	gradient_z += zij*pij*w;
      }
    }
    gradient_x *= a;
    gradient_y *= a;
    gradient_z *= a;
    AccelerationX[i]= (-1.0)*gradient_x/FluidDensity;
    AccelerationY[i]= (-1.0)*gradient_y/FluidDensity;
    AccelerationZ[i]= (-1.0)*gradient_z/FluidDensity;
  }
}


void moveParticleUsingPressureGradient( void ){
  int i;

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      VelocityX[i] +=AccelerationX[i]*DT;
      VelocityY[i] +=AccelerationY[i]*DT;
      VelocityZ[i] +=AccelerationZ[i]*DT;

      PositionX[i] +=AccelerationX[i]*DT*DT;
      PositionY[i] +=AccelerationY[i]*DT*DT;
      PositionZ[i] +=AccelerationZ[i]*DT*DT;
    }
    AccelerationX[i]=0.0;
    AccelerationY[i]=0.0;
    AccelerationZ[i]=0.0;
  }
}


void writeData_inProfFormat( void ){
  int i;
  FILE *fp;
  char fileName[256];

  sprintf(fileName, "output_%04d.prof",FileNumber);
  fp = fopen(fileName, "w");
  fprintf(fp,"%lf\n",Time);
  fprintf(fp,"%d\n",NumberOfParticles);
  for(i=0;i<NumberOfParticles;i++) {
    fprintf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n"
	    ,ParticleType[i], PositionX[i], PositionY[i], PositionZ[i]
	    ,VelocityX[i], VelocityY[i], VelocityZ[i], Pressure[i], ParticleNumberDensity[i]);
  }
  fclose(fp);
  FileNumber++;
}


void writeData_inVtuFormat( void ){
  int i;
  double absoluteValueOfVelocity;
  FILE *fp;
  char fileName[1024];

  sprintf(fileName, "particle_%04d.vtu", FileNumber);
  fp=fopen(fileName,"w");
  fprintf(fp,"<?xml version='1.0' encoding='UTF-8'?>\n");
  fprintf(fp,"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  fprintf(fp,"<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n",NumberOfParticles,NumberOfParticles);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%lf %lf %lf\n",PositionX[i],PositionY[i],PositionZ[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",ParticleType[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    absoluteValueOfVelocity=
      sqrt( VelocityX[i]*VelocityX[i] + VelocityY[i]*VelocityY[i] + VelocityZ[i]*VelocityZ[i] );
    fprintf(fp,"%f\n",(float)absoluteValueOfVelocity);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='Pressure' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%f\n",(float)Pressure[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='ParticleNumberDensity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%f\n",(float)ParticleNumberDensity[i]); 
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<Cells>\n");
  fprintf(fp,"<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",i);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type='Int32' Name='offsets' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",i+1);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type='UInt8' Name='types' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"1\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</UnstructuredGrid>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
} 
