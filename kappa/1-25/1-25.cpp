#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#define PI 3.141592653589793238462643383279

//Geometry of computation window and time parameters, in units of Dist*lex

#define N 81  //Size in X direction
#define M 81 //Size in Y direction

int n=N-1; //N-1
int m=M-1; //M-1
int l=(N-1)/2; //n/2
int L=(M-1)/2; //m/2
int R=(N-3)/2;  //Radius of the disc
double Dist=0.5; //Distance between divisions in exchange lengths
double invDist=2.0;//(1/0.5);
int Time=40000; //Limit time of RK4
int Time2=100000000;
double h=0.001; // Time step


//Physical parameters

double xi=1.68;//1.68;//1.290238; //DM parameter (xi)
double kappa=1.25;   //Anisotropy parameter (kappa)
double alpha=0.1; //Gilbert damping constant (alpha)
double alpha2=0;;// 1/(1+alpha^2)
double hex=0.0; //external applied field

//CUrrent parameters
double sigx=-1.0;
double sigy=0.0;
double sigz=0;

double dj=0.1;

//Temperature
double temp=0.0033;
double fact=sqrt(temp/(Dist*Dist*Dist*h));

//Definition of energys: E=Total E1=Exchange E2=DM E3=Anisotropy E4=Zeeman

double E=0;
double E1=0;
double E2=0;
double E3=0;
double E4=0;
double E0=0;
double E01=0;
double E02=0;
double E03=0;
double E04=0;
double Magnx, Magny, Magnz;
double mx[N][M], my[N][M], mz[N][M];
double mx0[N][M], my0[N][M], mz0[N][M];
double hx[N][M], hy[N][M], hz[N][M];
double phi[N][M];
double xixi[N][M];
double mdotsig, mdoth;
//Prints


double Abs(double KKK){
	if (KKK>=0){
		return KKK;
	}
	else{
		return -KKK;
	}
}


double Angle(){
	
	int i, j;
	double x, y;
		
	for(i=0; i<=n; i++){
		for(j=0; j<=m; j++){
		
			x=i-l;
			y=j-L;
			
			if(x==0){			
				phi[i][j]=2.0*atan(1.0)*((Abs(y))/(y));
			}
			else{
				phi[i][j]=atan(y/x);
			}
				
			if((x<0) & (y!=0)){	
				phi[i][j]=(4.0*atan(1.0)-Abs(phi[i][j]))*(-phi[i][j]/Abs(phi[i][j]));
			}
			else if (x<0){
				phi[i][j]=-4.0*atan(1.0);
			}
		
		}
	}
	
	phi[l][L]=0;
	return 0;
}


double AWGN_generator()
{/* Generates additive white Gaussian Noise samples with zero mean and a standard deviation of 1. */
  double temp1;
  double temp2;
  double result;
  int p;

  p = 1;
  while( p > 0 )
  {
	temp2 = ( rand() / ( (double)RAND_MAX ) ); 
    if ( temp2 == 0 )    {// temp2 is >= (RAND_MAX / 2)
	    p = 1;
    }// end if
    else
    {// temp2 is < (RAND_MAX / 2)
       p = -1;
    }// end else
  }// end while()
  temp1 = cos( ( 2.0 * (double)PI ) * rand() / ( (double)RAND_MAX ) );
  result = sqrt( -2.0 * log( temp2 ) ) * temp1;

  return result;	// return the generated random sample to the caller

}// end AWGN_generator()


void PrintM(double Matrix[N][M]){
	int i , j;
	for(j=M; j>=0; j--){
		for(i=0; i<N; i++){
			printf("%.5lf ", Matrix[i][j] );		
		}		
		printf("\n");	
	}	
	printf("\n");

}

//Derivatives

//First X

double dx(double v[N][M], int i, int j){

	if(i==n){
	return (v[0][j]-v[n-1][j])*0.5*invDist;
	}

	if(i==0){
	return (v[1][j]-v[n][j])*0.5*invDist;
	}
	else{
	return (v[i+1][j]-v[i-1][j])*0.5*invDist;
	
	}
}

//First Y

double dy(double v[N][M], int i, int j){

	if(j==m){
	return (v[i][0]-v[i][j-1])*0.5*invDist;
	}

	if(j==0){
	return (v[i][j+1]-v[i][m])*0.5*invDist;
	}
	else{
	return (v[i][j+1]-v[i][j-1])*0.5*invDist;
	}
}

//Second X


double dx2(double v[N][M], int i, int j){

	if(i==n){
	return (v[0][j]-2.0*v[i][j]+v[i-1][j])*invDist*invDist;
	}

	if(i==0){
	return (v[i+1][j]-2.0*v[i][j]+v[n][j])*invDist*invDist;
	}
	else{
	return (v[i+1][j]-2.0*v[i][j]+v[i-1][j])*invDist*invDist;
	}
}

//Second Y

double dy2(double v[N][M], int i, int j){

	if(j==m){
	return (v[i][0]-2.0*v[i][j]+v[i][j-1])*invDist*invDist;
	}

	if(j==0){
	return (v[i][j+1]-2.0*v[i][j]+v[i][m])*invDist*invDist;
	}
	else{
	return (v[i][j+1]-2.0*v[i][j]+v[i][j-1])*invDist*invDist;
	}
}

//Functions that compute the different terms of the effective field

double Exchange(double v[N][M], int i, int j){
	return (dx2(v, i , j)+dy2(v, i, j));
}

double Moriyax (int i, int j){
	return (-xixi[i][j]*dx(mz, i, j));
}

double Moriyay (int i, int j){
	return (-xixi[i][j]*dy(mz, i, j));
}

double Moriyaz (int i, int j){
	return xixi[i][j]*(dx(mx, i, j)+dy(my, i, j));
}

double Anisotropy (int i, int j){
	return kappa*mz[i][j];
}

double Zeeman (int i, int j){
	return hex;
}


double Thermal(int i,int j) {
	return fact*AWGN_generator();
}


//LLG Function

double LLGi(double mi, double mj, double mk, double hi, double hj, double hk, double sigi, double sigj, double sigk){

	return (-(mj*hk-mk*hj)-alpha*(mi*mdoth-hi)-dj*(mi*mdotsig-sigi)+dj*alpha*(mj*sigk-mk*sigj))/(1+alpha*alpha);
//return -(Pi+alpha*(vj*Pk-vk*Pj)+(1+alpha*beta)*ui)/(1+alpha*alpha);
}


//Effective field

void heff(int i, int j){

	hx[i][j]=Exchange(mx, i , j)+Moriyax(i, j)+Thermal(i,j);
	hy[i][j]=Exchange(my, i , j)+Moriyay(i, j)+Thermal(i,j);
	hz[i][j]=Exchange(mz, i , j)+Moriyaz(i, j)+Anisotropy(i, j)+Zeeman(i, j)+Thermal(i,j);

}

//Normalize Function

void Normalize( double vx, double vy, double vz){

	double Mod;
	
	Mod=sqrt(vx*vx+vy*vy+vz*vz);
	
	vx=vx/Mod;
	vy=vy/Mod;
	vz=vz/Mod;

}

void TimeEvolution(){

	double k1x[N][M], k1y[N][M], k1z[N][M];
	double k2x[N][M], k2y[N][M], k2z[N][M];
	double k3x[N][M], k3y[N][M], k3z[N][M];
	double k4x[N][M], k4y[N][M], k4z[N][M];
	double mxx, myy, mzz;
	double Px, Py, Pz, ux, uy, uz, Mod; // components of the vectorial product m x h
	char buffer [100];
	FILE * sFile;
	
	int i, tout, j, t, k=0;
	
	tout=Time/100;
	
	for(t=0; t<=Time; t++){
	
		for(i=0; i<=n; i++){
			for(j=0; j<=m; j++){
			
			
				//Compute the effective field
				heff(i, j);
				
				mxx=mx[i][j];
				myy=my[i][j];
				mzz=mz[i][j];
				
				mdoth=mxx*hx[i][j]+myy*hy[i][j]+mzz*hz[i][j];
				mdotsig=mxx*sigx+myy*sigy+mzz*sigz;
				//First constant
				
				
				k1x[i][j]=LLGi(mxx, myy, mzz, hx[i][j], hy[i][j], hz[i][j], sigx, sigy, sigz);
				k1y[i][j]=LLGi(myy, mzz, mxx, hy[i][j], hz[i][j], hx[i][j], sigy, sigz, sigx);
				k1z[i][j]=LLGi(mzz, mxx, myy, hz[i][j], hx[i][j], hy[i][j], sigz, sigx, sigy);
				
				mxx=mx[i][j]+0.5*h*k1x[i][j];
				myy=my[i][j]+0.5*h*k1y[i][j];
				mzz=mz[i][j]+0.5*h*k1z[i][j];
				
				
				//Second Constant
				mdoth=mxx*hx[i][j]+myy*hy[i][j]+mzz*hz[i][j];
				mdotsig=mxx*sigx+myy*sigy+mzz*sigz;
				
				k2x[i][j]=LLGi(mxx, myy, mzz, hx[i][j], hy[i][j], hz[i][j], sigx, sigy, sigz);
				k2y[i][j]=LLGi(myy, mzz, mxx, hy[i][j], hz[i][j], hx[i][j], sigy, sigz, sigx);
				k2z[i][j]=LLGi(mzz, mxx, myy, hz[i][j], hx[i][j], hy[i][j], sigz, sigx, sigy);
				
				mxx=mx[i][j]+0.5*h*k2x[i][j];
				myy=my[i][j]+0.5*h*k2y[i][j];
				mzz=mz[i][j]+0.5*h*k2z[i][j];
				//Third Constant
				
				mdoth=mxx*hx[i][j]+myy*hy[i][j]+mzz*hz[i][j];
				mdotsig=mxx*sigx+myy*sigy+mzz*sigz;
				
				k3x[i][j]=LLGi(mxx, myy, mzz, hx[i][j], hy[i][j], hz[i][j], sigx, sigy, sigz);
				k3y[i][j]=LLGi(myy, mzz, mxx, hy[i][j], hz[i][j], hx[i][j], sigy, sigz, sigx);
				k3z[i][j]=LLGi(mzz, mxx, myy, hz[i][j], hx[i][j], hy[i][j], sigz, sigx, sigy);
				
				mxx=mx[i][j]+h*k3x[i][j];
				myy=my[i][j]+h*k3y[i][j];
				mzz=mz[i][j]+h*k3z[i][j];
				
				//Fourth Constant
				
				mdoth=mxx*hx[i][j]+myy*hy[i][j]+mzz*hz[i][j];
				mdotsig=mxx*sigx+myy*sigy+mzz*sigz;
				
				k4x[i][j]=LLGi(mxx, myy, mzz, hx[i][j], hy[i][j], hz[i][j], sigx, sigy, sigz);
				k4y[i][j]=LLGi(myy, mzz, mxx, hy[i][j], hz[i][j], hx[i][j], sigy, sigz, sigx);
				k4z[i][j]=LLGi(mzz, mxx, myy, hz[i][j], hx[i][j], hy[i][j], sigz, sigx, sigy);
	
	}}
	
	//PrintM(k1x);
	//PrintM(k1y);
	//PrintM(k1z);
	
	//printf("%d\n", t);
	
	//Update m and normalize
	
		for(i=0; i<=n; i++){
			for(j=0; j<=n; j++){
			
			
				mx[i][j]+=(h/6.0)*(k1x[i][j]+2.0*k2x[i][j]+2.0*k3x[i][j]+k4x[i][j]);
				my[i][j]+=(h/6.0)*(k1y[i][j]+2.0*k2y[i][j]+2.0*k3y[i][j]+k4y[i][j]);
				mz[i][j]+=(h/6.0)*(k1z[i][j]+2.0*k2z[i][j]+2.0*k3z[i][j]+k4z[i][j]);
				
				Mod=mx[i][j]*mx[i][j]+my[i][j]*my[i][j]+mz[i][j]*mz[i][j];
				mx[i][j]=mx[i][j]/sqrt(Mod);
				my[i][j]=my[i][j]/sqrt(Mod);
				mz[i][j]=mz[i][j]/sqrt(Mod);
			
		
		}}
	
	
	
		if((t%tout)==0){
	
	
			sprintf(buffer,"C:/Users/anapascual/Desktop/s_skyrmions/kappa/1-25/kappa_%d.txt",k);
			k=k+1;
			sFile = fopen (buffer , "w" );
			printf("t= %d\n",k);
		
			for (j=0; j<=m; j++){
			for (i=0; i<=n; i++){
			
		
			fprintf(sFile,"%lf ",mz[i][j]);
			}
			fprintf(sFile,"\n");
			}
			fclose(sFile);
		}
	
	
	
	}


}

void Initialize(){

	int i, j;
	double x, y;
	double R=3.0;
	double Mod;
	
	//Condicions inicials
	
	
	
	for (i=0; i<=n; i++){
		for (j=0; j<=m; j++){
			x=Dist*(i-l);
			y=Dist*(j-L);
			
			mz[i][j]=(R*R-x*x-y*y)/(R*R+x*x+y*y);
			mx[i][j]=2.0*R*x/((R*R+x*x+y*y));
			my[i][j]=2.0*R*y/((R*R+x*x+y*y));
			
			hx[i][j]=0;
			hy[i][j]=0;
			hz[i][j]=0;
			xixi[i][j]=xi;
	}
	}
	
	
	
	for (i=0; i<=n; i++){
	
		mx[i][0]=0;
		mx[i][m]=0;
		
		my[i][0]=0;
		my[i][m]=0;
		
		mz[i][0]=-1;
		mz[i][m]=-1;
	}
	
	for (j=0; j<=m; j++){
	
		mx[n][j]=0;
		mx[0][j]=0;
		
		my[n][j]=0;
		my[0][j]=0;
		
		mz[n][j]=-1;
		mz[0][j]=-1;
	}
	
	
	mx[l][L]=0;
	my[l][L]=0;
	
	
	
	for (i=0; i<=n; i++){
		for (j=0; j<=m; j++){
			Mod=mx[i][j]*mx[i][j]+my[i][j]*my[i][j]+mz[i][j]*mz[i][j];
			mx[i][j]=mx[i][j]/sqrt(Mod);
			my[i][j]=my[i][j]/sqrt(Mod);
			mz[i][j]=mz[i][j]/sqrt(Mod);
	}
	}



}



void Check(){
	FILE* output;
	double res;
	output=fopen("check.dat","w+");
	for (int i=1;i<=1000;i++){
		res=AWGN_generator();
		fprintf(output,"%i %e\n",i,res);
	}
	fclose(output);
}


int main(){
	int i, j, k, x1, y1, x2, y2;
	double ux, uy, uz;
	
	//Check();
	//return 1;
	Initialize();
	TimeEvolution();
	
	
	return 0;

}

