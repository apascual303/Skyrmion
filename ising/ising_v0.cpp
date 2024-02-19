//MODEL D'ISING ------> 1D CLASSIC
//Versió 0

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define PI 3.141592653589793238462643383279
#define E 2.718281828459045235360
#define kb 1.380649 * pow(10, -23)
#define L 100

const double J = 1.0; //interaction constant
const double H = 1.0;
const double M = 1.0;
const double T = 293.0; //temperatura ambient
const double BETA =  1 / (kb * T);
const double MU = 4 * PI * pow(10, -7);


int Choose_Random_Position(double spin_matrix[L]){
	int rand_pos = rand() % L;	
	return rand_pos;
}

double Generator_New_Inclination(){
	double cos_min = -1.0;
	double cos_max = 1.0;
	return cos((((rand() / RAND_MAX()) * (cos_max - cos_min) + cos_min) * PI) / 180.0);
}


double Total_Energy(double energy_matrix[LX][LY], double total_energy){
	for(int i = 0; i < LX;i++){
		for(int j = 0; j < LY;j++){
			total_energy += energy_matrix[i][j];
		}
	}
	return total_energy;
}

double Dif_Energy(double spin_matrix[LX][LY], double energy_matrix[LX][LY], int x, int y){
	energy_matrix[x][y] = -MU*H*m*(1+spin_matrix[x][y]) + 4*J*pow(m,2)*(abs(spin_matrix[x][y] - spin_matrix[x-1][y]) + abs(spin_matrix[x][y] - spin_matrix[x+1][y]) + abs(spin_matrix[x][y] - spin_matrix[x][y-1]) + abs(spin_matrix[x][y] - spin_matrix[x][y+1]));
	return energy_matrix[x][y];
}

double Cost_Function(float temperature, int x, int y){
	return exp((-Dif_Energy(m[x][y], e[x][y], x, y)) / temperature);
}


double Temperature(float tmin, float tmax){
	double t = [];
	for((int i=0; i<=N; i++){
		t[i] =	tmin + i*0.01; //increasing 0.01 K each iteration
	}
	return t;
}




int main(){
	float t_min = 0.01, t_max = 5; //def minimum & maximum temperature in K
	double T = Temperature(t_min, t_max);
	
	m = [LX][LY]; //matrix with all magnetization inclination (theta) for each one of the N spins
	m = Initial_State(m);
	e = [LX][LY]; //matrix with the differences of energy of each 4-close neighbours spins for a certain spin in position [i][j]
	
	for(int z = 0; z < 100; z++){
		int pos[2];
		pos = Choose_Random_Position(m[LX][LY]);
		double spin = m[pos[0]][pos[1]];
		double change = Cost_Function(T[z], pos[0], pos[1]);
		if (change < 0){
			
		}			
	}
	
	
	
	
	
	
}



