//MODEL D'ISING ------> 1D QUÀNTIC
//Versió 0

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>

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

double Energia(double spin_matrix[L], int position){
	double e1 = 0;
	double e2 = 0;
	int c = 0;
	while(c < L){
		if(c != (L-1)){
			e1 += (-J)*(spin_matrix[c] * spin_matrix[c+1]) + (-M)*MU*(H/2)*(spin_matrix[c] + spin_matrix[c+1]);	
		}else{
			e1 += (-J)*(spin_matrix[c] * spin_matrix[0]) + (-M)*MU*(H/2)*(spin_matrix[c] + spin_matrix[0]);
		}
		c += 1;
	}
	spin_matrix[position] = -1.0;
	c = 0;
	while(c < L){
		if(c != (L-1)){
			e2 += (-J)*(spin_matrix[c] * spin_matrix[c+1]) + (-M)*MU*(H/2)*(spin_matrix[c] + spin_matrix[c+1]);	
		}else{
			e2 += (-J)*(spin_matrix[c] * spin_matrix[0]) + (-M)*MU*(H/2)*(spin_matrix[c] + spin_matrix[0]);
		}
		c += 1;
	}
	printf("\n%f\n", fabs(e1 - e2));
	return fabs(e1 - e2);
}



int main(){
	double spins[L];
	//DEFINIR ESTAT INICIAL
	for(int i=0;i<L;i++){
		spins[i] = 1.0;
		printf("%f\t", spins[i]);
	}
	printf("\n\n");
	
	srand(time(NULL));
	for(int t = 0; t < 50; t++){
		//s'escull un spin de posició random
		int s = Choose_Random_Position(spins);
		printf("\nRandom Position: %d\n", s);
		double dif_e = Energia(spins, s); // es calcula la diferencia entre l'energia del sistema tal i com està respecte a canviar de sentit el spin escollit
		if(dif_e <= 0){
			printf("dif_e <= 0\n");
			printf("Random spin before: %f\n", spins[s]);
			spins[s] = (-1.0) * spins[s];
			printf("Random spin after: %f\n", spins[s]);
		}else{
			printf("dif_e > 0\n");
			printf("Random spin before: %f\n", spins[s]);
			double p;
			p = exp(-BETA * dif_e);
			printf("Transition probability: %f\n", p); //REVISAR!!!!!!!!!!!!!!!!!! PQ SEMPRE RETORNA 0.0000?
			double r = rand();
			printf("r: %f\n", r);
			if(r <= p){
				spins[s] = (-1.0) * spins[s];
			}
			printf("Random spin after: %f\n", spins[s]);
		}
	}
	
	for(int i = 0; i < L; i++){
		printf("%f\t", spins[i]);
	}
	
}
