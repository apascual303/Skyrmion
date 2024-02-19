//MODEL D'ISING ------> 2D QUÀNTIC
//Versió 1

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

//FUNCIO PER ESTABLIR L'ESTAT INICIAL
/*double Initial_State(double spin_matrix[L]){ //PER ARA ESCULLO QUE L'ESTAT INICIAL SIGUI AMB TOTS ELS SPINS ANTI-PARAL·LELS A H
	for(int i=0;i<L;i++){
		spin_matrix[i] = 1.0;
		printf("%f\t", spin_matrix[i]);
	}
	return spin_matrix;
}*/

int Choose_Random_Position(double spin_matrix[L]){
	int rand_pos = rand() % L;	
	return rand_pos;
}

double Energia(double spin_matrix[L], int position1, int position2){
	double e1 = 0;
	double e2 = 0;
	
	for(int i = 0; i < L; i++){
		if(i != (L-1)){
			for(int j = 0; j < L; j++){
				e1 += (-J)*(spin_matrix[i][j] * spin_matrix[c+1]) + (-M)*MU*(H/2)*(spin_matrix[c] + spin_matrix[c+1]);
			}
		}else{
			
		}
		
	}
	
	
	
	
	while(c < L){
		if(c != (L-1)){
			for(int i = 0; i < L; )
			e1 += (-J)*(spin_matrix[c] * spin_matrix[c+1]) + (-M)*MU*(H/2)*(spin_matrix[c] * spin_matrix[c+1]);	
		}else{
			e1 += (-J)*(spin_matrix[0] * spin_matrix[1]) + (-M)*MU*(H/2)*(spin_matrix[0] * spin_matrix[1]);
		}
		c += 1;
	}
	spin_matrix[position] = -1.0;
	c = 0;
	while(c < L){
		if(c != (L-1)){
			e2 += (-J)*(spin_matrix[c] * spin_matrix[c+1]) + (-M)*MU*(H/2)*(spin_matrix[c] * spin_matrix[c+1]);	
		}else{
			e2 += (-J)*(spin_matrix[0] * spin_matrix[1]) + (-M)*MU*(H/2)*(spin_matrix[0] * spin_matrix[1]);
		}
		c += 1;
	}
	printf("\n%f\n", fabs(e1 - e2));
	return fabs(e1 - e2);
}



int main(){
	double spins[L][L];
	//DEFINIR ESTAT INICIAL
	for(int j=0;j<L;j++){
		for(int i = 0; i < L; i++){
			spins[i][j] = 1.0;
			printf("%f\t", spins[i][j]);
		}
		printf("\n");
	}
	printf("\n\n");
	
	srand(time(NULL));
	for(int t = 0; t < 100; t++){
		//s'escull un spin de posició random
		int s1 = Choose_Random_Position(spins);
		int s2 = Choose_Random_Position(spins);
		printf("\nRandom Position: %d, %d\n", s1, s2);
		double dif_e = Energia(spins, s1, s2); // es calcula la diferencia entre l'energia del sistema tal i com està respecte a canviar de sentit el spin escollit
		if(dif_e <= 0){
			printf("dif_e <= 0\n");
			printf("Random spin before: %f\n", spins[s1][s2]);
			spins[s] = -1.0;
			printf("Random spin after: %f\n", spins[s1][s2]);
		}else{
			printf("dif_e > 0\n");
			printf("Random spin before: %f\n", spins[s1][s2]);
			double p;
			p = pow(E, (-BETA * dif_e));
			printf("Transition probability: %f\n", p); //REVISAR!!!!!!!!!!!!!!!!!! PQ SEMPRE RETORNA 0.0000?
			double r = rand();
			printf("r: %f\n", r);
			if(r <= p){
				spins[s1][s2] = -1.0;
			}
			printf("Random spin after: %f\n", spins[s1][s2]);
		}
	}
	
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			printf("%f\t", spins[i][j]);			
		}
		printf("\n");
	}
	
}
