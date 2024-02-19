#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h> 


#define PI 3.141592653589793238462643383279
#define MU 4*PI*math.pow(10,-7)


const float LX = 100; //X size direction
const float LY = 100; //Y size direction
const float N = LX * LY; //number of dipoles
const float J = 1.0; //interaction constant 
const float m = 1.0, //magnetization module
const float H = -1.0; //external field module


//FUNCTION TO STABLISH INITIAL STATE OF THE SYSTEM
double Initial_State(double spin_matrix[LX][LY]){ //PER ARA ESCULLO QUE L'ESTAT INICIAL SIGUI AMB TOTS ELS SPINS ANTI-PARAL·LELS A H
	for(int i=0;i<LX;i++){
		for(int j=0;j<LY;j++){
			spin_matrix[i][j] = -1.0;
		}
	}
	return spin_matrix;
}


double Choose_Random_Position(double spin_matrix[LX][LY]){
	srand(time(NULL));
	int rand_pos[2];
	
	rand_pos[0] = rand() % LX;
	rand_pos[1] = rand() % LY;
	
	return rand_pos;
}
