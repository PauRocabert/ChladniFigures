#include <stdio.h>
#include <math.h>

#define Nt 20001
#define Nx 361

const double PI = 3.14159265359;
double da = 1./(Nx-1);
double dt = 2./(Nt-1);

double b = 0.02;
double w =  11;
//13.823
//4./3.0*sqrt(21)*PI;
double F0 = 5;
double u[3][Nx][Nx];


int ic = (Nx-1)/2;
int jc = (Nx-1)/6;

//FUNCIONS ####################################################################################################################

void cc(double u[3][Nx][Nx]){ //CONDICIONS DE CONTORN

    for (int i =2; i<ic; i++){
        //FRONTERA 1
        u[2][i][i] = 0.75*u[2][i+2][i-1]+0.25*u[2][i+1][i-1];
        //FRONTERA 2
        u[2][Nx-i][i] = 0.75*u[2][Nx-i-2][i-1]+0.25*u[2][Nx-i-1][i-1];
    }    

    //FRONTERA 3
    for (int i =2; i<Nx-2;i++){
        u[2][i][0] = u[2][i][1];  
    } 

    //VERTEX 3
    u[2][ic][ic] = u[2][ic][ic -1];

    //VERTEX 1
    u[2][0][0] = u[2][2][1];

    u[2][1][1] = u[2][2][1];
    u[2][1][0] = u[2][2][1];
    
    //VERTEX 2
    u[2][Nx-1][0] = u[2][Nx-3][1];

    u[2][Nx-2][1] = u[2][Nx-3][1];
    u[2][Nx-2][0] = u[2][Nx-3][1];

}

void actualitzar(double u[3][Nx][Nx]){ //EMPREM AQUET CODI PER EMPRAR MENYS MEMORIA

    for (int i =0; i<Nx; i++){
        for (int j =0; j<=ic; j++){
            u[0][i][j] = u[1][i][j];
            u[1][i][j] = u[2][i][j];
        }
    }

}

//PROGRAMA ####################################################################################################################
int main(){

    FILE *output = fopen("trig.txt","w");

    double gamma = pow(dt,2)/(pow(da,2)); 

    //CONDICIONS INICIALS
    for (int i = 0; i<Nx; i++){
        for (int j =0; j<=Nx; j++){
            u[0][i][j] = 0.;
            u[1][i][j] = 0.;
            u[2][i][j] = 0.;
        }
    }

    for (int n =2; n<Nt; n++){ //Euler explicit

        for (int i = 0; i<Nx; i++){

            if (i<ic){
                for (int j = 0; j<=i; j++){
                    u[2][i][j] = 1./(1+b*dt)*(gamma*(u[1][i+1][j] + u[1][i-1][j] +1./3.*(u[1][i][j+1] + u[1][i][j-1]) -8.0/3.*u[1][i][j])+2*u[1][i][j] - u[0][i][j]+b*u[0][i][j]*dt);
                }
            }

            else{
                for (int j = 0; j<=Nx-i; j++){
                    if ((i == ic) && (j == jc)){
                        u[2][ic][jc] = 1./(1+b*dt)*(F0*cos(w*(n-1)*dt)*pow(dt,2) + gamma*(u[1][i+1][j] + u[1][i-1][j] +1./3.*(u[1][i][j+1] + u[1][i][j-1]) -8.0/3.*u[1][i][j]) +2*u[1][i][j] - u[0][i][j] +b*u[0][i][j]*dt);
                    }
                    else{
                        u[2][i][j] = 1./(1+b*dt)*(gamma*(u[1][i+1][j] + u[1][i-1][j] +1./3.*(u[1][i][j+1] + u[1][i][j-1]) -8.0/3.*u[1][i][j]) +2*u[1][i][j] - u[0][i][j]+b*u[0][i][j]*dt);
                    }
                }
            }

        }

        cc(u);        

        if (n%100==0){
            for (int j =0; j<=ic; j++){
                for (int i =0; i<Nx-1; i++){
                    fprintf(output,"%lf ", u[2][i][j]);
                }
            fprintf(output,"%lf\n",u[2][Nx-1][j]);
            }
        }

        actualitzar(u);

    }

    return 0;

}