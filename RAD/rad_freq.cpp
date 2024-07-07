#include <stdio.h>
#include <math.h>

//VARIABLES I CONSTANTS ####################################################################################################################

#define N 201 //puntos en el mallado espacial
#define T 40001 //puntos en el mallado temporal

double dt = 10./(T-1);
double dr = 1./(N-1);

double u[3][N]; //MATRIU DE PERTORBACIONS NORMALITZADA

double b = 0.02;  
double f0=5.; //AMPLITUD NORMALITZADA

//FUNCIONS######################################################################################################################

void actualitzar(double u[3][N]){ //EMPREM AQUET CODI PER EMPRAR MENYS MEMORIA
    for (int j=0; j<N; j++){
        u[0][j]=u[1][j];
        u[1][j]=u[2][j];
    } 
}

//PROGRAMA ####################################################################################################################

int main(){
    FILE *output;
    output = fopen("rad_freq.txt", "w");
    double gamma = dt/dr;

    for (int dw = 1; dw<= 400; dw++){ //ITERAM PER A DIFERENTS W
        double omega = dw/20.0; 
        double umax = 0.;
        double emax = 0.;

        for (int j =0; j<N; j++){ //CONDICIONS INICIALS
            u[0][j] = 0.;
            u[1][j] = 0.;
        }

        for (int t=1; t<T; t++){ //EULER EXPLICIT
            
            u[2][1] = (f0*sin(omega*dt*(t-1))*pow(dt,2) +2*u[1][1]-u[0][1]+pow(gamma,2)*(u[1][2]-u[1][1]+(u[1][2]-u[1][1])/(2))+b*dt*u[0][1])/(1+b*dt);

            for (int j=2; j<N-1; j++){
                u[2][j] = (2*u[1][j]-u[0][j]+pow(gamma,2)*(u[1][j+1]-2*u[1][j]+u[1][j-1]+(u[1][j+1]-u[1][j-1])/(2*j))+b*dt*u[0][j])/(1+b*dt);
            }

            u[2][0]=u[2][1]; //CONDICIONS DE CONTORN
            u[2][N-1] = u[2][N-2];

            double v = (u[2][0] -u[1][0])/dt;
            double e= pow(v,2);

            if (fabs(u[2][0]) > umax){
                umax = fabs(u[2][0]);
            }

            if (e > emax){
                emax = e;
            }
            

            actualitzar(u);
        }

        fprintf(output,"%lf %lf %lf\n", omega, umax, emax);
    }
    
    return 0;
}