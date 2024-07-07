#include <stdio.h>
#include <math.h>

//VARIABLES I CONSTANTS ####################################################################################################################
#define T 2001 //PUNTS DEL MALLAT TEMPORAL
#define Nx 37 //PUNTS DEL MALLAT ESPAIAL
#define Ny 73

int ic = (Nx-1)/2; //Coordenades centre a la nova base
int jc = (Ny-1)/2;

double u[3][Nx][Ny]; 

const double PI = 3.14159265359;

double da = 1./(Nx-2);
double dt = 10./(T-1);

double b = 0.02; //ESMORTIMENT
double f0 = 5.; //AMPLITUD NORMALITZADA
//FUNCIONS ####################################################################################################################

void cc(double u[3][Nx][Ny]){ //CONDICIONS DE CONTORN

    //VERTEX 1(A)
    u[2][0][jc/2] = u[2][1][jc/2+1];
    //VERTEX 2(B)
    u[2][0][3*jc/2] = u[2][1][3*jc/2-1];
    //VERTEX 3(C)
    u[2][ic][Ny-1] = u[2][ic][Ny-2];
    //VERTEX 4 (D)
    u[2][Nx-1][3*jc/2] = u[2][Nx-2][3*jc/2-1];
    //VERTEX 5 (E)
    u[2][Nx-1][jc/2] = u[2][Nx-2][jc/2+1];
    //VERTEX 6 (F)
    u[2][ic][0] = u[2][ic][1];

    //ARESTES
    for (int i =1; i<ic; i++){
        //FRONTERA 1 [1 p]
        u[2][ic +i][Ny-i-1]   = 0.25*(u[2][ic +i-1][Ny-i-1] + 3*u[2][ic +i][Ny-i-2]);// 
        //FRONTERA 2 [1 p]
        u[2][i][3*jc/2+i]     = 0.25*(u[2][i+1][3*jc/2+i] + 3*u[2][i][3*jc/2+i-1]);//
        //FRONTERA 4
        u[2][i][jc/2-i]       = 0.25*(u[2][i+1][jc/2-i] + 3*u[2][i][jc/2-i+1]);
        //FRONTERA 5
        u[2][ic+i][i]         = 0.25*(u[2][ic+i-1][i] + 3*u[2][ic+i][i+1]);
    }
    
    for (int i =1; i<Nx-1; i++){
        //FRONTERA 3 
        u[2][0][jc/2 + i]     = u[2][1][jc/2 + i];
        //FRONTERA 6
        u[2][Nx-1][jc/2 + i]  = u[2][Nx-2][jc/2 + i];          
    }
}

void actualitzar(double u[3][Nx][Ny]){ //EMPREM AQUET CODI PER EMPRAR MENYS MEMORIA
for (int i =0; i<Nx; i++){ //ACTUALITZEM LA MATRIU u
        for (int j =0; j<Ny; j++){

            if(u[2][i][j]>1000){
                    printf("%d,%d\n",i,j);
                }

            u[0][i][j] = u[1][i][j];
            u[1][i][j] = u[2][i][j];
        }
    }
}
//PROGRAMA ####################################################################################################################
int main(){
    double gamma = (1.0/3.0)*pow(dt,2)/(pow(da,2)); 

    FILE *output = fopen("hex_freq.txt","w");

    // freq
    for (int f =1; f<=400; f++){
        double w = f/40.0*PI;
        double umax = 0.0;
        double emax = 0.00;

        //CONDICIONS INICIALS
        for (int i = 0; i<Nx; i++){
            for (int j =0; j<Ny; j++){
                u[0][i][j] = 0.;
                u[1][i][j] = 0.;
                u[2][i][j]=0.;
            }
        }

        //Euler explicit
        for (int n =1; n<T; n++){ //Euler explicit
        
            for (int i = 0; i<Nx; i++){
                if (i<ic){
                    for (int j = ic-i+1; j<=3*ic +i; j++){
                        u[2][i][j] = 1./(1+b*dt)*(gamma*((u[1][i+1][j] + u[1][i-1][j] -2*u[1][i][j]) + 3*(u[1][i][j+1] + u[1][i][j-1] -2*u[1][i][j]))+2*u[1][i][j] - u[0][i][j]+b*u[0][i][j]*dt);
                    }
                }

                else{
                    for (int j= i-ic; j< Ny-(i-ic); j++){
                        if ((i == ic) && (j == jc)){
                            u[2][ic][jc] = 1./(1+b*dt)*(f0*sin(w*(n-1)*dt)*pow(dt,2) + gamma*((u[1][i+1][j] + u[1][i-1][j] -2*u[1][i][j]) + 3*(u[1][i][j+1] + u[1][i][j-1] -2*u[1][i][j])) +2*u[1][i][j] - u[0][i][j] +b*u[0][i][j]*dt);
                        }
                        else{
                            u[2][i][j] = 1./(1+b*dt)*(gamma*((u[1][i+1][j] + u[1][i-1][j] -2*u[1][i][j]) + 3*(u[1][i][j+1] + u[1][i][j-1] -2*u[1][i][j]))+2*u[1][i][j] - u[0][i][j]+b*u[0][i][j]*dt);
                        }
                    }
                }
            }

            cc(u);

            
            double v = (u[2][ic][jc] -u[1][ic][jc])/dt;
            double e = pow(v,2);

            if (fabs(u[2][ic][jc]) > umax){
                umax = fabs(u[2][ic][jc]);
            }
            if (e > emax){
                emax = e;
            }
            
            
            actualitzar(u);
        }

        fprintf(output, "%lf %lf %lf\n", w, umax, emax);
    }
        
    return 0;
}