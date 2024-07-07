#include <stdio.h>
#include <math.h>

//VARIABLES I CONSTANTS ####################################################################################################################

#define N 101 //PUNTS DEL MALLAT ESPAIAL
#define T 10001 //PUNTS DEL MALLAT TEMPORAL

const double pi = 3.14159265358979323846;
const int ic = (N-1)/2; //PUNT CENTRAL

double dt = 4./(T-1);
double dx = 1./(N-1); 

double u[3][N][N]; //MATRIU DE PERTORBACIONS NORMALITZADA
double oldu[N][N];

double f0=5.; //AMPLITUD NORMALITZADA
double w= pi*sqrt(pow(1,2)+pow(2,2));
double b = 0.02;

//FUNCIONS ####################################################################################################################


double f(int i, int j,double t){ //FUNCIÓ PROPORCIONAL A LA FORÇA APLICADA A LA PLACA
    double ft;

    if (i==ic && j==ic){
        ft=f0*cos(w*t);
    }
    else{
        ft=0;
    }

    return(ft);
}

void cc(double u[3][N][N], int t){ //CONDICIONS DE CONTORN

    for (int i=1; i<N-1; i++){ //ARESTES
        u[t][i][0]=u[t][i][1];
        u[t][i][N-1]=u[t][i][N-2];
   
        u[t][0][i]=u[t][1][i];
        u[t][N-1][i]=u[t][N-2][i];      
    }

    u[t][0][0]=u[t][1][1]; //VERTEX
    u[t][N-1][N-1]=u[t][N-2][N-2];
    u[t][0][N-1]=u[t][1][N-2];
    u[t][N-1][0]=u[t][N-2][1];
}

void actualitzar(double u[3][N][N]){ //EMPREM AQUET CODI PER EMPRAR MENYS MEMORIA
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            u[0][i][j]=u[1][i][j];
            u[1][i][j]=u[2][i][j];
        } 
    }
}

//PROGRAMA ####################################################################################################################

int main(){

    double gamma= dt/dx; //CONSTANT GAMMA PER SIMPLIFICAR EXPRESSIONS

    FILE *output;
    output = fopen("quad_imp.txt", "w");

    for (int i=0; i<N; i++){ //CONDICIONS INICIALS (EN REPOS)
        for (int j=0; j<N; j++){
            u[0][i][j]=0.;
            u[1][i][j]=0.;
        }
    }
   
    for (int i=0; i<N; i++){
        fprintf(output,"%lf",u[0][i][0]);
        for (int j=1; j<N; j++){
            fprintf(output," %lf",u[0][i][j]);
        }
        fprintf(output,"\n");
    }

    for (int t=1; t<T; t++){ //EULER IMPLICIT

        for (int n=1; n<100000; n++){ //NUM ARBITRARIAMENT GRAN PER ASSEGURAR CONVERGENCIA

            double error = 0.; //DEFINIM UN CONTADOR PER L'ERROR

            for (int i=1; i<N-1; i++){
                for (int j=1; j<N-1; j++){

                    oldu[i][j]=u[2][i][j];

                    //EXTREMS
                    
                    if(i==1){

                        if(j==1){  //i-1= i, j-1=j

                            u[2][i][j]=
                                (2*u[1][i][j]-u[0][i][j]
                                +pow(gamma,2)*0.5*(u[2][i+1][j]+u[2][i][j+1]
                                -2*u[0][i][j]+u[0][i+1][j]+u[0][i][j+1])
                                +b*0.5*dt*(u[0][i][j]))/(1.+pow(gamma,2)+(b*dt)/2.);

                        }
                        else if (j==N-2){ //i-1= i, j+1=j

                            u[2][i][j]=
                                (2*u[1][i][j]-u[0][i][j]
                                +pow(gamma,2)*0.5*(u[2][i+1][j]+u[2][i][j-1]
                                -2*u[0][i][j]+u[0][i+1][j]+u[0][i][j-1])
                                +b*0.5*dt*(u[0][i][j]))/(1.+pow(gamma,2)+(b*dt)/2.);

                        }
                        else{ //i-1=i

                            u[2][i][j]=
                                (2*u[1][i][j]-u[0][i][j]
                                +pow(gamma,2)*0.5*(u[2][i+1][j]+u[2][i][j+1]+u[2][i][j-1]
                                -3*u[0][i][j]+u[0][i+1][j]+u[0][i][j+1]+u[0][i][j-1])
                                +b*0.5*dt*(u[0][i][j]))/(1.+3.*pow(gamma,2)/2.+(b*dt)/2.);

                        } 

                    }

                    else if(i==N-2){

                        if(j==1){ //i+1=i, i-1=j

                            u[2][i][j]=
                                (2*u[1][i][j]-u[0][i][j]
                                +pow(gamma,2)*0.5*(u[2][i-1][j]+u[2][i][j+1]
                                -2*u[0][i][j]+u[0][i-1][j]+u[0][i][j+1])
                                +b*0.5*dt*(u[0][i][j]))/(1.+pow(gamma,2)+(b*dt)/2.);

                        }
                        else if (j==N-2){ //i+1=i, j+1=j

                            u[2][i][j]=
                                (2*u[1][i][j]-u[0][i][j]
                                +pow(gamma,2)*0.5*(u[2][i-1][j]+u[2][i][j-1]
                                -2*u[0][i][j]+u[0][i-1][j]+u[0][i][j-1])
                                +b*0.5*dt*(u[0][i][j]))/(1.+pow(gamma,2)+(b*dt)/2.);

                        }
                        else{ //i+1=i

                            u[2][i][j]=
                                (2*u[1][i][j]-u[0][i][j]
                                +pow(gamma,2)*0.5*(u[2][i-1][j]+u[2][i][j+1]+u[2][i][j-1]
                                -3*u[0][i][j]+u[0][i-1][j]+u[0][i][j+1]+u[0][i][j-1])
                                +b*0.5*dt*(u[0][i][j]))/(1.+3.*pow(gamma,2)/2+(b*dt)/2.);

                        } 

                    }

                    else{

                        if (j==1){ //j-1=j

                            u[2][i][j]=
                                (2*u[1][i][j]-u[0][i][j]
                                +pow(gamma,2)*0.5*(u[2][i+1][j]+u[2][i-1][j]+u[2][i][j+1]
                                -3*u[0][i][j]+u[0][i+1][j]+u[0][i-1][j]+u[0][i][j+1])
                                +b*0.5*dt*(u[0][i][j]))/(1.+3.*pow(gamma,2)/2.+(b*dt)/2.);

                        }

                        else if (j==N-2){ //j+1=j

                            u[2][i][j]=
                                (2*u[1][i][j]-u[0][i][j]
                                +pow(gamma,2)*0.5*(u[2][i+1][j]+u[2][i-1][j]+u[2][i][j-1]
                                -3*u[0][i][j]+u[0][i+1][j]+u[0][i-1][j]+u[0][i][j-1])
                                +b*0.5*dt*(u[0][i][j]))/(1.+3.*pow(gamma,2)/2.+(b*dt)/2.);

                        }
                        
                        else { // PUNTS INTERIOR

                            u[2][i][j]=
                                (2*u[1][i][j]-u[0][i][j]
                                +pow(gamma,2)*0.5*(u[2][i+1][j]+u[2][i-1][j]+u[2][i][j+1]+u[2][i][j-1]
                                -4*u[0][i][j]+u[0][i+1][j]+u[0][i-1][j]+u[0][i][j+1]+u[0][i][j-1])
                                +f(i,j,(t-1)*dt)*pow(dt,2)
                                +b*0.5*dt*(u[0][i][j]))/(1.+2.*pow(gamma,2)+(b*dt)/2.);
                        }
                    }
                }
            }

            for (int i=1;i<N-1;i++){
                for (int j = 1; j<N-1; j++){
                    error += fabs(u[2][i][j]-oldu[i][j]); //VALOR ABSOLUT DE LA DIFERENCIA
                }
            }

            if (error <pow(10.,-6)){ //TRENQUEM EL FOR EN ARRIBAR A LA PRECICIÓ DESSITJADA (EN CAS CONTRARI CONTINUEM)
                break;
            }

        }

        cc(u,2); //CONDICIONS DE CONTORN
        actualitzar(u); //ACUALITZAR LA MATRIU u

        if(t%100==0){ //ACTIVAR O DESACTIVAR CONDICIÓ PER NO ESCRIURE/ESCRIURE TOTES LES MATRIUS
            for (int i=0; i<N; i++){
                fprintf(output,"%lf",u[2][i][0]);
                for (int j=1; j<N; j++){
                    fprintf(output," %lf", u[2][i][j]);
                }
                fprintf(output,"\n");   
            }
        }

    }

    return 0;
}