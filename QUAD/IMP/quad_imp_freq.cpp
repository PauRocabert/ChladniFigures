#include <stdio.h>
#include <math.h>
#include <iostream>

//VARIABLES I CONSTANTS ####################################################################################################################

#define N 37 //puntos en el mallado espacial
#define T 2001 //puntos en el mallado temporal

const double pi = 3.14159265358979323846;
const int ic = (N-1)/2; //PUNT CENTRAL

double dt = 10./(T-1);
double dx = 1./(N-1); 

double u[3][N][N]; //MATRIU DE PERTORBACIONS NORMALITZADA
double oldu[N][N];

double f0=5.; //AMPLITUD NORMALITZADA
double b = 0.02;

//FUNCIONS ####################################################################################################################


double f(int i, int j, double t, double w){ //FUNCIÓ PROPORCIONAL A LA FORÇA APLICADA A LA PLACA
    double ft;

    if (i==ic && j==ic){
        ft=f0*sin(w*t);
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

    double gamma= dt/dx;

    FILE *output;
    output = fopen("quad_imp_freq.txt", "w");

    for (int dw=1; dw<=250; dw++){ //ITERAM PER DIFERENTS FREQUENCIES

        double w = pi*dw/10;
        double umax = 0;
        double emax = 0;

        for (int i=0; i<N; i++){ //CONDICIONS INICIALS (EN REPOS)
            for (int j=0; j<N; j++){
                u[0][i][j]=0.;
                u[1][i][j]=0.;
            }
        }

        for (int t=1; t<T+1; t++){ //EULER IMPLICIT PER A CADA W

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
                            
                            else { //PUNTS INTERIORS

                                u[2][i][j]=
                                    (2*u[1][i][j]-u[0][i][j]
                                    +pow(gamma,2)*0.5*(u[2][i+1][j]+u[2][i-1][j]+u[2][i][j+1]+u[2][i][j-1]
                                    -4*u[0][i][j]+u[0][i+1][j]+u[0][i-1][j]+u[0][i][j+1]+u[0][i][j-1])
                                    +f(i,j,(t-1)*dt,w)*pow(dt,2)
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

            cc(u,2);

            double v = (u[2][ic][ic] -u[1][ic][ic])/dt;
            double e = pow(v,2)/dt; //relacio de proporcionalitat

            if (fabs(u[2][ic][ic])>umax){ //AMPLITUT MAXIMA

                umax=fabs(u[2][ic][ic]);
            }

            if (e>emax){ //POTENCIA MAXIMA
                emax=e;
            }
            
            actualitzar(u);

        }

        fprintf(output,"%lf %lf %lf\n", w, umax, emax);
    }

    return 0;
}