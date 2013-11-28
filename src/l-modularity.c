#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#define NPotts 30
#define NP 44
#define tinit .1
#define t_step .0005
#define nmatriz 10
#define V 10
#define nmonte 50000*V

double hamilton (int *S, double **J){
   int i, j;
   double H = 0;
   for (i = 1; i < V ; i++)
      for (j = 0; j < i ; j++)
         if(S[i]==S[j]) H-= J[i][j];
   return H;
}

int flipster (int *S, int x){

   int i, j, e_novo, e_velho;
   e_velho = S[x];
   e_novo = 1 + (int)((rand()%(NPotts)));
   S[x] = e_novo;
   return e_velho;
}

void annealing (double **J){

   int *S, i, j, k, l, ntot, N[NPotts+1], NE, e_velho, x;
   double a, t, H0, H1, rand01, sumH,  sumH2, nlocal, m=0, deltaH;

   S = (int *) malloc ((V)*sizeof(int)); //vetores de spins

   for (i = 0; i < V; i++)
      S[i] = 1 + (rand()%NPotts);

   /* Laco de Temperatura */
   t = tinit;
   H0 = hamilton (S, J);
   while(t>0.0){
      /* Inicializando variaveis */
      sumH = 0;
      sumH2 = 0;
      nlocal = 0;
      /* Fim de inicialização */

      /* Comeco Mesmo */
      for (l = 1; l < nmonte; l++){

         if(l%V==0){
            sumH += H0;
            sumH2 += H0*H0;
            nlocal++;
         }
         x = (int)((rand()%(V)));
         e_velho = flipster(S, x);
         H1 = hamilton(S, J);
         deltaH = H1-H0;
         if(deltaH>0){
            rand01 = (double)(rand()) / ((double)(RAND_MAX));
            if(rand01 < exp(-deltaH/t))
               {
                  H0=H1;
               }
            else
               {
                  S[x]=e_velho;
               }
         }
         else {H0 = H1;}

      }
      /*Final*/

      sumH /=nlocal;
      sumH2 /=nlocal;
      sumH2 = sumH2 - sumH*sumH;
      fprintf(Output, "%lf %.20lf %.20lf\n",t, sumH, sumH2);
      t-= t_step;
   } // Fim Laco de temperatura

}

void prepara_J (double **J){

   int i, j, k;
   double K[V], m = 0;

   for(i = 0; i < V; i++)
      for(j = 0; j < V; j++){
         K[i] = 0;
      }

   for (i = 0; i < V; i++)
      for( j = 0; j < V; j++)
         J[i][j] /= sqrt(J[i][i]*J[j][j]);

   for(i = 0; i < V; i++){
      for(j = 0; j < V; j++){
         J[i][i] = 0;
         J[j][j] = 0;
         K[i] += J[i][j];
      }
      m += K[i];
   }

   for( i = 0; i < V; i++){
      for( j = 0; j < V; j++){
         J[i][j] = J[i][j] - (K[i]*K[j])/m;
      }
   }
}
