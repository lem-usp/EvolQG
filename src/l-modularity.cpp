#include <Rcpp.h>
using namespace Rcpp;

#define tinit .1
#define NPotts 50
#define t_step .0005

double hamilton (NumericVector S, NumericMatrix J){
   int i, j;
   int V = J.nrow();
   double H = 0;
   for (i = 1; i < V ; i++)
      for (j = 0; j < i ; j++)
         if(S[i]==S[j]) H-= J(i, j);
   return H;
}

int flipster (NumericVector S, int x){
   int i, j, e_novo, e_velho;
   e_velho = S[x];
   e_novo = (int)(R::runif(1, NPotts));
   S[x] = e_novo;
   return e_velho;
}

NumericMatrix prepara_J (NumericMatrix x){

   int i, j, k, V;
   V = x.nrow();
   NumericMatrix J = clone(x);
   NumericVector K(V);
   double m = 0;
   for(i = 0; i < V; i++)
      for(j = 0; j < V; j++){
         K[i] = 0;
      }
   for (i = 0; i < V; i++)
      for( j = 0; j < V; j++)
         J(i,j) /= sqrt(J(i,i)*J(j,j));
   for(i = 0; i < V; i++){
      for(j = 0; j < V; j++){
         J(i,i) = 0;
         J(j,j) = 0;
         K[i] += J(i,j);
      }
      m += K[i];
   }
   for( i = 0; i < V; i++){
      for( j = 0; j < V; j++){
         J(i,j) = J(i,j) - (K[i]*K[j])/(2*m);
      }
   }
   return J;
}

// [[Rcpp::export]]
double annealing (NumericMatrix Corr, NumericVector S){

   NumericMatrix J = prepara_J (Corr);
   int V = Corr.nrow();
   int nmonte = 5000*V;
   int i, j, k, l, ntot, N[NPotts+1], NE, e_velho, x;
   double a, t, H0, H1, rand01, sumH,  sumH2, nlocal, m=0, deltaH;
   for (i = 0; i < V; i++)
      S[i] = (int)(R::runif(1, NPotts));
   t = tinit;
   H0 = hamilton (S, J);
   while(t>0.0){
      sumH = 0;
      sumH2 = 0;
      nlocal = 0;
      for (l = 1; l < nmonte; l++){
         if(l%V==0){
            sumH += H0;
            sumH2 += H0*H0;
            nlocal++;
         }
         x = (int)(R::runif(0, V));
         e_velho = flipster(S, x);
         H1 = hamilton(S, J);
         deltaH = H1-H0;
         if(deltaH>0){
            rand01 = (double)(R::runif(0,1));
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
      sumH /=nlocal;
      sumH2 /=nlocal;
      sumH2 = sumH2 - sumH*sumH;
      t-= t_step;
   }
   return sumH;
}
