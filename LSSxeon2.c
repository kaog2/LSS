#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#define MINNP 4 //!< Minimun number of minutiae in the computation of n_P
#define MAXNP 12 //!< Maximum number of minutiae in the computation of n_P
#define TAUP 0.4 //!< Sigmoid parameter 2 in the computation of n_P
#define MUP 20 //!< Sigmoid parameter 1 in the computation of n_P
#define MINME 0.6 //!< Minimum number of matching elements in two matchable cylinders
#define DELTAZETA 1.57079633 //!< Maximum global rotation allowed between two templates
#define REGTORAD 0.017453293 //!< Radians to degrees
#define PI 3.14159265 //!< Value of PI
#define PI2 6.283185307 //!< Value of PI*2
#define NUMCELL 384


/*Points*/
struct point {
  int x;
  int y;
};

typedef struct point Point;


__declspec (target(mic)) double dtime()
{
  double tseconds = 0.0;
  struct timeval mytime;
  gettimeofday(&mytime,(struct timezone*)0);
  tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
  return( tseconds );
}

__attribute__((target(mic))) float similarity(char *tq, char *tDb, int i, int j, int MINCELL);
__attribute__((target(mic))) float minimo(int n_A, int n_B);
__attribute__((target(mic))) float computeNP(int n_A, int n_B);
__attribute__((target(mic))) int roundInt(float r);
__attribute__((target(mic))) float psi(float v, float par1, float par2);
__attribute__((target(mic))) void consolidationLSS(Point *best, float *gamma, unsigned int nP, unsigned int nA, unsigned int nB);
__attribute__((target(mic))) float computeGlobalScore(Point *elements,float *gamma, int n_P, int nA);
__attribute__((target(mic))) float dFi(float ang1, float ang2);


/*
  offset = (row * NUMCOLS) + column

  matrix[ row ][ col ] = array[ row*m + j ].

  dim 2x3

  1 2 3
  4 5 6

  0*3 + 0 = 1 ; 0*3 +1 = 2 ; 0*3 + 2 = 3 ; 1* 3 + 0 = 4 ; 1*3 + 1 = 5 ; 1*3 + 2 = 6 
  */
  /*
  float minimo(int n_A, int n_B);
  float computeNP(int n_A, int n_B);
  int roundInt(float r);
  float psi(float v, float par1, float par2);
  void consolidationLSS(Point *best, float **gamma, unsigned int nP, unsigned int nA, unsigned int nB);
  float computeGlobalScore(Point *elements,float **gamma, int n_P);
  float dFi(float ang1, float ang2);
*/

int main(){//work good only with 200 fingerprints.. not with all the database . 

  int CantQuery = 27000; 
  int CantDatabase = 27000;

  double tinicio, tparada, ttotal = 0;
  int i,j;
  char num;
  char *database, *query;

  int *index_q; //idex of each cilinder of each finger
  int *index_db; // idex of each cilinder of each finger

  int *can_cilindros_query; //number of cilinder of a fingerprints
  int *can_cilindros_database; //number of cilinder of a fingerprints
  int ccq=0; //cilindros por huella
  int aux=0; //para avanzar en puntero cantidad de cilindros de cada huella 

  int pos_match=0;
  int huellas_con_cilindros = 5168704; //Cylinder for 100 fingerprints = 18422
  int MINCELL; //minimum of cells

  can_cilindros_query=(int *)malloc((CantQuery)*sizeof(int *));      
  can_cilindros_database=(int *)malloc((CantDatabase)*sizeof(int *));

  index_q=(int *)malloc((CantQuery)*sizeof(int *));      
  index_db=(int *)malloc((CantDatabase)*sizeof(int *));

  MINCELL = floor(MINME*NUMCELL);

  /*******************************para query************/
  FILE* tq = fopen("/home/pregrado/tesistas/kortega/measure/cylinders8NOhneGrad","r");
  //FILE* tq = fopen("/home/kevin/Desktop/basesDeDatos/100Cylinders8NFormato","r");
  char item[30];
  int c_cantCilindros;
  int ret = fscanf(tq, "%i" ,&c_cantCilindros);

  can_cilindros_query[aux]=c_cantCilindros; //number of cilinder of a fingerprints
  index_q[aux]=0; // first cilynder begint on the first place
  aux++;
  query=(char *)malloc((huellas_con_cilindros*NUMCELL)*sizeof(char *));

  for (i = 0; i < huellas_con_cilindros; ++i){

	    //if the fingerprint is save complete ... ready to read the next c_cantCilindros
	    if (ccq==c_cantCilindros)
	    { 
	      ret = fscanf(tq, "%i" ,&c_cantCilindros); 
	      can_cilindros_query[aux] = c_cantCilindros; // is saved the number of Cylinder of a fingerprints
	      index_q[aux] = i;// is saved index of the begin of a fingerprints
	      aux++;
	      ccq=0;
	    }
	
	    for (j = 0; j < NUMCELL; ++j)
	    {   
	               
	      ret = fscanf(tq, "%c %c" ,item,&num);
	      if(ret == 2){ 
	        if (num !=10)
	        { 
	          query[i* NUMCELL + j] = num;
	          //printf("%c\t",query[i* NUMCELL + j]);
	        }else{
	          printf("espacio");
	        }
	      }
	      else if(errno != 0) {
	         perror("scanf:");
	         break;
	      } else if(ret == EOF) {
	         break;
	      } else {
	         printf("Fin de lectura\n");
	      } 
	    }
	    ccq++;
  }
         
  //fclose(tq);
  printf("carga Q lista \n");
  /**********************************carga de minutas terminado************************************************/
     /*******************************para Database************/
  aux=0;//para avanzar en puntero cantidad de cilindros de cada huella 
  ccq=0;//cilindros por huella

  //FILE* tDb = fopen("/home/kevin/Desktop/basesDeDatos/huella_6_cilindros","r");
  //FILE* tDb = fopen("/home/kevin/Desktop/basesDeDatos/100Cylinders8NFormato","r");
  //ret = fscanf(tDb, "%i" ,&c_cantCilindros);
  ret = fscanf(tq, "%i" ,&c_cantCilindros);

  can_cilindros_database[aux] = c_cantCilindros; // first cilynder begint on the first place
  index_db[aux]=0; // first cilynder begint on the first place

  aux++;
  database=(char *)malloc((huellas_con_cilindros*NUMCELL)*sizeof(char *));

  for (i = 0; i < huellas_con_cilindros; ++i){

	    //if the fingerprint is save complete ... ready to read the next c_cantCilindros
	    if (ccq==c_cantCilindros)
	    {
	      ret = fscanf(tq, "%i" ,&c_cantCilindros); 
	      can_cilindros_database[aux] = c_cantCilindros; // is saved the number of Cylinder of a fingerprints
	      index_db[aux] = i;// is saved index of the begin of a fingerprints
	      aux++;
	      ccq=0;
	    }

	    for (j = 0; j < NUMCELL; ++j)
	    {   
	      int ret = fscanf(tq, "%c %c" ,item,&num);
	      if(ret == 2){ 
	        if (num !=10)
	        {
	          database[i* NUMCELL + j]=num;
	         
	        }else{
	          printf("espacio");
	        }        
	      }
	      else if(errno != 0) {
	        perror("scanf:");
	        break;
	      }else if(ret == EOF) {
	        break;
	      }else {
	        printf("Fin de lectura\n");
	      }                   
	    }

	    ccq++;
 	}

  //fclose(tDb);
  fclose(tq);
  printf("carga DB lista \n");
	/**********************************Fin carga DB************************************************/
	/**********************************Comparaciones************************************************/

		#pragma offload_transfer target(mic:0) mandatory\
    in(index_q:length(CantQuery) alloc_if(1) free_if(0))\
    in(index_db:length(CantDatabase) alloc_if(1) free_if(0))\
    in(can_cilindros_query:length(CantQuery) alloc_if(1) free_if(0))\
    in(can_cilindros_database:length(CantDatabase) alloc_if(1) free_if(0))\
    in(query:length(huellas_con_cilindros*NUMCELL) alloc_if(1) free_if(0))\
    in(database:length(huellas_con_cilindros*NUMCELL) alloc_if(1) free_if(0))\
    in(MINCELL:alloc_if(1) free_if(0))

    #pragma offload target(mic:0)\
    nocopy(index_q:length(CantQuery) alloc_if(0) free_if(0))\
    nocopy(index_db:length(CantDatabase) alloc_if(0) free_if(0))\
    nocopy(can_cilindros_query:length(CantQuery) alloc_if(0) free_if(0))\
    nocopy(can_cilindros_database:length(CantDatabase) alloc_if(0) free_if(0))\
    nocopy(query:length(huellas_con_cilindros*NUMCELL) alloc_if(0) free_if(0))\
    nocopy(database:length(huellas_con_cilindros*NUMCELL) alloc_if(0) free_if(0))\
    nocopy(MINCELL)
    {

	    omp_set_num_threads(20);
	    int co,db;//position of the Fingerabdruck
	    int pos_match;//position of the match fingerprint

	    #pragma omp parallel private(co,db,pos_match) /*shared(disVectores)*/ // un Hilo para cada cilindro .. hacer esto
  		{

		    int tid=omp_get_thread_num();//obtiene el hilo N° : XX

		    //block for 100 queryes
        tinicio = dtime();//captura el tiempo de inicio

		    for (co=tid ; co < CantQuery; co += omp_get_num_threads())
		    {
		        float  near =0;
		        double tPorMinutia=0;
		        float result=0;//similarity value

		        for ( db = 0; db < CantDatabase; ++db)
		        {
		          
		            unsigned int n_P = computeNP(can_cilindros_database[db], can_cilindros_query[co]);

		            if(n_P > can_cilindros_database[db] || n_P > can_cilindros_query[co])
		                printf("asdasdas\n");

		            float *gamma_matrix;
		            gamma_matrix=(float *)malloc((can_cilindros_query[co]*can_cilindros_database[db])*sizeof(float *));

		            Point *pairs = NULL;
		            pairs=(Point *)malloc((n_P)*sizeof(Point *));

		            int star_query = index_q[co]; //can_cilindros_query[co];
		            int star_db = index_db[db];// can_cilindros_database[db];
		            int end_query= index_q[co] + can_cilindros_query[co]; // **make a rest of -1
		            int end_db= index_db[db] + can_cilindros_database[db]; // **make a rest of -1
		            
		            int a,b;
		            int m,n; //for the matrix
		            for (m = 0 , a = star_query ; a < end_query  ;++a,++m){ 
		                for (n=0, b = star_db ; b < end_db; ++b,++n){ // a cilindro a comparar de la query 
		                  gamma_matrix[m*can_cilindros_database[db] +n] = similarity(query, database , b, a, MINCELL);

		                }
		            }

		            consolidationLSS(pairs, gamma_matrix, n_P, can_cilindros_database[db], can_cilindros_query[co]);
		            //l* ***********************************************Termina LSS ****************************************************************
		            result=computeGlobalScore(pairs, gamma_matrix,n_P,can_cilindros_database[db]);

		            if(result==1){

		              //tparada=dtime();  
		              //ttotal=(tparada-tinicio)+ttotal;
		              //tPorMinutia =  (tparada-tinicio) + tPorMinutia;
		              pos_match = db;  
		              near=result;
		              free(gamma_matrix);
		              //break;

		            }else{

		              if (result>=near){
		                //in the case of the result is not still the best, we go to the next Minutia query to comparate
		                near=result;
		                pos_match = db;
		              }

		              //tparada=dtime();
		              //ttotal=(tparada-tinicio)+ttotal;
		              //tPorMinutia =  (tparada-tinicio) + tPorMinutia;
		              //terminada la comparacion con el primer cilindro y se libera la memoria del cilindro ya comparado
		              free(gamma_matrix);
		              free(pairs);                     
		              //avanza a la siguiente minuta
		            } 
		        }//fin while de 54000

		        printf("Match %f\t with %d \t Query: %d \n",near,pos_match,co);
		        //printf("%f milisegundos \n",tPorMinutia );
		          
		    }//fin 100 queryes
		}//end of pragma

    tparada=dtime();
    ttotal=tparada-tinicio;
    printf("tiempo total es: %f \n",ttotal);
	}//fin offload
	#pragma offload_transfer target(mic:0) mandatory\
    out(index_q:length(CantQuery) alloc_if(0) free_if(1))\
    out(index_db:length(CantDatabase) alloc_if(0) free_if(1))\
    out(can_cilindros_query:length(CantQuery) alloc_if(0) free_if(1))\
    out(can_cilindros_database:length(CantDatabase) alloc_if(0) free_if(1))\
    out(query:length(huellas_con_cilindros*NUMCELL) alloc_if(0) free_if(1))\
    out(database:length(huellas_con_cilindros*NUMCELL) alloc_if(0) free_if(1))\
    out(MINCELL:alloc_if(0) free_if(1))\
    out(pos_match:alloc_if(0) free_if(1))\
    out(ttotal:alloc_if(0) free_if(1))

	return 0;
}

//tq Vector de queruy ... tDb vector de dababase
// i and j are the cilinders, so i is a cylinder and j another cylinder
float similarity(char *tq, char *tDb , int i, int j, int MINCELL){
  unsigned int count = 0;
  int counta_b = 0, countb_a = 0, count_diff = 0;

  //printf("(%d , %d )\n", i,j);

  unsigned int k;
  for (k=0; k<NUMCELL; ++k){ //recorrer cilindro i de DBk veces NUMCELL
    count++;
    if (tDb[i*NUMCELL +k] =='1'){
      counta_b++;
          
      if(tq[j*NUMCELL + k] =='0'){ //recorrer query i k veces NUMCELL
        count_diff++;
      }
    }
        
    if (tq[j*NUMCELL + k]=='1'){
      countb_a++;
          
      if(tDb[i*NUMCELL + k]=='0'){
        count_diff++;
      }
    }
  }   
  //Values of the Matching
    //printf("%f\n",(1 - (sqrt(count_diff)/(sqrt(counta_b)+sqrt(countb_a))))); 

    if (count >= MINCELL)
        return (1 - (sqrt(count_diff)/(sqrt(counta_b)+sqrt(countb_a))));
    else
        return 0;
}

float hamming(char **tq, char **tDb , int cant_cilindro, int i, int j){

    int hamming;
    const float P = 30; //constante seg´un paper

    int cilindrosDB;
    int k;
    hamming = 0;
        
        for (k = 0; k < NUMCELL; ++k)
        {
            if (tq[i][k] != tDb[j][k])
            {   
                hamming++;

            }
        }
    return pow(1.0-((float)hamming)/NUMCELL, P);
}


float computeNP(int n_A, int n_B){
  return MINNP + roundInt(psi(minimo(n_A,n_B),MUP,TAUP*(MAXNP-MINNP)));
}

float minimo(int n_A, int n_B){
  if (n_A < n_B) return n_A;
  else return n_B;
}

int roundInt(float r){
  
  if (r > 0.0)  return (int)(r + 0.5); //para positivos
  else return (int)(r - 0.5);//para negativos

}

float psi(float v, float par1, float par2){
    return 1.0 / (1.0 + exp((par2*(par1-v))));
}

void consolidationLSS(Point *best, float *gamma, unsigned int nP, unsigned int nA, unsigned int nB){
  unsigned int k;
  float scores[nP];
  unsigned int i;
  unsigned int j;
  unsigned int l;

  for(i = 0; i < nP; i++){
    scores[i] = -1.0;
  }

  for (i=0; i<nB; i++){

    for (j=0; j<nA; j++){

      for (k=0; k<nP && scores[k]>=gamma[i*nA + j]; k++);

      if (k < nP){

        for (l=nP-1; l>k; l--){

          scores[l] = scores[l-1];
          best[l] = best[l-1];
        }

        scores[k] = gamma[i*nA + j];
        best[k].x = i;
        best[k].y = j;
      }
    }
  }
}

float computeGlobalScore(Point *elements, float *gamma,int n_P , int nA){
  float sum = 0.0;
  int i;
  for (i = 0; i < n_P; ++i)
  {
      sum += gamma[elements[i].x*nA + elements[i].y]; //[elements[i].x][elements[i].y]; 
   //   printf("puntaje: %f\n", gamma[elements[i].x][elements[i].y]);
  }
  //printf("otra huella \n");

  return sum / n_P;
}
float dFi(float ang1, float ang2){
    float diff = ang1-ang2;

    if(diff >= -PI && diff < PI){
        return diff;
    } else if (diff < -PI) {
        return PI2 + diff;
    } else {
        return -PI2 + diff;
    }
}
