#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"

#define GB(i, j) g_block[i-1+(j-1)*n]
#define GI(i, j) g_inv[i-1+(j-1)*n]
#define PR(i, j) prop[i-1+(j-1)*n]
#define ARR_D(i) (double *) calloc(i, sizeof(double))
#define ARR_DC(i) (doublecomplex *) calloc(i, sizeof(doublecomplex))

int ind(int *list);
void spins(int index, int *s);
int ns;
integer n;

int main()
{
    int zgeevx_(char *balanc, char *doevl, char *doevr, char *rcnum,
	        integer *n, doublecomplex *prop, integer *ldprop,
	        doublecomplex *eval, doublecomplex *evecl, integer *ldevecl,
	        doublecomplex *evecr, integer *ldevecr, integer *ilo,
	        integer *ihi, doublereal *scale, doublereal *bnorm,
	        doublereal *rceval, doublereal *rcevecr,
	        doublecomplex *work, integer *dwork, doublereal *morework,
	        integer *info);
    int dsyev_(char *jobz, char *uplo, integer *n, double *prop,
                integer *lda, double *w, double *work,
                integer *lwork, integer *info);



    double *ham,*eval,*work,t,*dens,pr;
    doublecomplex *psi,*psit,*ov;
    integer lwork,info;

    FILE *data;
    int i,j,k,kp;
    int *list1,*list2;
    
    
    ns=13; // number of spins
    n=1 << ns; //matrix dimension
    double B=0.5; //magnetic field in z direction
    double Jz=0.0; //NN interaction in z
    double J=1.0;
    double alpha=0.5;
    double sumV,sign,sum;
    
    ham=ARR_D(n*n); //hamiltonian
    eval=ARR_D(n); dens=ARR_D(n);
    psi=ARR_DC(n); psit=ARR_DC(n); ov=ARR_DC(n);
    lwork=34*n;
    work = ARR_D(lwork);
    list1=calloc(ns,sizeof(int));
	list2=calloc(ns,sizeof(int));
    
    for (i=0;i<n*n;i++) ham[i]=0.0;
    
    for (i=0;i<n;i++){
    	spins(i,list1); // list of spins for basis state i
    	spins(i,list2);
    	for (k=0;k<ns;k++) {
    		list2[k]=1-list2[k]; //flip kth spin
    		j=ind(list2);
    		ham[i*n+j]=B;
    		list2[k]=1-list2[k];
		}
		for (k=0;k<ns-1;k++) {
    		list2[k]=1-list2[k]; //flip kth spin
    		list2[k+1]=1-list2[k+1]; //flip k+1th spin
    		j=ind(list2);
    		ham[i*n+j]=Jz;
    		list2[k]=1-list2[k];
    		list2[k+1]=1-list2[k+1];
		}
		sumV=0.0;
		for(k=0;k<ns;k++) for (kp=0;kp<k;kp++){
			if (list1[k]==list1[kp]) sign=1; else sign=-1;
			sumV+=sign*J/pow(k-kp,alpha);
		}
		ham[i*n+i]=sumV;
	}
    
   dsyev_("V","U",&n,ham,&n,eval,work,&lwork,&info); //diagonalize ham, on output contains eigenstates
   printf("info %d %lf\n",info,work[0]/n);
   
   //for (i=0;i<n;i++) printf("%lf\n",eval[i]);
   
   for (i=0;i<n;i++) psi[i].r=psi[i].i=0.;
   for (i=0;i<ns;i++) list1[i]=1;
   list1[(ns-1)/2]=0;
   psi[ind(list1)].r=1;    // initialize psi, central spin left
   //for(i=0;i<n;i++) printf("%lf %lf\n",psi[i].r,psi[i].i);
   
   for(i=0;i<n;i++) {
   		ov[i].r=ov[i].i=0.0; // overlap of initial state with i-th eigenstate
		for (j=0;j<n;j++) ov[i].r+=psi[j].r*ham[i*n+j]; //here psi and eigenstates are real
		//printf("%d %lf\n",i,ov[i].r);
	}		
	
//	for (t=0;t<9.9;t+=0.5) { // start time loop
 	for (t=0;t<3e9;t+=0.3e9) { // start time loop
		for(j=0;j<n;j++) {  //sum over basis states
			psit[j].r=psit[j].i=0.0;
			for(i=0;i<n;i++) { //sum over eigenstates
				psit[j].r+=ov[i].r*cos(eval[i]*t)*ham[i*n+j];
				psit[j].i-=ov[i].r*sin(eval[i]*t)*ham[i*n+j];
			}
		}
		//sum=0.0;
		//for(j=0;j<n;j++) sum+=pow(psit[j].r,2)+pow(psit[j].i,2);
		//printf("%lf %lf\n",t,sum);
		
		// now let's calculate spin densitties at time t
		for (k=0;k<ns;k++) dens[k]=0.; //initialize densities to 0
		for (i=0;i<n;i++) { //loop over basis states
			pr=pow(psit[i].r,2)+pow(psit[i].i,2); // prob to be in that basis state
			spins(i,list1); // get list of spins
			for (k=0;k<ns;k++) {
				if (list1[k]) dens[k]+=pr; else dens[k]-=pr;
			}
		}
		printf("%lf",t);
		for (k=0;k<ns;k++) printf(" %lf",dens[k]);
		printf("\n");
	}	
    
/*    for (i=0;i<n;i++) {
		printf("%d: ",i);
		spins(i,list1);
		for (j=0;j<ns;j++) printf("%d ",list1[j]);
		printf("\n");
		printf("%d\n",/*ind(list1));
	}
*/

    
//   time(&n); 


    return 0;
}

int ind(int *list)   // find basis element conrresponding to spin chain "list"
{
	int i,j=0;
	for(i=0;i<ns;i++)
		j+=list[i]*pow(2,i);
	return j;
}

void spins(int index, int *s) // convert spin chain "s" to integer basis element index

{
	int i;
	for (i=0;i<ns;i++)
		s[i]=(1<<i & index)>0;
	return;	
}


