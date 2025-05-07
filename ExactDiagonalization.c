#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"

#define GB(i, j) g_block[i-1+(j-1)*n]
#define GI(i, j) g_inv[i-1+(j-1)*n]
#define PR(i, j) prop[i-1+(j-1)*n]
#define ARR_D(i) (double *) calloc(i, sizeof(double))
#define ARR_DC(i) (doublecomplex *) calloc(i, sizeof(doublecomplex))

int ind(int* list);
void spins(int index, int* s);
int ns;
int lmid;
integer n;

int main()
{
	int zgeevx_(char* balanc, char* doevl, char* doevr, char* rcnum,
		integer * n, doublecomplex * prop, integer * ldprop,
		doublecomplex * eval, doublecomplex * evecl, integer * ldevecl,
		doublecomplex * evecr, integer * ldevecr, integer * ilo,
		integer * ihi, doublereal * scale, doublereal * bnorm,
		doublereal * rceval, doublereal * rcevecr,
		doublecomplex * work, integer * dwork, doublereal * morework,
		integer * info);
	int dsyev_(char* jobz, char* uplo, integer * n, double* prop,
		integer * lda, double* w, double* work,
		integer * lwork, integer * info);



	double* ham, * eval, * work, t, * dens, * probSpinFlip, pr;
	doublecomplex* psi, * psit, * ov;
	integer lwork, info;

	
	int i, j, k, kp, d;
	int* list1, * list2;

	//Writing to file to another directory
	char* file = "D:\\OneDrive - Tulane University\\RESEARCH\\Quantum transport\\Code\\Plots\\NewMeasures\\Test_a.txt";
	FILE* fh_output;
	fh_output = fopen(file, "w");

	FILE* fh_output2;
	fh_output2 = fopen("D:\\OneDrive - Tulane University\\RESEARCH\\Quantum transport\\Code\\Plots\\NewMeasures\\Test1_a.txt", "w");

	ns = 11; // number of spins
	// 1<<ns equivalent to 2^{ns}
	n = 1 << ns; //matrix dimension
	lmid = (ns - 1) / 2;
	d = 1;//site distance from center
	double B = 0.5; //magnetic field in z direction
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double J = 0.0;
	double Jz = 10.0; //NN interaction in z
	double alpha = 20.0;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double sumV, sign, sum;

	ham = ARR_D(n * n); //hamiltonian
	eval = ARR_D(n); dens = ARR_D(n); probSpinFlip = ARR_D(n);
	psi = ARR_DC(n); psit = ARR_DC(n); ov = ARR_DC(n);
	lwork = 34 * n;
	work = ARR_D(lwork);
	list1 = calloc(ns, sizeof(int));
	list2 = calloc(ns, sizeof(int));

	//printf("hello\n");

	for (i = 0; i < n * n; i++) ham[i] = 0.0;

	for (i = 0; i < n; i++) {
		spins(i, list1); // list of spins for basis state i
		spins(i, list2);
		for (k = 0; k < ns; k++) {
			list2[k] = 1 - list2[k]; //flip kth spin
			j = ind(list2);
			ham[i * n + j] = B;
			list2[k] = 1 - list2[k];
		}
		for (k = 0; k < ns - 1; k++) {
			list2[k] = 1 - list2[k]; //flip kth spin
			list2[k + 1] = 1 - list2[k + 1]; //flip k+1th spin
			j = ind(list2);
			ham[i * n + j] = Jz;
			list2[k] = 1 - list2[k];
			list2[k + 1] = 1 - list2[k + 1];
		}
		sumV = 0.0;
		//printf("%s\n", "-------k,kp-------");
		//for (k = 0; k < ns; k++) printf("%d\n", list1[k]);
		//printf("%s\n", "---------");
		for (k = 0; k < ns; k++) for (kp = 0; kp < k; kp++) {
			//if k-kp==1
			int dist = 1; //NN only
			if ((k - kp == 1) || abs(k-kp) == ns-1)//////////////////////////////////////////////////// Open boundary condition: NN only
			{
				//printf("%s%d   %s%d\n", "k=", k, "kp=", kp);
				if (list1[k] == list1[kp])
				{
					sign = 1;
					//printf("%f\n", sign);
				}
				else
				{
					sign = -1;
					//printf("%f\n", sign);
				}
				//if (list1[k] == list1[kp] && k-kp==1) sign = 1; else sign = -1;
				sumV += sign * J / pow(dist, alpha);
				//printf("%s%d  \n", "k-kp=", k - kp);
			}
		}
		//printf("\n");
		ham[i * n + i] = sumV;
	}
	
	//print ham : matches mathematica code
	printf("%s\n", "-----------------HAMILTONIAN H------------------");
	for (int i = 0; i < n * n; i++)
	{
		printf("%.2f   ", ham[i]);
		if ((i + 1) % n == 0) { printf("\n"); }
	}
	
	//'V':  Compute eigenvalues and eigenvectors
	//'U':  Upper triangle of ham is stored
	//On entry, ham - the symmetric matrix. Since 'U' the leading N - by - N upper triangular part of A contains the upper triangular part of the matrix
	// On exit, if 'V', then if INFO = 0, ham contains the orthonormal eigenvectors of the matrix A
	//e val- Eigenvalue (the eigenvalues in ascending order)
		
	dsyev_("V", "U", &n, ham, &n, eval, work, &lwork, &info); //diagonalize ham, on output contains eigenstates
	//printf("info %d %f\n",info,work[0]/n);
	//for (i = 0; i < n; i++) eval[i] = 1;

	//eigenvalue matches Mathematica
	printf("%s\n", "-----------------EIGENVALUES En------------------");
	for (i=0;i<n;i++) printf("%lf\n",eval[i]);


	for (i = 0; i < n; i++) psi[i].r = psi[i].i = 0.;
	for (i = 0; i < ns; i++) list1[i] = 1;
	//list1[lmid] = 0;	
	//list1[0] = 0; // perturbation at site 1
	psi[ind(list1)].r = 1;    // initialize psi, central spin left
	//printf("%d\n", ind(list1));
	//printf("%d\n", (ns - 1) / 2);
	
	//psi(0) matches Mathematica
	printf("%s\n", "-----------------PSI(0) |psi(0)>------------------");
	for(i=0;i<n;i++) printf("%0.1f %0.1f\n",psi[i].r,psi[i].i);
	
	//print eigenvectors : matches mathematica code
	printf("%s\n", "-----------------EIGENVECTORS |n>------------------");
	for (int i = 0; i < n * n; i++)
	{
		printf("%f   ", ham[i]);
		if ((i + 1) % n == 0) { printf("\n"); }
	}

	printf("%s\n", "-----------------OVERLAP <n|psi(0)>------------------");
	for (i = 0; i < n; i++) {
		ov[i].r = ov[i].i = 0.0; // overlap of initial state with i-th eigenstate
		for (j = 0; j < n; j++) ov[i].r += psi[j].r * ham[i * n + j]; //here psi and eigenstates are real
		
		printf("%d %lf\n",i,ov[i].r);
	}


	printf("%s\n", "-----------------<sigma_x>------------------");
	//-----------------------------------------------------------------------------------------------------------TIME LOOP-------------------------------------------------
	for (t=0;t<10;t+=0.01) { // start time loop
	//for (t = 0; t < 3e7; t += 3e5) { // start time loop
	//for (t = 0; t < 5; t += 0.05) { // start time loop 
	//for (t = 0; t < 1e9; t += 1e7) { // start time loop 
		for (j = 0; j < n; j++) {  //sum over basis states
			psit[j].r = psit[j].i = 0.0;
			for (i = 0; i < n; i++) { //sum over eigenstates
				psit[j].r += ov[i].r * cos(eval[i] * t) * ham[i * n + j];
				psit[j].i -= ov[i].r * sin(eval[i] * t) * ham[i * n + j];
			}
		}
		//printf("%s\n", "-----------------|psi(t).r>------------------");
		//for (j = 0; j < n; j++) printf("%d %lf\n", j, psit[j].r);

		//printf("%s\n", "-----------------|psi(t).i>------------------");
		//for (j = 0; j < n; j++) printf("%d %lf\n", j, psit[j].i);
		//sum=0.0;
		//for(j=0;j<n;j++) sum+=pow(psit[j].r,2)+pow(psit[j].i,2);
		//printf("%lf %lf\n",t,sum);

		// now let's calculate spin densities at time t
		for (k = 0; k < ns; k++)
		{
			//initialize densities to 0 //also initialize probability of spin flip to 0
			dens[k] = 0.; 
			probSpinFlip[k] = 0.; 
		}
		
		for (i = 0; i < n; i++) { //loop over basis states
			pr = pow(psit[i].r, 2) + pow(psit[i].i, 2); // prob to be in that basis state
			//printf("%s\n", "-----------------<psi(t)|psi(t)>------------------");
			//printf("%f\n", pr);
			//printf("%s\n", "-----------------list1------------------");
			//for (k = 0; k < ns; k++) printf("%d   ", list1[k]);
			spins(i, list1); // get list of spins
			//printf("\n");
			//printf("%s\n", "-----------------list1 after spins()------------------");
			//for (k = 0; k < ns; k++) printf("%d   ", list1[k]);
			//printf("\n");

			//add probability of state occuring if 
			for (k = 0; k < ns; k++) {
				if (list1[k]) 
				{
					//printf("%s","true");
					dens[k] += pr; 
					probSpinFlip[k] += pr;
				}
				else
				{
					//printf("%s", "false");
					dens[k] -= pr;
				}
				
			}

			/*//Spin is initially 1 but flips to 0
			//if Lmid+d flips to 0, then add pr to probSpinFlip, else subtract
			if (list1[lmid + d] == 0)
			{
				probSpinFlip[k] += pr;
			}
			else
			{
				probSpinFlip[k] -= pr;
			}*/
		}
		/*
		printf("%1.0e", t);
		for (k = 0; k < ns; k++) printf(" %lf", dens[k]);
		printf("\n");
		*/

		printf("%1.0e", t);
		for (k = 0; k < ns; k++) printf(" %lf", probSpinFlip[k]);
		printf("\n");

		fprintf(fh_output, "%1.0e", t);
		for (k = 0; k < ns; k++)
			fprintf(fh_output, " %f", dens[k]);
		fprintf(fh_output, "\n");

		fprintf(fh_output2, "%1.0e", t);
		for (k = 0; k < ns; k++)
			fprintf(fh_output2, " %f", probSpinFlip[k]);
		fprintf(fh_output2, "\n");
		
	}
	
	/*printf("%s\n", "-----------------spins------------------");
	    for (i=0;i<n;i++) {
			printf("%d: ",i);
			spins(i,list1);
			for (j=0;j<ns;j++) printf("%d ",list1[j]);
			printf("\n");
			printf("%d\n",ind(list1));
		}*/
	


	//   time(&n); 
	

	//close the file
	fclose(fh_output);
	fclose(fh_output2);
	return 0;
}

//converts binary to decimal
int ind(int* list)   // find basis element corresponding to spin chain "list"
{
	int i, j = 0;
	for (i = 0; i < ns; i++)
	{
		j += list[i] * pow(2, i);
		//printf("%d\n", j);
	}
	return j;
}

//all possible configurations of the spin chains
void spins(int index, int* s) // convert spin chain "s" to integer basis element index

{
	int i;
	for (i = 0; i < ns; i++)
		s[i] = (1 << i & index) > 0;
	return;
}
