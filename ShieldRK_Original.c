#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include <time.h>

#define ARR_D(i) (double *) calloc(i, sizeof(double))
#define ARR_DC(i) (doublecomplex *) calloc(i, sizeof(doublecomplex))
#define ARR_I(i) (int *) calloc(i, sizeof(int))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))



int ind(int* list);
void spins(int index, int* s);
void analyzetime(double* ts, double* ds, int nt);
void Hpsi(double* ham, doublecomplex* psi, doublecomplex* psi1);
void Hpsi2(doublecomplex* psi, doublecomplex* psi1);
void rk(double* ham, doublecomplex* psi, doublecomplex* psi1, double dt);
void multvec(double fact, doublecomplex* v, doublecomplex* v1);
void multveci(double fact, doublecomplex* v, doublecomplex* v1);
void addmultvec(doublecomplex* v0, double fact, doublecomplex* v, doublecomplex* v1);

int ns;
integer n;

doublecomplex* t1, * t2, * t3, * k1, * k2, * k3, * k4;

double pi_value = 3.14159265358979323846;
int nham;
double* ham2;
int* hind1, * hind2;

int main()
{
	srand(time(NULL));
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



	double* ham, * eval, * work, t, * dens, * probSpinFlip0, * probSpinFlip1, pr, rem, e0, mx, * mxsave, * crossCorrelation, * Correlation2;
	double* rho11a, * rho12ra, * rho12ia;
	double* rho11na, * rho12rna, * rho12ina;
	double d0, d1, d2, d3, p, maxd, maxd2, * ts, * ds, * d2s, * ds0, * ds2, savenorm;
	doublecomplex* psi, * psit, * psitdt, * ov;
	integer lwork, info;


	int i, j, k, kp, it, initb1, initb2; //removed initb
	int* list1, * list2;
	int dist;
	int centerspin;
	int varyVariable;
	FILE* data, * survf, * fdens, * frhoa, * frhona, * fCorrelationFunc, * fCorrelation2;

	int dyn = 1; //0 means use diagonalization, 1 means use runge kutta
	int flip = 0;
	int per = 1;	//////

	FILE* fh_output;
	FILE* fh_output2;

	//mentions base file location
	char* base = "D:\\OneDrive - Tulane University\\RESEARCH\\Quantum transport\\Code\\Plots\\B0LFIMwithLR\\N15B0LFIMEmergentHfullb1lc";
	//char* base = "TestIPHJlongVariedN13Saturation";//// code to send to cypress
	FILE* energy = fopen("D:\\OneDrive - Tulane University\\RESEARCH\\Quantum transport\\Code\\Plots\\B0LFIMwithLR\\Energy.txt", "w");
	//FILE* energy = fopen("Energy.txt", "w");

	//Writing to file to another directory
	for (varyVariable = 10; varyVariable < 11; varyVariable++)
	{



		char filenameNA[256];
		char filenameA[256];
		char filenameFrhoNA[256];
		char filenameFrhoA[256];
		char filenameCorrelationFunction[256];
		char filenameCorrelation2[256];

		// dynamically changes filename while in varyVariable loop
		sprintf(filenameNA, "%s%d%s.txt", base, varyVariable, "na");
		sprintf(filenameA, "%s%d%s.txt", base, varyVariable, "a");
		sprintf(filenameFrhoA, "%s%d%s.txt", base, varyVariable, "frhoA");
		sprintf(filenameFrhoNA, "%s%d%s.txt", base, varyVariable, "frhoNA");
		sprintf(filenameCorrelationFunction, "%s%d%s.txt", base, varyVariable, "CorrelationFunction");
		sprintf(filenameCorrelation2, "%s%d%s.txt", base, varyVariable, "Correlation2");

		fh_output = fopen(filenameNA, "w");
		fh_output2 = fopen(filenameA, "w");
		frhoa = fopen(filenameFrhoA, "w");
		frhona = fopen(filenameFrhoNA, "w");
		fCorrelationFunc = fopen(filenameCorrelationFunction, "w");
		fCorrelation2 = fopen(filenameCorrelation2, "w");
		fdens = fopen("D:\\OneDrive - Tulane University\\RESEARCH\\Quantum transport\\Code\\Plots\\Manifold\\fdens.txt", "w");
		//fdens = fopen("fdens.txt", "w");




		//single deviation: centerspin = 0; // psia 
		//ferromagnetic: centerspin = 1 //psib
		for (centerspin = 0; centerspin < 2; centerspin++) //Automating the subtraction process
		{
			int nsmax = 15;
			//int nsmax = varyVariable;
			ns = nsmax; // number of spins
			n = 1 << nsmax; //matrix dimension, left shifting: a<<b = a*2^b. In this case 1<<ns = 2^ns


			////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			double B = 0; //magnetic field in z direction
			//double B = (double)varyVariable/2;
			double W = 0.0; // disorder strength
			double J = 0; //NN interaction in x
			//double J = (double)varyVariable / 10;
			//double Jlong = 0; //long range in x direction
			double Jlong = (double)varyVariable / 30;
			//double Jz = (double)varyVariable / 10; //NN interaction in z
			double Jz = 1;
			double alpha = 20.0; //short range interaction in x
			double alphalong = 0.0; //long range interaction in x
			//double alphalong = (double)varyVariable/10;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//superposition state parametrization
			//double epsilon = (double)varyVariable / 100;
			double epsilon = 0.01;
			double thetaa = pi_value / 2;
			double thetab = pi_value / 2 + epsilon;
			double phia = 0;
			double phib = 0;

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			double sumV, sumVlong, sign, sum, dt;
			////h = B+W = Magnetic field + random disorder
			double Hexpect, Hexpectnew, h[100];



			srand(time(NULL));
			//generates sequence of random numbers and stores them in the array `h`, where each number is calculated as `B + W` times a random value within the range [-0.5, 0.5].
			for (k = 0; k < ns; k++)
			{
				h[k] = B + W * (rand() / RAND_MAX - 0.5);
				//printf("h[k] = %lf", h[k]);
			}



			for (ns = nsmax; ns <= nsmax; ns += 2) {  //loop over number of spins

				n = 1 << ns; //matrix dimension
				nham = 2 * ns * n;   // number of nonzero entries in Hamiltonian

				//if (alpha<1) J=J*pow(ns,alpha-1); // asymptotic Kac prescription
				sumV = 0; // rescale J using exact Kac prescription 
				for (k = 0; k < ns; k++) for (kp = 0; kp < k; kp++)
					sumV += 1 / pow(k - kp, alpha);
				//J = J * ns / sumV;
				printf("Kac J %lf %lf \n", J, J * pow(ns, 1 - alpha));


				sumVlong = 0; // rescale Jlong using exact Kac prescription 
				for (k = 0; k < ns; k++) for (kp = 0; kp < k; kp++)
					sumVlong += 1 / pow(k - kp, alphalong);
				Jlong = Jlong * ns / sumVlong;
				//Jlong = 1;
				//check
				double sumVlongcheck = sumVlong;
				printf("Kac Jlong %lf %lf %lf \n", sumVlongcheck, Jlong, Jlong * pow(ns, 1 - alphalong));


				/*
				if (!dyn) {
					ham = ARR_D(n * n); //full hamiltonian matrix , not needed if dynamics is done with ham2 (nonzero entries only)
					eval = ARR_D(n);
					lwork = 34 * n;
					work = ARR_D(lwork);
					ov = ARR_DC(n);
					for (i = 0; i < n * n; i++) ham[i] = 0.0;
				}
				*/


				if (dyn) {
					t1 = ARR_DC(n); t2 = ARR_DC(n); t3 = ARR_DC(n);
					k1 = ARR_DC(n); k2 = ARR_DC(n); k3 = ARR_DC(n); k4 = ARR_DC(n);
				}

				dens = ARR_D(n); mxsave = ARR_D(n); probSpinFlip0 = ARR_D(n); probSpinFlip1 = ARR_D(n); crossCorrelation = ARR_D(n); Correlation2 = ARR_D(n);
				rho11a = ARR_D(n); rho12ra = ARR_D(n); rho12ia = ARR_D(n);
				rho11na = ARR_D(n); rho12rna = ARR_D(n); rho12ina = ARR_D(n);


				nham = 2 * ns * n + n;   // number of nonzero entries in Hamiltonian
				ham2 = ARR_D(nham); hind1 = ARR_I(nham); hind2 = ARR_I(nham);  //ham2 stores nonzero entrries only 
				list1 = ARR_I(ns);
				list2 = ARR_I(ns);
				psi = ARR_DC(n); psit = ARR_DC(n); psitdt = ARR_DC(n);
				ts = ARR_D(20000); ds = ARR_D(20000); d2s = ARR_D(20000);
				ds0 = ARR_D(20000); ds2 = ARR_D(20000);
				nham = 0;


				//for (i = 0; i < ns; i++) list1[i] = 1;
				//if (flip) list1[(ns-1)/2]=0;
				//if (flip == 1) list1[0] = 0;   // flip leftmost spin
				//if (flip == -1) list1[(ns - 1) / 2] = 0; //flip center spin?
				//initb = ind(list1);

				//for (i = 0; i < n; i++) psi[i].r = psi[i].i = 0.;
				//for (i = 0; i < ns; i++) list1[i] = 1;
				//list1[2] = centerspin;
				////list1[(ns - 1)/2] = centerspin;/////////////////////////////////////////////////////////////////////////////////////////////////////center spin
				////printf("centerspin in code %d\n", centerspin);
				//initb = ind(list1);
				//psi[initb].r = 1;    // initialize psi, central spin left
				////for(i=0;i<n;i++) printf("%lf %lf\n",psi[i].r,psi[i].i);

				//change variables from a to b on the second run
				if (centerspin == 1)
				{
					thetaa = thetab;
					phia = phib;
				}

				for (i = 0; i < n; i++) psi[i].r = psi[i].i = 0.;
				// set up list1 corresponding to basis state psi1
				for (i = 0; i < ns; i++) list1[i] = 1;

				int perturbedSite = 1; //perturbation on site 2

				//psi1
				list1[perturbedSite] = 0; //perturbation on site 2 //not a full spin flip
				initb1 = ind(list1);
				//real part of psi
				psi[initb1].r = cos(thetaa / 2);

				//psi2
				// set up list2 corresponding to basis state psi2
				list1[perturbedSite] = 1; // reverting site 2 perturbation
				//list1[nsmax-2] = 0; // perturbation on site ns-1
				//list1[varyVariable-1] = 0; // perturbation on site varyVariable
				list1[3] = 0; // perturbation on site 4
				initb2 = ind(list1);
				psi[initb2].r += cos(phia) * sin(thetaa / 2);
				psi[initb2].i = sin(phia) * sin(thetaa / 2);

				//psi1 + psi2 (perturbation at sites 2 and ns-1)
				for (i = 0; i < ns; i++) list1[i] = 1;
				list1[perturbedSite] = 0; //perturbation on site 2
				//list1[nsmax-2] = 0; // perturbation on site ns-1
				//list1[varyVariable - 1] = 0; // perturbation on site varyVariable
				list1[3] = 0; // perturbation on site 4*
				initb1 = ind(list1);

				for (i = 0; i < n; i++) {
					spins(i, list1); // list of spins for basis state i
					spins(i, list2);
					for (k = 0; k < ns; k++) {
						list2[k] = 1 - list2[k]; //flip kth spin
						j = ind(list2);
						//if (!dyn) ham[i * n + j] = B;
						ham2[nham] = h[k]; hind1[nham] = i; hind2[nham] = j; nham++;
						list2[k] = 1 - list2[k];
					}


					//for (k = 0; k < ns; k++) {/////////////////////////////////////////////////////////////////////////////////////////////PBC
					for (k = 0; k < ns - 1; k++) {////////////////////////////////////////////////////////////////////////////////////////OBC

						//Calculate the index of the next spin with periodic boundary conditions
						//int next_spin_index = (k + 1) % (ns);/////////////////////////////////////////////////////////////////////////////PBC
						int next_spin_index = k + 1;////////////////////////////////////////////////////////////////////////////////////OBC

						// Flip the spins at indices k and next_spin_index
						list2[k] = 1 - list2[k];
						list2[next_spin_index] = 1 - list2[next_spin_index];

						// Calculate the indices for the Hamiltonian
						int j = ind(list2);

						// Update the Hamiltonian
						ham2[nham] = Jz;
						hind1[nham] = i;
						hind2[nham] = j;
						nham++;

						// Restore the original spin configuration
						list2[k] = 1 - list2[k];
						list2[next_spin_index] = 1 - list2[next_spin_index];
					}


					sumV = 0.0;
					for (k = 0; k < ns; k++) for (kp = 0; kp < k; kp++) {
						dist = k - kp;/////////////////////////////////////////////////////////////////////////////////////////////OBC
						//dist = MIN(k - kp, ns + kp - k); ///////////////////////////////////////////////////////////////////////////////////////PBC
						if (abs(k - kp) == 1)/////////////////////////////////////////////////////////////////////////////////////// |n-m|=1 NN coupling only, OBC
							//if ((k - kp) == 1 || abs(k - kp) == ns - 1)//////////////////////////////////////////////////////////////////////////// Boundary condition: NN only, PBC
						{
							//printf("%s %d %s %d\n","n=", k, "m=", kp);
							if (list1[k] == list1[kp]) sign = 1; else sign = -1;
							sumV += sign * J / pow(dist, alpha);
							//if (per && kp == 0 && k == ns - 1) sumV += sign * J;	////
						}

					}


					sumVlong = 0.0;
					for (k = 0; k < ns; k++) for (kp = 0; kp < k; kp++) {
						dist = k - kp;//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////OBC
						//dist = MIN(k - kp, ns + kp - k); /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////PBC
						//{
							//printf("%s %d %s %d\n","n=", k, "m=", kp);
						if (list1[k] == list1[kp]) sign = 1; else sign = -1;
						sumVlong += sign * Jlong / pow(dist, alphalong);
						//printf("sumVlong = %.5e\n", sumVlong);  // prints in scientific notation with 5 decimal places
						//}
					}
					//if (!dyn) ham[i * n + i] = sumV;
					if (i == initb1) printf("Q %lf\n", sumV);
					ham2[nham] = sumV + sumVlong; hind1[nham] = i; hind2[nham] = i; nham++;
				}

				//printf("nham %d\n", nham);  //nonzero entries in Hamiltonian
				//for(i=0;i<nham;i++) printf("ham %d %d %lf\n",hind1[i],hind2[i],ham2[i]);
				//return 0;

				//for (i=0;i<n;i++) for (j=0;j<n;j++) printf("ham %d %d %lf\n",i,j,ham[i*n+j]);

				/*for (i = 0; i < n; i++) psi[i].r = psi[i].i = 0.;
				for (i = 0; i < ns; i++) list1[i] = 1;
				list1[(ns - 1) / 2] = 0;
				initb = ind(list1);*/

				for (i = 0; i < n; i++) {  // calculate total x spin Mx for each basis state
					spins(i, list1);
					mx = 0;
					for (k = 0; k < ns; k++) mx += list1[k] * 2 - 1;
					mxsave[i] = mx;
				}

				if (0) for (i = 0; i < n; i++) for (j = 0; j < n; j++) //for testing purposes
					if (i != j && (mxsave[i] < ns - 4 || mxsave[j] < ns - 4)) ham[i * n + j] = 0;

				/*
				if (!dyn) { //diagonalize if needed
					dsyev_("V", "U", &n, ham, &n, eval, work, &lwork, &info); //diagonalize ham, on output contains eigenstates
					printf("info %d %lf\n", info, work[0] / n);
				}

				*/
				//for (i=0;i<n;i++) printf("%lf\n",eval[i]); //eigenvalues



				/*
				if (!dyn) for (i = 0; i < n; i++) {
					ov[i].r = ov[i].i = 0.0; // overlap of initial state with i-th eigenstate
					for (j = 0; j < n; j++) ov[i].r += psi[j].r * ham[i * n + j]; //here psi and eigenstates are real
					//printf("%d %lf\n",i,ov[i].r);
				}*/

				/*e0=ham[initb*n+initb]; //testing only
				for(i=0;i<n;i++) {
					 rem=ham[initb*n+i];
					 if (i!=initb && fabs(rem)>1e-11) {
						 printf("%d %lf %lf\n",i,rem,ham[i*n+i]-e0);
						 spins(i,list1);
						 for (k=0;k<ns;k++) printf("%d ",list1[k]*2-1); printf("\n");
					 }
				 }
				 return 0;*/

				Hpsi2(psi, psit);
				pr = 0.; for (i = 0; i < n; i++) pr += psi[i].r * psit[i].r + psi[i].i * psit[i].i;
				Hexpect = pr;
				//printf("H expect initial %lf\n", Hexpect); //expectation value of energy in initial state

				//survf = fopen("surv.out", "w");

				it = 0;
				dt = 0.01;  //time step
				for (i = 0; i < n; i++) psit[i] = psi[i];  // psi at time 0
				for (t = 0; t < 10.001; t += dt) { // start time loop//////////////////////////////////////////////////////////////////////////time

					/*
					if (dyn == 0) for (j = 0; j < n; j++) {  //sum over basis states if using diagonalization method
						psit[j].r = psit[j].i = 0.0;
						for (i = 0; i < n; i++) { //sum over eigenstates
							psit[j].r += ov[i].r * cos(eval[i] * t) * ham[i * n + j];
							psit[j].i -= ov[i].r * sin(eval[i] * t) * ham[i * n + j];
						}
					}
					/*
					*
					/*printf("t %lf psi ",t);
					for (i=0;i<n;i++) printf("%lf +i %lf ,",psit[i].r,psit[i].i);
					printf("\n");*/

					// check normalization
					pr = 0;
					for (i = 0; i < n; i++) pr += pow(psit[i].r, 2) + pow(psit[i].i, 2);
					//printf("norm %lf %lf %lg\n", t, pr, 1 - pr);
					savenorm = pr;

					//H|psi(t)>
					Hpsi2(psit, psitdt);
					pr = 0.; for (i = 0; i < n; i++) pr += psitdt[i].r * psit[i].r + psitdt[i].i * psit[i].i;
					Hexpectnew = pr;
					//printf("H expect %lf %lf %lf\n", t, Hexpect, Hexpectnew); //expectation value of energy at time t compared with initial

					// now let's calculate spin densities at time t
					if (fabs(t * 100 - round(t * 100)) < 0.001) {

						for (k = 0; k < ns; k++)
						{
							dens[k] = 0.; //initialize densities to 0

							if (centerspin == 0)
							{
								crossCorrelation[k] = 0;//initialize crosscorrelation term to 0
								Correlation2[k] = 0;
								probSpinFlip0[k] = 0.;
								rho11na[k] = rho12rna[k] = rho12ina[k] = 0.;
							}
							else
							{
								probSpinFlip1[k] = 0.;
								rho11a[k] = rho12ra[k] = rho12ia[k] = 0.;
							}
						}


						for (i = 0; i < n; i++) { //loop over basis states |b>, ex: |11>, |10>, |01>, |00>
							pr = pow(psit[i].r, 2) + pow(psit[i].i, 2); // prob to be in that basis state |<psi|b>|^2
							spins(i, list1); // get list of spins
							for (k = 0; k < ns; k++) { //iterates over spins in the list.
								//In Cij, fix i=perturbed site, if k==perturbed site
								if (list1[k] == list1[perturbedSite])
								{
									crossCorrelation[k] += pr;
								}
								else
								{
									crossCorrelation[k] -= pr;
								}

								//<psi|sigma_2|psi>
								if (k == perturbedSite)
								{
									Correlation2[k] += pr;
								}
								else
								{
									Correlation2[k] -= pr;
								}

								if (list1[k])// if spin is 1
								{

									dens[k] += pr; //dens doesn't have a/na suffix since I'm checking dens for one initial state at a time
									//rho11[k] += pr;
									list1[k] = 1 - list1[k];//pauli x flips 
									j = ind(list1);
									list1[k] = 1 - list1[k];//flipping it back
									//rho12r[k] += psit[i].r * psit[j].r + psit[i].i * psit[j].i;
									//rho12i[k] += psit[i].i * psit[j].r - psit[i].r * psit[j].i;
									if (centerspin == 0) //i.e spins are not aligned
									{

										dens[k] += pr;
										//probSpinFlip0[k] += pr;
										rho11na[k] += pr;
										rho12rna[k] += psit[i].r * psit[j].r + psit[i].i * psit[j].i;
										rho12ina[k] += psit[i].i * psit[j].r - psit[i].r * psit[j].i;
									}
									else //spins are aligned
									{
										dens[k] += pr;
										//probSpinFlip1[k] += pr;
										rho11a[k] += pr;
										rho12ra[k] += psit[i].r * psit[j].r + psit[i].i * psit[j].i;
										rho12ia[k] += psit[i].i * psit[j].r - psit[i].r * psit[j].i;
									}
								}
								else
								{
									////printf("%s", "false");
									////probSpinFlip[k] += pr;
									//if (centerspin == 0)
									//{
									//	probSpinFlip0[k] += pr;
									//	dens[k] -= pr;
									//	rho11na[k] += pr; // takes the same value that probSpinFlip0 does
									//	//rho11na[k] = probSpinFlip0[k];
									//	rho12rna[k] += psit[i].r * psit[j].r + psit[i].i * psit[j].i;
									//	rho12ina[k] += psit[i].i * psit[j].r - psit[i].r * psit[j].i;
									//	
									//}
									//else
									//{
									//	probSpinFlip1[k] += pr;
									//	dens[k] -= pr;
									//	rho11a[k] += pr; //rho11a = probSpinFlip1
									//	//rho11a[k] = probSpinFlip1[k];
									//	rho12ra[k] += psit[i].r * psit[j].r + psit[i].i * psit[j].i;
									//	rho12ia[k] += psit[i].i * psit[j].r - psit[i].r * psit[j].i;
									//	
									//}
									////dens[k] -= pr;
									////rho11[k] -= pr;
								}
							}
						}


					}



					if (1) {   //print out densities
						//printf("%lf", t);
						//for (k = 0; k < ns; k++) printf(" %lf", dens[k]);
						//printf("\n");
						//printf("Kac J %lf %lf \n", J, J * pow(ns, 1 - alpha));
						//printf("Kac Jlong %lf %lf %lf \n", sumVlongcheck, Jlong, Jlong * pow(ns, 1 - alphalong));
					}

					if (0 && it == 2351) {  //testing only
						for (mx = ns; mx >= -ns; mx -= 2) {
							d0 = 0;
							for (i = 0; i < n; i++) if (mxsave[i] == mx) d0 += pow(psit[i].r, 2) + pow(psit[i].i, 2);
							//printf("2351 mx %d %lf\n", mx, d0);
						}
					}

					//what happens here
					rem = pow(psit[initb1].r, 2) + pow(psit[initb1].i, 2); //prob to be in initial basis state
					//printf("rem %lf\n",rem);
					d0 = d1 = d2 = d3 = 0;
					for (i = 0; i < n; i++) {
						p = pow(psit[i].r, 2) + pow(psit[i].i, 2);
						if (mxsave[i] == ns) d0 += p;  //prob for all spins right
						else if (mxsave[i] == ns - 2) d1 += p;  //prob for one spin left (same as initial state)
						else if (mxsave[i] == ns - 4) d2 += p;  //prob for two spins left (same as initial state)
						else if (mxsave[i] == ns - 6) d3 += p;
					}


					d1 -= rem;
					//printf("prob %lf %lf %lf %lf %lf %lf %lf\n", t, 1 - rem, d0, d1, d2, d3, 1 - rem - d0 - d1 - d2 - d3);
					//fprintf(survf, "%lf %lf %lf %lf %lf %lf %lf %lf\n", t, 1 - rem, d1 + rem, d0, d1, d2, d3, 1 - rem - d0 - d1 - d2 - d3);
					ts[it] = t; ds[it] = 1 - rem; d2s[it] = 1 - (d1 + rem); it++; //save probabilities at time t
					ds0[it] = d0; ds2[it] = d2;
					//fprintf(survf,"%lf %lf\n",t,1-rem);

					if (dyn == 1) {  //do evolution using runge kutta
						rk(ham2, psit, psitdt, dt);
						for (i = 0; i < n; i++) psit[i] = psitdt[i];
					}


					//fprintf(fdens, "%1.2e", t);
					for (k = 0; k < ns; k++)
						fprintf(fdens, " %0.15f", dens[k]);
					fprintf(fdens, "\n");

					////print out densities, should be (1-probspin flip)/2
					//fprintf(fdens, "%lf ", t);
					//printf("%lf ", t);
					//for (k = 0; k < ns; k++) fprintf(fdens, " %20.14lg", dens[k]);
					//fprintf(fdens, "\n");

					//center deviation
					if (centerspin == 0)
					{
						printf("centerspin %d", centerspin);
						printf("\n");
						fprintf(fh_output, "%1.2e", t);
						fprintf(frhona, "%1.2e", t);
						fprintf(fCorrelationFunc, "%1.2e", t);
						fprintf(fCorrelation2, "%1.2e", t);
						printf("%1.2f ", t);
						for (k = 0; k < ns; k++) {
							fprintf(fh_output, " %0.15f", rho11na[k]);
							fprintf(fCorrelationFunc, " %0.15f", crossCorrelation[k]);
							fprintf(fCorrelation2, " %0.15f", Correlation2[k]);
							fprintf(frhona, " %0.15f %0.15f %0.15f", rho11na[k], rho12rna[k], rho12ina[k]);
							//printf(" %0.15f %0.15f %0.15f", rho11na[k], rho12rna[k], rho12ina[k]);
							//printf("%0.15f", probSpinFlip0[k]);
						}
						//printf("\n");
						fprintf(fh_output, "\n");
						fprintf(frhona, "\n");
						fprintf(fCorrelationFunc, "\n");
						fprintf(fCorrelation2, "\n");
						//printf("centerspin %d", centerspin);
						//printf("\n");
						//printf("\n");
						//printf("\n");
						printf("varyVariable %d", varyVariable);
					}

					//aligned
					if (centerspin == 1)
					{

						printf("centerspin %d", centerspin);
						printf("\n");
						fprintf(fh_output2, "%1.2e", t);
						fprintf(frhoa, "%1.2e", t);
						printf("%1.2e ", t);
						for (k = 0; k < ns; k++)
						{
							//printf("%s", "Code came here!!!!!!!!!!!!!!!!!");
							//printf("\n");
							fprintf(fh_output2, " %0.15f", rho11a[k]);
							fprintf(frhoa, " %0.15f %0.15f %0.15f", rho11a[k], rho12ra[k], rho12ia[k]);
							//printf(" %0.15f %0.15f %0.15f", rho11a[k], rho12ra[k], rho12ia[k]);
							//printf(" %0.15f", probSpinFlip1[k]);
						}
						//printf("\n");
						fprintf(fh_output2, "\n");
						fprintf(frhoa, "\n");

						//printf("\n");
						//printf("\n");
						printf("\n");
						printf("varyVariable %d", varyVariable);

						//printf("centerspin %d", centerspin);						
					}

					//printf("%s", "Code came here!!!!!!!!!!!!!!!!!");
				}

				data = fopen("shieldrks.check", "a");
				fprintf(data, "%d %lf %lf %lf\n", ns, dt, savenorm, Hexpectnew / Hexpect);
				fclose(data);

				//printf("analyze  ");
				//analyzetime(ts,ds,it);
				//printf("d2  ");
				analyzetime(ts, d2s, it);

			}
			fprintf(energy, "%0.15f\n", Hexpect);
			//printf("normalization to check for convergence %lf %e\n", t, 1 - savenorm);
			//printf("Kac Jlong %lf %lf %lf \n", sumVlongcheck, Jlong, Jlong * pow(ns, 1 - alphalong));
			//close the file
			// 
			//Free all dynamically allocated variables before starting the next loop
			/*if (dyn) {	free(t1); free(t2);	free(t3); free(k1); free(k2);	free(k3);	free(k4);}

			free(dens);	free(mxsave); free(probSpinFlip0); free(probSpinFlip1);
			free(rho11); free(rho12r); free(rho12i); free(ham2); free(hind1);
			free(hind2); free(list1); free(list2); free(psi);
			free(psit); free(psitdt); free(ts); free(ds); free(d2s); free(ds0);
			free(ds2);*/
		}//loop over hamiltonian variables


		fclose(fh_output);
		fclose(fdens);
		fclose(frhoa);
		fclose(frhona);
		fclose(fh_output2);
		fclose(fCorrelationFunc);
		fclose(fCorrelation2);
	}

	fclose(energy);
	return 0;
}


int ind(int* list)   // find basis element conrresponding to spin chain "list"
{
	int i, j = 0;
	for (i = 0; i < ns; i++)
		j += list[i] * pow(2, i);
	return j;
}

void spins(int index, int* s) // convert spin chain "s" to integer basis element index

{
	int i;
	for (i = 0; i < ns; i++)
		s[i] = (1 << i & index) > 0;
	return;
}

void analyzetime(double* ts, double* ds, int nt)  //analyze survival probabilities over time
{
	double dmax, tmax, t90, davg;
	int navg;
	int it;
	FILE* data;
	dmax = tmax = t90 = 0;
	for (it = 0; it < nt; it++)
		if (ds[it] > dmax) { dmax = ds[it]; tmax = ts[it]; }
	int found = 0;
	it = 0;
	while (!found) {
		if (ds[it] > 0.9 * dmax) { t90 = ts[it]; found = 1; }
		it++;
	}
	davg = 0; navg = 0;
	for (it = 0; it < nt; it++)
		if (ts[it] > 1.0) { davg += ds[it]; navg++; }
	//if(ts[it]>10.0) {davg+=ds[it]; navg++;}
	davg /= navg;
	double q1 = ds[1] / pow(ts[1], 2);
	double q2 = ds[2] / pow(ts[2], 2);
	printf("t90 %lf tmax %lf max prob %lf\n", t90, tmax, dmax);
	printf("davg %lf\n", davg);
	printf("quad %lf acc %lf\n", q1, q2 / q1);
	data = fopen("shieldrks.out", "a");
	fprintf(data, "%d %lf\n", ns, davg);
	fclose(data);
}


//dense matrix
void Hpsi(double* ham, doublecomplex* psi, doublecomplex* psi1) { // psi1 = H*psi
	for (int i = 0; i < n; i++) {
		psi1[i].r = psi1[i].i = 0.0;
		for (int j = 0; j < n; j++) {
			psi1[i].r += ham[i * n + j] * psi[j].r;
			psi1[i].i += ham[i * n + j] * psi[j].i;
		}
	}
}

//sparse matrix
void Hpsi2(doublecomplex* psi, doublecomplex* psi1) {  //psi1 = H*psi but using nonzero elements of H (stored in ham2)
	for (int i = 0; i < n; i++)
		psi1[i].r = psi1[i].i = 0.0;
	for (int j = 0; j < nham; j++) {
		psi1[hind1[j]].r += ham2[j] * psi[hind2[j]].r;
		psi1[hind1[j]].i += ham2[j] * psi[hind2[j]].i;
	}
}


void rk(double* ham, doublecomplex* psi, doublecomplex* psi1, double dt) {  //RK4
	//doublecomplex t1[n],t2[n],t3[n],k1[n],k2[n],k3[n],k4[n];
	Hpsi2(psi, psi1); multveci(dt, psi1, k1); addmultvec(psi, 0.5, k1, t1);

	/*printf("psi1 ");
for (int i=0;i<n;i++) printf("%lf +i %lf ,",psi1[i].r,psi1[i].i);
printf("\n");

			printf("k11 ");
for (int i=0;i<n;i++) printf("%lf +i %lf ,",k1[i].r,k1[i].i);
printf("\n");*/

	Hpsi2(t1, psi1); multveci(dt, psi1, k2); addmultvec(psi, 0.5, k2, t2);
	Hpsi2(t2, psi1); multveci(dt, psi1, k3); addmultvec(psi, 1.0, k3, t3);
	Hpsi2(t3, psi1); multveci(dt, psi1, k4);
	for (int i = 0; i < n; i++) psi1[i].r = psi[i].r + (k1[i].r + 2 * k2[i].r + 2 * k3[i].r + k4[i].r) / 6.0;
	for (int i = 0; i < n; i++) psi1[i].i = psi[i].i + (k1[i].i + 2 * k2[i].i + 2 * k3[i].i + k4[i].i) / 6.0;

	//addmultvec(psi,1.0,k1,psi1);

		/*		printf("psi1 ");
		for (int i=0;i<n;i++) printf("%lf +i %lf ,",psi1[i].r,psi1[i].i);
		printf("\n");*/
}

void multvec(double fact, doublecomplex* v, doublecomplex* v1) {   // v1 = fact*v
	for (int i = 0; i < n; i++) {
		v1[i].r = fact * v[i].r; v1[i].i = fact * v[i].i;
	}
}

void multveci(double fact, doublecomplex* v, doublecomplex* v1) {  // v1 = -i fact*v
	for (int i = 0; i < n; i++) {
		v1[i].r = fact * v[i].i; v1[i].i = -fact * v[i].r;
	}
}


void addmultvec(doublecomplex* v0, double fact, doublecomplex* v, doublecomplex* v1) {  // v1 = v0 + fact*v
	for (int i = 0; i < n; i++) {
		v1[i].r = v0[i].r + fact * v[i].r; v1[i].i = v0[i].i + fact * v[i].i;
	}
}
