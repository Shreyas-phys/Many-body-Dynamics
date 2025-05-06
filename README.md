Here is a detailed documentation of the `ShieldRK-Original.c` file in markdown format. The documentation includes explanations of the code and relevant code snippets.

---

# Documentation for `ShieldRK-Original.c`

This file is part of the `Quantum-Dynamics` repository and is written in C. The code implements computations related to quantum shielding and transport phenomena, using concepts like spin chains, Hamiltonians, and Runge-Kutta methods for solving differential equations.

## Table of Contents
1. [File Overview](#file-overview)
2. [Key Functions and Macros](#key-functions-and-macros)
3. [Global Variables](#global-variables)
4. [Main Function](#main-function)
5. [Helper Functions](#helper-functions)
6. [Runge-Kutta Implementation](#runge-kutta-implementation)
7. [Sparse and Dense Matrix Operations](#sparse-and-dense-matrix-operations)

---

## File Overview

This file performs simulations and analysis of spin chains using a combination of diagonalization and Runge-Kutta numerical methods. The main goal is to model and study quantum shielding effects in spin systems.

---

## Key Functions and Macros

### Defined Macros
The following macros are defined for memory allocation and mathematical convenience:

```c
#define ARR_D(i) (double *) calloc(i, sizeof(double))
#define ARR_DC(i) (doublecomplex *) calloc(i, sizeof(doublecomplex))
#define ARR_I(i) (int *) calloc(i, sizeof(int))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
```

- `ARR_D(i)`: Allocates an array of `double` with size `i`.
- `ARR_DC(i)`: Allocates an array of `doublecomplex` with size `i`.
- `ARR_I(i)`: Allocates an array of `int` with size `i`.
- `MIN(x, y)`: Returns the smaller of `x` and `y`.

### Function Declarations
The file declares several helper functions used throughout the simulation:

```c
int ind(int* list);
void spins(int index, int* s);
void analyzetime(double* ts, double* ds, int nt);
void Hpsi(double* ham, doublecomplex* psi, doublecomplex* psi1);
void Hpsi2(doublecomplex* psi, doublecomplex* psi1);
void rk(double* ham, doublecomplex* psi, doublecomplex* psi1, double dt);
void multvec(double fact, doublecomplex* v, doublecomplex* v1);
void multveci(double fact, doublecomplex* v, doublecomplex* v1);
void addmultvec(doublecomplex* v0, double fact, doublecomplex* v, doublecomplex* v1);
```

#### Key Functions:
- `ind`: Converts a spin chain to an integer basis state index.
- `spins`: Converts an integer index back to a spin chain.
- `analyzetime`: Analyzes survival probabilities over time.
- `Hpsi` and `Hpsi2`: Perform dense and sparse matrix-vector multiplications for Hamiltonians.
- `rk`: Implements the Runge-Kutta 4th order method for time evolution.

---

## Global Variables

Several global variables are declared to store simulation parameters, Hamiltonian data, and intermediate results:

```c
int ns; // Number of spins
integer n; // Matrix dimension (2^ns)
doublecomplex* t1, * t2, * t3, * k1, * k2, * k3, * k4; // Temporary vectors for RK4
int nham; // Number of non-zero Hamiltonian entries
double* ham2; // Array storing non-zero Hamiltonian entries
int* hind1, * hind2; // Indices for non-zero Hamiltonian elements
```

These variables are used across different functions for memory efficiency and reusability.

---

## Main Function

The `main` function initializes the simulation parameters, allocates memory, constructs the Hamiltonian, and performs time evolution.

### Key Steps in `main`

#### 1. Initialization
The program initializes global variables and simulation parameters, such as the number of spins (`ns`), the magnetic field (`B`), and interaction strengths (`J`, `Jlong`).

```c
ns = 10;  // Number of spins
n = 1 << ns;  // Matrix dimension, 2^ns
double B = 0.5;  // Magnetic field in z direction
double J = 2.0;  // Nearest-neighbor interaction in x direction
```

#### 2. Hamiltonian Construction
The Hamiltonian is constructed by iterating over all spin configurations. Sparse matrix representation is used for efficiency.

```c
for (i = 0; i < n; i++) {
    spins(i, list1); // Get spin configuration for state i
    for (k = 0; k < ns - 1; k++) {
        list2[k] = 1 - list2[k]; // Flip k-th spin
        list2[k + 1] = 1 - list2[k + 1]; // Flip (k+1)-th spin
        j = ind(list2);
        ham2[nham] = J; hind1[nham] = i; hind2[nham] = j; nham++;
    }
}
```

#### 3. Time Evolution
Time evolution is performed using the Runge-Kutta 4th order method (`rk` function).

```c
for (t = 0; t < 30.001; t += dt) {
    rk(ham2, psit, psitdt, dt); // Perform RK4 step
    for (i = 0; i < n; i++) psit[i] = psitdt[i];
}
```

---

## Helper Functions

### 1. `ind` and `spins`
These functions convert between spin configurations and integer indices.

```c
int ind(int* list) {
    int i, j = 0;
    for (i = 0; i < ns; i++)
        j += list[i] * pow(2, i);
    return j;
}

void spins(int index, int* s) {
    for (int i = 0; i < ns; i++)
        s[i] = (1 << i & index) > 0;
}
```

### 2. `analyzetime`
Analyzes survival probabilities over time and calculates metrics like maximum probability (`dmax`) and average probability (`davg`).

```c
void analyzetime(double* ts, double* ds, int nt) {
    double dmax = 0, davg = 0;
    for (int it = 0; it < nt; it++) {
        if (ds[it] > dmax) dmax = ds[it];
        davg += ds[it];
    }
    davg /= nt;
    printf("Max probability: %lf, Avg probability: %lf\n", dmax, davg);
}
```

---

## Runge-Kutta Implementation

The `rk` function implements the 4th-order Runge-Kutta method for time-dependent evolution.

```c
void rk(double* ham, doublecomplex* psi, doublecomplex* psi1, double dt) {
    Hpsi2(psi, psi1); multveci(dt, psi1, k1);
    addmultvec(psi, 0.5, k1, t1);
    Hpsi2(t1, psi1); multveci(dt, psi1, k2);
    addmultvec(psi, 0.5, k2, t2);
    Hpsi2(t2, psi1); multveci(dt, psi1, k3);
    addmultvec(psi, 1.0, k3, t3);
    Hpsi2(t3, psi1); multveci(dt, psi1, k4);
    for (int i = 0; i < n; i++) {
        psi1[i].r = psi[i].r + (k1[i].r + 2 * k2[i].r + 2 * k3[i].r + k4[i].r) / 6.0;
        psi1[i].i = psi[i].i + (k1[i].i + 2 * k2[i].i + 2 * k3[i].i + k4[i].i) / 6.0;
    }
}
```

This function uses temporary vectors (`t1`, `t2`, `k1`, etc.) to calculate intermediate steps.

---

## Sparse and Dense Matrix Operations

The file provides two functions for matrix-vector multiplication:
- `Hpsi` for dense matrices.
- `Hpsi2` for sparse matrices.

```c
void Hpsi(double* ham, doublecomplex* psi, doublecomplex* psi1) {
    for (int i = 0; i < n; i++) {
        psi1[i].r = psi1[i].i = 0.0;
        for (int j = 0; j < n; j++) {
            psi1[i].r += ham[i * n + j] * psi[j].r;
            psi1[i].i += ham[i * n + j] * psi[j].i;
        }
    }
}

void Hpsi2(doublecomplex* psi, doublecomplex* psi1) {
    for (int i = 0; i < n; i++)
        psi1[i].r = psi1[i].i = 0.0;
    for (int j = 0; j < nham; j++) {
        psi1[hind1[j]].r += ham2[j] * psi[hind2[j]].r;
        psi1[hind1[j]].i += ham2[j] * psi[hind2[j]].i;
    }
}
```

`Hpsi2` is optimized for sparse matrices by using precomputed non-zero elements (`ham2`).

---

## Initial States

This is how we set up our initial states:  
(*IPH* – Inside manifold with Polar perturbation from the Horizontal plane).  
Here, $\psi_1$ and $\psi_2$ are two states from band 1.

Set angles as:
- $\theta_a = \frac{\pi}{2}$
- $\theta_b = \frac{\pi}{2} + \epsilon$
- $\phi_a = \phi_b = 0$

The states are defined as:
$$
\psi_a = \frac{1}{\sqrt{2}} \psi_1 + \frac{1}{\sqrt{2}} \psi_2
$$

and

$$
\psi_b = \cos\left(\frac{\pi}{4} + \epsilon\right) \psi_1 + \sin\left(\frac{\pi}{4} + \epsilon\right) \psi_2
$$

band 1:

```c
// set up list1 corresponding to basis state psi1
for (i = 0; i < ns; i++) list1[i] = 1;

int perturbedSite = 1; //perturbation on site 2

//psi1
list1[perturbedSite] = 0; //perturbation on site 2 
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
```

The following below is the code for when the psi1 and psi2 are in the second band
```c
// set up list1 corresponding to basis state psi1
for (i = 0; i < ns; i++) list1[i] = 1;

int perturbedSite = 1; //perturbation on site 2

//psi1
list1[perturbedSite] = list1[3] = 0; //perturbation on site 2 //not a full spin flip
initb1 = cutind(list1);
//real part of psi
psi[initb1].r = cos(thetaa / 2);

//psi2
// set up list2 corresponding to basis state psi2
list1[perturbedSite] = list1[3] = 1; // reverting site 2 perturbation
//list1[nsmax-2] = 0; // perturbation on site ns-1
//list1[varyVariable-1] = 0; // perturbation on site varyVariable
list1[2] = list1[4] = 0; // perturbation on site 4
initb2 = cutind(list1);
psi[initb2].r += cos(phia) * sin(thetaa / 2);
psi[initb2].i = sin(phia) * sin(thetaa / 2);

//psi1 + psi2 (perturbation at sites 2 and ns-1)
for (i = 0; i < ns; i++) list1[i] = 1;
list1[perturbedSite] = list1[3] = 0; //perturbation on site 2
//list1[nsmax-2] = 0; // perturbation on site ns-1
//list1[varyVariable - 1] = 0; // perturbation on site varyVariable
list1[2] = list1[4] = 0; // perturbation on site 4*
initb1 = cutind(list1);
```


## Setting W to zero:
```c
for (k = 0; k < ns; k++) for (kp = 0; kp < k; kp++) {
	dist = min(k - kp, ns + kp - k); //dist=1 for NN
	//to set W=0, skip contributions where both k and kp belong to D
	if (list1[k] == 0 && list1[kp] == 0)// sign = 3 or 4 or -3 or -4 to set W to zero continue;  // Set W to zero
	if (list1[k] == list1[kp]) sign = 1; else sign = -1;
	sumVlong += sign * Jlong / pow(dist, alphalong);
	//}
}
ham2_effW[nham] = sumV + sumVlong;

```
This documentation provides a comprehensive breakdown of the `ShieldRK-Original.c` file, explaining its purpose, structure, and key components with code snippets.
