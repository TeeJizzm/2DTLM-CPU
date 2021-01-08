// Includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>

// Definitions
#define M_PI 3.14276
#define c 299792458
#define mu0 M_PI*4e-7
#define eta0 c*mu0



// CPU function for source calculation
void stageSource(double* V1, double* V2, double* V3, double* V4, int x, int y, double E0, int NY) {
    /* Stage 1: Source */

    // Adapted to be 1D
    V1[x * NY + y] = V1[x * NY + y] + E0;
    V2[x * NY + y] = V2[x * NY + y] - E0;
    V3[x * NY + y] = V3[x * NY + y] - E0;
    V4[x * NY + y] = V4[x * NY + y] + E0;
    // Using 1 dimensional arrays is more obvious to work with when porting to GPU

} // end func

// CPU function
void stageScatter(double* V1, double* V2, double* V3, double* V4, int NX, int NY, double Z) {
    /* Stage 2: Scatter */
    // Variables 
    double I = 0, V = 0;

    // Parallelisable code

    // for int i = 0; i < NX*NY; i++
    for (int x = 0; x < NX; x++) {
        for (int y = 0; y < NY; y++) {
            I = (V1[(x * NY) + y] + V4[(x * NY) + y] - V2[(x * NY) + y] - V3[(x * NY) + y]) / (2);
            // factorized by 2 for unnecessary mathematics

            V = 2 * V1[x * NY + y] - I;         //port1
            V1[x * NY + y] = V - V1[x * NY + y];

            V = 2 * V2[x * NY + y] + I;         //port2
            V2[x * NY + y] = V - V2[x * NY + y];

            V = 2 * V3[x * NY + y] + I;         //port3
            V3[x * NY + y] = V - V3[x * NY + y];

            V = 2 * V4[x * NY + y] - I;         //port4
            V4[x * NY + y] = V - V4[x * NY + y];
        }
    }
} // end func



// CPU Function
void stageConnect(double* V1, double* V2, double* V3, double* V4, // Arrays
    int NX, int NY, // Array arguments
    double rXmin, double rXmax, double rYmin, double rYmax) { // Boundary conditions
/* Stage 3: Connect */
// Variables
    double tempV = 0;

    // Connect internals
    for (int x = 1; x < NX; x++) {
        for (int y = 0; y < NY; y++) {
            tempV = V2[x * NY + y];
            V2[x * NY + y] = V4[(x - 1) * NY + y];
            V4[(x - 1) * NY + y] = tempV;
        }
    }
    for (int x = 0; x < NX; x++) {
        for (int y = 1; y < NY; y++) {
            tempV = V1[x * NY + y];
            V1[x * NY + y] = V3[x * NY + y - 1];
            V3[x * NY + y - 1] = tempV;
        }
    }

    // Connect boundaries
    for (int x = 0; x < NX; x++) {
        V3[x * NY + NY - 1] = rYmax * V3[x * NY + NY - 1];
        V1[x * NY] = rYmin * V1[x * NY]; // V1[x * NY + 0] = rYmin * V1[x * NY + 0];
    }
    for (int y = 0; y < NY; y++) {
        V4[(NX - 1) * NY + y] = rXmax * V4[(NX - 1) * NY + y];
        V2[y] = rXmin * V2[y]; // V2[0 * NY + y] = rXmin * V2[0 * NY + y];
    }
} // end func



int main() {

    // Start timer
    std::clock_t start = std::clock();

    /* Variables */
    // Changable variables
    int NX = 200; // number of X
    int NY = 200; // number of Y
    int NT = 8192; // number of Times/Iterations
    double dl = 1;

    // Retrieval from GPU
    //*/
    double* V1 = new double[int(NX * NY)](); // new double[int(NX*NY)](); // Sets all values to 0 
    double* V2 = new double[int(NX * NY)]();
    double* V3 = new double[int(NX * NY)]();
    double* V4 = new double[int(NX * NY)]();
    //*/


    // Variables and coefficients
    // Scatter coefficient
    double Z = eta0 / sqrt(2.);

    // Boundary connect Coefficiants
    double rXmin = -1;
    double rXmax = -1;
    double rYmin = -1;
    double rYmax = -1;

    // impulse parameters
    double dt = dl / (sqrt(2.) * c);
    double width = 20 * dt * sqrt(2.);
    double delay = 100 * dt * sqrt(2.);

    // input position
    int Ein[] = { 10,10 };
    // output/reading position
    int Eout[] = { 15,15 };

    // file output
    std::ofstream output("CPU.csv");


    for (int n = 0; n < NT; n++) {
        // Variables dependant on n
        double E0 = (1 / sqrt(2.)) * exp(-(n * dt - delay) * (n * dt - delay) / (width * width));


        /* Stage 1: Source */
        stageSource(V1, V2, V3, V4, Ein[0], Ein[1], E0, NY);


        /* Stage 2: Scatter */
        stageScatter(V1, V2, V3, V4, NX, NY, Z);
        

        /* Stage 3: Connect */
        stageConnect(V1, V2, V3, V4, NX, NY, rXmin, rXmax, rYmin, rYmax);
      

        output << n * dt << "," << V2[Eout[0] * NY + Eout[1]] + V4[Eout[0] * NY + Eout[1]] << std::endl;
        if (n % 100 == 0)
            std::cout << n << std::endl;

    } // End of loop


    output.close();
    std::cout << "Done: " << ((std::clock() - start) / (double)CLOCKS_PER_SEC) << std::endl;
    std::cin.get();


} // end main

// EOF