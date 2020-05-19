// * 2020/04/20 Programed by Sota Goto
#include <iostream>
#include <fstream>
#include <omp.h>
#include <vector>
#include <math.h>
#include <algorithm> // for sort_function, max_element_function
#include <random>
#include <iomanip>
#include <string>
#include <time.h> // calculation time
#include <cstdio>

#include "global.h" // define arrays and physical const
#include "input_file_reader.cpp" // input mesh data
#include "geometry.cpp" // assign arrays and construct geometry info
#include "display_array.cpp" // for confirming arrays
#include "Calculate_matrix.cpp" // calculate matrices and vectors
#include "set_bc.cpp" // set boundary condtions
#include "visualize.cpp" // visualize the obtained results
#include "EigenValueProblem.cpp" // solve an eigenvalue problem
#include "utilities.cpp" // solve equations

// The program to solve Mu'+Ku=F. ' is a time-differential operator.

int main(void){
    clock_t start=clock();
    int i;
    printf("節点数 : %d\n要素数 : %d\n", n_p, n_e);

    Assign_array(); // define arrays
    
    input_file_reader(); // input node and element data
    
    InitialAdjacency(); // construct geometry info
    
    Calculate_matrix(); // construct matrices

    Calculate_Covmatrix(); // construct covariance matrix

    Power_method(); // solve an eigenvalue problem of covariance matrix

    Normalize_Eigenfunction();
    
    // set_bc(); // setting boundary condition
    
    // Jacobi_static(); // solve the linear equation by Jacobi scheme
    
    // Dynamics(); // dynamical analysis
    
    //Calculation time
    clock_t end=clock();
    const double time = static_cast <double> (end-start)/CLOCKS_PER_SEC*1000;
    printf("Calculation time : %.3f [ms]\n",time);
    return 0;
}
