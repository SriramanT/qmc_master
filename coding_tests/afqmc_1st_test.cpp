//
//  afqmc_1st_test.cpp
//  
//
//  Created by Francisco Brito on 09/03/2018.
//

#include <stdio.h>
#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <typeinfo>
#include <cmath>
#include <random>
#include <time.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#define N_SITES 8

using namespace Eigen;
using namespace std;

int main()
{
    int i;
    int j;
    int l;
    int step;
    int i_chosen;
    int l_chosen;
    
    int seed = 12345;                       // set a seed
    
    //random number to make decision
    double decisionMaker;
    
    mt19937 gen(seed);
    uniform_real_distribution<> dis(0.0, 1.0);
    
    int totalMCSteps = 30;
    int L = 10;                             // number of imaginary time subintervals
    double beta = 1.;                       // imaginary time interval or equivalent maximum temperature of the (d+1) classical system
    double dt = beta/L;                     // time subinterval width
    double t = 1.;                          // hopping parameters
    double U = 1.;                          // interaction energy
    int N = 5;                              // number of sites
    double nu = acosh( exp( U * dt / 2 ) );
    
    
    // Initialize the HS field at all time slices in [0, beta] with +1 and -1 randomly
    
    MatrixXd h = MatrixXd::Random(L,N);     // HS field h_{l, i}
    
    // Hopping matrix 1 D - 5 sites PBC
    
    MatrixXd K(5,5) ;
    
    K << 0, 1, 0, 0, 1,
         1, 0, 1, 0, 0,
         0, 1, 0, 1, 0,
         0, 0, 1, 0, 1,
         1, 0, 0, 1, 0;
    
    cout <<"\n\n K = \n\n" << K << "\n\n";
    
    // Compute matrix 'prefactor' of the B-matrices e^{t\delta\tau K}
    
    MatrixXd B_preFactor = t*dt*K.exp();
    
    cout <<"\n\n exp (t * dt * K) = \n\n" << B_preFactor << "\n\n";
    
    // Fill the HS field matrix
    
    for (l = 0; l < L; l++)
    {
        
        for (i = 0; i < N; i++)
        {
            if ( h(l, i) < 0 )
            {
                h(l, i) = -1;
            }
            else
            {
                h(l, i) = 1;
            }
        }
    }
    
    MatrixXd h_new = h;     // New HS field h_{l, i} - to make determinant updates
    
    cout << "Here's the Hubbard Stratonovich Field\n\n h = \n\n" << h << endl ;
    
    MatrixXd BpOld[L]; // B plus
    MatrixXd BmOld[L]; // B minus
    
    MatrixXd BpNew[L]; // B plus
    MatrixXd BmNew[L]; // B minus
    
    // Build the B-matrices
    
    for (l = 0; l < L; l++)
    {
        BpOld[l] =  MatrixXd::Identity(N,N);
        
        BpOld[l].diagonal() = nu * h.row(l);
        
        BpOld[l] = B_preFactor * BpOld[l].exp();
        
//        cout << "\n\nHere is a B_plus matrix\n\n" << BpOld[l] << endl;
        
        BmOld[l] =  MatrixXd::Identity(N,N);
        
        BmOld[l].diagonal() = (-nu) * h.row(l);
        
        BmOld[l] = B_preFactor * BmOld[l].exp();
        
//        cout << "\n\nHere is a B_minus matrix\n\n" << BmOld[l] << endl;
        
        // The new matrices are initially set as equal to the old ones
        
        BpNew[l] = BpOld[l];
        
        BmNew[l] = BmOld[l];
        
    }
    
    // Build the M-matrices
    
    MatrixXd MPlusOld = MatrixXd::Identity(N,N);
    MatrixXd MMinusOld = MatrixXd::Identity(N,N);
    MatrixXd MPlusNew;
    MatrixXd MMinusNew;
    
    for (l = 0; l < L; l++)
    {
        MPlusOld *= BpOld[L - 1 - l];
        MMinusOld *= BmOld[L - 1 - l];
    }
    
    MPlusOld += MatrixXd::Identity(N,N);
    MMinusOld += MatrixXd::Identity(N,N);
    
//    cout << "\n\nHere is a M_plus matrix\n\n" << MPlusOld << endl;
//    cout << "\n\nHere is a M_minus matrix\n\n" << MMinusOld << endl;
//
//    cout << "\n\nHere is the M_plus determinant\n\n" << MPlusOld.determinant() << endl;
//    cout << "\n\nHere is the M_minus determinant\n\n" << MMinusOld.determinant() << endl;
    
    double detOldPlus = MPlusOld.determinant();
    double detOldMinus = MMinusOld.determinant();
    double detsProdOld = detOldPlus * detOldMinus;
    
    double detsProdNew;
    
    double detNewPlus;
    double detNewMinus;
    
    // TODO: When you accept a step, you must update the Bp and Bm matrices corresponding to the chosen l - > Figure out a way to update this efficiently
    
    double acceptanceRatio;
    
    //inititialize chosen entry of HS field matrix
    l_chosen = 0;
    i_chosen = 0;
    
    for (step = 0; step < totalMCSteps; step++)
    {

//        cout << "Chosen l: " << l_chosen << endl;
//        cout << "Chosen i: " << i_chosen << endl;
        
        // flip
        h_new(l_chosen, i_chosen) *= -1;
        
        // rebuild B-matrices

        BpNew[l_chosen] =  MatrixXd::Identity(N,N);
    
        BpNew[l_chosen].diagonal() = nu * h_new.row(l_chosen);
        
        BpNew[l_chosen] = B_preFactor * BpNew[l_chosen].exp();
        
//        cout << "\n\nHere is a B_plus matrix\n\n" << BpOld[l_chosen] << endl;
//
//        cout << "\n\nHere is a B_plus matrix\n\n" << BpNew[l_chosen] << endl;

        BmNew[l_chosen] =  MatrixXd::Identity(N,N);

        BmNew[l_chosen].diagonal() = (-nu) * h_new.row(l_chosen);

        BmNew[l_chosen] = B_preFactor * BmNew[l_chosen].exp();
        
//        cout << "\n\nHere is a B_minus matrix\n\n" << BmOld[l_chosen] << endl;
//
//        cout << "\n\nHere is a B_minus matrix\n\n" << BmNew[l_chosen] << endl;
        
        // rebuild M-matrices
        
        MPlusNew = MatrixXd::Identity(N,N);
        MMinusNew = MatrixXd::Identity(N,N);
        
        for (l = 0; l < L; l++)
        {
            MPlusNew *= BpNew[L - 1 - l];
            MMinusNew *= BmNew[L - 1 - l];
        }

        MPlusNew += MatrixXd::Identity(N,N);
        MMinusNew += MatrixXd::Identity(N,N);

//        cout << "\n\nHere is a M_plus matrix\n\n" << MPlusNew << endl;
//        cout << "\n\nHere is a M_minus matrix\n\n" << MMinusNew << endl;
//
//        cout << "\n\nHere is the M_plus determinant\n\n" << MPlusNew.determinant() << endl;
//        cout << "\n\nHere is the M_minus determinant\n\n" << MMinusNew.determinant() << endl;

        detNewPlus = MPlusNew.determinant();
        detNewMinus = MMinusNew.determinant();
        detsProdNew = detNewPlus * detNewMinus;

        acceptanceRatio = detsProdNew / detsProdOld;

        cout << "\n\nacceptance ratio: " << acceptanceRatio << endl;
        
        if (acceptanceRatio >= 1 )
        {
            // revert changes
            BpNew[l_chosen] = BpOld[l_chosen];
            BmNew[l_chosen] = BmOld[l_chosen];
            MPlusNew = MPlusOld;
            MMinusNew = MMinusOld;
            h_new(l_chosen, i_chosen) = h(l_chosen, i_chosen);
            detsProdNew = detsProdOld;
        }
        else
        {
            // draw random number to decide
            
            decisionMaker = dis(gen);
            
            cout << "Random number decides whether or not to accept: " << decisionMaker << endl;
            
            if (decisionMaker <= acceptanceRatio) // else do nothing
            {
                //update everything
                BpOld[l_chosen] = BpNew[l_chosen];
                BmOld[l_chosen] = BmNew[l_chosen];
                MPlusOld = MPlusNew;
                MMinusOld = MMinusNew;
                h(l_chosen, i_chosen) = h_new(l_chosen, i_chosen);
                detsProdOld = detsProdNew;
            }
            else
            {
                // revert changes
                BpNew[l_chosen] = BpOld[l_chosen];
                BmNew[l_chosen] = BmOld[l_chosen];
                MPlusNew = MPlusOld;
                MMinusNew = MMinusOld;
                h_new(l_chosen, i_chosen) = h(l_chosen, i_chosen);
                detsProdNew = detsProdOld;
            }
            
        }
        
        // 'pick' a site in (d+1)-dimensional space
        
        if (i_chosen < N - 1)
        {
            i_chosen += 1;
            continue;
        }
        
        else
        {
            if (l_chosen < L - 1)
            {
                l_chosen += 1;
                i_chosen = 1;
            }
            else
            {
                l_chosen = 1;
                i_chosen = 1;
            }
        }
        
    
    }
    
    
    return 0;
}
