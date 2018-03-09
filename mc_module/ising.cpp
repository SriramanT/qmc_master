//
//  ising.cpp
//  
//
//  Created by Francisco Brito on 03/03/2018.
//

#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <typeinfo>
#include <cmath>
#include <random>
#include <time.h>
#include <Eigen/Dense>

#define N_SITES 8

using namespace Eigen;
using namespace std;
int main()
{
    MatrixXd m = MatrixXd::Random(3,3);
    //m = (m + MatrixXd::Constant(3,3,1.2)) * 50;
    cout << "m =" << endl << m << endl;
    //VectorXd v(3);
    //v << 1, 2, 3;
    //cout << "m * v =" << endl << m * v << endl;
}
