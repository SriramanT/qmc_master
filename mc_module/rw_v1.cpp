//
//  rw.cpp
//  
//
//  Created by Francisco Brito on 28/02/2018.
//

#include <iostream>
#include <armadillo>
#include <stdio.h>
#include <cstdlib>
#include <typeinfo>
#include <cmath>
#include <random>
#include <fstream>

#define X0 3

using namespace std;
using namespace arma;

// Armadillo documentation is available at:
// http://arma.sourceforge.net/docs.html

int main(){
    
    ofstream x_data;
    x_data.open("x_data.txt");
    
    int i;
    // probability to go to the right
    double p_r = 0.4;
    // probability to go to the left
    double p_l = 0.4;
    // probability of sticking to the wall
    double p_wall = 0.2; // 1 - pr - pl would be the "normal" scenario
    // box size
    int L = 10;
    // initialize position and sums to compute averages. write to file
    int x = X0;
    int n_steps = 1000;
    double x_sum = x*1.;
    double x2_sum = pow(x, 2);
    double decision_maker;
    
    // first entry of txt file. x, < x > , < x^2 >
    x_data << i + 1 << " " << x << " " << x_sum/(i+1) << " " << x2_sum/(i+1) << endl;
    
    for (i = 1; i < n_steps; i++){
        
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0.0, 1.0);
        decision_maker = dis(gen);
        
        if (x == 0)
        {
            if (decision_maker >= p_wall)
            {
                x += 1;
            }
            x_sum += x*1.;
            x2_sum += pow(x, 2);
            // corr
            x_data << i + 1 << " " << x << " " << x_sum/(i+1) << " " << x2_sum/(i+1) << endl;
            continue;
        }
        if (x == L)
        {
            if (decision_maker >= p_wall)
            {
                x -= 1;
            }
            x_sum += x*1.;
            x2_sum += pow(x, 2);
            x_data << i + 1 << " " << x << " " << x_sum/(i+1) << " " << x2_sum/(i+1) << endl;
            continue;
        }
        if (decision_maker <= p_r)
        {
            x += 1;
        }
        if (decision_maker >= 1 - p_l)
        {
            x -= 1;
        }
        x_sum += x*1.;
        x2_sum += pow(x, 2);
        x_data << i + 1 << " " << x << " " << x_sum/(i+1) << " " << x2_sum/(i+1) << endl;
    }
    
    x_data.close();
    
    return 0;
}
