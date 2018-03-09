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
    int k;
    // probability to go to the right
    double p_r = 0.4;
    // probability to go to the left
    double p_l = 0.4;
    // self loop
    double p_0 = 1 - p_r - p_l;
    // probability of sticking to the wall
    double p_wall = p_0; //0.99; // p_0 = 1 - p_r - p_l would be the "normal" scenario
    // box size
    int L = 100;

    int n_steps = 1000;
    int n_avgs = 100;
    double decision_maker;
    
    //for (k = 1; k < n_avgs + 1; k++){
    
        // initialize position and sums to compute averages. write to file
        mat x(n_steps, 2);
        x(0) = X0;
        double x_sum = x(0)*1.;
        double x2_sum = pow(x(0), 2);
    
        for (i = 1; i < n_steps; i++){
            
            // Draw a random number to decide whether of not to walk and if so which way
            random_device rd;
            mt19937 gen(rd());
            uniform_real_distribution<> dis(0.0, 1.0);
            decision_maker = dis(gen);
            
            x(i, 1) = i;
            
            if (x(i-1, 0) == 0)
            {
                if (decision_maker >= p_wall)
                {
                    x(i, 0) = x(i-1, 0) + 1;
                }
                x_sum += x(i, 0)*1.;
                x2_sum += pow(x(i, 0), 2);
                continue;
            }
            if (x(i-1, 0) == L)
            {
                if (decision_maker >= p_wall)
                {
                    x(i, 0) = x(i-1, 0) - 1;
                }
                x_sum += x(i, 0)*1.;
                x2_sum += pow(x(i, 0), 2);
                continue;
            }
            if (decision_maker <= p_r)
            {
                x(i, 0) = x(i-1, 0) + 1;
            }
            if (decision_maker >= 1 - p_l)
            {
                x(i, 0) = x(i-1, 0) - 1;
            }
            x_sum += x(i, 0)*1.;
            x2_sum += pow(x(i, 0 ), 2);
        }
    
    x.save("x.txt", raw_ascii);
    
    //x_data << k << " " << x_sum/n_steps << endl;
        
    //}
    
    //x_data.close();
    
    return 0;
}
