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

using namespace std;
using namespace arma;

// Armadillo documentation is available at:
// http://arma.sourceforge.net/docs.html

int main(){
    
    int i;
    double p_r = 0.4; // probability to go to the right
    double p_l = 0.4; // probability to go to the left
    int L = 10; // box size
    int x = 3; // initial position
    int N_steps = 10;
    double decision_maker;
    
    for (i = 0; i < N_steps; i++){
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        decision_maker = dis(gen);
        
        if ( (x != 0) or (x != L) ){
    
        
            if (decision_maker <= p_r){
                x += 1;
            }
        
            if (decision_maker >= 1 - p_l){
                x -= 1;
            }
        cout << x << " " << i + 1 << endl;
            
        }
    }
    
    
    
    return 0;
}
