//
//  matrices.cpp
//  
//
//  Created by Francisco Brito on 23/04/2018.
//

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

Eigen::MatrixXd build_Bmatrix(bool sign, double nu, int nSites, Eigen::RowVectorXd h_l, Eigen::MatrixXd B_preFactor, double dt, double mu)
{
    Eigen::MatrixXd B;
    Eigen::MatrixXd preFactor = B_preFactor;
    B = Eigen::MatrixXd::Identity(nSites, nSites);
    int i;
    
    if (sign == true)
    {
        for (i = 0; i < nSites; i++)
        {
            B(i, i) = nu * h_l(i) + dt * mu ;
        }
    }
    else
    {
        for (i = 0; i < nSites; i++)
        {
            B(i, i) = -1 * nu * h_l(i) + dt * mu ;
        }
    }
    B = preFactor * B.exp();
    return B;
}

Eigen::VectorXd uSigma(int nSites, Eigen::MatrixXd Green, int i_chosen)
{
    Eigen::MatrixXd e_i;
    e_i = Eigen::MatrixXd::Zero(nSites , 1) ;
    e_i(i_chosen, 0) = 1;
    return (Eigen::MatrixXd::Identity(nSites, nSites) - Green) * e_i;
}

Eigen::VectorXd wSigma(int nSites, Eigen::MatrixXd Green, int i_chosen)
{
    Eigen::MatrixXd e_i;
    e_i = Eigen::MatrixXd::Zero(nSites , 1) ;
    e_i(i_chosen, 0) = 1;
    return Green.transpose() * e_i;
}