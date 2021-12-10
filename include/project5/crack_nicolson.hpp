#pragma once

#include <armadillo>

class CrackSystem {
    public:
        CrackSystem(int M, double h, double dt);
    private:
        // void initialize_A_B(int M, arma::cx_mat & A, arma::cx_mat & B, double dt, double h, arma::cx_mat v);
        int ij_to_k(int i, int j);
        int M, M_star, M_star_square;
        double dt, h;
        arma::sp_cx_mat A, B; // TODO: Remove this line if it aint woking, and change to arma::cx_mat
        void initialize_A_B();
        arma::cx_vec solve_for_u_next(arma::cx_vec u);
        arma::cx_vec initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);

};
