#pragma once

#include <armadillo>

// TODO: Rename class name
class CrackSystem {
    public:
        CrackSystem(double dt, double h, double T, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double v_0);
    private:
        // void initialize_A_B(int M, arma::cx_mat & A, arma::cx_mat & B, double dt, double h, arma::cx_mat v);
        int ij_to_k(int i, int j);
        double j_to_x(int j);
        double i_to_y(int i);
        int M, M_star, M_star_square;
        double h, dt, T;
        arma::sp_cx_mat A, B;
        void initialize_A_B();
        arma::cx_vec solve_for_u_next(arma::cx_vec u);
        arma::cx_vec initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);
        arma::sp_mat initialize_V_double_slit(double v_0);
};
