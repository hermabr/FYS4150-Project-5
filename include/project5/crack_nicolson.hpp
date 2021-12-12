#pragma once

#include <armadillo>

enum struct slits{
    one, two, three
};

// TODO: Rename class name
class CrackSystem {
    public:
        CrackSystem(double dt, double h, double T, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double v_0);
    private:
        // void initialize_A_B(int M, arma::cx_mat & A, arma::cx_mat & B, double dt, double h, arma::cx_mat v);
        int ij_to_k(int i, int j);
        double j_to_x(int j);
        double i_to_y(int i);
        arma::superlu_opts opts;
        int M, M_star, M_star_square;
        double h, dt, T;
        arma::sp_cx_mat A, B;
        arma::sp_mat V;
        arma::cx_vec u;
        void initialize_A_B();
        arma::cx_vec solve_for_u_next(arma::cx_vec u);
        arma::cx_vec initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);
        /**
         * @brief Initialize V for a double slit set up. Sets values at grid points where there is a wall to v_0, else 0
         * 
         * @param v_0 
         * @return arma::sp_mat 
         */
        arma::sp_mat initialize_V(double v_0, slits slits);
        /**
         * @brief Probability at grid point (i, j).
         * Calculated as |U(i, j)|²
         * 
         * @param i 
         * @param j 
         * @return double The probabilty at the grid point
         */
        double probability_at(int i, int j);
        /**
        * @brief The sum of the probabilities at each point in the grid. 
        * Since our boundary condition is 0, we expect this to be 1 at all times
        * 
        * @return double The total probabilty
        */
        double total_probability();
};
