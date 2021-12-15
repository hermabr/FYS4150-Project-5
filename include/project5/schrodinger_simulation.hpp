#pragma once

#include <armadillo>

enum Slits{
    one, two, three
};

class SchrodingerSimulation {
    public:
        SchrodingerSimulation(Config config, Slits Slits);
        void simulate(std::string outfile);
    private:
        // void initialize_A_B(int M, arma::cx_mat & A, arma::cx_mat & B, double dt, double h, arma::cx_mat v);
        /**
         * @brief Maps values of i and j to k, raising error if values for i and j are out of range
         * 
         * @param i 
         * @param j 
         * @return int 
         */
        int ij_to_k(int i, int j);
        /**
         * @brief Mapping indices j to values x
         * 
         * @param j 
         * @return double 
         */
        double j_to_x(int j);
        /**
         * @brief Mapping indices i to values y
         * 
         * @param i 
         * @return double 
         */
        double i_to_y(int i);
        /**
         * @brief Options for the SuperLU solver used in our Crank-Nicolson implementation
         * 
         */
        arma::superlu_opts opts;
        int M, M_star, M_star_square;
        double h, dt, T;
        /**
         * @brief The matricies required for the Crank-Nicolson approach
         * 
         */
        arma::sp_cx_mat A, B;
        /**
         * @brief The potential matrix
         * 
         */
        arma::sp_mat V;
        /**
         * @brief Holding the current state of the wave function
         * 
         */
        arma::cx_vec u;
        /**
         * @brief Setting up the matrices A and B required for the Crank-Nicolson approach
         * 
         */
        void initialize_A_B();
        arma::cx_vec solve_for_u_next(arma::cx_vec u);
        /**
         * @brief Initialize the vector u as according to a unnormalized gaussian wave packet, 
         * before also normalizing so that the total probability (summing ecah element multiplied 
         * with it's complex conjugate) is 1
         * 
         * @param x_c x-coordinate of the center of the initial wave packet
         * @param y_c y-coordinate of the center of the initial wave packet
         * @param sigma_x width in x-direction of the initial wave packet
         * @param sigma_y width in y-direction of the initial wave packet
         * @param p_x wave packet momenta in x-direction
         * @param p_y wave packet momenta in y-direction
         * @return arma::cx_vec u_0
         */
        arma::cx_vec initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);
        /**
         * @brief Initialize V for a double slit set up. Sets values at grid points where there is a wall to v_0, else 0
         * 
         * @param v_0 
         * @return arma::sp_mat 
         */
        arma::sp_mat initialize_V(double v_0, Slits Slits);
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
