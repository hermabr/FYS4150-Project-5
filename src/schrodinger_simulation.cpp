#include <iostream>
#include <iomanip>
#include <assert.h> 
#include <armadillo>
#include "project5/config.hpp"
#include "project5/schrodinger_simulation.hpp"

using namespace std;
using namespace std::complex_literals;

// IMPORTANT NOTE: Documentation is found in the header file.

SchrodingerSimulation::SchrodingerSimulation(Config config, Slits Slits) :  h(config.h), dt(config.dt), T(config.T) {
    if (1.0/h != floor(1.0/h)) throw invalid_argument("1/h must be a whole number");

    opts.symmetric = true;

    // define M, M_star and M_star_square, where M_star=M-2
    M = (int) (1.0 / h) + 1;
    M_star = M - 2; M_star_square = M_star * M_star;
    
    // initialize the matrix V
    V = initialize_V(config.v_0, Slits);
    // initialize the matrices A and B
    initialize_A_B();
    // initialize the matrix u
    u = initialize_u(config.x_c, config.y_c, config.s_x, config.s_y, config.p_x, config.p_y);
}

int SchrodingerSimulation::ij_to_k(int i, int j){
    // verify that i and j are in the correct range
    if (j >= M_star_square or i >= M_star_square or j < 0 or i < 0) throw out_of_range("Requires 0 <= i, j < M - 2");
    return i * M_star + j;
}

double SchrodingerSimulation::i_to_y(int i){
    // get the y coordinate of the ith point
    return 1 - (double)(i + 1) / M_star;
}

double SchrodingerSimulation::j_to_x(int j){
    // get the x coordinate of the jth point
    return (double)(j + 1) / M_star;
}



void SchrodingerSimulation::initialize_A_B() {
    // initialize the matrices A and B
    A = arma::sp_cx_mat(M_star_square, M_star_square);
    B = arma::sp_cx_mat(M_star_square, M_star_square);

    // initialize empty vectors a and b
    arma::cx_vec a(M_star_square, arma::fill::zeros);
    arma::cx_vec b(M_star_square, arma::fill::zeros);

    // set the value of r
    complex<double> r = 1i * dt / (2. * h * h);

    // set the values of a and b
    for (int k = 0; k < M_star_square; k++) {
        int i = k / M_star;
        int j = k % M_star;
        double v_i_j = V(i, j);
        a(k) = 1. + 4. * r + 1i * dt / 2. * v_i_j;
        b(k) = 1. - 4. * r - 1i * dt / 2. * v_i_j;
    }

    // set the values of A and B
    for (int k = 0; k < M_star_square; k++) {
        int i = k / M_star;
        int j = k % M_star;

        A(k, k) = a(k);
        B(k, k) = b(k);

        if (i != 0) {
            A(k, (i-1) * M_star + j) = -r;
            B(k, (i-1) * M_star + j) = r;
        }

        if (i != M_star - 1) {
            A(k, (i + 1) * M_star + j) = -r;
            B(k, (i + 1) * M_star + j) = r;
        }

        if (j != 0) {
            A(k, k - 1) = -r;
            B(k, k - 1) = r;
        }

        if (j != M_star - 1) {
            A(k, k + 1) = -r;
            B(k, k + 1) = r;
        }
    }
}

arma::cx_vec SchrodingerSimulation::solve_for_u_next(arma::cx_vec u) {
    // multiply the matrix B by the vector u
    arma::cx_vec b = B * u;
    // solve the system using superlu
    arma::cx_vec new_u = arma::spsolve(A, b, "superlu", opts);
    return new_u;
}

arma::cx_vec SchrodingerSimulation::initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y) {
    // initialize the vector u
    arma::cx_vec u(M_star_square);
    double s = 0;
    for (int i = 0; i < M_star; i++){
        double y = i_to_y(i);
        for (int j = 0; j < M_star; j++){
            double x = j_to_x(j);
            // set the value of v to the equation from the report
            complex<double> v = exp(
                -(x - x_c) * (x - x_c) / (2 * sigma_x * sigma_x) 
                - (y - y_c) * (y - y_c) / (2 * sigma_y * sigma_y)
                + 1.i * p_x * (x - x_c) + 1.i * p_y * (y - y_c)
                );
            // add the sum of the probability to s, for the normalization
            s += real(v) * real(v) + imag(v) * imag(v);
            u(ij_to_k(i, j)) = v;
        }
    }
    u /= sqrt(s);
    return u;
}


arma::sp_mat SchrodingerSimulation::initialize_V(double v_0, Slits slits){ 
    arma::sp_mat V(M_star, M_star);
    arma::vec slit_tops;

    // set the correct number of slits according to the slit parameter
    switch(slits) {
        case Slits::one: 
            slit_tops = {.525};
            break;
        case Slits::two: 
            slit_tops = {.575, .475};
            break;
        case Slits::three:
            slit_tops = {.625, .525, .425};
            break;
    };

    // define is_slit to be if the current y coordinate is in the slit
    auto is_slit = [slit_tops](double y){
        float slit_width = .05;
        for (auto i = 0; (long long unsigned int) i < slit_tops.size(); i++){
            if (y <= slit_tops(i) and y > slit_tops(i) - slit_width)
                return true;
        }
        return false;
    };

    // set the value of V
    for (int i = 0; i < M_star; i++){
        double y = i_to_y(i);
        if (is_slit(y))
            continue;
        int j = (.49 - h) / h;
        while (j_to_x(j) < .51){
            V(i, j) = v_0;
            j++;
        }
    }
    return V;
}

void SchrodingerSimulation::simulate(string outfile){
    // set the time step to be the total time T divided by the size of dt
    int timesteps = (int) (T / dt);

    // initialize the matrix U, to be the simulation for all the time steps
    arma::cx_mat U = arma::cx_mat(timesteps + 1, M_star_square);
    for (int i = 0; i < M_star_square; i ++){
            U(0, i) = u(i);
    }
    
    for (int t = 1; t <= timesteps; t ++){
        cerr << "Simulating time step " << t << "/" << timesteps << "\r" << flush;
        // solve for the next u
        u = solve_for_u_next(u);
        // write the current u to the matrix U
        for (int i = 0; i < M_star_square; i ++){
            U(t, i) = u(i);
        }
    }
    // save the matrix U to a binary armadillo file
    U.save(outfile);
}
