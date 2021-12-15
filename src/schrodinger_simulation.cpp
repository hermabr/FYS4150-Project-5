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

    M = (int) (1.0 / h) + 1;
    M_star = M - 2; M_star_square = M_star * M_star;
    
    V = initialize_V(config.v_0, Slits);
    initialize_A_B();
    u = initialize_u(config.x_c, config.y_c, config.s_x, config.s_y, config.p_x, config.p_y);
    cout << scientific << setprecision(15);
}

// TODO: MIGHT THIS BE (i-1) and (j-1), not i and j?
int SchrodingerSimulation::ij_to_k(int i, int j){
    if (j >= M_star_square or i >= M_star_square or j < 0 or i < 0) throw out_of_range("Requires 0 <= i, j < M - 2");
    return i * M_star + j;
}

double SchrodingerSimulation::i_to_y(int i){
    return 1 - (double)(i + 1) / M_star;
}

double SchrodingerSimulation::j_to_x(int j){
    return (double)(j + 1) / M_star;
}



void SchrodingerSimulation::initialize_A_B() {
    A = arma::sp_cx_mat(M_star_square, M_star_square);
    B = arma::sp_cx_mat(M_star_square, M_star_square);

    arma::cx_vec a(M_star_square, arma::fill::zeros);
    arma::cx_vec b(M_star_square, arma::fill::zeros);

    complex<double> r = 1i * dt / (2. * h * h);

    for (int k = 0; k < M_star_square; k++) {
        int i = k / M_star;
        int j = k % M_star;
        double v_i_j = V(i, j);
        a(k) = 1. + 4. * r + 1i * dt / 2. * v_i_j;
        b(k) = 1. - 4. * r - 1i * dt / 2. * v_i_j;
    }

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
    arma::cx_vec b = B * u;
    arma::cx_vec new_u = arma::spsolve(A, b, "superlu", opts);
    return new_u;
}

arma::cx_vec SchrodingerSimulation::initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y) {
    arma::cx_vec u(M_star_square);
    double s = 0;
    for (int i = 0; i < M_star; i++){
        double y = i_to_y(i);
        for (int j = 0; j < M_star; j++){
            double x = j_to_x(j);
            complex<double> v = exp(
                -(x - x_c) * (x - x_c) / (2 * sigma_x * sigma_x) 
                - (y - y_c) * (y - y_c) / (2 * sigma_y * sigma_y)
                + 1.i * p_x * (x - x_c) + 1.i * p_y * (y - y_c)
                );
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

    auto is_slit = [slit_tops](double y){
        float slit_width = .05;
        for (int i = 0; i < slit_tops.size(); i++){
            if (y <= slit_tops(i) and y > slit_tops(i) - slit_width)
                return true;
        }
        return false;
    };

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
    int timesteps = (int) (T / dt);

    arma::cx_mat U = arma::cx_mat(timesteps + 1, M_star_square);
    for (int i = 0; i < M_star_square; i ++){
            U(0, i) = u(i);
    }
    
    for (int t = 1; t <= timesteps; t ++){
        u = solve_for_u_next(u);
        for (int i = 0; i < M_star_square; i ++){
            U(t, i) = u(i);
        }
    }
    U.save(outfile);
}
