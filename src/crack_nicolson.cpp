#include <iostream>
#include <iomanip>
#include <assert.h> 
#include <armadillo>
#include "project5/config.hpp"
#include "project5/crack_nicolson.hpp"

using namespace std;
using namespace std::complex_literals;

// TODO: Remove these and replace with the non-macros
#define cd complex<double>
#define cmat arma::cx_mat

// TODO: Rename hehe
CrackSystem::CrackSystem(Config config, slits slits) :  h(config.h), dt(config.dt), T(config.T) {
    if (1.0/h != floor(1.0/h)) throw invalid_argument("1/h must be a whole number");
    // TODO: Check that T/dt is a whole number?

    arma::superlu_opts opts;
    opts.symmetric = true;

    M = (int) (1.0 / h) + 1;
    M_star = M - 2; M_star_square = M_star * M_star;
    
    // TODO: DO SOMETHING WITH v
    cmat v = cmat(M, M);
    v.fill(cd(1, 0));
    V = initialize_V(config.v_0, slits);
    initialize_A_B();
    u = initialize_u(config.x_c, config.y_c, config.s_x, config.s_y, config.p_x, config.p_y);
    u.save("output/data/u.bin");
    cout << scientific << setprecision(15);
    cout << "deviation from 1 of total probability after initialization: " << 1 - total_probability() << endl;
    
}

// TODO: MIGHT THIS BE (i-1) and (j-1), not i and j?
int CrackSystem::ij_to_k(int i, int j){
    // TODO: is this assert correct?
    if (j >= M_star_square or i >= M_star_square) throw exception(); // TODO: What exception?
    // return i * M_star_square + j;
    return i * M_star + j;
}

double CrackSystem::i_to_y(int i){
    return 1 - (double)(i + 1) / M_star;
}

double CrackSystem::j_to_x(int j){
    return (double)(j + 1) / M_star;
}



// TODO: Possible to initialize the diagonal elements in one liners?
// void CrackSystem::initialize_A_B(int M, cmat & A, cmat & B, double dt, double h, cmat v) {
void CrackSystem::initialize_A_B() {
    A = arma::sp_cx_mat(M_star_square, M_star_square);
    B = arma::sp_cx_mat(M_star_square, M_star_square);

    arma::cx_vec a(M_star_square, arma::fill::zeros);
    arma::cx_vec b(M_star_square, arma::fill::zeros);

    cd r = 1i * dt / (2. * h * h);

    // TODO: I assume this one is categorically better (except the fact that it doesn't work yet)
    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         int k = ij_to_k(i, j);
    //
    //         a(k) = 1. + 4. * r + 1i * dt / 2. * v_i_j;
    //         b(k) = 1. - 4. * r - 1i * dt / 2. * v_i_j;
    //     }
    // }
    
    for (int k = 0; k < M_star_square; k++) {
        int i = k / M_star;
        int j = k % M_star;
        double v_i_j = V(i, j);
        a(k) = 1. + 4. * r + 1i * dt / 2. * v_i_j;
        b(k) = 1. - 4. * r - 1i * dt / 2. * v_i_j;
    }

    for (int k = 0; k < M_star_square; k++) {
        int i = 1 + k / M_star;
        int j = 1 + k % M_star;

        A(k, k) = a(k);
        B(k, k) = b(k);

        if (i != 1) {
            A(k, (i-2) * M_star + j - 1) = -r;
            B(k, (i-2) * M_star + j - 1) = r;
        }

        if (i != M_star) {
            A(k, i * M_star + j - 1) = -r;
            B(k, i * M_star + j - 1) = r;
        }

        if (j != 1) {
            A(k, k - 1) = -r;
            B(k, k - 1) = r;
        }

        if (j != M_star) {
            A(k, k + 1) = -r;
            B(k, k + 1) = r;
        }
    }
}

arma::cx_vec CrackSystem::solve_for_u_next(arma::cx_vec u) {
    arma::cx_vec b = B * u;
    arma::cx_vec new_u = arma::spsolve(A, b, "superlu", opts);
    return new_u;
}

arma::cx_vec CrackSystem::initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y) {
    arma::cx_vec u(M_star_square);
    double s = 0;
    for (int i = 0; i < M_star; i++){
        double y = i_to_y(i);
        for (int j = 0; j < M_star; j++){
            double x = j_to_x(j);
            cd v = exp(
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


arma::sp_mat CrackSystem::initialize_V(double v_0, slits slits){ 
    //arma::sp_mat V(M_star, M_star, arma::fill::zeros);
    arma::sp_mat V(M_star, M_star);
    arma::vec slit_tops;

    switch(slits) {
        case slits::one: 
            slit_tops = {.525};
            break;
        case slits::two: 
            slit_tops = {.575, .475};
            break;
        case slits::three:
            slit_tops = {.625, .525, .425};
            break;
    };

    auto is_slit = [slit_tops](double y){
        float slit_width = .05;
        for (int i = 0; i < slit_tops.size(); i++){
            if (y <= slit_tops(i) and y >= slit_tops(i) - slit_width) // is inclusive correct?
                return true;
        }
        return false;
    };

    for (int i = 0; i < M_star; i++){
        double y = i_to_y(i);
        if (is_slit(y))
            continue;
        int j = (.49 - h) / h; // is this correct?
        while (j_to_x(j) <= .51){
            V(i, j) = v_0;
            j++;
        }
    }
    return V;
}

double CrackSystem::probability_at(int i, int j){
    cd uij = u(ij_to_k(i, j));
    return real(uij) * real(uij) + imag(uij) * imag(uij);
}

double CrackSystem::total_probability(){
    double p = 0;
    for (int i = 0; i < M_star; i++)
        for (int j = 0; j < M_star; j++){
            p += probability_at(i, j);
        }
            
    return p;
}

void CrackSystem::simulate(){
    // TODO: CHANGE THIS PLZ
    int timesteps = (int) (T / dt);

    arma::cx_mat U = arma::cx_mat(timesteps, M_star_square);
    // TODO: Make sure we have correct time steps 
    for (int t = 1; t <= timesteps; t ++){
        u = solve_for_u_next(u);
        for (int i = 0; i < M_star_square; i ++){
            U(t - 1, i) = u(i);
        }
        cout << "deviation from 1 of total probability after " << t <<" step(s): " << 1 - total_probability() << endl;
    }
    U.save("output/data/UBER.bin");
}
