#include <iostream>
#include <assert.h> 
#include <armadillo>
#include "project5/config.hpp"
#include "project5/crack_nicolson.hpp"

using namespace std;
using namespace std::complex_literals;

// TODO: Remove these and replace with the non-macros
#define cd complex<double>
#define cmat arma::cx_mat

CrackSystem::CrackSystem(double h, double dt, double T, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double v_0) :  h(h), dt(dt), T(T) {
    if (1.0/h != floor(1.0/h)) throw invalid_argument("1/h must be a whole number");

    M = (int) (1.0 / h) + 1;
    M_star = M - 2; M_star_square = M_star * M_star;
    
    // TODO: DO SOMETHING WITH v
    cmat v = cmat(M, M);
    v.fill(cd(1, 0));
    cerr << "C" << endl;
    initialize_V_double_slit(v_0);
    initialize_A_B();
    arma::cx_vec u = initialize_u(x_c, y_c, sigma_x, sigma_y, p_x, p_y);

    u = solve_for_u_next(u);
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

    cerr << "I am not sure how to initialize v(i,j)" << endl;
    // TODO: I assume this one is categorically better (except the fact that it doesn't work yet)
    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         int k = ij_to_k(i, j);
    //
    //         a(k) = 1. + 4. * r + 1i * dt / 2. * v_i_j;
    //         b(k) = 1. - 4. * r - 1i * dt / 2. * v_i_j;
    //     }
    // }
    
    cd v_i_j = 1.;
    for (int k = 0; k < M_star_square; k++) {
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
    arma::cx_vec new_u = arma::spsolve(A, b);
    return new_u;
}

arma::cx_vec CrackSystem::initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y) {
    arma::cx_vec u(M_star * M_star);
    cd s = 0;
    for (int i = 0; i < M_star; i++){
        double y = i_to_y(i);
        for (int j = 0; j < M_star; j++){
            double x = j_to_x(j);
            cd v = exp(
                -(x - x_c) * (x - x_c) / (2 * sigma_x * sigma_x) 
                - (y - y_c) * (y - y_c) / (2 * sigma_y * sigma_y)
                + 1i * p_x * (x - x_c) + 1i * p_y * (y - y_c)
                );
            s += v * v;
            u(ij_to_k(i, j)) = v;
        }
    }
    u /= s;
    return u;
}


arma::sp_mat CrackSystem::initialize_V_double_slit(double v_0){ 
    //arma::sp_mat V(M_star, M_star, arma::fill::zeros);
    arma::sp_mat V(M_star, M_star);
    for (int i = 0; i < M_star; i++){
        double y = i_to_y(i);
        cout << "y: "<< y << endl;
        if ((y <= .575 and y >= .525) or (y <= .475 and y >= .425))
            continue;
        int j = (.49 - h) / h; // is this correct?
        while (j_to_x(j) <= .51){
            cout << "x: "<< j_to_x(j) << endl;
            V(i, j) = v_0;
            j++;
        }
    }
    V.print();
    return V;
}