#include <iostream>
#include <assert.h> 
#include <armadillo>
#include "project5/config.hpp"
#include "project5/crack_nicolson.hpp"

using namespace std;

// TODO: Remove these and replace with the non-macros
#define cd complex<double>
#define cmat arma::cx_mat

CrackSystem::CrackSystem(int M, double dt, double h) :  M(M), dt(dt), h(h) {
    // std::cout << M << " " << dt << " " << h << std::endl;
    M_star = M - 2; M_star_square = M_star * M_star;
    
    // TODO: DO SOMETHING WITH v
    cmat v = cmat(M, M);
    v.fill(cd(1, 0));

    initialize_A_B();
    initialize_u(1, 1, 1, 1, 1, 1);
}

// TODO: MIGHT THIS BE (i-1) and (j-1), not i and j?
int CrackSystem::ij_to_k(int i, int j){
    // TODO: is this assert correct?
    assert (j <= M_star && i <= M_star);
    // if (j >= M_star_square or i >= M_star_square) throw exception();
    // return i * M_star_square + j;
    return i * M_star + j;
}

// TODO: Possible to initialize the diagonal elements in one liners?
// void CrackSystem::initialize_A_B(int M, cmat & A, cmat & B, double dt, double h, cmat v) {
void CrackSystem::initialize_A_B() {
    A = cmat(M_star_square, M_star_square, arma::fill::zeros);
    B = cmat(M_star_square, M_star_square, arma::fill::zeros);

    arma::cx_vec a(M_star_square, arma::fill::zeros);
    arma::cx_vec b(M_star_square, arma::fill::zeros);

    cd r = e_i * dt / (2. * h * h);

    cerr << "I am not sure how to initialize v(i,j )" << endl;
    // TODO: I assume this one is categorically better (except the fact that it doesn't work yet)
    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         int k = ij_to_k(i, j);
    //
    //         a(k) = 1. + 4. * r + e_i * dt / 2. * v_i_j;
    //         b(k) = 1. - 4. * r - e_i * dt / 2. * v_i_j;
    //     }
    // }
    
    cd v_i_j = 1.;
    for (int k = 0; k < M_star_square; k++) {
        a(k) = 1. + 4. * r + e_i * dt / 2. * v_i_j;
        b(k) = 1. - 4. * r - e_i * dt / 2. * v_i_j;
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

void CrackSystem::solve_for_u_next(arma::cx_vec u_next, arma::cx_vec u) {
    cerr << "If this doesnt work: Fix in the header file, by changing A and B from sp_cx_mat to cx_mat" << endl;
    arma::cx_vec b = B * u;
    // solve(u_next, A, b);
    u_next = arma::spsolve(A, b);
    // arma::cx_vec new_u = arma::spsolve(A, b, "superlu");
    // new_u.print();
}


// arma::cx_vec SystemMatrix::calc_new_u(const arma::cx_vec & u) const {
//     arma::cx_vec b = B * u;
//     arma::cx_vec u_new = arma::spsolve(A, b, "superlu");
//     return u_new;
// }

arma::cx_vec CrackSystem::initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y) {
    arma::cx_vec u(M_star * M_star);
    double h = 1. / M;
    cd s = 0;
    for (int i = 0; i < M_star; i++){
        double y = (i + 1) * h;
        for (int j = 0; j < M_star; j++){
            double x = (j + 1) * h;
            cd v = exp(
                -(x - x_c) * (x - x_c) / (2 * sigma_x * sigma_x) 
                - (y - y_c) * (y - y_c) / (2 * sigma_y * sigma_y)
                + e_i * p_x * (x - x_c) + e_i * p_y * (y - y_c)
                );
            s += v * v;
            u(ij_to_k(i, j)) = v;
        }
    }
    u /= s;
    return u;
}
