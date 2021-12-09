#include <armadillo>
#include <complex>
#include <iostream>
#include <config.hpp>

using namespace std;

#define cd complex<double>
#define cmat arma::cx_mat

int idx(int i, int j, int M){
    if (j >= M - 2 or i >= M - 2) throw exception();
    return i * (M - 2) + j;
}

bool has_flag(const std::string& option, char** begin, char** end){
    return std::find(begin, end, option) != end;
}

void print_help_message() {
    cout << "Usage" << endl;
    cout << "\t./runner [flags]" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "\t-h\tShow this help message" << endl;
}

void initialize_A_B(int M, cmat & A, cmat & B, cd dt, cd h) {
    int N = M - 2;
    int N_square = N * N;

    A = cmat(N_square, N_square, arma::fill::zeros);
    B = cmat(N_square, N_square, arma::fill::zeros);

    arma::cx_vec a(N_square, arma::fill::zeros);
    arma::cx_vec b(N_square, arma::fill::zeros);

    cerr << "I am not sure how to initialize i" << endl;

    cd r = e_i * dt / (2. * h * h);

    cerr << "I am not sure how to initialize v(i,j )" << endl;
    cd v_i_j = 1.;
    for (int k = 0; k < N_square; k++) {
        a(k) = 1. + 4. * r + e_i * dt / 2. * v_i_j;
        b(k) = 1. - 4. * r - e_i * dt / 2. * v_i_j;
    }

    for (int k = 0; k < N_square; k++) {
        int i = 1 + k / N;
        int j = 1 + k % N;

        A(k, k) = a(k);
        B(k, k) = b(k);

        if (i != 1) {
            A(k, (i-2) * N + j - 1) = -r;
            B(k, (i-2) * N + j - 1) = r;
        }

        if (i != N) {
            A(k, i * N + j - 1) = -r;
            B(k, i * N + j - 1) = r;
        }

        if (j != 1) {
            A(k, k - 1) = -r;
            B(k, k - 1) = r;
        }

        if (j != N) {
            A(k, k + 1) = -r;
            B(k, k + 1) = r;
        }
    }
}

void solve_for_u_next(cmat & A, arma::cx_vec u_next, cmat & B, arma::cx_vec u){
    arma::cx_vec b = B * u;
    solve(u_next, A, b);
}

arma::cx_vec initialize_u(int M, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y){
    arma::cx_vec u((M - 2) * (M - 2));
    double h = 1. / M;
    cd s = 0;
    for (int i = 0; i < M - 2; i++){
        double y = (i + 1) * h;
        for (int j = 0; j < M - 2; j++){
            double x = (j + 1) * h;
            cd v = exp(
                -(x - x_c) * (x - x_c) / (2 * sigma_x * sigma_x) 
                - (y - y_c) * (y - y_c) / (2 * sigma_y * sigma_y)
                + e_i * p_x * (x - x_c) + e_i * p_y * (y - y_c)
                );
            s += v * v;
            u(idx(i, j, M)) = v;
        }
    }
    u /= s;
    return u;
}

int main() {
    int M = 5;
    cmat A, B;
    cd dt = 0.1, h = 0.1;

    initialize_A_B(M, A, B, dt, h);

    A.print();

    return 0;
}
