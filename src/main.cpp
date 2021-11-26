#include <complex>
#include <iostream>
#include <armadillo>

using namespace std;

#define cd complex<double>
#define cmat arma::cx_mat

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
    const cd i(0, 1);

    cd r = i * dt / (2. * h * h);

    cerr << "I am not sure how to initialize v(i,j )" << endl;
    cd v_i_j = 1.;
    for (int k = 0; k < N_square; k++) {
        a(k) = 1. + 4. * r + i * dt / 2. * v_i_j;
        b(k) = 1. - 4. * r - i * dt / 2. * v_i_j;
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

int main() {
    int M = 5;
    cmat A, B;
    cd dt = 0.1, h = 0.1;

    initialize_A_B(M, A, B, dt, h);

    A.print();

    return 0;
}
