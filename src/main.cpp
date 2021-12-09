#include <armadillo>
#include <complex>
#include <iostream>
#include <project5/config.hpp>
#include <project5/crack_nicolson.hpp>

using namespace std;

#define cd complex<double>
#define cmat arma::cx_mat


// TODO: MIGHT THIS BE (i-1) and (j-1), not i and j?
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

int main() {
    // TODO: How should M be initialized?
    int M = 5;
    double dt = 0.1, h = 0.1;

    CrackSystem cs(M, dt, h);
    // int M = 5;
    // cmat A, B;
    // cd dt = 0.1, h = 0.1;
    //
    // initialize_A_B(M, A, B, dt, h);

    // A.print();

    return 0;
}
