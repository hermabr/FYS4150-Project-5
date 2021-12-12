// #include <iostream>
// #include <complex>
//
// using namespace std::complex_literals;
//
// int main() {
//     auto drugs = 1.0 + 3.0i;
//     std::cout << drugs << '\n';
// }

#include <complex>
#include <fstream>
#include <assert.h> 
#include <iostream>
#include <armadillo>
#include <project5/config.hpp>
#include <project5/crack_nicolson.hpp>

#define cd complex<double>
#define cmat arma::cx_mat

using namespace std;
using namespace std::complex_literals;

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

Config parse_config(string filename) {
    ifstream config_file(filename);

    string variable_name;
    string variable_names[10] = {"h", "dt", "T", "x_c", "y_c", "s_x", "s_y", "p_x", "p_y", "v_0"};

    double values[10];

    for (int i = 0; i < 10; i++) {
        config_file >> variable_name;
        assert (variable_name == variable_names[i]);
        config_file >> values[i];
    }

    Config config {values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9]};

    return config;
}

int main() {
    Config config = parse_config("config.in");
    cerr << "Using config "; config.print();
    CrackSystem cs(config, slits::two);
    cs.simulate();

    return 0;
}
