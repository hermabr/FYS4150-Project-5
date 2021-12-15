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

/**
 * @brief Parse the content of the infile to a setup for the system
 * 
 * @param filename name of the infile
 * @return Config the setup of the system
 */
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

/**
 * @brief Produces filename to correspond with the filename
 * 
 * @param config The setup of the system
 * @param slits The number of slits
 * @return string A name for the outfile
 */
string to_string(Slits slits){
    string s = "";
    switch(slits){
        case Slits::one: s += "simple"; break;
        case Slits::two: s += "double"; break;
        case Slits::three: s += "triple"; break;
    }
    return s;
}

void setup_and_run_system(Slits slits, int config_nr, Config config){
    CrackSystem cs(config, slits);
    // cs.simulate("output/data/" + to_string(slits) + "_" + "config" + to_string(config_nr) + "dt=" + to_string(config.dt) + ".bin");
    cs.simulate("output/data/" + to_string(slits) + "_" + "config" + to_string(config_nr) + ".bin");
}

int main() {
    for (int i = 1; i <= 3; i++){
        string config_file_name = "config" + to_string(i) + ".in";
        Config config = parse_config(config_file_name);
        cerr << "Using config "; config.print();
        Slits slits;
        if (i == 3){
            for (int nslits = Slits::one; nslits <= Slits::three; nslits++){
                slits = static_cast<Slits>(nslits);
                setup_and_run_system(slits, i, config);
            }
        }
        else {
            slits = Slits::two;
            setup_and_run_system(slits, i, config);
        }   
    }
    return 0;
}
