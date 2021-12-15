#include <complex>
#include <fstream>
#include <assert.h> 
#include <iostream>
#include <armadillo>
#include <project5/config.hpp>
#include <project5/schrodinger_simulation.hpp>

using namespace std;
using namespace std::complex_literals;

/**
 * @brief Parse the content of the infile to a setup for the system
 * 
 * @param filename name of the infile
 * @return Config the setup of the system
 */
Config parse_config(string filename) {
    // read the file
    ifstream config_file(filename);

    string variable_name;
    string variable_names[10] = {"h", "dt", "T", "x_c", "y_c", "s_x", "s_y", "p_x", "p_y", "v_0"};

    double values[10];

    // check that the file contains the correct variables and read them
    for (int i = 0; i < 10; i++) {
        config_file >> variable_name;
        // check that the variable is correct
        if (variable_name != variable_names[i]) throw invalid_argument("Invalid config file. The variable " + variable_name + " is missing.");
        
        config_file >> values[i];
    }

    // create the config object using the values read from the file
    Config config {values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9]};

    return config;
}

/**
 * @brief Parses number of slits to 'simple', 'double' or 'triple'
 * 
 * @param slits The number of slits
 * @return string 'simple', 'double' or 'triple'
 */
string to_string(Slits slits) {
    switch (slits) {
        case Slits::one:   return "simple";
        case Slits::two:   return "double";
        case Slits::three: return "triple";
        default:           throw invalid_argument("Invalid number of slits");
    }
}

/** 
 * @brief Setup the system
 * 
 * @param slits The number of slits
 * @param config The setup of the system
 * @return Schrodinger_simulation The system
 */
void setup_and_run_system(Slits slits, int config_nr, Config config){
    SchrodingerSimulation ss(config, slits);
    ss.simulate("output/data/" + to_string(slits) + "_" + "config" + to_string(config_nr) + ".bin");
}

int main() {
    for (int i = 1; i <= 3; i++){
        // set the filename of the config file and parse it
        string config_file_name = "config" + to_string(i) + ".in";
        Config config = parse_config(config_file_name);
        cerr << "Running using config " + to_string(i) +": "; config.print();
        Slits slits;
        // run for all slits if config 3, else run for just two-slit
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
