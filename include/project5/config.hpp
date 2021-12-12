#pragma once

#include <armadillo>
#include <complex>

struct Config {
    double h, dt, T, x_c, y_c, s_x, s_y, p_x, p_y, v_0;

    // initialize
    // void init(std::string filename) {
    //     string variable_names[10] = {"h", "dt", "T", "x_c", "y_c", "s_x", "s_y", "p_x", "p_y", "v_0"}

    void print() {
        std::cerr << "h=" << h << ", ";
        std::cerr << "dt=" << dt << ", ";
        std::cerr << "T=" << T << ", ";
        std::cerr << "x_c=" << x_c << ", ";
        std::cerr << "y_c=" << y_c << ", ";
        std::cerr << "s_x=" << s_x << ", ";
        std::cerr << "s_y=" << s_y << ", ";
        std::cerr << "p_x=" << p_x << ", ";
        std::cerr << "p_y=" << p_y << ", ";
        std::cerr << "v_0=" << v_0 << std::endl;
    }
};

// Config config1 = {
//     .h = 0.005,
//     .dt = 2.5*1e-5,
//     .T = 0.008,
//     .x_c = 0.25,
//     .s_x = 0.05,
//     .p_x = 200.,
//     .y_c = 0.5,
//     .s_y = 0.05,
//     .p_y = 0,
//     .v_0 = 0
// };
//
// Config config2 = {
//     .h = 0.005,
//     .dt = 2.5*1e-5,
//     .T = 0.002,
//     .x_c = 0.25,
//     .s_x = 0.05,
//     .p_x = 200.,
//     .y_c = 0.5,
//     .s_y = 0.20,
//     .p_y = 0,
//     .v_0 = 1e10
// };
