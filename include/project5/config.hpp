#pragma once

#include <armadillo>
#include <complex>

struct Config {
    double h, dt, T, x_c, s_x, p_x, y_c, s_y, p_y, v_0;

    // initialize
    // void init(std::string filename) {
    void init(double h_, double dt_, double T_, double x_c_, double s_x_, double p_x_, double y_c_, double s_y_, double p_y_, double v_0_) {
        h = h_;
        dt = dt_;
        T = T_;
        x_c = x_c_;
        s_x = s_x_;
        p_x = p_x_;
        y_c = y_c_;
        s_y = s_y_;
        p_y = p_y_;
        v_0 = v_0_;
    }

    void print() {
        std::cout << "h: " << h << ", ";
        std::cout << "dt: " << dt << ", ";
        std::cout << "T: " << T << ", ";
        std::cout << "x_c: " << x_c << ", ";
        std::cout << "s_x: " << s_x << ", ";
        std::cout << "p_x: " << p_x << ", ";
        std::cout << "y_c: " << y_c << ", ";
        std::cout << "s_y: " << s_y << ", ";
        std::cout << "p_y: " << p_y << ", ";
        std::cout << "v_0: " << v_0 << std::endl;
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
