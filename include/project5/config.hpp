#pragma once

#include <armadillo>
#include <complex>

struct Config {
    /** Contains the config variables */
    double h, dt, T, x_c, y_c, s_x, s_y, p_x, p_y, v_0;

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
