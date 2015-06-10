
#include "SystemParams.h"

#include <limits>

// slices of the initial circle
int SystemParams::circle_init_slices = 4; // 6 hexagon
// radius of the initial circle
float SystemParams::circle_radius = 15.0f;


// max length of each line segment
float SystemParams::D = 1.0f;       // 1.0
// split threshold
float SystemParams::k_max = 0.3f;  // 1.2
// merge threshold
float SystemParams::k_min = 0.3f;  // 0.25

// delta constant
float SystemParams::delta_const = 1.0f;

// brownian constant
float SystemParams::f_b = 0.08f;   // 0.085f

// fairing constant
float SystemParams::f_f = 0.06f; // 0.01f

// attraction-repulsion constant
float SystemParams::f_a = 0.0040f;

// standard deviation of the random distribution (brownian motion)
int SystemParams::dist_std_dev = 628; // 628 = 3.14 * 2 * 100

// Lennard-Jones constant
float SystemParams::delta_l_j = 1.0f; // 1.25

// search radius
float SystemParams::radius_1 = 20.0f;   // 100

float SystemParams::kdtree_radius = SystemParams::radius_1 + SystemParams::D;

int SystemParams::n_min = 0;

// max iteration
int SystemParams::max_iter = std::numeric_limits<int>::max();
