#ifndef SYSTEMPARAMS_H
#define SYSTEMPARAMS_H

class SystemParams
{
public:
    SystemParams();
    ~SystemParams();

public:
    static int   circle_init_slices;
    static float circle_radius;

    static float D;
    static float d_split;
    //static float k_min;

    static float delta_const;

    static float f_b;       // brownian
    static float f_f;       // fairing
    static float f_a;       // attraction-repulsion

    static float ar_clamp;

    static int dist_std_dev;

    static float sigma_l_j;

    static float search_radius;

    static float kdtree_radius;

    //static int n_min;

    static int max_iter;    // debug, fix me
};

#endif // SYSTEMPARAMS_H
