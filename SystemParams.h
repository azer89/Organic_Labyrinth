#ifndef SYSTEMPARAMS_H
#define SYSTEMPARAMS_H

class SystemParams
{
public:
    SystemParams();
    ~SystemParams();

public:
    static int init_slices;
    static float radius;

    static float D;

    static float f_b;       // brownian
    static float f_f;       // fairing
    static float f_a;       // attraction-repulsion

    static int dist_std_dev;

    static float delta_l_j;
    static int search_a_r;
    static float radius_a_r;

    static int max_iter;    // debug, fix me
};

#endif // SYSTEMPARAMS_H
