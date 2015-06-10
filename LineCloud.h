
/**
 * Fast query of line by using KD-Tree
 *
 * Author: Reza Adhitya Saputra (reza.adhitya.saputra@gmail.com)
 * Version: 2014
 *
 */

#ifndef __Line_Cloud__
#define __Line_Cloud__

#include  <vector>


template <typename T>
struct LineCloud
{
    struct LineData
    {
        //// you can choose whether x,y are first, second, or mid point.
        T x;
        T y;

        //// first point
        int index0;

        //// second point
        int index1;

        //// Something 1
        int info1;

        //// Something 2
        int info2;
    };

    std::vector<LineData> lines;

    // for nanoflann internal use
    inline size_t kdtree_get_point_count() const { return lines.size(); }

    // for nanoflann internal use
    inline T kdtree_distance(const T *p1, const size_t idx_p2, size_t size) const
    {
        const T d0 = p1[0] - lines[idx_p2].x;
        const T d1 = p1[1] - lines[idx_p2].y;
        return sqrt(d0 * d0 + d1 * d1);
    }

    // for nanoflann internal use
    inline T kdtree_get_pt(const size_t idx, int dim) const
    {
        if (dim == 0) return lines[idx].x;
        else return lines[idx].y;
    }

    // for nanoflann internal use
    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const { return false; }
};

#endif
