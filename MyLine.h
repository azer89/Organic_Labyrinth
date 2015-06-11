
/**
 *
 * Line representation with two points, start and end points
 * Each points doesn't hold the real position
 *
 * Author: Reza Adhitya Saputra (reza.adhitya.saputra@gmail.com)
 * Version: 2014
 *
 *
 */

#ifndef __MYLINE_H__
#define __MYLINE_H__

#include "MyPoint.h"

//namespace CVSystem
//{

    struct MyLine
    {
    public:
        float XA;	float YA;	// start
        float XB;	float YB;	// end

        // custom
        int index1;
        int index2;

        // Constructor #1
        MyLine()
        {
            this->XA = -1;	this->YA = -1;
            this->XB = -1;	this->YB = -1;

            this->index1 = -1;
            this->index2 = -1;
        }

        // Constructor #2
        MyLine(float XA, float YA, float XB, float YB)
        {
            this->XA = XA;	this->YA = YA;
            this->XB = XB;	this->YB = YB;

            this->index1 = -1;
            this->index2 = -1;
        }

        MyLine(MyPoint pt1, MyPoint pt2)
        {
            this->XA = pt1.x;	this->YA = pt1.y;
            this->XB = pt2.x;	this->YB = pt2.y;

            this->index1 = -1;
            this->index2 = -1;
        }


        MyLine Resize(float val)
        {
            MyLine newL;

            newL.XA = this->XA * val;
            newL.YA = this->YA * val;

            newL.XB = this->XB * val;
            newL.YB = this->YB * val;

            return newL;
        }

        // Uninitialized ?
        bool Invalid()
        {
            if(((int)XA) == -1 && ((int)YA) == -1 && ((int)XB) == -1 && ((int)YB) == -1) return true;
            return false;
        }

        // Start point
        MyPoint GetPointA() { return MyPoint(XA, YA); }

        // End point
        MyPoint GetPointB() { return MyPoint(XB, YB); }

        // Direction of the vector
        MyPoint Direction() { return MyPoint(XB - XA, YB - YA);}

        float Magnitude()
        {
            MyPoint dir = Direction();
            return sqrt(dir.x * dir.x + dir.y * dir.y);
        }

        bool LiesHere(MyPoint pt)
        {
            float det = (XB - XA) * (pt.y - YA) - (YB - YA) * (pt.x - XA);
            if(det > -1e-8 && det < 1e-8) return true;
            return false;
        }

        // returns:
        //	 1 : same direction
        //  -1 : opposite direction
        //   0 : other else
        int HasSameDirection(MyLine otherLine)
        {
            float mag1 = Magnitude();
            float mag2 = otherLine.Magnitude();

            MyPoint dir1 = Direction();
            MyPoint dir2 = otherLine.Direction();

            float a_dot_b = dir1.Dot(dir2);
            float a_b_mag = mag1 *  mag2;

            float addValue = a_dot_b + a_b_mag;
            if(addValue > -1e-8 && addValue < 1e-8 ) { return -1; }

            float subsValue = a_dot_b - a_b_mag;

            if(subsValue > -1e-8 && subsValue < 1e-8 ) { return 1; }

            return 0;
        }
    };
//}

#endif
