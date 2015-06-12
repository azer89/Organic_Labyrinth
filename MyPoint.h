#ifndef MYPOINT_H
#define MYPOINT_H

#include <cmath>

struct MyPoint
{
public:
    // x
    float x;

    // y
    float y;

    // custom
    int index;

    int age;

    // Default constructor
    MyPoint()
    {
        this->x = -1;
        this->y = -1;
        this->index = -1;
        this->age = 0;
    }

    // Constructor
    MyPoint(float x, float y)
    {
        this->x = x;
        this->y = y;
        this->index = -1;
        this->age = 0;
    }

    // Scale a point
    MyPoint Resize(float val)
    {
        MyPoint newP;
        newP.x = this->x * val;
        newP.y = this->y * val;
        return newP;
    }

    // if a point is (-1, -1)
    bool Invalid()
    {
        if(((int)x) == -1 && ((int)y) == -1)
            { return true; }
        return false;
    }

    // Normalize
    MyPoint Norm()
    {
        float vlength = std::sqrt( x * x + y * y );
        return MyPoint(this->x / vlength, this->y / vlength);
    }

    // Euclidean distance
    float Distance(MyPoint other)
    {
        float xDist = x - other.x;
        float yDist = y - other.y;
        return sqrt(xDist * xDist + yDist * yDist);
    }

    // Euclidean distance
    float Distance(float otherX, float otherY)
    {
        float xDist = x - otherX;
        float yDist = y - otherY;
        return sqrt(xDist * xDist + yDist * yDist);
    }

    // squared euclidean distance
    float DistanceSquared(MyPoint other)
    {
        float xDist = x - other.x;
        float yDist = y - other.y;
        return (xDist * xDist + yDist * yDist);
    }

    // squared euclidean distance
    float DistanceSquared(float otherX, float otherY)
    {
        float xDist = x - otherX;
        float yDist = y - otherY;
        return (xDist * xDist + yDist * yDist);
    }

    // operator overloading
    MyPoint operator+ (const MyPoint& other) { return MyPoint(x + other.x, y + other.y); }

    // operator overloading
    MyPoint operator- (const MyPoint& other) { return MyPoint(x - other.x, y - other.y); }
    bool operator== (const MyPoint& other)
    { return (abs(this->x - other.x) < 1e-8 && abs(this->y - other.y) < 1e-8); }

    // operator overloading
    bool operator!= (const MyPoint& other)
    { return (abs(this->x - other.x) >= 1e-8 || abs(this->y - other.y) >= 1e-8); }

    // operator overloading
    MyPoint operator+= (const MyPoint& other)
    {
        x += other.x;
        y += other.y;
        return *this;
    }

    // operator overloading
    MyPoint operator-= (const MyPoint& other)
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    // operator overloading
    MyPoint operator/ (const float& val) { return MyPoint(x / val, y / val); }

    // operator overloading
    MyPoint operator* (const float& val) { return MyPoint(x * val, y * val); }

    // operator overloading
    MyPoint operator*= (const float& val)
    {
        x *= val;
        y *= val;
        return *this;
    }

    // operator overloading
    MyPoint operator/= (const float& val)
    {
        x /= val;
        y /= val;
        return *this;
    }

    // length of a vector
    float Length() { return sqrt(x * x + y * y); }

    // squared length of a vector
    float LengthSquared() { return x * x + y * y; }

    // dot product
    float Dot(MyPoint otherPt) { return x * otherPt.x + y * otherPt.y; }

    // cross product
    MyPoint Cross(MyPoint otherPt)
    {
        //U x V = Ux*Vy-Uy*Vx
        return MyPoint(x * otherPt.y, y * otherPt.x);
    }

    // linear dependency test
    bool IsLinearDependent(MyPoint otherPoint)
    {
        float det = (this->x * otherPoint.y) - (this->y * otherPoint.x);
        if(det > -1e-8 && det < 1e-8) { return true; }
        return false;
    }

    // angle direction
    MyPoint DirectionTo(MyPoint otherPt)
    {
        return MyPoint(otherPt.x - this->x, otherPt.y - this->y);
    }
};

#endif
