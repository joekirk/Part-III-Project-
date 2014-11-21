//
//  core.cpp
//  Part III Project - version1
//
//  Created by Joe Kirk on 17/01/2013.
//  Copyright (c) 2013 Joe Kirk. All rights reserved.
//

#include "core.h"

namespace core {

    
string enumtostring(atmtype A)
    {
        switch (A) {
            case Fe:
                return "Fe";
            case Ge:
                return "Ge";
            case Sb:
                return "Sb";
            case Te:
                return "Te";
        }
        
    }
    
    string convertdouble(double number)
    {
        if (number == 3.0)
            return "3.0";
        if (number == 3.1)
            return "3.1";
        if (number == 3.2)
            return "3.2";
        if (number == 3.3)
            return "3.3";
        if (number == 3.4)
            return "3.4";
        if (number == 3.5)
            return "3.5";
        if (number == 3.22)
            return "3.22";
        
        else {
            cout << "Undefined cutoff";
            return "undefined";}
                
    }
        
    
    
void PrintAtomcoord (Atomcoord a)
{
    cout << '(' << a.x << ',' << a.y << ',' << a.z << ')' << endl;
}

Atomcoord operator+(Atomcoord a, Atomcoord b)
{
    return Atomcoord(a.x +b.x, a.y + b.y, a.z + b.z);
}

Atomcoord operator-(Atomcoord a, Atomcoord b)
{
    return Atomcoord(a.x - b.x, a.y - b.y, a.z - b.z);
}


//Determine the absolute size of the vector of atomic co-ordinates
double mod(Atomcoord a)
{
    return (sqrt((a.x*a.x) + (a.y*a.y) + (a.z*a.z)));
}

//returns the x component on the dot product of two Atomcoords
double dotproductx(Atomcoord a, Atomcoord b)
{
    return (a.x*b.x);
}

//returns the y component on the dot product of two Atomcoords
double dotproducty(Atomcoord a, Atomcoord b)
{
    return (a.y*b.y);
}

//returns the z component on the dot product of two Atomcoords
double dotproductz(Atomcoord a, Atomcoord b)
{
    return (a.z*b.z);
}

//returns the dot product of two Atomcoords
double dotproduct(Atomcoord a, Atomcoord b)
{
    return (a.x*b.x + a.y*b.y + a.z*b.z);
}

//returns the cross product of two Atomcoords
Atomcoord crossproduct(Atomcoord a, Atomcoord b)
{
    return Atomcoord (a.y*b.z -a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
    

double angledegrees(Atomcoord a, Atomcoord b)
{
    return (acos(dotproduct(a, b) / (mod(a) * mod(b))))*180.0/M_PI;
}
    
}
