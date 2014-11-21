//
//  core.h
//  Part III Project - version1
//
//  Created by Joe Kirk on 17/01/2013.
//  Copyright (c) 2013 Joe Kirk. All rights reserved.
//

#ifndef __Part_III_Project___version1__core__
#define __Part_III_Project___version1__core__

#include <iostream>
#include <vector>
#include "math.h"
#include "Matrix.h"


using namespace std;

namespace core {
    
    enum atmtype
    {
        Fe = 26,
        Ge = 32,
        Sb = 51,
        Te = 52
    };
    

    struct Atomcoord {          //create a struct to hold all the atomic positions for each frame
        double x, y, z;
        int type;
        Atomcoord (double a, double b, double c) : x(a), y(b), z(c){}
    };
    
    struct DOS {          //create a struct to hold all the density of states both up and down
        double up, down;
        int type;
        DOS (double a, double b) : up(a), down(b) {}
    };

    
    string enumtostring(atmtype A);
    
    string convertdouble (double number);
   
    void PrintAtomcoord (Atomcoord a);
    
    Atomcoord operator+(Atomcoord a, Atomcoord b);
    
    Atomcoord operator-(Atomcoord a, Atomcoord b);
    
    double mod(Atomcoord a);
    
    double dotproductx(Atomcoord a, Atomcoord b);
    
    double dotproducty(Atomcoord a, Atomcoord b);
    
    double dotproductz(Atomcoord a, Atomcoord b);
    
    double dotproduct(Atomcoord a, Atomcoord b);
    
    Atomcoord crossproduct(Atomcoord a, Atomcoord b);
    
    double angledegrees(Atomcoord a, Atomcoord b);

}

#endif /* defined(__Part_III_Project___version1__core__) */
