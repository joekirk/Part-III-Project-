//
//  main.cpp
//  Part III Project - version1
//
//  Created by Joe Kirk on 17/01/2013.
//  Copyright (c) 2013 Joe Kirk. All rights reserved.
//

#include <iostream>
#include "structuralanalysis.h"
#include "core.h"
#include "input.h"
#include <vector>
#include <complex>
#include "MatrixI:O.h"
#include <iomanip>
#include "fftw3.h"


using namespace core;
using namespace std;

int main()
try{
  /*  Fileread file;
    file.name = "Song_FeGST_7_200 - Quench.axsf";
    
    file.readdynamicanalysis(file.name, 1000);
    dynamicanalysis::nearestneighbourevolution(file); */

     //program to read a CONTCAR file and then produce it for VMD for the relaxed models
    
    DOSCARread file;
    file.readfile(2, "Fe_GST_1_Amorphous_200");
    
    return 0;
    
    
}
catch (exception& e) {
    cerr << e.what() << endl;
    return 1;
    }
    catch (...) {
        cerr << "exception \n";
        return 2;
    }
    


