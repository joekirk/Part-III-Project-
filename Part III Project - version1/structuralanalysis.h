//
//  structural analysis.h
//  Part III Project - version1
//
//  Created by Joe Kirk on 17/01/2013.
//  Copyright (c) 2013 Joe Kirk. All rights reserved.
//

#ifndef __Part_III_Project___version1__structural_analysis__
#define __Part_III_Project___version1__structural_analysis__

#include "core.h"
#include "input.h"
#include "fftw3.h"



struct staticanalysis {
     
    //radial distribution function for atom A to all other atoms in the model
    static void  g(int firstframe,  int lastframe, int datapoints,  atmtype A,  Fileread& file);
    
    //total radial distribution function for the model
    static void  totalg(int firstframe,  int lastframe, int datapoints, Fileread& file);
    
    //outputs the radial distribution functions for all atom pairs in the system
    static void rdf (int firstframe,  int lastframe, int datapoints, Fileread& file);

    
    static void bondangledistributionsingleatom (int firstframe, int lastframe, int datapoints, atmtype A, Fileread& file, double cutoff);

    static void badf (int firstframe, int lastframe, int datapoints, Fileread& file, double cutoff);
    
    static void bondstatistics(Fileread file, int datapoints, int firstframe, int lastframe, double cutoff);
};

struct dynamicanalysis {
    static void msdforatom(atmtype A,int datapoints, Fileread& file);
    static void averageMSD(Fileread& file, int datapoints);
    static void msd(Fileread& file, int datapoints);
    
    //very expensive pieces of code, better to use readfordynamic analyis then the two below these two
    static void wrongbonds(Fileread& file, int datapoints, int firstframe, int lastframe, double cutoff);
    static void fourfoldrings(Fileread& file, int datapoints, int firstframe, int lastframe);
    
    //use when used readfordynamic analysis
    static void wrongbonds(Fileread& file, double cutoff);
    static void fourfoldrings(Fileread& file);
    
    static void VMDviewer(Fileread& file, int atomtocentre);
    
    static void maxfourieranalysis(Fileread& file);
    
    static void nearestneighbourevolution(Fileread& file);
};

#endif  /*defined(__Part_III_Project___version1__structural_analysis__)*/
