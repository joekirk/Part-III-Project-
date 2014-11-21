//
//  input.h
//  Part III Project - version1
//
//  Created by Joe Kirk on 17/01/2013.
//  Copyright (c) 2013 Joe Kirk. All rights reserved.
//

#ifndef __Part_III_Project___version1__input__
#define __Part_III_Project___version1__input__

#include "core.h"
#include <fstream>

using namespace core;

class Fileread {             //create a class to read all of the frames
private:
    int dimn;
    double x, y, z;          //co-ordinates (x, y, z)
    int type;
    
    void readcoords(istream& ist);
    void readprimvec(istream& ist);
    void readframes(istream& ist);
    void readframes(istream& ist, int datapoints);
    
    vector<Atomcoord> data;
    vector<Atomcoord> rijoneatom;                      //interatomic distances for atom i
    vector< vector<Atomcoord> > rijallatoms;           //interatomic distances for all atoms in a single frame
    
    void boudaryconditions(int firstframe, int lastframe);
    void boundaryconditions();
    
    
    Numeric_lib::Matrix <double,2>  neighbourtable (int,int);
    void generateneighbourtable();
    void generateneighbourtable(int firstframe, int lastframe);
    
    
   

public:
    int framesinfile;
    int modelsize;
    vector < vector<Atomcoord> > frames;
    vector <Atomcoord> cellvectors;
    vector<int> atomtypes;
    string name;
    vector<vector< vector<Atomcoord> > > rijallframes; //all frames all interatomic distances
    vector < Numeric_lib::Matrix<double, 2> > neighbourtables;
    void read(string name, int firstframe, int lastframe);
    void read(string name);
    void readformsd(string name);
    void readdynamicanalysis(string name, int datapoints);
    
    
};

class DOSCARread {             //create a class to read all of the frames
private:
  
public:
    double sup, sdown, pxup, pxdown, pyup, pydown, pzup, pzdown, dxyup, dxydown, dyzup, dyzdown, dxzup,dxzdown, dx2y2up, dx2y2down, dz2up, dz2down;
    double E, DOSup, DOSdown;
    double highestE, lowestE, Ef;
    int datapoints;
    vector<double> energy;
    vector<DOS> densityofstates;
    vector<DOS> s, p , d, dxy, dxz, dyz, dx2y2, dz2, t2g, eg;
    vector< vector <DOS> > Fes, Fep, Fed, Fet2g, Feeg, Fedxy, Fedxz, Fedyz, Fedx2y2, Fedz2;
    vector <DOS> totalFes, totalFep, totalFed;
    void readtotalDOS(istream& ist);
    void readFeDOS(istream& ist, int NFe);
    void readfile(int NFe, string outputfolder);

};





#endif /* defined(__Part_III_Project___version1__input__) */
