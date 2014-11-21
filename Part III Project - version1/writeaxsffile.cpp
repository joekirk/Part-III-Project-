//
//  writeaxsffile.cpp
//  Part III Project - version1
//
//  Created by Joe Kirk on 30/01/2013.
//  Copyright (c) 2013 Joe Kirk. All rights reserved.
//

#include "writeaxsffile.h"


void write(Fileread&file) {
string name = "/Users/joekirk/Desktop/test.axsf";
ofstream ost(name.c_str());
if (!ost) cout << "cannot open file " << name;

ost << "AMNISTEPS " << file.frames.size() << endl;
ost << endl;
ost << "CRYSTAL" << endl;
ost << endl;
ost << "PRIMVEC" << endl;
for (int i = 0; i < file.cellvectors.size(); ++i)
{
    ost << file.cellvectors[i].x << "\t" << file.cellvectors[i].y << "\t" << file.cellvectors[i].z << endl;
}
for (int i = 0; i <file.frames.size(); ++i)
{
    ost << endl;
    ost << "PRIMCOORD " << i+1 << endl;
    ost << file.modelsize << " " << 1 << endl;;
    for (int j=0; j < file.modelsize; ++j)
        ost << " " << file.atomtypes[j] << "\t" << file.frames[i][j].x << "\t" <<file.frames[i][j].y <<"\t" << file.frames[i][j].z << endl;
}

}