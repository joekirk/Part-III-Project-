//
//  input.cpp
//  Part III Project - version1
//
//  Created by Joe Kirk on 17/01/2013.
//  Copyright (c) 2013 Joe Kirk. All rights reserved.
//

#include "input.h"
#include <sstream>

/*--Reads a frame of atomic co-ordinates into a vector called frames -- the files must be of format
 *  atom type   x-co-ordinate   y-co-ordinate   z-co-ordinate
 *  26          0.05            0.87            0.95
 *  32          0.76            0.34            0.12
 *  32          0.54            0.39            0.46
 *  .           .               .               .
 *  .           .               .               .
 *  .           .               .               .
 */

void DOSCARread::readtotalDOS(istream& ist)
{
        double discard = 0;
       
        //read in the header file
        ist >> highestE >> lowestE >> datapoints  >> Ef >>  discard;
            
            //read in the total density of states and place in the vector
            for (int i = 0; i < datapoints; ++i)
            {
                ist >> E >> DOSup >> DOSdown >> discard >> discard;
                energy.push_back(E);
                densityofstates.push_back(DOS(DOSup, DOSdown));
            }
}

void DOSCARread::readFeDOS(istream& ist, int NFe)
{
    double discard = 0;
    
    for (int j = 0; j < NFe; ++j) {
        
        ist >> highestE >> lowestE >> datapoints  >> Ef >>  discard;
        
        for (int i = 0; i < datapoints; ++i)
        {
            //ist >> discard >> sup >> sdown >> pyup >> pydown >> pzup >> pzdown >> pxup >> pxdown
                           //>> dxyup >> dxydown >> dyzup >> dyzdown >> dxzup >> dxzdown >> dx2y2up >> dx2y2down >> dz2up >> dz2down;
            
            ist >> discard >> sup >> sdown >> pyup >> pydown >> pzup >> pzdown >> pxup >> pxdown
            >> dxyup >> dxydown >> dyzup >> dyzdown >> dz2up >> dz2down >> dxzup >> dxzdown >> dx2y2up >> dx2y2down;
            
            s.push_back(DOS(sup, sdown));
            p.push_back(DOS(pyup+pzup+pxup, pydown +pzdown + pxdown));
            d.push_back(DOS(dxyup+dyzup+dxzup+ dx2y2up+dz2up, dxydown + dyzdown + dxzdown + dx2y2down + dz2down));
            dxy.push_back(DOS(dxyup,dxydown));
            dxz.push_back(DOS(dxzup,dxzdown));
            dyz.push_back(DOS(dyzup,dyzdown));
            dx2y2.push_back(DOS(dx2y2up,dx2y2down));
            dz2.push_back(DOS(dz2up,dz2down));
            t2g.push_back(DOS(dxyup + dxzup+ dyzup, dxydown + dyzdown + dxzdown));
            eg.push_back(DOS(dx2y2up + dz2up, dx2y2down + dz2down));
            
        }
        
        Fes.push_back(s);
        Fep.push_back(p);
        Fed.push_back(d);
        Fedxy.push_back(dxy);
        Fedxz.push_back(dxz);
        Fedyz.push_back(dyz);
        Fedx2y2.push_back(dx2y2);
        Fedz2.push_back(dz2);
        Fet2g.push_back(t2g);
        Feeg.push_back(eg);


        s.clear();
        p.clear();
        d.clear();
        dxy.clear();
        dxz.clear();
        dyz.clear();
        dx2y2.clear();
        dz2.clear();
        t2g.clear();
        eg.clear();
        
    }
    
    //initialise the vectors to 0
    for (int i = 0; i < datapoints; ++i)
    {
        totalFes.push_back(DOS(0,0));
        totalFep.push_back(DOS(0,0));
        totalFed.push_back(DOS(0,0));
        
    }
 
    for (int i = 0; i < datapoints; ++i)
    {
        for (int j = 0; j < NFe; ++j) {
            totalFes[i].up += Fes[j][i].up;
            totalFes[i].down += Fes[j][i].down;
            totalFep[i].up += Fep[j][i].up;
            totalFep[i].down += Fep[j][i].down;
            totalFed[i].up += Fed[j][i].up;
            totalFed[i].down += Fed[j][i].down;
            
        }
        
    }
    
    
}

  
    


void DOSCARread::readfile(int NFe, string outputfolder)
{
    
    string nameaccess= "/Users/joekirk/Desktop/DOSCAR.txt";
    ifstream ist (nameaccess.c_str());
    if (!ist) cout << "cannot open file " << "DOSCAR";
    
    string FeDoped_GST = " Fe-Doped GST";
    string CAR = "  CAR";

    while (!ist.eof()) {
        
        string line;
        getline(ist, line);
        
        if (line[1] == FeDoped_GST[1] && line[3] == FeDoped_GST[3] && line[5] == FeDoped_GST[5] && line[7] == FeDoped_GST[7]){
            readtotalDOS(ist);
            readFeDOS(ist, NFe);
            }
}
   
    string name = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_DOS.csv";
    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    ost << "Model: " << outputfolder << "\n";
    ost << "*CHECK* #Fe: ," << NFe << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
    ost << "E, E - Ef (UNITS?), TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down\n";

    for (int i = 0; i < datapoints; ++i){
        ost << energy[i] << "," << energy[i] - Ef << "," << densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "\n";
    }
    ost.close();
        
    if (NFe == 2) {
       
        string name0 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe0DOS.csv";
        ofstream ost0(name0.c_str());
        if (!ost0) cout << "cannot open file " << name0;
        
        ost0 << "Model: " << outputfolder << "\n";
        ost0 << "Fe index ," << 0 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
        ost0 << "E, E - Ef (UNITS?), Fe0 DOS s - up, Fe0 DOS s - down, Fe0 DOS p - up, Fe0 DOS p - down, Fe0 DOS d - up, Fe0 DOS d - down, (+)Fe0 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe0 dxy up, Fe0 dxy down, Fe0 dxz up, Fe0 dxz down, Fe0 dyz uo, Fe0 dyz down, Fe0 dx2-y2 up, Fe0 dx2-y2 down, Fe0 dz2 up, Fe0 dz2 down, Fe0 t2g up, Fe0 t2g down, Fe0 eg up, Fe0 eg down\n";
        
        for (int i = 0; i < datapoints; ++i){
            ost0 << energy[i] << "," << energy[i] - Ef << "," << Fes[0][i].up <<  "," << -Fes[0][i].down << "," << Fep[0][i].up <<  "," << -Fep[0][i].down << "," << Fed[0][i].up << "," << -Fed[0][i].down << "," <<  Fed[0][i].down << "," <<  densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[0][i].up << "," << -Fedxy[0][i].down << "," << Fedxz[0][i].up << "," << -Fedxz[0][i].down << "," << Fedyz[0][i].up << "," << -Fedyz[0][i].down << "," << Fedx2y2[0][i].up << "," << -Fedx2y2[0][i].down << "," << Fedz2[0][i].up << "," << -Fedz2[0][i].down << "," << Fet2g[0][i].up << "," << -Fet2g[0][i].down << "," << Feeg[0][i].up << "," << -Feeg[0][i].down <<  "\n";
        }
        
        ost0.close();
        
        
        string name1 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe1DOS.csv";
        ofstream ost1(name1.c_str());
        if (!ost1) cout << "cannot open file " << name1;
        
        ost1 << "Model: " << outputfolder << "\n";
        ost1 << "Fe index ," << 1 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
       ost1 << "E, E - Ef (UNITS?), Fe1 DOS s - up, Fe1 DOS s - down, Fe1 DOS p - up, Fe1DOS p - down, Fe1 DOS d - up, Fe1 DOS d - down, (+)Fe1 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe1 dxy up, Fe1 dxy down, Fe1 dxz up, Fe1 dxz down, Fe1 dyz uo, Fe1 dyz down, Fe1 dx2-y2 up, Fe1 dx2-y2 down, Fe1 dz2 up, Fe1 dz2 down, Fe1 t2g up, Fe1 t2g down, Fe1 eg up, Fe1 eg down\n";
        for (int i = 0; i < datapoints; ++i){
            ost1 << energy[i] << "," << energy[i] - Ef << "," << Fes[1][i].up <<  "," << -Fes[1][i].down << "," << Fep[1][i].up <<  "," << -Fep[1][i].down << "," << Fed[1][i].up << "," << -Fed[1][i].down << ","  << Fed[1][i].down << "," << densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[1][i].up << "," << -Fedxy[1][i].down << "," << Fedxz[1][i].up << "," << -Fedxz[1][i].down << "," << Fedyz[1][i].up << "," << -Fedyz[1][i].down << "," << Fedx2y2[1][i].up << "," << -Fedx2y2[1][i].down << "," << Fedz2[1][i].up << "," << -Fedz2[1][i].down << "," << Fet2g[1][i].up << "," << -Fet2g[1][i].down << "," << Feeg[1][i].up << "," << -Feeg[1][i].down <<  "\n";
        
        }
        
         ost1.close();
    }
    
    if (NFe == 8) {
        {
            
            string name0 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe0DOS.csv";
            ofstream ost0(name0.c_str());
            if (!ost0) cout << "cannot open file " << name0;
            
            ost0 << "Model: " << outputfolder << "\n";
            ost0 << "Fe index ," << 0 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
            ost0 << "E, E - Ef (UNITS?), Fe0 DOS s - up, Fe0 DOS s - down, Fe0 DOS p - up, Fe0 DOS p - down, Fe0 DOS d - up, Fe0 DOS d - down, (+)Fe0 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe0 dxy up, Fe0 dxy down, Fe0 dxz up, Fe0 dxz down, Fe0 dyz uo, Fe0 dyz down, Fe0 dx2-y2 up, Fe0 dx2-y2 down,Fe0 dz2 up, Fe0 dz2 down, Fe0 t2g up, Fe0 t2g down, Fe0 eg up, Fe0 eg down\n";
            
            for (int i = 0; i < datapoints; ++i){
                ost0 << energy[i] << "," << energy[i] - Ef << "," << Fes[0][i].up <<  "," << -Fes[0][i].down << "," << Fep[0][i].up <<  "," << -Fep[0][i].down << "," << Fed[0][i].up << "," << -Fed[0][i].down << "," <<  Fed[0][i].down << "," <<  densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[0][i].up << "," << -Fedxy[0][i].down << "," << Fedxz[0][i].up << "," << -Fedxz[0][i].down << "," << Fedyz[0][i].up << "," << -Fedyz[0][i].down << "," << Fedx2y2[0][i].up << "," << -Fedx2y2[0][i].down << "," << Fedz2[0][i].up << "," << -Fedz2[0][i].down << "," << Fet2g[0][i].up << "," << -Fet2g[0][i].down << "," << Feeg[0][i].up << "," << -Feeg[0][i].down <<  "\n";
            }
            
            ost0.close();
            
            
            string name1 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe1DOS.csv";
            ofstream ost1(name1.c_str());
            if (!ost1) cout << "cannot open file " << name1;
            
            ost1 << "Model: " << outputfolder << "\n";
            ost1 << "Fe index ," << 1 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
            ost1 << "E, E - Ef (UNITS?), Fe1 DOS s - up, Fe1 DOS s - down, Fe1 DOS p - up, Fe1DOS p - down, Fe1 DOS d - up, Fe1 DOS d - down, (+)Fe1 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe1 dxy up, Fe1 dxy down, Fe1 dxz up, Fe1 dxz down, Fe1 dyz uo, Fe1 dyz down, Fe1 dx2-y2 up, Fe1 dx2-y2 down,Fe1 dz2 up, Fe1 dz2 down, Fe1 t2g up, Fe1 t2g down, Fe1 eg up, Fe1 eg down\n";
            for (int i = 0; i < datapoints; ++i){
                ost1 << energy[i] << "," << energy[i] - Ef << "," << Fes[1][i].up <<  "," << -Fes[1][i].down << "," << Fep[1][i].up <<  "," << -Fep[1][i].down << "," << Fed[1][i].up << "," << -Fed[1][i].down << ","  << Fed[1][i].down << "," << densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[1][i].up << "," << -Fedxy[1][i].down << "," << Fedxz[1][i].up << "," << -Fedxz[1][i].down << "," << Fedyz[1][i].up << "," << -Fedyz[1][i].down << "," << Fedx2y2[1][i].up << "," << -Fedx2y2[1][i].down << "," << Fedz2[1][i].up << "," << -Fedz2[1][i].down << "," << Fet2g[1][i].up << "," << -Fet2g[1][i].down << "," << Feeg[1][i].up << "," << -Feeg[1][i].down <<  "\n";
                
            }
            
            ost1.close();
            
            
            string name2 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe2DOS.csv";
            ofstream ost2(name2.c_str());
            if (!ost2) cout << "cannot open file " << name2;
            
            ost2 << "Model: " << outputfolder << "\n";
            ost2 << "Fe index ," << 2 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
            ost2 << "E, E - Ef (UNITS?), Fe2 DOS s - up, Fe2 DOS s - down, Fe2 DOS p - up, Fe2DOS p - down, Fe2 DOS d - up, Fe2 DOS d - down, (+)Fe2 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe2 dxy up, Fe2 dxy down, Fe2 dxz up, Fe2 dxz down, Fe2 dyz up, Fe2dyz down, Fe2 dx2-y2 up, Fe2 dx2-y2 down, Fe2 dz2 up, Fe2 dz2 down, Fe2 t2g up, Fe2 t2g down, Fe2 eg up, Fe2 eg down\n";
            for (int i = 0; i < datapoints; ++i){
                ost2 << energy[i] << "," << energy[i] - Ef << "," << Fes[2][i].up <<  "," << -Fes[2][i].down << "," << Fep[2][i].up <<  "," << -Fep[2][i].down << "," << Fed[2][i].up << "," << -Fed[2][i].down << "," << Fed[2][i].down << "," <<   densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[2][i].up << "," << -Fedxy[2][i].down << "," << Fedxz[2][i].up << "," << -Fedxz[2][i].down << "," << Fedyz[2][i].up << "," << -Fedyz[2][i].down << "," << Fedx2y2[2][i].up << "," << -Fedx2y2[2][i].down << "," << Fedz2[2][i].up << "," << -Fedz2[2][i].down << "," << Fet2g[2][i].up << "," << -Fet2g[2][i].down << "," << Feeg[2][i].up << "," << -Feeg[2][i].down <<  "\n";
            }
             ost2.close();
            
            string name3 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe3DOS.csv";
            ofstream ost3(name3.c_str());
            if (!ost3) cout << "cannot open file " << name3;
            
            ost3 << "Model: " << outputfolder << "\n";
            ost3 << "Fe index ," << 3 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
           ost3 << "E, E - Ef (UNITS?), Fe3 DOS s - up, Fe3 DOS s - down, Fe3 DOS p - up, Fe3DOS p - down, Fe3 DOS d - up, Fe3 DOS d - down, (+)Fe3 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe3 dxy up, Fe3 dxy down, Fe3 dxz up, Fe3 dxz down, Fe3 dyz up, Fe3 dyz down, Fe3 dx2-y2 up, Fe3 dx2-y2 down, Fe3 dz2 up, Fe3 dz2 down,Fe3 t2g up, Fe3 t2g down, Fe3 eg up, Fe3 eg down\n";
            
            for (int i = 0; i < datapoints; ++i){
                ost3 << energy[i] << "," << energy[i] - Ef << "," << Fes[3][i].up <<  "," << -Fes[3][i].down << "," << Fep[3][i].up <<  "," << -Fep[3][i].down << "," << Fed[3][i].up << "," << -Fed[3][i].down << "," << Fed[3][i].down <<  "," <<  densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[3][i].up << "," << -Fedxy[3][i].down << "," << Fedxz[3][i].up << "," << -Fedxz[3][i].down << "," << Fedyz[3][i].up << "," << -Fedyz[3][i].down << "," << Fedx2y2[3][i].up << "," << -Fedx2y2[3][i].down << "," << Fedz2[3][i].up << "," << -Fedz2[3][i].down << "," << Fet2g[3][i].up << "," << -Fet2g[3][i].down << "," << Feeg[3][i].up << "," << -Feeg[3][i].down <<  "\n";
            }
             ost3.close();
            
            
            string name4 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe4DOS.csv";
            ofstream ost4(name4.c_str());
            if (!ost4) cout << "cannot open file " << name4;
            
            ost4 << "Model: " << outputfolder << "\n";
            ost4 << "Fe index ," << 4 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
            ost4 << "E, E - Ef (UNITS?), Fe4 DOS s - up, Fe4 DOS s - down, Fe4 DOS p - up, Fe4 DOS p - down, Fe4 DOS d - up, Fe4 DOS d - down, (+)Fe4 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe4 dxy up, Fe4 dxy down, Fe4 dxz up, Fe4 dxz down, Fe4 dyz up, Fe4 dyz down, Fe4 dx2-y2 up, Fe4 dx2-y2 down, Fe4 dz2 up, Fe4 dz2 down,Fe4 t2g up, Fe4 t2g down, Fe4 eg up, Fe4 eg down\n";
            
            for (int i = 0; i < datapoints; ++i){
                ost4 << energy[i] << "," << energy[i] - Ef << "," << Fes[4][i].up <<  "," << -Fes[4][i].down << "," << Fep[4][i].up <<  "," << -Fep[4][i].down << "," << Fed[4][i].up << "," << -Fed[4][i].down << "," << Fed[4][i].down <<  "," <<   densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down <<  "," << Fedxy[4][i].up << "," << -Fedxy[4][i].down << "," << Fedxz[4][i].up << "," << -Fedxz[4][i].down << "," << Fedyz[4][i].up << "," << -Fedyz[4][i].down << "," << Fedx2y2[4][i].up << "," << -Fedx2y2[4][i].down << "," << Fedz2[4][i].up << "," << -Fedz2[4][i].down << "," << Fet2g[4][i].up << "," << -Fet2g[4][i].down << "," << Feeg[4][i].up << "," << -Feeg[4][i].down <<  "\n";

            }
             ost4.close();
            
            
            string name5 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe5DOS.csv";
            ofstream ost5(name5.c_str());
            if (!ost5) cout << "cannot open file " << name5;
            
            ost5 << "Model: " << outputfolder << "\n";
            ost5 << "Fe index ," << 5 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
            ost5 << "E, E - Ef (UNITS?), Fe5 DOS s - up, Fe5 DOS s - down, Fe5 DOS p - up, Fe5 DOS p - down, Fe5 DOS d - up, Fe5 DOS d - down, (+)Fe5 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe5 dxy up, Fe5 dxy down, Fe5 dxz up, Fe5 dxz down, Fe5 dyz up, Fe5 dyz down, Fe5 dx2-y2 up, Fe5 dx2-y2 down, Fe5 dz2 up, Fe5 dz2 down,Fe5 t2g up, Fe5 t2g down, Fe5 eg up, Fe5 eg down\n";

            
            for (int i = 0; i < datapoints; ++i){
                ost5 << energy[i] << "," << energy[i] - Ef << "," << Fes[5][i].up <<  "," << -Fes[5][i].down << "," << Fep[5][i].up <<  "," << -Fep[5][i].down << "," << Fed[5][i].up << "," << -Fed[5][i].down << "," << Fed[5][i].down <<  "," <<   densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down <<  "," << Fedxy[5][i].up << "," << -Fedxy[5][i].down << "," << Fedxz[5][i].up << "," << -Fedxz[5][i].down << "," << Fedyz[5][i].up << "," << -Fedyz[5][i].down << "," << Fedx2y2[5][i].up << "," << -Fedx2y2[5][i].down << "," << Fedz2[5][i].up << "," << -Fedz2[5][i].down << "," << Fet2g[5][i].up << "," << -Fet2g[5][i].down << "," << Feeg[5][i].up << "," << -Feeg[5][i].down <<  "\n";
            }
            
             ost5.close();
            
            string name6 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe6DOS.csv";
            ofstream ost6(name6.c_str());
            if (!ost6) cout << "cannot open file " << name6;
            
            ost6 << "Model: " << outputfolder << "\n";
            ost6 << "Fe index ," << 6 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
            ost6 << "E, E - Ef (UNITS?), Fe6 DOS s - up, Fe6 DOS s - down, Fe6 DOS p - up, Fe6 DOS p - down, Fe6 DOS d - up, Fe6 DOS d - down, (+)Fe6 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe6 dxy up, Fe6 dxy down, Fe6 dxz up, Fe6 dxz down, Fe6 dyz up, Fe6 dyz down, Fe6 dx2-y2 up, Fe6 dx2-y2 down, Fe6 dz2 up, Fe6 dz2 down, Fe6 t2g up, Fe6 t2g down, Fe6 eg up, Fe6 eg down\n";
            
            for (int i = 0; i < datapoints; ++i){
                ost6 << energy[i] << "," << energy[i] - Ef << "," << Fes[6][i].up <<  "," << -Fes[6][i].down << "," << Fep[6][i].up <<  "," << -Fep[6][i].down << "," << Fed[6][i].up << "," << -Fed[6][i].down << "," << Fed[6][i].down << "," <<   densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down <<  "," << Fedxy[6][i].up << "," << -Fedxy[6][i].down << "," << Fedxz[6][i].up << "," << -Fedxz[6][i].down << "," << Fedyz[6][i].up << "," << -Fedyz[6][i].down << "," << Fedx2y2[6][i].up << "," << -Fedx2y2[6][i].down << "," << Fedz2[6][i].up << "," << -Fedz2[6][i].down << "," << Fet2g[6][i].up << "," << -Fet2g[6][i].down << "," << Feeg[6][i].up << "," << -Feeg[6][i].down <<  "\n";
            }
            
             ost6.close();
            string name7 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe7DOS.csv";
            ofstream ost7(name7.c_str());
            if (!ost7) cout << "cannot open file " << name7;
            
            ost7 << "Model: " << outputfolder << "\n";
            ost7 << "Fe index ," << 7 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
            ost7 << "E, E - Ef (UNITS?), Fe7 DOS s - up, Fe7 DOS s - down, Fe7 DOS p - up, Fe7 DOS p - down, Fe7 DOS d - up, Fe7 DOS d - down, (+)Fe7 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe7 dxy up, Fe7 dxy down, Fe7 dxz up, Fe7 dxz down, Fe7 dyz up, Fe7 dyz down, Fe7 dx2-y2 up, Fe7 dx2-y2 down, Fe7 dz2 up, Fe7 dz2 down,Fe7 t2g up, Fe7 t2g down, Fe7 eg up, Fe7 eg down\n";
            
            for (int i = 0; i < datapoints; ++i){
                ost7 << energy[i] << "," << energy[i] - Ef << "," << Fes[7][i].up <<  "," << -Fes[7][i].down << "," << Fep[7][i].up <<  "," << -Fep[7][i].down << "," << Fed[7][i].up << "," << -Fed[7][i].down << "," << Fed[7][i].down << ","  <<  densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[7][i].up << "," << -Fedxy[7][i].down << "," << Fedxz[7][i].up << "," << -Fedxz[7][i].down << "," << Fedyz[7][i].up << "," << -Fedyz[7][i].down << "," << Fedx2y2[7][i].up << "," << -Fedx2y2[7][i].down << "," << Fedz2[7][i].up << "," << -Fedz2[7][i].down << "," << Fet2g[7][i].up << "," << -Fet2g[7][i].down << "," << Feeg[7][i].up << "," << -Feeg[7][i].down <<  "\n";
            }
            
             ost7.close();
            
            
            
            
            
            
        }

    }
    
    if (NFe == 14)
    {
        
        string name0 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe0DOS.csv";
        ofstream ost0(name0.c_str());
        if (!ost0) cout << "cannot open file " << name0;
        
        ost0 << "Model: " << outputfolder << "\n";
        ost0 << "Fe index ," << 0 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
        ost0 << "E, E - Ef (UNITS?), Fe0 DOS s - up, Fe0 DOS s - down, Fe0 DOS p - up, Fe0 DOS p - down, Fe0 DOS d - up, Fe0 DOS d - down, (+)Fe0 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe0 dxy up, Fe0 dxy down, Fe0 dxz up, Fe0 dxz down, Fe0 dyz uo, Fe0 dyz down, Fe0 dx2-y2 up, Fe0 dx2-y2 down, Fe0 dz2 up, Fe0 dz2 down, Fe0 t2g up, Fe0 t2g down, Fe0 eg up, Fe0 eg down\n";
        
        for (int i = 0; i < datapoints; ++i){
            ost0 << energy[i] << "," << energy[i] - Ef << "," << Fes[0][i].up <<  "," << -Fes[0][i].down << "," << Fep[0][i].up <<  "," << -Fep[0][i].down << "," << Fed[0][i].up << "," << -Fed[0][i].down << "," <<  Fed[0][i].down << "," <<  densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[0][i].up << "," << -Fedxy[0][i].down << "," << Fedxz[0][i].up << "," << -Fedxz[0][i].down << "," << Fedyz[0][i].up << "," << -Fedyz[0][i].down << "," << Fedx2y2[0][i].up << "," << -Fedx2y2[0][i].down << "," << Fedz2[0][i].up << "," << -Fedz2[0][i].down << "," << Fet2g[0][i].up << "," << -Fet2g[0][i].down << "," << Feeg[0][i].up << "," << -Feeg[0][i].down <<  "\n";
        }
        
        ost0.close();
        
        
        string name1 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe1DOS.csv";
        ofstream ost1(name1.c_str());
        if (!ost1) cout << "cannot open file " << name1;
        
        ost1 << "Model: " << outputfolder << "\n";
        ost1 << "Fe index ," << 1 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
        ost1 << "E, E - Ef (UNITS?), Fe1 DOS s - up, Fe1 DOS s - down, Fe1 DOS p - up, Fe1DOS p - down, Fe1 DOS d - up, Fe1 DOS d - down, (+)Fe1 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe1 dxy up, Fe1 dxy down, Fe1 dxz up, Fe1 dxz down, Fe1 dyz uo, Fe1 dyz down, Fe1 dx2-y2 up, Fe1 dx2-y2 down, Fe1 dz2 up, Fe1 dz2 down, Fe1 t2g up, Fe1 t2g down, Fe1 eg up, Fe1 eg down\n";
        for (int i = 0; i < datapoints; ++i){
            ost1 << energy[i] << "," << energy[i] - Ef << "," << Fes[1][i].up <<  "," << -Fes[1][i].down << "," << Fep[1][i].up <<  "," << -Fep[1][i].down << "," << Fed[1][i].up << "," << -Fed[1][i].down << ","  << Fed[1][i].down << "," << densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[1][i].up << "," << -Fedxy[1][i].down << "," << Fedxz[1][i].up << "," << -Fedxz[1][i].down << "," << Fedyz[1][i].up << "," << -Fedyz[1][i].down << "," << Fedx2y2[1][i].up << "," << -Fedx2y2[1][i].down << "," << Fedz2[1][i].up << "," << -Fedz2[1][i].down << "," << Fet2g[1][i].up << "," << -Fet2g[1][i].down << "," << Feeg[1][i].up << "," << -Feeg[1][i].down <<  "\n";
            
        }
        
        ost1.close();
        
        
        string name2 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe2DOS.csv";
        ofstream ost2(name2.c_str());
        if (!ost2) cout << "cannot open file " << name2;
        
        ost2 << "Model: " << outputfolder << "\n";
        ost2 << "Fe index ," << 2 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
        ost2 << "E, E - Ef (UNITS?), Fe2 DOS s - up, Fe2 DOS s - down, Fe2 DOS p - up, Fe2DOS p - down, Fe2 DOS d - up, Fe2 DOS d - down, (+)Fe2 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe2 dxy up, Fe2 dxy down, Fe2 dxz up, Fe2 dxz down, Fe2 dyz up, Fe2dyz down, Fe2 dx2-y2 up, Fe2 dx2-y2 down, Fe2 dz2 up, Fe2 dz2 down,Fe2 t2g up, Fe2 t2g down, Fe2 eg up, Fe2 eg down\n";
        for (int i = 0; i < datapoints; ++i){
            ost2 << energy[i] << "," << energy[i] - Ef << "," << Fes[2][i].up <<  "," << -Fes[2][i].down << "," << Fep[2][i].up <<  "," << -Fep[2][i].down << "," << Fed[2][i].up << "," << -Fed[2][i].down << "," << Fed[2][i].down << "," <<   densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[2][i].up << "," << -Fedxy[2][i].down << "," << Fedxz[2][i].up << "," << -Fedxz[2][i].down << "," << Fedyz[2][i].up << "," << -Fedyz[2][i].down << "," << Fedx2y2[2][i].up << "," << -Fedx2y2[2][i].down << "," << Fedz2[2][i].up << "," << -Fedz2[2][i].down << "," << Fet2g[2][i].up << "," << -Fet2g[2][i].down << "," << Feeg[2][i].up << "," << -Feeg[2][i].down <<  "\n";
        }
        ost2.close();
        
        string name3 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe3DOS.csv";
        ofstream ost3(name3.c_str());
        if (!ost3) cout << "cannot open file " << name3;
        
        ost3 << "Model: " << outputfolder << "\n";
        ost3 << "Fe index ," << 3 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
        ost3 << "E, E - Ef (UNITS?), Fe3 DOS s - up, Fe3 DOS s - down, Fe3 DOS p - up, Fe3DOS p - down, Fe3 DOS d - up, Fe3 DOS d - down, (+)Fe3 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe3 dxy up, Fe3 dxy down, Fe3 dxz up, Fe3 dxz down, Fe3 dyz up, Fe3 dyz down, Fe3 dx2-y2 up, Fe3 dx2-y2 down, Fe3 dz2 up, Fe3 dz2 down,Fe3 t2g up, Fe3 t2g down, Fe3 eg up, Fe3 eg down\n";
        
        for (int i = 0; i < datapoints; ++i){
            ost3 << energy[i] << "," << energy[i] - Ef << "," << Fes[3][i].up <<  "," << -Fes[3][i].down << "," << Fep[3][i].up <<  "," << -Fep[3][i].down << "," << Fed[3][i].up << "," << -Fed[3][i].down << "," << Fed[3][i].down <<  "," <<  densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[3][i].up << "," << -Fedxy[3][i].down << "," << Fedxz[3][i].up << "," << -Fedxz[3][i].down << "," << Fedyz[3][i].up << "," << -Fedyz[3][i].down << "," << Fedx2y2[3][i].up << "," << -Fedx2y2[3][i].down << "," << Fedz2[3][i].up << "," << -Fedz2[3][i].down << "," << Fet2g[3][i].up << "," << -Fet2g[3][i].down << "," << Feeg[3][i].up << "," << -Feeg[3][i].down <<  "\n";
        }
        ost3.close();
        
        
        string name4 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe4DOS.csv";
        ofstream ost4(name4.c_str());
        if (!ost4) cout << "cannot open file " << name4;
        
        ost4 << "Model: " << outputfolder << "\n";
        ost4 << "Fe index ," << 4 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
        ost4 << "E, E - Ef (UNITS?), Fe4 DOS s - up, Fe4 DOS s - down, Fe4 DOS p - up, Fe4 DOS p - down, Fe4 DOS d - up, Fe4 DOS d - down, (+)Fe4 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe4 dxy up, Fe4 dxy down, Fe4 dxz up, Fe4 dxz down, Fe4 dyz up, Fe4 dyz down, Fe4 dx2-y2 up, Fe4 dx2-y2 down, Fe4 dz2 up, Fe4 dz2 down, Fe4 t2g up, Fe4 t2g down, Fe4 eg up, Fe4 eg down\n";
        
        for (int i = 0; i < datapoints; ++i){
            ost4 << energy[i] << "," << energy[i] - Ef << "," << Fes[4][i].up <<  "," << -Fes[4][i].down << "," << Fep[4][i].up <<  "," << -Fep[4][i].down << "," << Fed[4][i].up << "," << -Fed[4][i].down << "," << Fed[4][i].down <<  "," <<   densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down <<  "," << Fedxy[4][i].up << "," << -Fedxy[4][i].down << "," << Fedxz[4][i].up << "," << -Fedxz[4][i].down << "," << Fedyz[4][i].up << "," << -Fedyz[4][i].down << "," << Fedx2y2[4][i].up << "," << -Fedx2y2[4][i].down << "," << Fedz2[4][i].up << "," << -Fedz2[4][i].down << "," << Fet2g[4][i].up << "," << -Fet2g[4][i].down << "," << Feeg[4][i].up << "," << -Feeg[4][i].down <<  "\n";
            
        }
        ost4.close();
        
        
        string name5 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe5DOS.csv";
        ofstream ost5(name5.c_str());
        if (!ost5) cout << "cannot open file " << name5;
        
        ost5 << "Model: " << outputfolder << "\n";
        ost5 << "Fe index ," << 5 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
        ost5 << "E, E - Ef (UNITS?), Fe5 DOS s - up, Fe5 DOS s - down, Fe5 DOS p - up, Fe5 DOS p - down, Fe5 DOS d - up, Fe5 DOS d - down, (+)Fe5 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe5 dxy up, Fe5 dxy down, Fe5 dxz up, Fe5 dxz down, Fe5 dyz up, Fe5 dyz down, Fe5 dx2-y2 up, Fe5 dx2-y2 down, Fe5 dz2 up, Fe5 dz2 down,Fe5 t2g up, Fe5 t2g down, Fe5 eg up, Fe5 eg down\n";
        
        
        for (int i = 0; i < datapoints; ++i){
            ost5 << energy[i] << "," << energy[i] - Ef << "," << Fes[5][i].up <<  "," << -Fes[5][i].down << "," << Fep[5][i].up <<  "," << -Fep[5][i].down << "," << Fed[5][i].up << "," << -Fed[5][i].down << "," << Fed[5][i].down <<  "," <<   densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down <<  "," << Fedxy[5][i].up << "," << -Fedxy[5][i].down << "," << Fedxz[5][i].up << "," << -Fedxz[5][i].down << "," << Fedyz[5][i].up << "," << -Fedyz[5][i].down << "," << Fedx2y2[5][i].up << "," << -Fedx2y2[5][i].down << "," << Fedz2[5][i].up << "," << -Fedz2[5][i].down << "," << Fet2g[5][i].up << "," << -Fet2g[5][i].down << "," << Feeg[5][i].up << "," << -Feeg[5][i].down <<  "\n";
        }
        
        ost5.close();
        
        string name6 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe6DOS.csv";
        ofstream ost6(name6.c_str());
        if (!ost6) cout << "cannot open file " << name6;
        
        ost6 << "Model: " << outputfolder << "\n";
        ost6 << "Fe index ," << 6 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
        ost6 << "E, E - Ef (UNITS?), Fe6 DOS s - up, Fe6 DOS s - down, Fe6 DOS p - up, Fe6 DOS p - down, Fe6 DOS d - up, Fe6 DOS d - down, (+)Fe6 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe6 dxy up, Fe6 dxy down, Fe6 dxz up, Fe6 dxz down, Fe6 dyz up, Fe6 dyz down, Fe6 dx2-y2 up, Fe6 dx2-y2 down, Fe6 dz2 up, Fe6 dz2 down,Fe6 t2g up, Fe6 t2g down, Fe6 eg up, Fe6 eg down\n";
        
        for (int i = 0; i < datapoints; ++i){
            ost6 << energy[i] << "," << energy[i] - Ef << "," << Fes[6][i].up <<  "," << -Fes[6][i].down << "," << Fep[6][i].up <<  "," << -Fep[6][i].down << "," << Fed[6][i].up << "," << -Fed[6][i].down << "," << Fed[6][i].down << "," <<   densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down <<  "," << Fedxy[6][i].up << "," << -Fedxy[6][i].down << "," << Fedxz[6][i].up << "," << -Fedxz[6][i].down << "," << Fedyz[6][i].up << "," << -Fedyz[6][i].down << "," << Fedx2y2[6][i].up << "," << -Fedx2y2[6][i].down << "," << Fedz2[6][i].up << "," << -Fedz2[6][i].down << "," << Fet2g[6][i].up << "," << -Fet2g[6][i].down << "," << Feeg[6][i].up << "," << -Feeg[6][i].down <<  "\n";
        }
        
        ost6.close();
        string name7 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe7DOS.csv";
        ofstream ost7(name7.c_str());
        if (!ost7) cout << "cannot open file " << name7;
        
        ost7 << "Model: " << outputfolder << "\n";
        ost7 << "Fe index ," << 7 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
        ost7 << "E, E - Ef (UNITS?), Fe7 DOS s - up, Fe7 DOS s - down, Fe7 DOS p - up, Fe7 DOS p - down, Fe7 DOS d - up, Fe7 DOS d - down, (+)Fe7 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe7 dxy up, Fe7 dxy down, Fe7 dxz up, Fe7 dxz down, Fe7 dyz up, Fe7 dyz down, Fe7 dx2-y2 up, Fe7 dx2-y2 down, Fe7 dz2 up, Fe7 dz2 down, Fe7 t2g up, Fe7 t2g down, Fe7 eg up, Fe7 eg down\n";
        
        for (int i = 0; i < datapoints; ++i){
            ost7 << energy[i] << "," << energy[i] - Ef << "," << Fes[7][i].up <<  "," << -Fes[7][i].down << "," << Fep[7][i].up <<  "," << -Fep[7][i].down << "," << Fed[7][i].up << "," << -Fed[7][i].down << "," << Fed[7][i].down << ","  <<  densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[7][i].up << "," << -Fedxy[7][i].down << "," << Fedxz[7][i].up << "," << -Fedxz[7][i].down << "," << Fedyz[7][i].up << "," << -Fedyz[7][i].down << "," << Fedx2y2[7][i].up << "," << -Fedx2y2[7][i].down << "," << Fedz2[7][i].up << "," << -Fedz2[7][i].down << "," << Fet2g[7][i].up << "," << -Fet2g[7][i].down << "," << Feeg[7][i].up << "," << -Feeg[7][i].down <<  "\n";
        }
        
        ost7.close();
            
            string name8 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe8DOS.csv";
            ofstream ost8(name8.c_str());
            if (!ost8) cout << "cannot open file " << name8;
            
            ost8 << "Model: " << outputfolder << "\n";
            ost8 << "Fe index ," << 8 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
            ost8 << "E, E - Ef (UNITS?), Fe8 DOS s - up, Fe8 DOS s - down, Fe8 DOS p - up, Fe8 DOS p - down, Fe8 DOS d - up, Fe8 DOS d - down, (+)Fe8 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe8 dxy up, Fe8 dxy down, Fe8 dxz up, Fe8 dxz down, Fe8 dyz up, Fe8 dyz down, Fe8 dx2-y2 up, Fe8 dx2-y2 down, Fe8 dz2 up, Fe8 dz2 down, Fe8 t2g up, Fe8 t2g down, Fe8 eg up, Fe8 eg down\n";
            
            for (int i = 0; i < datapoints; ++i){
                ost8 << energy[i] << "," << energy[i] - Ef << "," << Fes[8][i].up <<  "," << -Fes[8][i].down << "," << Fep[8][i].up <<  "," << -Fep[8][i].down << "," << Fed[8][i].up << "," << -Fed[8][i].down <<  "," <<  Fed[8][i].down << "," << densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[8][i].up << "," << -Fedxy[8][i].down << "," << Fedxz[8][i].up << "," << -Fedxz[8][i].down << "," << Fedyz[8][i].up << "," << -Fedyz[8][i].down << "," << Fedx2y2[8][i].up << "," << -Fedx2y2[8][i].down << "," << Fedz2[8][i].up << "," << -Fedz2[8][i].down << "," << Fet2g[8][i].up << "," << -Fet2g[8][i].down << "," << Feeg[8][i].up << "," << -Feeg[8][i].down <<  "\n";

            }
        
         ost8.close();
            
            
            string name9 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe9DOS.csv";
            ofstream ost9(name9.c_str());
            if (!ost9) cout << "cannot open file " << name9;
            
            ost9 << "Model: " << outputfolder << "\n";
            ost9 << "Fe index ," << 9 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
            ost9 << "E, E - Ef (UNITS?), Fe9 DOS s - up, Fe9 DOS s - down, Fe9 DOS p - up, Fe9 DOS p - down, Fe9 DOS d - up, Fe9 DOS d - down, (+)Fe9 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe9 dxy up, Fe9 dxy down, Fe9 dxz up, Fe9 dxz down, Fe9 dyz up, Fe9 dyz down, Fe9 dx2-y2 up, Fe9 dx2-y2 down, Fe9 dz2 up, Fe9 dz2 down,Fe9 t2g up, Fe9 t2g down, Fe9 eg up, Fe9 eg down\n";
            
            for (int i = 0; i < datapoints; ++i){
                ost9 << energy[i] << "," << energy[i] - Ef << "," << Fes[9][i].up <<  "," << -Fes[9][i].down << "," << Fep[9][i].up <<  "," << -Fep[9][i].down << "," << Fed[9][i].up << "," << -Fed[9][i].down <<  "," << Fed[9][i].down << "," << densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[9][i].up << "," << -Fedxy[9][i].down << "," << Fedxz[9][i].up << "," << -Fedxz[9][i].down << "," << Fedyz[9][i].up << "," << -Fedyz[9][i].down << "," << Fedx2y2[9][i].up << "," << -Fedx2y2[9][i].down << "," << Fedz2[9][i].up << "," << -Fedz2[9][i].down << "," << Fet2g[9][i].up << "," << -Fet2g[9][i].down << "," << Feeg[9][i].up << "," << -Feeg[9][i].down <<  "\n";
            }
            
         ost9.close();
        
            string name10 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe10DOS.csv";
            ofstream ost10(name10.c_str());
            if (!ost10) cout << "cannot open file " << name10;
            
            ost10 << "Model: " << outputfolder << "\n";
            ost10 << "Fe index ," << 10 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
            ost10 << "E, E - Ef (UNITS?), Fe10 DOS s - up, Fe10 DOS s - down, Fe10 DOS p - up, Fe10 DOS p - down, Fe10 DOS d - up, Fe10 DOS d - down, (+)Fe10 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe10 dxy up, Fe10 dxy down, Fe10 dxz up, Fe10 dxz down, Fe10 dyz up, Fe10 dyz down, Fe10 dx2-y2 up, Fe10 dx2-y2 down, Fe10 dz2 up, Fe10 dz2 down, Fe10 t2g up, Fe10 t2g down, Fe10 eg up, Fe10 eg down\n";
            
            for (int i = 0; i < datapoints; ++i){
                ost10 << energy[i] << "," << energy[i] - Ef << "," << Fes[10][i].up <<  "," << -Fes[10][i].down << "," << Fep[10][i].up <<  "," << -Fep[10][i].down << "," << Fed[10][i].up << "," << -Fed[10][i].down <<  "," << Fed[10][i].down << "," << densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[10][i].up << "," << -Fedxy[10][i].down << "," << Fedxz[10][i].up << "," << -Fedxz[10][i].down << "," << Fedyz[10][i].up << "," << -Fedyz[10][i].down << "," << Fedx2y2[10][i].up << "," << -Fedx2y2[10][i].down << "," << Fedz2[10][i].up << "," << -Fedz2[10][i].down << "," << Fet2g[10][i].up << "," << -Fet2g[10][i].down << "," << Feeg[10][i].up << "," << -Feeg[10][i].down <<  "\n";
            }
             ost10.close();
            
            string name11 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe11DOS.csv";
            ofstream ost11(name11.c_str());
            if (!ost11) cout << "cannot open file " << name11;
            
            ost11 << "Model: " << outputfolder << "\n";
            ost11 << "Fe index ," << 11 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
            ost11 << "E, E - Ef (UNITS?), Fe11 DOS s - up, Fe11 DOS s - down, Fe11 DOS p - up, Fe11 DOS p - down, Fe11 DOS d - up, Fe11 DOS d - down, (+)Fe11 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe11 dxy up, Fe11 dxy down, Fe11 dxz up, Fe11 dxz down, Fe11 dyz up, Fe11 dyz down, Fe11 dx2-y2 up, Fe11 dx2-y2 down,Fe11 dz2 up, Fe11 dz2 down, Fe11 t2g up, Fe11 t2g down, Fe11 eg up, Fe11 eg down\n";
            
            for (int i = 0; i < datapoints; ++i){
                ost11 << energy[i] << "," << energy[i] - Ef << "," << Fes[11][i].up <<  "," << -Fes[11][i].down << "," << Fep[11][i].up <<  "," << -Fep[11][i].down << "," << Fed[11][i].up << "," << -Fed[11][i].down <<  "," << Fed[11][i].down << "," << densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[11][i].up << "," << -Fedxy[11][i].down << "," << Fedxz[11][i].up << "," << -Fedxz[11][i].down << "," << Fedyz[11][i].up << "," << -Fedyz[11][i].down << "," << Fedx2y2[11][i].up << "," << -Fedx2y2[11][i].down << "," << Fedz2[11][i].up << "," << -Fedz2[11][i].down << "," << Fet2g[11][i].up << "," << -Fet2g[11][i].down << "," << Feeg[11][i].up << "," << -Feeg[11][i].down <<  "\n";
            }
         ost11.close();
            
            
            
            string name12 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe12DOS.csv";
            ofstream ost12(name12.c_str());
            if (!ost12) cout << "cannot open file " << name12;
            
            ost12 << "Model: " << outputfolder << "\n";
            ost12 << "Fe index ," << 12 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
        ost12 << "E, E - Ef (UNITS?), Fe12 DOS s - up, Fe12 DOS s - down, Fe12 DOS p - up, Fe12 DOS p - down, Fe12 DOS d - up, Fe12 DOS d - down, (+)Fe12 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe12 dxy up, Fe12 dxy down, Fe12 dxz up, Fe12 dxz down, Fe12 dyz up, Fe12 dyz down, Fe12 dx2-y2 up, Fe12 dx2-y2 down, Fe12 dz2 up, Fe12 dz2 down,Fe12 t2g up, Fe12 t2g down, Fe12 eg up, Fe12 eg down\n";

            
            for (int i = 0; i < datapoints; ++i){
                ost12 << energy[i] << "," << energy[i] - Ef << "," << Fes[12][i].up <<  "," << -Fes[12][i].down << "," << Fep[12][i].up <<  "," << -Fep[12][i].down << "," << Fed[12][i].up << "," << -Fed[12][i].down <<  "," << Fed[12][i].down << "," << densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[12][i].up << "," << -Fedxy[12][i].down << "," << Fedxz[12][i].up << "," << -Fedxz[12][i].down << "," << Fedyz[12][i].up << "," << -Fedyz[12][i].down << "," << Fedx2y2[12][i].up << "," << -Fedx2y2[12][i].down << "," << Fedz2[12][i].up << "," << -Fedz2[12][i].down << "," << Fet2g[12][i].up << "," << -Fet2g[12][i].down << "," << Feeg[12][i].up << "," << -Feeg[12][i].down <<  "\n";
            }
            
         ost12.close();
        
            string name13 = "/Users/joekirk/Desktop/DOSCAR/" + outputfolder +"/" + outputfolder+"_Fe13DOS.csv";
            ofstream ost13(name13.c_str());
            if (!ost13)  cout << "cannot open file " << name13;
            
            ost13 << "Model: " << outputfolder << "\n";
            ost13 << "Fe index ," << 13 << ",Highest E: ," << highestE << ",Lowest E: ," << lowestE << ",#datapoints: ," << datapoints << ",Ef: ," << Ef << "\n";
        ost13 << "E, E - Ef (UNITS?), Fe13 DOS s - up, Fe13 DOS s - down, Fe13 DOS p - up, Fe13 DOS p - down, Fe13 DOS d - up, Fe13 DOS d - down, (+)Fe13 DOS d - down, TotalDOS up, TotalDOS down, Fe DOS s - up, Fe DOS s - down, Fe DOS p - up, Fe DOS p - down, Fe DOS d - up, Fe DOS d - down, Fe13 dxy up, Fe13 dxy down, Fe13 dxz up, Fe13 dxz down, Fe13 dyz up, Fe13 dyz down, Fe13 dx2-y2 up, Fe13 dx2-y2 down, Fe13 dz2 up, Fe13 dz2 down, Fe13 t2g up, Fe13 t2g down, Fe13 eg up, Fe13 eg down\n";

            
            for (int i = 0; i < datapoints; ++i){
                ost13 << energy[i] << "," << energy[i] - Ef << "," << Fes[13][i].up <<  "," << -Fes[13][i].down << "," << Fep[13][i].up <<  "," << -Fep[13][i].down << "," << Fed[13][i].up << "," << -Fed[13][i].down <<  "," <<  Fed[13][i].down << "," << densityofstates[i].up << "," << -densityofstates[i].down << "," << totalFes[i].up <<  "," << -totalFes[i].down << "," << totalFep[i].up <<  "," << -totalFep[i].down << "," << totalFed[i].up << "," << -totalFed[i].down << "," << Fedxy[13][i].up << "," << -Fedxy[13][i].down << "," << Fedxz[13][i].up << "," << -Fedxz[13][i].down << "," << Fedyz[13][i].up << "," << -Fedyz[13][i].down << "," << Fedx2y2[13][i].up << "," << -Fedx2y2[13][i].down << "," << Fedz2[13][i].up << "," << -Fedz2[13][i].down << "," << Fet2g[13][i].up << "," << -Fet2g[13][i].down << "," << Feeg[13][i].up << "," << -Feeg[13][i].down <<  "\n";
            }
            
         ost13.close();
    }
}

void Fileread::readcoords(istream& ist)
{
    for (int i = 0; i < modelsize ; ++i) {
        ist >> type >> x >> y >> z;
        data.push_back(Atomcoord(x, y, z));
        atomtypes.push_back(type);                           }
    frames.push_back(data);
    data.clear();
}

// -- Reads in the cell vectors into a vector of Atomcoords --

void Fileread::readprimvec (istream& ist)
{
    
    for (int i = 0; i < 3; ++i) {
        ist >> x >> y >> z;
        cellvectors.push_back(Atomcoord(x, y, z));
    }
}

// -- Reads all the frames within the file with the layout --

/*   ANIMSTEPS 4000
 *
 *   CRYSTAL
 *
 *   PRIMVEC
 *   14.67000   00.00000   00.00000
 *   00.00000   14.67000   00.00000
 *   00.00000   00.00000   14.67000
 *
 *   PRIMCOORD 1
 *   100 1
 *   26   11.00250   01.22250   03.66750
 *   32   01.22250   11.00250   06.11250
 *   32   11.00250   01.22250   06.11250
 */

void Fileread::readframes(istream& ist)
{
    string PRIMCOORD = "PRIMCOORD";
    string PRIMVEC = "PRIMVEC";
    string AMNISTEPS = "AMNISTEPS";
    
    int number1;
    int framenumber;
    
    
    while (!ist.eof()) {
        
        string line;
        getline(ist, line);
       
        if (line[0] == AMNISTEPS[0] && line[8] == AMNISTEPS[8]){
            istringstream ( line.substr(10) ) >> framesinfile;
        }
        if (line[0] == PRIMVEC[0] && line[6] == PRIMVEC[6])
            readprimvec(ist);
        if (line[0] == PRIMCOORD[0] && line[8]==PRIMCOORD[8]) {
            istringstream ( line.substr(9) ) >> framenumber;
            ist >> modelsize >> number1;
            readcoords(ist);
        }
            
    }

    
}

void Fileread::readframes(istream& ist, int datapoints)
{
    string PRIMCOORD = "PRIMCOORD";
    string PRIMVEC = "PRIMVEC";
    string AMNISTEPS = "AMNISTEPS";
    
    int number1;
    int framenumber;
    
    int increment= 0;
    int i = 1;
    
    while (!ist.eof()) {
        
        string line;
        getline(ist, line);
        
        if (line[0] == AMNISTEPS[0] && line[8] == AMNISTEPS[8]){
            istringstream ( line.substr(10) ) >> framesinfile;
            increment = floor(framesinfile/datapoints);
        }
        if (line[0] == PRIMVEC[0] && line[6] == PRIMVEC[6])
            readprimvec(ist);
        if (line[0] == PRIMCOORD[0] && line[8]==PRIMCOORD[8]) {
            istringstream ( line.substr(9) ) >> framenumber;
            if (framenumber != (increment * i)){}
            else {
            ist >> modelsize >> number1;
            readcoords(ist);
                ++i;}
        }
        
        
    }
    
}
void Fileread::boudaryconditions(int firstframe, int lastframe) {
    
    //fuction works out the fractional component of the interatomic vector along the x, y and z axis
    //respectively, if greater than 0.5 then we must subract the x,y or z cellvector
    
    
    for (int k = firstframe; k < lastframe; ++k) {
        for(int j = 0; j<modelsize; ++j) {
            for (int i=0; i<modelsize; ++i)
            {
                Atomcoord interatomicvector = frames[k][j] - frames[k][i];
                double p = dotproduct(interatomicvector, cellvectors[0]);
                
                double q = dotproduct(interatomicvector, cellvectors[1]);
                
                double r = dotproduct(interatomicvector, cellvectors[2]);
                
                double s, t, u;
                
                if (((p/mod(cellvectors[0]))/mod(cellvectors[0])) > 0.5)
                    s = interatomicvector.x - cellvectors[0].x;
                
                else if (((p/mod(cellvectors[0]))/mod(cellvectors[0])) < -0.5)
                    s = interatomicvector.x + cellvectors[0].x;
                
                else s = interatomicvector.x;
                
                if ((q/mod(cellvectors[1]))/mod(cellvectors[1]) > 0.5)
                    t = interatomicvector.y - cellvectors[1].y;
                
                else if ((q/mod(cellvectors[1]))/mod(cellvectors[1]) < -0.5)
                    t = interatomicvector.y + cellvectors[1].y;
                
                else t = interatomicvector.y;
                
                
                if ((r/mod(cellvectors[2]))/mod(cellvectors[2]) > 0.5)
                    u = interatomicvector.z - cellvectors[2].z;
                
                else if ((r/mod(cellvectors[2]))/mod(cellvectors[2]) < -0.5)
                    u = interatomicvector.z + cellvectors[2].z;
                
                else u = interatomicvector.z;
                
                rijoneatom.push_back(Atomcoord(s, t, u));
            }
            
            rijallatoms.push_back(rijoneatom);
            rijoneatom.clear();
            
            
        }
        rijallframes.push_back(rijallatoms);
        rijallatoms.clear();
        
    }
    
    //rijallframes[k][j][i]
    // k is the frame
    // j is the atom we are looking at the distance between
    // i is the origin atom
}


void Fileread::boundaryconditions() {
    
    //fuction works out the fractional component of the interatomic vector along the x, y and z axis
    //respectively, if greater than 0.5 then we must subract the x,y or z cellvector
    
    
    for (int k = 0; k < frames.size(); ++k) {
        for(int j = 0; j<modelsize; ++j) {
            for (int i=0; i<modelsize; ++i)
            {
                Atomcoord interatomicvector = frames[k][j] - frames[k][i];
                
                double p = dotproduct(interatomicvector, cellvectors[0]);
                
                double q = dotproduct(interatomicvector, cellvectors[1]);
                
                double r = dotproduct(interatomicvector, cellvectors[2]);
                
                double s, t, u;
                
                if (((p/mod(cellvectors[0]))/mod(cellvectors[0])) > 0.5)
                    s = interatomicvector.x - cellvectors[0].x;
                
                else if (((p/mod(cellvectors[0]))/mod(cellvectors[0])) < -0.5)
                    s = interatomicvector.x + cellvectors[0].x;
                
                else s = interatomicvector.x;
                
                if ((q/mod(cellvectors[1]))/mod(cellvectors[1]) > 0.5)
                    t = interatomicvector.y - cellvectors[1].y;
                
                else if ((q/mod(cellvectors[1]))/mod(cellvectors[1]) < -0.5)
                    t = interatomicvector.y + cellvectors[1].y;
                
                else t = interatomicvector.y;
                
                
                if ((r/mod(cellvectors[2]))/mod(cellvectors[2]) > 0.5)
                    u = interatomicvector.z - cellvectors[2].z;
                
                else if ((r/mod(cellvectors[2]))/mod(cellvectors[2]) < -0.5)
                    u = interatomicvector.z + cellvectors[2].z;
                
                else u = interatomicvector.z;
                
                rijoneatom.push_back(Atomcoord(s, t, u));
            }
            
            rijallatoms.push_back(rijoneatom);
            rijoneatom.clear();
            
            
        }
        rijallframes.push_back(rijallatoms);
        rijallatoms.clear();
        
    }
    
    //rijallframes[k][j][i]
    // k is the frame
    // j is the atom we are looking at the distance between
    // i is the origin atom
}


void Fileread::generateneighbourtable()
{
    
    
    Numeric_lib::Matrix <double,2>  neighbourtable (modelsize, modelsize);
    
    for (int k=0; k < frames.size(); ++k) {
        
        for (int j=0; j< modelsize; ++j) {
            
            for (int i=0; i< modelsize; ++i) {
                
                //neighbourtable[j][i] => j is the row i is the column
                neighbourtable[j][i] = mod(rijallframes[k][i][j]);
            }
        }
        neighbourtables.push_back(neighbourtable);
    }
}

void Fileread::generateneighbourtable(int firstframe, int lastframe)
{
    
    Numeric_lib::Matrix <double,2>  neighbourtable (modelsize, modelsize);
    
    for (int k=0; k < lastframe - firstframe; ++k) {
        
        for (int j=0; j< modelsize; ++j) {
            
            for (int i=0; i< modelsize; ++i) {
                
                //neighbourtable[j][i] => j is the row i is the column
                neighbourtable[j][i] = mod(rijallframes[k][i][j]);
            }
        }
        neighbourtables.push_back(neighbourtable);
    }
}

//reads the entire file but only imposes boundary conditions generating a neighbour table bewteen specified frames
//the name entered should be the file name

void Fileread::read(string name, int firstframe, int lastframe)
{
    string nameaccess= "/Users/joekirk/Desktop/Models/" + name;
    ifstream ist (nameaccess.c_str());
    if (!ist) cout << "cannot open file " << name;
    
    readframes(ist);
    boudaryconditions(firstframe, lastframe);
    generateneighbourtable(firstframe, lastframe);
    
}

// only reads in the necessary, it will parse the whole file but only read in the frames that we need **--DO NOT USE FOR MSD--**
/*void Fileread::read(string name, int firstframe, int lastframe, int datapoints)
{
    string nameaccess= "/Users/joekirk/Desktop/Models/" + name;
    ifstream ist (nameaccess.c_str());
    if (!ist) cout << "cannot open file " << name;
    
    readframes(ist);
    boudaryconditions(firstframe, lastframe);
    generateneighbourtable(firstframe, lastframe);
    
}*/

void Fileread::read(string name)
{

    string nameaccess= "/Users/joekirk/Desktop/Models/" + name;
    ifstream ist (nameaccess.c_str());
    if (!ist) cout << "cannot open file " << name;
    
    readframes(ist);
    boundaryconditions();
    generateneighbourtable();
    
}

void Fileread::readformsd(string name)
{
    string nameaccess= "/Users/joekirk/Desktop/Models/" + name;
    ifstream ist (nameaccess.c_str());
    if (!ist) cout << "cannot open file " << name;
    
    readframes(ist);
}

void Fileread::readdynamicanalysis(string name, int datapoints)
{
    
    
    string nameaccess= "/Users/joekirk/Desktop/Models/" + name;
    ifstream ist (nameaccess.c_str());
    if (!ist) cout << "cannot open file " << name;
    
    readframes(ist, datapoints);
    boundaryconditions();
    generateneighbourtable();
    

}


