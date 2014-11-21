//
//  structural analysis.cpp
//  Part III Project - version1
//
//  Created by Joe Kirk on 17/01/2013.
//  Copyright (c) 2013 Joe Kirk. All rights reserved.
//

#include "structuralanalysis.h"


#define REAL 0
#define IMAG 1

//generates radial distribution functions for atom A to all other atoms; Fe, Ge, Sb and Te and X
void staticanalysis::g(int firstframe, int lastframe, int datapoints,  atmtype A,  Fileread& file)
{
    
    int numberofframes = lastframe - firstframe;
    
    double rmin = 0.0;
    double rmax = 7.0;
    double deltaR = (rmax - rmin) / datapoints;
    vector <double> Rvalues;
    
    for (double R = rmin; R < rmax; R += deltaR)
    {
        Rvalues.push_back(R);
    }
    
    vector <double> numberinbinAX  (datapoints);
    vector <double> numberinbinAFe (datapoints);
    vector <double> numberinbinAGe (datapoints);
    vector <double> numberinbinASb (datapoints);
    vector <double> numberinbinATe (datapoints);
    
    for (int k = 0; k < numberofframes; ++k) { // k=framenumber => cycle through all relevant frames
        
        for (int j = 0; j <file.modelsize; ++j){  // j=index for atom B
            
            for (int i = 0; i <file.modelsize; ++i) { // i=index for atom A
                
                if (file.atomtypes[i] == A && file.neighbourtables[k][i][j] < rmax)  {
                    int h = floor((file.neighbourtables[k][i][j]-rmin)/deltaR);
                    numberinbinAX[h] += 1.0;
                }
                
                if (file.atomtypes[i] == A && file.atomtypes[j] == Fe && file.neighbourtables[k][i][j] <= rmax)
                {
                    int a =  floor((file.neighbourtables[k][i][j] - rmin) / deltaR);
                    numberinbinAFe[a]+=1.0;
                }
                
                if (file.atomtypes[i] == A && file.atomtypes[j] == Ge && file.neighbourtables[k][i][j] <= rmax)
                {
                    int b =  floor((file.neighbourtables[k][i][j] - rmin) / deltaR);
                    numberinbinAGe[b]+=1.0;
                }
                
                if (file.atomtypes[i] == A && file.atomtypes[j] == Sb && file.neighbourtables[k][i][j] <= rmax)
                {
                    int c =  floor((file.neighbourtables[k][i][j] - rmin) / deltaR);
                    numberinbinASb[c]+=1.0;
                }
                if (file.atomtypes[i] == A && file.atomtypes[j] == Te && file.neighbourtables[k][i][j] <= rmax)
                {
                    int d =  floor((file.neighbourtables[k][i][j] - rmin) / deltaR);
                    numberinbinATe[d]+=1.0;
                }
            } 
        }
    }
    
    
    int NFe =0, NGe = 0, NSb = 0, NTe = 0, Na = 0;
    
    for (int i = 0; i < file.modelsize; ++i)
    {
        if (file.atomtypes[i] == Fe)
            ++NFe;
        if (file.atomtypes[i] == Ge)
            ++NGe;
        if (file.atomtypes[i] == Sb)
            ++NSb;
        if (file.atomtypes[i] == Te)
            ++NTe;
        if (file.atomtypes[i] == A)
            ++Na;
    }
    
    vector<double> binAX;
    vector<double> binAFe;
    vector<double> binAGe;
    vector<double> binASb;
    vector<double> binATe;
    
    double coordinationnox = 0.0;
    double coordinationnoFe = 0.0;
    double coordinationnoGe = 0.0;
    double coordinationnoSb = 0.0;
    double coordinationnoTe = 0.0;
    
    vector <double> totaln; //vector containing running totals of the co-ordination number (i.e intg(r)) to all atoms
    vector <double> totalnFe; //vector containing running totals of the co-ordination number (i.e intg(r)) to  Fe
    vector <double> totalnGe;
    vector <double> totalnSb;
    vector <double> totalnTe;
    
    
    for (int i = 1; i < datapoints; ++i) {
        
       double gx = numberinbinAX[i]/(4*M_PI*Na*(file.modelsize - 1)*Rvalues[i]*Rvalues[i]*deltaR*numberofframes/(dotproduct(file.cellvectors[0], crossproduct(file.cellvectors[1], file.cellvectors[2]))));
        
        binAX.push_back(gx);
        
        double nx = gx*Rvalues[i]*Rvalues[i] * deltaR * 4 * M_PI * (file.modelsize - 1) /(dotproduct(file.cellvectors[0], crossproduct(file.cellvectors[1], file.cellvectors[2])));
        
        coordinationnox += nx;
        totaln.push_back(coordinationnox);
        
        
        double gFe = numberinbinAFe[i]/(4*M_PI*Na*NFe*Rvalues[i]*Rvalues[i] * deltaR * numberofframes/(dotproduct(file.cellvectors[0], crossproduct(file.cellvectors[1], file.cellvectors[2]))));
        
        binAFe.push_back(gFe);
        
        double nFe = gFe*Rvalues[i]*Rvalues[i] * deltaR * 4 * M_PI * NFe /(dotproduct(file.cellvectors[0], crossproduct(file.cellvectors[1], file.cellvectors[2])));
        
        coordinationnoFe += nFe;
        totalnFe.push_back(coordinationnoFe);
        
        
        double gGe = numberinbinAGe[i]/(4*M_PI*Na*NGe*Rvalues[i]*Rvalues[i] * deltaR * numberofframes/(dotproduct(file.cellvectors[0], crossproduct(file.cellvectors[1], file.cellvectors[2]))));
        
        binAGe.push_back(gGe);
        
        double nGe = gGe*Rvalues[i]*Rvalues[i] * deltaR * 4 * M_PI * NGe /(dotproduct(file.cellvectors[0], crossproduct(file.cellvectors[1], file.cellvectors[2])));
        
        coordinationnoGe += nGe;
        totalnGe.push_back(coordinationnoGe);
        
        
        double gSb = numberinbinASb[i]/(4*M_PI*Na*NSb*Rvalues[i]*Rvalues[i] * deltaR * numberofframes/(dotproduct(file.cellvectors[0], crossproduct(file.cellvectors[1], file.cellvectors[2]))));
        
        binASb.push_back(gSb);
        
        double nSb = gSb*Rvalues[i]*Rvalues[i] * deltaR * 4 * M_PI * NSb /(dotproduct(file.cellvectors[0], crossproduct(file.cellvectors[1], file.cellvectors[2])));
        
        coordinationnoSb += nSb;
        totalnSb.push_back(coordinationnoSb);
        
        double gTe = numberinbinATe[i]/(4*M_PI*Na*NTe*Rvalues[i]*Rvalues[i] * deltaR * numberofframes/(dotproduct(file.cellvectors[0], crossproduct(file.cellvectors[1], file.cellvectors[2]))));
        
        binATe.push_back(gTe);
        
        double nTe = gTe*Rvalues[i]*Rvalues[i] * deltaR * 4 * M_PI * NTe /(dotproduct(file.cellvectors[0], crossproduct(file.cellvectors[1], file.cellvectors[2])));
        
        coordinationnoTe += nTe;
        totalnTe.push_back(coordinationnoTe);
        
    }
    
    int window = 10;
    int completewindow = datapoints - (window-1);
    
    vector<double> extendeddataAFe(datapoints + (datapoints-completewindow));
    vector<double> extendeddataAGe(datapoints + (datapoints-completewindow));
    vector<double> extendeddataASb(datapoints + (datapoints-completewindow));
    vector<double> extendeddataATe(datapoints + (datapoints-completewindow));
    vector<double> extendeddataAX(datapoints + (datapoints-completewindow));
    
    int edges = extendeddataAFe.size() - datapoints;
    
    if (edges%2 == 0) {
        
        if (edges == 0) {}
        
        else {
            for (int i=0; i < edges/2; ++i)
            {
                extendeddataAFe[i] = binAFe[0];
                extendeddataAGe[i] = binAGe[0];
                extendeddataASb[i] = binASb[0];
                extendeddataATe[i] = binATe[0];
                extendeddataAX[i] = binAX[0];
            }}
        
        if (edges/2+datapoints == extendeddataAFe.size()) {}
        
        else {
            for (int i=edges/2+datapoints; i < extendeddataAFe.size(); ++i)
            {
                extendeddataAFe[i] = binAFe[datapoints-1];
                extendeddataAGe[i] = binAGe[datapoints-1];
                extendeddataASb[i] = binASb[datapoints-1];
                extendeddataATe[i] = binATe[datapoints-1];
                extendeddataAX[i] = binAX[datapoints-1];
            }}
        
        for (int i = edges/2; i < datapoints+edges/2; ++i)
        {
            extendeddataAFe[i] = binAFe[i-edges/2];
            extendeddataAGe[i] = binAGe[i-edges/2];
            extendeddataASb[i] = binASb[i-edges/2];
            extendeddataATe[i] = binATe[i-edges/2];
            extendeddataAX[i] = binAX[i-edges/2];

        }
    }
    
    else if(edges%2 ==1) {
        
        for (int i=0; i < floor(edges/2)+1; ++i)
        {
            extendeddataAFe[i] = binAFe[0];
            extendeddataAGe[i] = binAGe[0];
            extendeddataASb[i] = binASb[0];
            extendeddataATe[i] = binATe[0];
            extendeddataAX[i] = binAX[0];
        }
        
        if (floor(edges/2) + 1 + datapoints ==extendeddataAFe.size())
        {}
        
        else {
            for (int i=floor(edges/2) + 1 + datapoints; i < extendeddataAFe.size(); ++i)
            {
                extendeddataAFe[i] = binAFe[datapoints-1];
                extendeddataAGe[i] = binAGe[datapoints-1];
                extendeddataASb[i] = binASb[datapoints-1];
                extendeddataATe[i] = binATe[datapoints-1];
                extendeddataAX[i] = binAX[datapoints-1];
                
            }
        }
        
        for (int i = floor(edges/2)+1; i < datapoints+floor(edges/2)+1; ++i)
        {
            extendeddataAFe[i] = binAFe[i-(floor(edges/2)+1)];
            extendeddataAGe[i] = binAGe[i-(floor(edges/2)+1)];
            extendeddataASb[i] = binASb[i-(floor(edges/2)+1)];
            extendeddataATe[i] = binATe[i-(floor(edges/2)+1)];
            extendeddataAX[i] = binAX[i-(floor(edges/2)+1)];
        }
    }
    
    vector<double> tobesortedAFe;
    vector<double> tobesortedAGe;
    vector<double> tobesortedASb;
    vector<double> tobesortedATe;
    vector<double> tobesortedAX;
    
    vector<double> filtervalueAFe (datapoints);
    vector<double> filtervalueAGe (datapoints);
    vector<double> filtervalueASb (datapoints);
    vector<double> filtervalueATe (datapoints);
    vector<double> filtervalueAX (datapoints);
    
    for(int i =0; i < extendeddataAFe.size() - (window-1); ++i)
    {
        
        tobesortedAFe = extendeddataAFe;
        tobesortedAGe = extendeddataAGe;
        tobesortedASb = extendeddataASb;
        tobesortedATe = extendeddataATe;
        tobesortedAX = extendeddataAX;
    
        sort (tobesortedAFe.begin()+i, tobesortedAFe.begin()+(i+(window)));
        
        if (window%2 == 1)
            filtervalueAFe[i] = tobesortedAFe[(i + floor((window)/2))];
        
        else if (window%2 == 0)
            filtervalueAFe[i] = (tobesortedAFe[i+((window/2)-1)] + tobesortedAFe[i+(window/2)])/2;
        
        tobesortedAFe.clear();
        
        sort (tobesortedAGe.begin()+i, tobesortedAGe.begin()+(i+(window)));
        
        if (window%2 == 1)
            filtervalueAGe[i] = tobesortedAGe[(i + floor((window)/2))];
        
        else if (window%2 == 0)
            filtervalueAGe[i] = (tobesortedAGe[i+((window/2)-1)] + tobesortedAGe[i+(window/2)])/2;
        
        tobesortedAGe.clear();
        
        sort (tobesortedASb.begin()+i, tobesortedASb.begin()+(i+(window)));
        
        if (window%2 == 1)
            filtervalueASb[i] = tobesortedASb[(i + floor((window)/2))];
        
        else if (window%2 == 0)
            filtervalueASb[i] = (tobesortedASb[i+((window/2)-1)] + tobesortedASb[i+(window/2)])/2;
        
        tobesortedASb.clear();
        
        sort (tobesortedATe.begin()+i, tobesortedATe.begin()+(i+(window)));
        
        if (window%2 == 1)
            filtervalueATe[i] = tobesortedATe[(i + floor((window)/2))];
        
        else if (window%2 == 0)
            filtervalueATe[i] = (tobesortedATe[i+((window/2)-1)] + tobesortedATe[i+(window/2)])/2;
        
        tobesortedATe.clear();
        
        sort (tobesortedAX.begin()+i, tobesortedAX.begin()+(i+(window)));
        
        if (window%2 == 1)
            filtervalueAX[i] = tobesortedAX[(i + floor((window)/2))];
        
        else if (window%2 == 0)
            filtervalueAX[i] = (tobesortedAX[i+((window/2)-1)] + tobesortedAX[i+(window/2)])/2;
        
        tobesortedAX.clear();
    
    }

    
    
    double maximumAFe =0.0, maximumAGe =0.0, maximumASb =0.0, maximumATe =0.0, maximumAX = 0.0;
    double maxrAFe = 0.0, maxrAGe = 0.0, maxrASb = 0.0, maxrATe = 0.0, maxrAX = 0.0;
    double cutoff = 4.00;
    
    int iFemax =0, iGemax=0, iSbmax=0, iTemax=0, iXmax=0;
    
    for (int i = 1; i < datapoints-1; ++i)
    {
        if (filtervalueAFe[i] > maximumAFe && Rvalues[i] < cutoff){
            maximumAFe = filtervalueAFe[i];
            maxrAFe = Rvalues[i];
            iFemax = i;
            
        }
        if (filtervalueAGe[i] > maximumAGe && Rvalues[i] < cutoff){
            maximumAGe = filtervalueAGe[i];
            maxrAGe = Rvalues[i];
            iGemax = i;
        }
        if (filtervalueASb[i] > maximumASb && Rvalues[i] < cutoff){
            maximumASb = filtervalueASb[i];
            maxrASb = Rvalues[i];
            iSbmax=i;
        }
        if (filtervalueATe[i] > maximumATe && Rvalues[i] < cutoff){
            maximumATe = filtervalueATe[i];
            maxrATe = Rvalues[i];
            iTemax = i;
        }
        if (filtervalueAX[i] > maximumAX && Rvalues[i] < cutoff) {
            maximumAX = filtervalueAX[i];
            maxrAX = Rvalues[i];
            iXmax = i;
        }
    }
    
    double minimumAFe =maximumAFe, minimumAGe = maximumAGe, minimumASb =maximumASb, minimumATe =maximumATe, minimumAX = maximumAX;
    double minrAFe = 0.0, minrAGe = 0.0, minrASb = 0.0, minrATe = 0.0, minrAX = 0.0;
    
    int iFemin =0, iGemin=0, iSbmin=0, iTemin=0, iXmin=0;
  
    for (int i = 1; i < datapoints-1; ++i)
    {
        if (filtervalueAFe[i] < minimumAFe && Rvalues[i] < cutoff && Rvalues[i] > maxrAFe){
            minimumAFe = filtervalueAFe[i];
            minrAFe = Rvalues[i];
            iFemin = i;
        }
        if (filtervalueAGe[i] < minimumAGe && Rvalues[i] < cutoff && Rvalues[i] > maxrAGe){
            minimumAGe = filtervalueAGe[i];
            minrAGe = Rvalues[i];
            iGemin = i;
        }

        if (filtervalueASb[i] < minimumASb && Rvalues[i] < cutoff && Rvalues[i] > maxrASb){
            minimumASb = filtervalueASb[i];
            minrASb = Rvalues[i];
            iSbmin = i;
        }

        if (filtervalueATe[i] < minimumATe && Rvalues[i] < cutoff && Rvalues[i] > maxrATe){
            minimumATe = filtervalueATe[i];
            minrATe = Rvalues[i];
            iTemin = i;
        }

        if (filtervalueAX[i] < minimumAX && Rvalues[i] < cutoff && Rvalues[i] > maxrAX){
            minimumAX = filtervalueAX[i];
            minrAX = Rvalues[i];
            iXmin = i;
        }

        }
    
    
    
    //open the relevant csv file and output g(r), along with intg(r)
    string output = "/Users/joekirk/Desktop/Part III Project - Initial Results/" + file.name + "/RDF/RDF - " + enumtostring(A) + ".csv";
    ofstream ost(output.c_str());
    if (!ost) cout << "cannot open file " << output;
    
    ost << "Firstframe: " << firstframe << ", Lastframe: " << lastframe << ", Numberofframes: " << numberofframes << "\n";
    
    ost << "Data extracted from file: " << file.name << "\n";
   
    ost << enumtostring(A) << "-Fe maximum: ," << maximumAFe << " , rmax value: ," << maxrAFe << "," << "intg(rmax): ," << totalnFe[iFemax] << ","
        << enumtostring(A) << "-Fe minimum: ," << minimumAFe << " , rmin value: ," << minrAFe << "," << "intg(rmin): ," << totalnFe[iFemin] << "\n"
        << enumtostring(A) << "-Ge maximum: ," << maximumAGe << " , rmax value: ," << maxrAGe << "," << "intg(rmax): ," << totalnGe[iGemax] << ","
        << enumtostring(A) << "-Ge minimum: ," << minimumAGe << " , rmin value: ," << minrAGe << "," << "intg(rmin): ," << totalnGe[iGemin] << "\n"
        << enumtostring(A) << "-Sb maximum: ," << maximumASb << " , rmax value: ," << maxrASb << "," << "intg(rmax): ," << totalnSb[iSbmax] << ","
        << enumtostring(A) << "-Sb minimum: ," << minimumASb << " , rmin value: ," << minrASb << "," << "intg(rmin): ," << totalnSb[iSbmin] << "\n"
        << enumtostring(A) << "-Te maximum: ," << maximumATe << " , rmax value: ," << maxrATe << "," << "intg(rmax): ," << totalnTe[iTemax] << ","
        << enumtostring(A) << "-Te minimum: ," << minimumATe << " , rmin value: ," << minrATe << "," << "intg(rmin): ," << totalnTe[iTemin] << "\n"
        << enumtostring(A) << "-X maximum: ," << maximumAX << " , rmax value: ," << maxrAX << "," << "intg(rmax): ," << totaln[iXmax] << ","
        << enumtostring(A) << "-X minimum: ," << minimumAX << " , rmin value: ," << minrAX << "," << "intg(rmin): ," << totaln[iXmin] << "\n\n";
    
    ost <<  "r / A , g(r) " << enumtostring(A) << "- X / AU , Median Filtered " << enumtostring(A) << "- X / AU ,intg(r) " + enumtostring(A) << "- X/ AU , "
    <<  " g(r) " << enumtostring(A) << "- Fe / AU , Median Filtered " << enumtostring(A) << "- Fe / AU , intg(r) " + enumtostring(A) << "- Fe / AU , "
    <<  " g(r) " << enumtostring(A) << "- Ge / AU , Median Filtered " << enumtostring(A) << "- Ge / AU , intg(r) " + enumtostring(A) << "- Ge / AU , "
    <<  " g(r) " << enumtostring(A) << "- Sb / AU , Median Filtered " << enumtostring(A) << "- Sb / AU , intg(r) " + enumtostring(A) << "- Sb / AU , "
    <<  " g(r) " << enumtostring(A) << "- Te / AU , Median Filtered " << enumtostring(A) << "- Te / AU , intg(r) " + enumtostring(A) << "- Te / AU \n ";
    
    for (int i=1; i<datapoints-1; ++i)
    {
        ost <<  Rvalues[i] << " , " << binAX[i] << " , " << filtervalueAX[i] << " , "<< totaln[i] <<  ","
                                    << binAFe[i] << " , " << filtervalueAFe[i] << " , "<< totalnFe[i] <<  ","
                                    << binAGe[i] << " , " << filtervalueAGe[i] << " , "<< totalnGe[i] <<  ","
                                    << binASb[i] << " , " << filtervalueASb[i] << " , "<< totalnSb[i] <<  ","
                                    << binATe[i] << " , " << filtervalueATe[i] << " , " << totalnTe[i] <<  "\n";
    }
}

void staticanalysis::totalg(int firstframe,  int lastframe, int datapoints, Fileread& file)
{
    vector <double> bins;
    vector <double> totaln;  //vector to determine co-ordination number
    double coordinationno = 0.0;
    
    int numberofframes = lastframe - firstframe;
    
    double rmin = 0.0;
    double rmax = 7.0;
    double deltaR = (rmax - rmin) / datapoints;
    vector <double> Rvalues;
    
    for (double R = rmin; R < rmax; R+=deltaR)
    {
        Rvalues.push_back(R);
    }
    
    
    vector <double> numberinbin (datapoints);
    
    for (int k = 0; k < numberofframes; ++k) {
        
        for (int j = 0; j <file.modelsize; ++j){
            
            for (int i = 0; i <file.modelsize; ++i) {
                
                if (file.neighbourtables[k][i][j] < rmax) {
                    int h = floor((file.neighbourtables[k][i][j]-rmin)/deltaR);
                    numberinbin[h]+=1.0;
                }
            }
        }
    }
    
    
    
    int NFe =0, NGe = 0, NSb = 0, NTe = 0;

    for (int i = 0; i < file.modelsize; ++i)
    {
        if (file.atomtypes[i] == Fe)
            ++NFe;
        if (file.atomtypes[i] == Ge)
            ++NGe;
        if (file.atomtypes[i] == Sb)
            ++NSb;
        if (file.atomtypes[i] == Te)
            ++NTe;
    }
    
    
    for (int i = 1; i < numberinbin.size(); ++i) {
        double g = numberinbin[i]/(4*M_PI*(NFe+NGe+NSb+NTe)*(NFe+NGe+NSb+NTe-1)*Rvalues[i]*Rvalues[i]*deltaR*numberofframes/(dotproduct(file.cellvectors[0], crossproduct(file.cellvectors[1], file.cellvectors[2]))));
        
        bins.push_back(g);
        
        double n = g*Rvalues[i]*Rvalues[i] * deltaR * 4 * M_PI * (NFe+NGe+NSb+NTe) /(dotproduct(file.cellvectors[0], crossproduct(file.cellvectors[1], file.cellvectors[2])));
        
        coordinationno += n;
        totaln.push_back(coordinationno);
        
    }
    
    int window = 10;
    int completewindow = datapoints - (window-1);
    
    vector<double> extendeddata(datapoints + (datapoints-completewindow));
    
    
    int edges = extendeddata.size() - datapoints;
    
    if (edges%2 == 0) {
        
        if (edges == 0) {}
        
        else {
            for (int i=0; i < edges/2; ++i)
            {
                extendeddata[i] = bins[0];
              ;
            }}
        
        if (edges/2+datapoints == extendeddata.size()) {}
        
        else {
            for (int i=edges/2+datapoints; i < extendeddata.size(); ++i)
            {
                extendeddata[i] = bins[datapoints-1];
            
            }}
        
        for (int i = edges/2; i < datapoints+edges/2; ++i)
        {
            extendeddata[i] = bins[i-edges/2];
            
        }
    }
    
    else if(edges%2 ==1) {
        
        for (int i=0; i < floor(edges/2)+1; ++i)
        {
            extendeddata[i] = bins[0];
            
        }
        
        if (floor(edges/2) + 1 + datapoints == extendeddata.size())
        {}
        
        else {
            for (int i=floor(edges/2) + 1 + datapoints; i < extendeddata.size(); ++i)
            {
                extendeddata[i] = bins[datapoints-1];
                
            }
        }
        
        for (int i = floor(edges/2)+1; i < datapoints+floor(edges/2)+1; ++i)
        {
            extendeddata[i] = bins[i-(floor(edges/2)+1)];
           
        }
    }
    
    vector<double> tobesorted;
   
    
    vector<double> filtervalue (datapoints);
   
    
    for(int i =0; i < extendeddata.size() - (window-1); ++i)
    {
        
        tobesorted = extendeddata;
        
        
        sort (tobesorted.begin()+i, tobesorted.begin()+(i+(window)));
        
        if (window%2 == 1)
            filtervalue[i] = tobesorted[(i + floor((window)/2))];
        
        else if (window%2 == 0)
            filtervalue[i] = (tobesorted[i+((window/2)-1)] + tobesorted[i+(window/2)])/2;
        
        tobesorted.clear();

    }

    
    
    double maximumXX =0.0;
    double maxrXX = 0.0;
    double cutoff = 4.00;
    int imax=0;
    for (int i = 1; i < datapoints-1; ++i)
    {
        if (filtervalue[i] > maximumXX && Rvalues[i] < cutoff){
            maximumXX = bins[i];
            maxrXX = Rvalues[i];
            imax = i;
        }
        
    }
    
    double minimumXX =maximumXX;
    double minrXX = 0.0;
    int imin=0;
    
    for (int i = 1; i < datapoints-1; ++i)
    {
        if (bins[i] < minimumXX && Rvalues[i] < cutoff && Rvalues[i] > maxrXX){
            minimumXX = bins[i];
            minrXX = Rvalues[i];
            imin = i;
        }
       
    }
    
    string name = "/Users/joekirk/Desktop/Part III Project - Initial Results/" + file.name + "/RDF/RDF - x-x.csv";
    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    ost << "Firstframe: " << firstframe << ", Lastframe: " << lastframe << ", Numberofframes: " << numberofframes << "\n";
    
    ost << "Data extracted from file: " << file.name << "\n";
    
    ost << "X-X maximum: ," << maximumXX << " , rmax value: ," << maxrXX << ", intg(rmax): ," << totaln[imax]
        << ", X-X minimum: ," << minimumXX << " , rmin value: ," << minrXX << ", intg(rmin): ," << totaln[imin] << "\n\n";

    
    ost <<  "r / A , g(r) X-X / AU , intg(r) x-x / AU\n";
    
    for (int i=1; i<bins.size(); ++i)
    {
        ost <<  Rvalues[i] << " , " << bins[i] << " , " << totaln[i] << "\n";
        
    }
}

void staticanalysis::rdf (int firstframe,  int lastframe, int datapoints, Fileread& file)
{
    g(firstframe, lastframe, datapoints,  Fe,  file);
    g(firstframe, lastframe, datapoints,  Ge,  file);
    g(firstframe, lastframe, datapoints,  Sb,  file);
    g(firstframe, lastframe, datapoints,  Te,  file);
    
}

// -- need to define the cut off within the function ie. the first minimum of the radial distribution function -- done by inspection -- not a great way to do it -- can always alter later.

void staticanalysis::badf (int firstframe, int lastframe, int datapoints, Fileread& file, double cutoff)
{
    double deltatheta =0.0;
    double thetamin = 40.0;
    double thetamax = 180.0;
    
    int numberofframes = lastframe - firstframe;
    
    deltatheta = (thetamax - thetamin) / datapoints;
    vector <double> thetavalues;
    
    for (double theta = thetamin; theta < thetamax; theta += deltatheta)
    {
        thetavalues.push_back(theta);
    }
    
    
    vector <double> numberinbinFe (datapoints);
    vector <double> numberinbinGe (datapoints);
    vector <double> numberinbinSb (datapoints);
    vector <double> numberinbinTe (datapoints);
    
    int NFe =0, NGe = 0, NSb = 0, NTe = 0;

    for (int i = 0; i < file.modelsize; ++i)
    {
        if (file.atomtypes[i] == Fe)
            ++NFe;
        if (file.atomtypes[i] == Ge)
            ++NGe;
        if (file.atomtypes[i] == Sb)
            ++NSb;
        if (file.atomtypes[i] == Te)
            ++NTe;
    }

    
    
    
    for (int k = 0; k < numberofframes; ++ k) {
        
        for (int l = 0; l < file.modelsize; ++l) {
            
            for (int j = 0; j < file.modelsize; ++j) {
                
                for (int i =0; i < file.modelsize; ++i) {
                    
                    if (file.atomtypes[i] == Fe && file.neighbourtables[k][i][j] < cutoff
                        && file.neighbourtables[k][i][l] < cutoff && j != l && j != i && l != i)
                        
                    {
                        double dot = dotproduct(file.rijallframes[k][j][i],file.rijallframes[k][l][i])/(file.neighbourtables[k][i][j]*file.neighbourtables[k][i][l]);
                        
                        double thetai = acos(dot);
                        thetai *= 180.0/M_PI;
                        
                        if (thetai < thetamax && thetai > thetamin) {
                            int a = floor ((thetai - thetamin) / deltatheta);
                            numberinbinFe[a] += 1.0;
                        }
                        
                        
                        
                    }
                    if (file.atomtypes[i] == Ge && file.neighbourtables[k][i][j] < cutoff
                        && file.neighbourtables[k][i][l] < cutoff && j != l && j != i && l != i)
                        
                    {
                        double dot = dotproduct(file.rijallframes[k][j][i],file.rijallframes[k][l][i])/(file.neighbourtables[k][i][j]*file.neighbourtables[k][i][l]);
                        
                        double thetai = acos(dot);
                        thetai *= 180.0/M_PI;
                        
                        if (thetai < thetamax && thetai > thetamin) {
                            int b = floor ((thetai - thetamin) / deltatheta);
                            numberinbinGe[b] += 1.0;
                        }
                        
                        
                        
                    }
                    if (file.atomtypes[i] == Sb && file.neighbourtables[k][i][j] < cutoff
                        && file.neighbourtables[k][i][l] < cutoff && j != l && j != i && l != i)
                        
                    {
                        double dot = dotproduct(file.rijallframes[k][j][i],file.rijallframes[k][l][i])/(file.neighbourtables[k][i][j]*file.neighbourtables[k][i][l]);
                        
                        double thetai = acos(dot);
                        thetai *= 180.0/M_PI;
                        
                        if (thetai < thetamax && thetai > thetamin) {
                            int c = floor ((thetai - thetamin) / deltatheta);
                            numberinbinSb[c] += 1.0;
                        }
                        
                        
                        
                    }
                    if (file.atomtypes[i] == Te && file.neighbourtables[k][i][j] < cutoff
                        && file.neighbourtables[k][i][l] < cutoff && j != l && j != i && l != i)
                        
                    {
                        double dot = dotproduct(file.rijallframes[k][j][i],file.rijallframes[k][l][i])/(file.neighbourtables[k][i][j]*file.neighbourtables[k][i][l]);
                        
                        double thetai = acos(dot);
                        thetai *= 180.0/M_PI;
                        
                        if (thetai < thetamax && thetai > thetamin) {
                            int d = floor ((thetai - thetamin) / deltatheta);
                            numberinbinTe[d] += 1.0;
                        }
                        
                        
                        
                    }
                    
                }
            }
        }
    }
    
    
    for (int i = 0; i < datapoints; ++i)
    {
        numberinbinFe[i] /= (2*numberofframes*NFe);
        numberinbinGe[i] /= (2*numberofframes*NGe);
        numberinbinSb[i] /= (2*numberofframes*NSb);
        numberinbinTe[i] /= (2*numberofframes*NTe);
    }
    
    int window = 50;
    int completewindow = datapoints - (window-1);
    
    vector<double> extendeddataFe(datapoints + (datapoints-completewindow));
    vector<double> extendeddataGe(datapoints + (datapoints-completewindow));
    vector<double> extendeddataSb(datapoints + (datapoints-completewindow));
    vector<double> extendeddataTe(datapoints + (datapoints-completewindow));
     
    int edges = extendeddataFe.size() - datapoints;
    
    if (edges%2 == 0) {
    
        if (edges == 0) {}
        
        else {
            for (int i=0; i < edges/2; ++i)
            {
                extendeddataFe[i] = numberinbinFe[0];
                extendeddataGe[i] = numberinbinGe[0];
                extendeddataSb[i] = numberinbinSb[0];
                extendeddataTe[i] = numberinbinTe[0];
            }}
        
        if (edges/2+datapoints == extendeddataFe.size()) {}
       
        else {
            for (int i=edges/2+datapoints; i < extendeddataFe.size(); ++i)
            {
                extendeddataFe[i] = numberinbinFe[datapoints-1];
                extendeddataGe[i] = numberinbinGe[datapoints-1];
                extendeddataSb[i] = numberinbinSb[datapoints-1];
                extendeddataTe[i] = numberinbinTe[datapoints-1];
            }}
        
        for (int i = edges/2; i < datapoints+edges/2; ++i)
        {
            extendeddataFe[i] = numberinbinFe[i-edges/2];
            extendeddataGe[i] = numberinbinGe[i-edges/2];
            extendeddataSb[i] = numberinbinSb[i-edges/2];
            extendeddataTe[i] = numberinbinTe[i-edges/2];
        }
    }
    
    else if(edges%2 ==1) {
        
        for (int i=0; i < floor(edges/2)+1; ++i)
        {
            extendeddataFe[i] = numberinbinFe[0];
            extendeddataGe[i] = numberinbinGe[0];
            extendeddataSb[i] = numberinbinSb[0];
            extendeddataTe[i] = numberinbinTe[0];
        }
        
        if (floor(edges/2) + 1 + datapoints ==extendeddataFe.size())
        {}
        
        else {
            for (int i=floor(edges/2) + 1 + datapoints; i < extendeddataFe.size(); ++i)
            {
                extendeddataFe[i] = numberinbinFe[datapoints-1];
                extendeddataGe[i] = numberinbinGe[datapoints-1];
                extendeddataSb[i] = numberinbinSb[datapoints-1];
                extendeddataTe[i] = numberinbinTe[datapoints-1];
                
            }
        }
        
        for (int i = floor(edges/2)+1; i < datapoints+floor(edges/2)+1; ++i)
        {
            extendeddataFe[i] = numberinbinFe[i-(floor(edges/2)+1)];
            extendeddataGe[i] = numberinbinGe[i-(floor(edges/2)+1)];
            extendeddataSb[i] = numberinbinSb[i-(floor(edges/2)+1)];
            extendeddataTe[i] = numberinbinTe[i-(floor(edges/2)+1)];
        }
    }
    
    vector<double> tobesortedFe;
    vector<double> tobesortedGe;
    vector<double> tobesortedSb;
    vector<double> tobesortedTe;
    
    vector<double> filtervalueFe (datapoints);
    vector<double> filtervalueGe (datapoints);
    vector<double> filtervalueSb (datapoints);
    vector<double> filtervalueTe (datapoints);
    
    for(int i =0; i < extendeddataFe.size() - (window-1); ++i)
    {
        
        tobesortedFe = extendeddataFe;
        tobesortedGe = extendeddataGe;
        tobesortedSb = extendeddataSb;
        tobesortedTe = extendeddataTe;

        
        sort (tobesortedFe.begin()+i, tobesortedFe.begin()+(i+(window)));
        
        if (window%2 == 1)
            filtervalueFe[i] = tobesortedFe[(i + floor((window)/2))];
        
        else if (window%2 == 0)
            filtervalueFe[i] = (tobesortedFe[i+((window/2)-1)] + tobesortedFe[i+(window/2)])/2;
        
        tobesortedFe.clear();
        
        sort (tobesortedGe.begin()+i, tobesortedGe.begin()+(i+(window)));
        
        if (window%2 == 1)
            filtervalueGe[i] = tobesortedGe[(i + floor((window)/2))];
        
        else if (window%2 == 0)
            filtervalueGe[i] = (tobesortedGe[i+((window/2)-1)] + tobesortedGe[i+(window/2)])/2;
        
        tobesortedGe.clear();
        
        sort (tobesortedSb.begin()+i, tobesortedSb.begin()+(i+(window)));
        
        if (window%2 == 1)
            filtervalueSb[i] = tobesortedSb[(i + floor((window)/2))];
        
        else if (window%2 == 0)
            filtervalueSb[i] = (tobesortedSb[i+((window/2)-1)] + tobesortedSb[i+(window/2)])/2;
        
        tobesortedSb.clear();
        
        sort (tobesortedTe.begin()+i, tobesortedTe.begin()+(i+(window)));
        
        if (window%2 == 1)
            filtervalueTe[i] = tobesortedTe[(i + floor((window)/2))];
        
        else if (window%2 == 0)
            filtervalueTe[i] = (tobesortedTe[i+((window/2)-1)] + tobesortedTe[i+(window/2)])/2;
        
        tobesortedTe.clear();
    }
    
    double maximumFe =0.0, maximumGe =0.0, maximumSb =0.0, maximumTe =0.0;
    double maxthetaFe = 0.0, maxthetaGe = 0.0, maxthetaSb = 0.0, maxthetaTe = 0.0;
    
    for (int i = 0; i < filtervalueFe.size(); ++i)
    {
        
        if (filtervalueFe[i] > maximumFe){
            maximumFe = filtervalueFe[i];
        maxthetaFe = thetavalues[i];
        }
        if (filtervalueGe[i] > maximumGe){
            maximumGe = filtervalueGe[i];
        maxthetaGe = thetavalues[i];
        }
        if (filtervalueSb[i] > maximumSb){
            maximumSb = filtervalueSb[i];
        maxthetaSb = thetavalues[i];
        }
        if (filtervalueTe[i] > maximumTe){
            maximumTe = filtervalueTe[i];
        maxthetaTe = thetavalues[i];
        }
    }

    
    
    string name = "/Users/joekirk/Desktop/Part III Project - Initial Results/" + file.name + "/BADF/cutoff" + convertdouble(cutoff) + "BADF.csv";
    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    ost << "cutoff: " << cutoff << ", Firstframe: " << firstframe << ", Lastframe: " << lastframe << ", Numberofframes: " << numberofframes << "\n";
    ost << "Data extracted from file: " << file.name << "\n";
    ost << "Maximum Fe , " << maximumFe << " , Theta max (Fe) , " << maxthetaFe << "\n"
        << "Maximum Ge , " << maximumGe << " , Theta max (Ge) , " << maxthetaGe << "\n"
        << "Maximum Sb , " << maximumSb << " , Theta max (Sb) , " << maxthetaSb << "\n"
        << "Maximum Te , " << maximumTe << " , Theta max (Te) , " << maxthetaTe << "\n";
    ost << "theta , " << "Fe - x Average number of bonds" << ", Median filtered (Fe) ,"
                      << "Ge - x Average number of bonds" << ", Median filtered (Ge) ,"
                      << "Sb - x Average number of bonds" << ", Median filtered (Sb) ,"
                      << "Te - x Average number of bonds" << ", Median filtered (Te)\n";
    
    
    for (int i=0; i<datapoints; ++i)
    {
        ost <<  thetavalues[i] << " , " << numberinbinFe[i] << "," << filtervalueFe[i]
                               << " , " << numberinbinGe[i] << "," << filtervalueGe[i]
                               << " , " << numberinbinSb[i] << "," << filtervalueSb[i]
                               << " , " << numberinbinTe[i] << "," << filtervalueTe[i] << "\n";
    }
    
}

void staticanalysis::bondangledistributionsingleatom (int firstframe, int lastframe, int datapoints, atmtype A, Fileread& file, double cutoff)
{
    double deltatheta =0.0;
    double thetamin = 40.0;
    double thetamax = 180.0;
    
    int numberofframes = lastframe - firstframe;
    
    deltatheta = (thetamax - thetamin) / datapoints;
    vector <double> thetavalues;
    
    for (double theta = thetamin; theta < thetamax; theta += deltatheta)
    {
        thetavalues.push_back(theta);
    }
    
    vector <double> numberinbin (datapoints);

    
    int Na =0;
    
    for (int i = 0; i < file.modelsize; ++i)
    {
        if (file.atomtypes[i] == A)
            ++Na;
    }
    
    
    
    for (int k = 0; k < numberofframes; ++ k) {
        
        for (int l = 0; l < file.modelsize; ++l) {
            
            for (int j = 0; j < file.modelsize; ++j) {
                
                for (int i =0; i < file.modelsize; ++i) {
                    
                    if (file.atomtypes[i] == A && file.neighbourtables[k][i][j] < cutoff
                        && file.neighbourtables[k][i][l] < cutoff && j != l && j != i && l != i)
                        
                    {
                        double dot = dotproduct(file.rijallframes[k][j][i],file.rijallframes[k][l][i])/(file.neighbourtables[k][i][j]*file.neighbourtables[k][i][l]);
                        
                        double thetai = acos(dot);
                        thetai *= 180.0/M_PI;
                        
                        if (thetai < thetamax && thetai > thetamin) {
                            int h = floor ((thetai - thetamin) / deltatheta);
                            numberinbin[h] += 1.0;
                        }
                        
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < numberinbin.size(); ++i)
    {
        numberinbin[i] /= (2*numberofframes);
    }
    
    int window = 50;
    vector<double> filtervalue (datapoints);
    int completewindow = datapoints - (window-1);
    
    vector<double> extendeddata(datapoints + (datapoints-completewindow));
    int edges = extendeddata.size() - datapoints;
    if (edges%2 == 0) {
        if (edges == 0) {}
        else {
            for (int i=0; i < edges/2; ++i)
            {
                extendeddata[i] = numberinbin[0];
            }}
        
        if (edges/2+datapoints == extendeddata.size()) {}
        else {
            for (int i=edges/2+datapoints; i < extendeddata.size(); ++i)
            {
                extendeddata[i] = numberinbin[datapoints-1];
            }}
        for (int i = edges/2; i < datapoints+edges/2; ++i)
        {
            extendeddata[i] = numberinbin[i-edges/2];
        }
    }
    
    else if(edges%2 ==1) {
        
        for (int i=0; i < floor(edges/2)+1; ++i)
        {
            extendeddata[i] = numberinbin[0];
        }
        
        if (floor(edges/2) + 1 + datapoints ==extendeddata.size())
        {}
        else {
            for (int i=floor(edges/2) + 1 + datapoints; i < extendeddata.size(); ++i)
            {
                extendeddata[i] = numberinbin[datapoints-1];
            }
        }
        
        for (int i = floor(edges/2)+1; i < datapoints+floor(edges/2)+1; ++i)
        {
            extendeddata[i] = numberinbin[i-(floor(edges/2)+1)];
        }
    }
    
    vector<double> tobesorted;
    
    for(int i =0; i < extendeddata.size() - (window-1); ++i)
    {
        
        tobesorted = extendeddata;
        
        sort (tobesorted.begin()+i, tobesorted.begin()+(i+(window)));
        if (window%2 == 1)
            filtervalue[i] = tobesorted[(i + floor((window)/2))];
        else if (window%2 == 0)
            filtervalue[i] = (tobesorted[i+((window/2)-1)] + tobesorted[i+(window/2)])/2;
        
        tobesorted.clear();
    }
    
    double maximum =0.0;
    double maxtheta = 0.0;
    
    for (int i = 0; i < filtervalue.size(); ++i)
    {
        if (filtervalue[i] > maximum){
            maximum = filtervalue[i];
            maxtheta = thetavalues[i];
        }
    }
        
    
    string name = "/Users/joekirk/Desktop/Part III Project - Initial Results/" + file.name +"/BADF/BADF - " + enumtostring(A) + "-x.csv";
    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    ost << "Cutoff: " << cutoff << ", Firstframe: " << firstframe << ", Lastframe: " << lastframe << ", Numberofframes: " << numberofframes << "\n";
    ost << "Data extracted from file: " << file.name << "\n"
        << "Maximum , " << maximum << "Theta max," << maxtheta <<  "\n";
    ost << "theta , " << enumtostring(A) <<  " - x Average number of bonds" << ", Median filtered\n";
    
    for (int i=0; i<numberinbin.size(); ++i)
    {
        ost <<  thetavalues[i] << " , " << numberinbin[i] << "," << filtervalue[i] << "\n";
    }
    
}

void dynamicanalysis::msdforatom(atmtype A, int datapoints, Fileread& file)
{
    
    double timestep= 0.005;  //i ( 0 -> 1 ) is a time step of 5 fs, i.e 0.005 ps
    double rsquaredsum = 0.0;
    int Na = 0;
    
    for (int i = 0; i < file.modelsize; ++i){
        if (file.atomtypes[i] == A)
            ++Na;}
    
    
    vector <int> boundariescrossedx (file.modelsize);
    vector <int> boundariescrossedy (file.modelsize);
    vector <int> boundariescrossedz (file.modelsize);
    vector<double> rsquared;
    
    for (int j = 0; j < file.frames.size(); ++j) {
        
        for (int i = 0; i < file.modelsize; ++i) {
            
            if (file.atomtypes[i] == A) {
                
                Atomcoord ri = Atomcoord(0.0, 0.0, 0.0);
                
                if (j > 1 && (file.frames[j][i].x - file.frames[j-1][i].x) < ((-file.cellvectors[0].x)/2))
                    ++ boundariescrossedx[i];
                else if (j > 1 && (file.frames[j][i].x - file.frames[j-1][i].x) > ((file.cellvectors[0].x)/2))
                    --boundariescrossedx[i];
                else {};
                
                double rx = file.frames[j][i].x + (boundariescrossedx[i]*file.cellvectors[0].x)- file.frames[0][i].x;
                
                if (j > 1 && (file.frames[j][i].y - file.frames[j-1][i].y) < ((-file.cellvectors[1].y)/2))
                    ++boundariescrossedy[i];
                else if (j > 1 && (file.frames[j][i].y - file.frames[j-1][i].y) > ((file.cellvectors[1].y)/2))
                    --boundariescrossedy[i];
                else {};
                
                double ry = file.frames[j][i].y + (boundariescrossedy[i]*file.cellvectors[1].y)- file.frames[0][i].y;
                
                if (j > 1 && (file.frames[j][i].z - file.frames[j-1][i].z) < ((-file.cellvectors[2].z)/2))
                    ++boundariescrossedz[i];
                else if (j > 1 && (file.frames[j][i].z - file.frames[j-1][i].z) > ((file.cellvectors[2].z)/2))
                    --boundariescrossedz[i];
                else{};
                
                double rz = file.frames[j][i].z + (boundariescrossedz[i]*file.cellvectors[2].z)- file.frames[0][i].z;
                
                ri = Atomcoord(rx, ry, rz);
                rsquaredsum += (mod(ri)*mod(ri));
                
            }
            
            
        }
        
        
        rsquaredsum /= Na;
        rsquared.push_back(rsquaredsum);
        rsquaredsum = 0.0; //set rsquaredsum back equal to zero for the next frame
        
    }
    
    
    int increment = floor (rsquared.size()/datapoints);
    
    string name = "/Users/joekirk/Desktop/Part III Project - Initial Results/" + file.name +"/MSD/MSD - " + enumtostring(A) + ".csv";
    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    ost << "Data extracted from file: " << file.name << "\n";
    ost <<  "i , t / ps , " << enumtostring(A) << "-MSD / A2\n";
    
    for (int i=0; i<rsquared.size(); i+=increment)
    {
        ost <<  i << " , " << (i*timestep) + timestep << " , " << rsquared[i] << "\n";
        
    }
    
}


void dynamicanalysis::averageMSD(Fileread& file, int datapoints)
{
    
    double timestep= 0.005;
    double rsquaredsum = 0.0;
    
    vector <int> boundariescrossedx (file.modelsize);
    vector <int> boundariescrossedy (file.modelsize);
    vector <int> boundariescrossedz (file.modelsize);
    vector<double> rsquared;
    
    for (int j = 0; j < file.frames.size(); ++j) {
        
        for (int i = 0; i < file.modelsize; ++i) {
            
            
            
            Atomcoord ri = Atomcoord(0.0, 0.0, 0.0);
            
            if (j > 1 && (file.frames[j][i].x - file.frames[j-1][i].x) < ((-file.cellvectors[0].x)/2))
                ++ boundariescrossedx[i];
            else if (j > 1 && (file.frames[j][i].x - file.frames[j-1][i].x) > ((file.cellvectors[0].x)/2))
                --boundariescrossedx[i];
            else {};
            
            double rx = file.frames[j][i].x + (boundariescrossedx[i]*file.cellvectors[0].x)- file.frames[0][i].x;
            
            if (j > 1 && (file.frames[j][i].y - file.frames[j-1][i].y) < ((-file.cellvectors[1].y)/2))
                ++boundariescrossedy[i];
            else if (j > 1 && (file.frames[j][i].y - file.frames[j-1][i].y) > ((file.cellvectors[1].y)/2))
                --boundariescrossedy[i];
            else {};
            
            double ry = file.frames[j][i].y + (boundariescrossedy[i]*file.cellvectors[1].y)- file.frames[0][i].y;
            
            if (j > 1 && (file.frames[j][i].z - file.frames[j-1][i].z) < ((-file.cellvectors[2].z)/2))
                ++boundariescrossedz[i];
            else if (j > 1 && (file.frames[j][i].z - file.frames[j-1][i].z) > ((file.cellvectors[2].z)/2))
                --boundariescrossedz[i];
            else{};
            
            double rz = file.frames[j][i].z + (boundariescrossedz[i]*file.cellvectors[2].z)- file.frames[0][i].z;
            
            ri = Atomcoord(rx, ry, rz);
            rsquaredsum += (mod(ri)*mod(ri));
            
            
            
            
        }
        
        
        rsquaredsum /= file.modelsize;
        rsquared.push_back(rsquaredsum);
        rsquaredsum = 0.0; //set rsquaredsum back equal to zero for the next frame
        
    }
    
    int increment = floor (rsquared.size()/datapoints);
    
    string name = "/Users/joekirk/Desktop/Part III Project - Initial Results/" + file.name + "/MSD/MSD - Average.csv";
    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    
    ost << "Data extracted from file: " << file.name << "\n";
    ost <<  "i , t / ps , Average MSD / A2\n";
    
    for (int i=0; i<rsquared.size(); i+=increment)
    {
        ost <<  i << " , " << (i*timestep) + timestep  << " , " << rsquared[i] << "\n";
        
    }
    
}

void dynamicanalysis::msd(Fileread& file, int datapoints)
{
    msdforatom(Fe, datapoints, file);
    msdforatom(Ge, datapoints, file);
    msdforatom(Sb, datapoints, file);
    msdforatom(Te, datapoints, file);
    averageMSD(file, datapoints);
}


void dynamicanalysis::fourfoldrings(Fileread& file) {
    
    vector < Numeric_lib::Matrix<double, 2> > neighbourtablecopy = file.neighbourtables;
    

    double cutoff = 3.5;
    double anglemax = 110.0;
    double anglemin = 70.0;
    int numberoffourfoldrings = 0.0;
    vector<int> fourfoldrings;
    int feinfourfoldrings = 0.0;
    vector<int> fefourfoldrings;
    
    vector<vector<int>> vmdouts;
    
    for (int m = 0; m < neighbourtablecopy.size(); ++m) { // k=framenumber => cycle through all relevant frames
        vector<int> vmdout (file.modelsize);
        for (int a = 0; a <file.modelsize; ++a) {  // a=index for atom A
            
            for (int b = 0; b <file.modelsize; ++b) { // b=index for atom B
                
                for (int c = 0; c <file.modelsize; ++c){  // c=index for atom C
                    
                    for (int d = 0; d <file.modelsize; ++d) { // d=index for atom D
                        
                        if (neighbourtablecopy[m][d][c] < cutoff
                            && neighbourtablecopy[m][c][b] < cutoff
                            && neighbourtablecopy[m][b][a] < cutoff
                            && neighbourtablecopy[m][a][d] < cutoff
                            && d != c && d != b && d != a && c !=b && c != a && b != a) {
                            
                            double angleDAB = angledegrees(file.rijallframes[m][d][a], file.rijallframes[m][b][a]);
                            
                            if (angleDAB < anglemax && angleDAB> anglemin) {
                                
                                double angleABC = angledegrees(file.rijallframes[m][a][b], file.rijallframes[m][c][b]);
                                
                                if (angleABC < anglemax && angleABC > anglemin) {
                                    
                                    double angleBCD = angledegrees(file.rijallframes[m][b][c], file.rijallframes[m][d][c]);
                                    
                                    if (angleBCD < anglemax && angleBCD > anglemin) {
                                        
                                        double angleCDA = angledegrees(file.rijallframes[m][c][d], file.rijallframes[m][a][d]);
                                        
                                        if (angleCDA < anglemax && angleCDA > anglemin) {
                                            
                                            //see if the two normals have a small angle between them
                                            double anglenormalABAD = angledegrees(crossproduct(file.rijallframes[m][c][d], file.rijallframes[m][c][b]), crossproduct(file.rijallframes[m][a][b], file.rijallframes[m][a][d]));
                                            
                                            if (anglenormalABAD < 20) {
                                                ++numberoffourfoldrings;
                                                neighbourtablecopy[m][a][b] = 10;
                                                neighbourtablecopy[m][b][a] = 10;
                                                
                                                vmdout[a] += 1;
                                                vmdout[b] += 1;
                                                vmdout[c] += 1;
                                                vmdout[d] += 1;
                                                
                                                if (file.atomtypes[a] == Fe) {
                                                    
                                                    ++feinfourfoldrings;
                                                    
                                                }
                                                
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        fourfoldrings.push_back(numberoffourfoldrings);
        numberoffourfoldrings=0;
        fefourfoldrings.push_back(feinfourfoldrings);
        feinfourfoldrings=0;
        vmdouts.push_back(vmdout);
        vmdout.clear();
        
    }
    
    int numberofframes = file.frames.size();
    double a = numberofframes;
    double b = file.framesinfile;
    
    int increment = floor (b/a);
    double timestep = 0.005;
    
    
    string name = "/Users/joekirk/Desktop/Part III Project - Initial Results/" + file.name + "/Ring_analysis/Fourfoldrings.csv";
    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    
    ost << "Data extracted from file: " << file.name << "\n";
    ost <<  "i , t / ps , Average Number of fourfold rings\n";
    
    for (int i=0; i<fourfoldrings.size(); ++i)
    {
        ost <<  i << " , " << (i*timestep*increment) + timestep  << " , " << fourfoldrings[i] << "\n";
        
    }
    
    ost.close();
    
    string namevmd =  "/Users/joekirk/Desktop/Part III Project - Initial Results/" + file.name + "/Ring_analysis/VMD_fourfold_" + file.name;
    ofstream ost1(namevmd.c_str());
    if (!ost1) cout << "cannot open file " << namevmd;
    ost1 << "AMNISTEPS " << numberofframes << endl;
    ost1 << endl;
    ost1 << "CRYSTAL" << endl;
    ost1 << endl;
    ost1 << "PRIMVEC" << endl;
    for (int i = 0; i < file.cellvectors.size(); ++i)
    {
        ost1 << file.cellvectors[i].x << "\t" << file.cellvectors[i].y << "\t" << file.cellvectors[i].z << endl;
    }
    for (int i = 0; i < vmdouts.size(); ++i)
    {
        ost1 << endl;
        ost1 << "PRIMCOORD " << i+1 << endl;
        ost1 << file.modelsize << " " << 1 << endl;;
        for (int j=0; j < file.modelsize; ++j)
            if (vmdouts[i][j] == 0) {
                ost1 << " " << file.atomtypes[j] << "\t" << 0.0000 << "\t" << 0.0000 <<"\t" << 0.0000 << endl;}
            else {
                ost1 << " " << file.atomtypes[j] << "\t" << file.frames[i][j].x << "\t" <<file.frames[i][j].y <<"\t" << file.frames[i][j].z << endl;}
    }
    
    ost1.close();
    
    string namevmdbackground =  "/Users/joekirk/Desktop/Part III Project - Initial Results/" + file.name + "/Ring_analysis/VMD_background_" + file.name;
    ofstream ost2(namevmdbackground.c_str());
    if (!ost2) cout << "cannot open file " << namevmd;
    ost2 << "AMNISTEPS " << numberofframes << endl;
    ost2 << endl;
    ost2 << "CRYSTAL" << endl;
    ost2 << endl;
    ost2 << "PRIMVEC" << endl;
    for (int i = 0; i < file.cellvectors.size(); ++i)
    {
        ost2 << file.cellvectors[i].x << "\t" << file.cellvectors[i].y << "\t" << file.cellvectors[i].z << endl;
    }
    for (int i = 0; i < vmdouts.size(); ++i)
    {
        ost2 << endl;
        ost2 << "PRIMCOORD " << i+1 << endl;
        ost2 << file.modelsize << " " << 1 << endl;;
        for (int j=0; j < file.modelsize; ++j){
                ost2 << " " << file.atomtypes[j] << "\t" << file.frames[i][j].x << "\t" <<file.frames[i][j].y <<"\t" << file.frames[i][j].z << endl;}
    }
    

    ost2.close();
    
    
}

void dynamicanalysis::wrongbonds(Fileread& file, int datapoints, int firstframe, int lastframe, double cutoff){
    
    int NFe =0, NGe = 0, NSb = 0, NTe = 0;
    
    for (int i = 0; i < file.modelsize; ++i)
    {
        if (file.atomtypes[i] == Fe)
            ++NFe;
        if (file.atomtypes[i] == Ge)
            ++NGe;
        if (file.atomtypes[i] == Sb)
            ++NSb;
        if (file.atomtypes[i] == Te)
            ++NTe;
    }
    
    
    
    int numberofframes =  lastframe - firstframe;
    double timestep = 0.005;
    int increment = floor (numberofframes/datapoints);
    
    int FeFe=0, FeGe=0, FeSb=0, FeTe=0, GeFe = 0, GeGe=0,GeSb=0, GeTe=0, SbFe = 0, SbGe = 0, SbSb=0, SbTe=0, TeFe = 0, TeGe = 0, TeSb = 0, TeTe=0;
    vector<double> totalFeFe, totalFeGe, totalFeSb, totalFeTe ,totalGeFe, totalGeGe, totalGeSb, totalGeTe, totalSbFe, totalSbGe, totalSbSb, totalSbTe, totalTeFe, totalTeGe, totalTeSb, totalTeTe;
    
    vector<double> fractionofwrongbonds;
    
    for (int m = 0; m < numberofframes; m+=increment){
        
        for (int j = 0; j < file.modelsize; ++j){
            
            for (int i = 0; i < file.modelsize; ++i) {
                
                if (file.atomtypes[j] == Fe && file.atomtypes[i] == Fe && i != j && file.neighbourtables[m][j][i] < cutoff) ++FeFe;
                if (file.atomtypes[j] == Fe && file.atomtypes[i] == Ge && i != j && file.neighbourtables[m][j][i] < cutoff) ++FeGe;
                if (file.atomtypes[j] == Fe && file.atomtypes[i] == Sb && i != j && file.neighbourtables[m][j][i] < cutoff) ++FeSb;
                if (file.atomtypes[j] == Fe && file.atomtypes[i] == Te && i != j && file.neighbourtables[m][j][i] < cutoff) ++FeTe;
                
                if (file.atomtypes[j] == Ge && file.atomtypes[i] == Fe && i != j && file.neighbourtables[m][j][i] < cutoff) ++GeFe;
                if (file.atomtypes[j] == Ge && file.atomtypes[i] == Ge && i != j && file.neighbourtables[m][j][i] < cutoff) ++GeGe;
                if (file.atomtypes[j] == Ge && file.atomtypes[i] == Sb && i != j && file.neighbourtables[m][j][i] < cutoff) ++GeSb;
                if (file.atomtypes[j] == Ge && file.atomtypes[i] == Te && i != j && file.neighbourtables[m][j][i] < cutoff) ++GeTe;
                
                if (file.atomtypes[j] == Sb && file.atomtypes[i] == Fe && i != j && file.neighbourtables[m][j][i] < cutoff) ++SbFe;
                if (file.atomtypes[j] == Sb && file.atomtypes[i] == Ge && i != j && file.neighbourtables[m][j][i] < cutoff) ++SbGe;
                if (file.atomtypes[j] == Sb && file.atomtypes[i] == Sb && i != j && file.neighbourtables[m][j][i] < cutoff) ++SbSb;
                if (file.atomtypes[j] == Sb && file.atomtypes[i] == Te && i != j && file.neighbourtables[m][j][i] < cutoff) ++SbTe;
                
                
                if (file.atomtypes[j] == Te && file.atomtypes[i] == Fe && i != j && file.neighbourtables[m][j][i] < cutoff) ++TeFe;
                if (file.atomtypes[j] == Te && file.atomtypes[i] == Ge && i != j && file.neighbourtables[m][j][i] < cutoff) ++TeGe;
                if (file.atomtypes[j] == Te && file.atomtypes[i] == Sb && i != j && file.neighbourtables[m][j][i] < cutoff) ++TeSb;
                if (file.atomtypes[j] == Te && file.atomtypes[i] == Te && i != j && file.neighbourtables[m][j][i] < cutoff) ++TeTe;
            }
        }
        
        double numberofwrongbonds = FeGe +FeSb +GeSb + GeGe + SbSb + TeTe;
        double totalnumberofbonds = FeGe +FeSb +GeSb + GeGe + SbSb + TeTe + FeFe + FeTe+ GeTe + SbTe;
        double fraction = numberofwrongbonds/totalnumberofbonds;
        
        totalFeFe.push_back((double)FeFe/NFe);
        totalFeGe.push_back((double)FeGe/NFe);
        totalFeSb.push_back((double)FeSb/NFe);
        totalFeTe.push_back((double)FeTe/NFe);
        totalGeFe.push_back((double)GeFe/NGe);
        totalGeGe.push_back((double)GeGe/NGe);
        totalGeSb.push_back((double)GeSb/NGe);
        totalGeTe.push_back((double)GeTe/NGe);
        totalSbFe.push_back((double)SbFe/NSb);
        totalSbGe.push_back((double)SbGe/NSb);
        totalSbSb.push_back((double)SbSb/NSb);
        totalSbTe.push_back((double)SbTe/NSb);
        totalTeFe.push_back((double)TeFe/NTe);
        totalTeGe.push_back((double)TeGe/NTe);
        totalTeSb.push_back((double)TeSb/NTe);
        totalTeTe.push_back((double)TeTe/NTe);
        fractionofwrongbonds.push_back(fraction);
        
        
        
        FeFe=0, FeGe=0, FeSb=0, FeTe=0, GeFe = 0, GeGe=0,GeSb=0, GeTe=0, SbFe=0, SbGe=0, SbSb=0, SbTe=0, TeFe = 0, TeGe = 0, TeSb = 0, TeTe=0, numberofwrongbonds = 0, totalnumberofbonds=0;
    }
    
    //convert to doubles before the division of the two, otherwise will get an int
    
    string name = "/Users/joekirk/Desktop/Part III Project - Initial Results/" +file.name + "/Bond_Statistics/Dynamic/Wrongbonds " + file.name + ".csv";
    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    ost << "Data extracted from file: " << file.name << "\n";
    
    
    ost <<  "i ,time / s, FeFe, FeGe, FeSb, FeTe, GeFe, GeGe, GeSb, GeTe,SbFe, SbGe, SbSb, SbTe,TeFe, TeGe, TeSb, TeTe, Fraction of wrong bonds\n";
    
    for (int i=0; i< totalFeFe.size(); ++i)
    {
        ost <<  i << " , " << (i*timestep*increment) + timestep << " , "
        <<  totalFeFe[i] << " , " << totalFeGe[i] << " , " << totalFeSb[i] << " , " << totalFeTe[i] << " , "
        <<  totalGeFe[i] << " , " << totalGeGe[i] << " , " << totalGeSb[i] << " , " << totalGeTe[i] << " , "
        <<  totalSbFe[i] << " , " << totalSbGe[i] << " , " << totalSbSb[i] << " , " << totalSbTe[i] << " , "
        <<  totalTeFe[i] << " , " << totalTeGe[i] << " , " << totalTeSb[i] << " , " << totalTeTe[i] << " , "
        <<  fractionofwrongbonds[i] << "\n";
        
    }
}

void dynamicanalysis::wrongbonds(Fileread& file, double cutoff){
    
    int NFe =0, NGe = 0, NSb = 0, NTe = 0;
    
    for (int i = 0; i < file.modelsize; ++i)
    {
        if (file.atomtypes[i] == Fe)
            ++NFe;
        if (file.atomtypes[i] == Ge)
            ++NGe;
        if (file.atomtypes[i] == Sb)
            ++NSb;
        if (file.atomtypes[i] == Te)
            ++NTe;
    }

    
    
    int numberofframes =  file.frames.size();
    double timestep = 0.005;
    
    int FeFe=0, FeGe=0, FeSb=0, FeTe=0, GeFe = 0, GeGe=0,GeSb=0, GeTe=0, SbFe = 0, SbGe = 0, SbSb=0, SbTe=0, TeFe = 0, TeGe = 0, TeSb = 0, TeTe=0;
    vector<double> totalFeFe, totalFeGe, totalFeSb, totalFeTe ,totalGeFe, totalGeGe, totalGeSb, totalGeTe, totalSbFe, totalSbGe, totalSbSb, totalSbTe, totalTeFe, totalTeGe, totalTeSb, totalTeTe;
    
    vector<double> fractionofwrongbonds;
    
    for (int m = 0; m < numberofframes; ++m){
        
        for (int j = 0; j < file.modelsize; ++j){
            
            for (int i = 0; i < file.modelsize; ++i) {
                
                if (file.atomtypes[j] == Fe && file.atomtypes[i] == Fe && i != j && file.neighbourtables[m][j][i] < cutoff) ++FeFe;
                if (file.atomtypes[j] == Fe && file.atomtypes[i] == Ge && i != j && file.neighbourtables[m][j][i] < cutoff) ++FeGe;
                if (file.atomtypes[j] == Fe && file.atomtypes[i] == Sb && i != j && file.neighbourtables[m][j][i] < cutoff) ++FeSb;
                if (file.atomtypes[j] == Fe && file.atomtypes[i] == Te && i != j && file.neighbourtables[m][j][i] < cutoff) ++FeTe;
                
                if (file.atomtypes[j] == Ge && file.atomtypes[i] == Fe && i != j && file.neighbourtables[m][j][i] < cutoff) ++GeFe;
                if (file.atomtypes[j] == Ge && file.atomtypes[i] == Ge && i != j && file.neighbourtables[m][j][i] < cutoff) ++GeGe;
                if (file.atomtypes[j] == Ge && file.atomtypes[i] == Sb && i != j && file.neighbourtables[m][j][i] < cutoff) ++GeSb;
                if (file.atomtypes[j] == Ge && file.atomtypes[i] == Te && i != j && file.neighbourtables[m][j][i] < cutoff) ++GeTe;
                
                if (file.atomtypes[j] == Sb && file.atomtypes[i] == Fe && i != j && file.neighbourtables[m][j][i] < cutoff) ++SbFe;
                if (file.atomtypes[j] == Sb && file.atomtypes[i] == Ge && i != j && file.neighbourtables[m][j][i] < cutoff) ++SbGe;
                if (file.atomtypes[j] == Sb && file.atomtypes[i] == Sb && i != j && file.neighbourtables[m][j][i] < cutoff) ++SbSb;
                if (file.atomtypes[j] == Sb && file.atomtypes[i] == Te && i != j && file.neighbourtables[m][j][i] < cutoff) ++SbTe;
                
                
                if (file.atomtypes[j] == Te && file.atomtypes[i] == Fe && i != j && file.neighbourtables[m][j][i] < cutoff) ++TeFe;
                if (file.atomtypes[j] == Te && file.atomtypes[i] == Ge && i != j && file.neighbourtables[m][j][i] < cutoff) ++TeGe;
                if (file.atomtypes[j] == Te && file.atomtypes[i] == Sb && i != j && file.neighbourtables[m][j][i] < cutoff) ++TeSb;
                if (file.atomtypes[j] == Te && file.atomtypes[i] == Te && i != j && file.neighbourtables[m][j][i] < cutoff) ++TeTe;
            }
        }
        
        double numberofwrongbonds = FeGe +FeSb +GeSb + GeGe + SbSb + TeTe;
        double totalnumberofbonds = FeGe +FeSb +GeSb + GeGe + SbSb + TeTe + FeFe + FeTe+ GeTe + SbTe;
        double fraction = numberofwrongbonds/totalnumberofbonds;
        
        totalFeFe.push_back((double)FeFe/NFe);
        totalFeGe.push_back((double)FeGe/NFe);
        totalFeSb.push_back((double)FeSb/NFe);
        totalFeTe.push_back((double)FeTe/NFe);
        totalGeFe.push_back((double)GeFe/NGe);
        totalGeGe.push_back((double)GeGe/NGe);
        totalGeSb.push_back((double)GeSb/NGe);
        totalGeTe.push_back((double)GeTe/NGe);
        totalSbFe.push_back((double)SbFe/NSb);
        totalSbGe.push_back((double)SbGe/NSb);
        totalSbSb.push_back((double)SbSb/NSb);
        totalSbTe.push_back((double)SbTe/NSb);
        totalTeFe.push_back((double)TeFe/NTe);
        totalTeGe.push_back((double)TeGe/NTe);
        totalTeSb.push_back((double)TeSb/NTe);
        totalTeTe.push_back((double)TeTe/NTe);
        fractionofwrongbonds.push_back(fraction);
        
        
        
        FeFe=0, FeGe=0, FeSb=0, FeTe=0, GeFe = 0, GeGe=0,GeSb=0, GeTe=0, SbFe=0, SbGe=0, SbSb=0, SbTe=0, TeFe = 0, TeGe = 0, TeSb = 0, TeTe=0, numberofwrongbonds = 0, totalnumberofbonds=0;
    }
    
    //convert to doubles before the division of the two, otherwise will get an int
    double a = numberofframes;
    double b = file.framesinfile;
    
    int increment = floor (b/a);
    
    string name = "/Users/joekirk/Desktop/Part III Project - Initial Results/" +file.name + "/Bond_Statistics/Dynamic/Wrongbonds " + file.name + ".csv";
    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    ost << "Data extracted from file: " << file.name << "\n";
    
    
    ost <<  "i ,time / s, FeFe, FeGe, FeSb, FeTe, GeFe, GeGe, GeSb, GeTe,SbFe, SbGe, SbSb, SbTe,TeFe, TeGe, TeSb, TeTe, Fraction of wrong bonds\n";
    
    for (int i=0; i< numberofframes; ++i)
    {
        ost <<  i << " , " << (i*timestep*increment) + timestep << " , "
        <<  totalFeFe[i] << " , " << totalFeGe[i] << " , " << totalFeSb[i] << " , " << totalFeTe[i] << " , "
        <<  totalGeFe[i] << " , " << totalGeGe[i] << " , " << totalGeSb[i] << " , " << totalGeTe[i] << " , "
        <<  totalSbFe[i] << " , " << totalSbGe[i] << " , " << totalSbSb[i] << " , " << totalSbTe[i] << " , "
        <<  totalTeFe[i] << " , " << totalTeGe[i] << " , " << totalTeSb[i] << " , " << totalTeTe[i] << " , "
        <<  fractionofwrongbonds[i] << "\n";
        
    }
}

void staticanalysis::bondstatistics(Fileread file, int datapoints,  int firstframe, int lastframe, double cutoff)
{
    int numberofframes =  lastframe - firstframe;
    double timestep = 0.005;
    int increment = floor (numberofframes/datapoints);

    int NFe =0, NGe = 0, NSb = 0, NTe = 0;
    
    for (int i = 0; i < file.modelsize; ++i)
    {
        if (file.atomtypes[i] == Fe)
            ++NFe;
        if (file.atomtypes[i] == Ge)
            ++NGe;
        if (file.atomtypes[i] == Sb)
            ++NSb;
        if (file.atomtypes[i] == Te)
            ++NTe;
    }
    
    int FeFe=0, FeGe=0, FeSb=0, FeTe=0, GeFe = 0, GeGe=0,GeSb=0, GeTe=0, SbFe = 0, SbGe = 0, SbSb=0, SbTe=0, TeFe = 0, TeGe = 0, TeSb = 0, TeTe=0;
    vector<double> totalFeFe, totalFeGe, totalFeSb, totalFeTe ,totalGeFe, totalGeGe, totalGeSb, totalGeTe, totalSbFe, totalSbGe, totalSbSb, totalSbTe, totalTeFe, totalTeGe, totalTeSb, totalTeTe;
    
    vector<double> fractionofwrongbonds;
    
    for (int m = 0; m < numberofframes; m+=increment){
        
        for (int j = 0; j < file.modelsize; ++j){
            
            for (int i = 0; i < file.modelsize; ++i) {
                
                if (file.atomtypes[j] == Fe && file.atomtypes[i] == Fe && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++FeFe;
                if (file.atomtypes[j] == Fe && file.atomtypes[i] == Ge && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++FeGe;
                if (file.atomtypes[j] == Fe && file.atomtypes[i] == Sb && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++FeSb;
                if (file.atomtypes[j] == Fe && file.atomtypes[i] == Te && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++FeTe;
                
                if (file.atomtypes[j] == Ge && file.atomtypes[i] == Fe && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++GeFe;
                if (file.atomtypes[j] == Ge && file.atomtypes[i] == Ge && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++GeGe;
                if (file.atomtypes[j] == Ge && file.atomtypes[i] == Sb && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++GeSb;
                if (file.atomtypes[j] == Ge && file.atomtypes[i] == Te && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++GeTe;
                
                if (file.atomtypes[j] == Sb && file.atomtypes[i] == Fe && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++SbFe;
                if (file.atomtypes[j] == Sb && file.atomtypes[i] == Ge && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++SbGe;
                if (file.atomtypes[j] == Sb && file.atomtypes[i] == Sb && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++SbSb;
                if (file.atomtypes[j] == Sb && file.atomtypes[i] == Te && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++SbTe;
                
                
                if (file.atomtypes[j] == Te && file.atomtypes[i] == Fe && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++TeFe;
                if (file.atomtypes[j] == Te && file.atomtypes[i] == Ge && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++TeGe;
                if (file.atomtypes[j] == Te && file.atomtypes[i] == Sb && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++TeSb;
                if (file.atomtypes[j] == Te && file.atomtypes[i] == Te && i != j &&  file.neighbourtables[m][j][i] < cutoff) ++TeTe;
            }
        }
        
        double numberofwrongbonds = FeGe +FeSb +GeSb + GeGe + SbSb + TeTe;
        double totalnumberofbonds = FeGe +FeSb +GeSb + GeGe + SbSb + TeTe + FeFe + FeTe+ GeTe + SbTe;
        double fraction = numberofwrongbonds/totalnumberofbonds;
        
        totalFeFe.push_back((double)FeFe/NFe);
        totalFeGe.push_back((double)FeGe/NFe);
        totalFeSb.push_back((double)FeSb/NFe);
        totalFeTe.push_back((double)FeTe/NFe);
        totalGeFe.push_back((double)GeFe/NGe);
        totalGeGe.push_back((double)GeGe/NGe);
        totalGeSb.push_back((double)GeSb/NGe);
        totalGeTe.push_back((double)GeTe/NGe);
        totalSbFe.push_back((double)SbFe/NSb);
        totalSbGe.push_back((double)SbGe/NSb);
        totalSbSb.push_back((double)SbSb/NSb);
        totalSbTe.push_back((double)SbTe/NSb);
        totalTeFe.push_back((double)TeFe/NTe);
        totalTeGe.push_back((double)TeGe/NTe);
        totalTeSb.push_back((double)TeSb/NTe);
        totalTeTe.push_back((double)TeTe/NTe);
        fractionofwrongbonds.push_back(fraction);
        
        
        
        FeFe=0, FeGe=0, FeSb=0, FeTe=0, GeFe = 0, GeGe=0,GeSb=0, GeTe=0, SbFe=0, SbGe=0, SbSb=0, SbTe=0, TeFe = 0, TeGe = 0, TeSb = 0, TeTe=0, numberofwrongbonds = 0, totalnumberofbonds=0;
    }
    
;
    
    string name = "/Users/joekirk/Desktop/Part III Project - Initial Results/" +file.name + "/Bond_Statistics/Static/Bondstatistics " + file.name + ".csv";
    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    ost << "Data extracted from file: " << file.name << "\n";
    ost << "Cutoff: " << cutoff << "\n";
    
    
    ost <<  "i ,time / s, FeFe, FeGe, FeSb, FeTe, GeFe, GeGe, GeSb, GeTe,SbFe, SbGe, SbSb, SbTe,TeFe, TeGe, TeSb, TeTe, Fraction of wrong bonds\n";
    
    for (int i=0; i < totalFeFe.size() ; ++i)
    {
        ost <<  i << " , " << (i*timestep*increment) + timestep << " , "
        <<  totalFeFe[i] << " , " << totalFeGe[i] << " , " << totalFeSb[i] << " , " << totalFeTe[i] << " , "
        <<  totalGeFe[i] << " , " << totalGeGe[i] << " , " << totalGeSb[i] << " , " << totalGeTe[i] << " , "
        <<  totalSbFe[i] << " , " << totalSbGe[i] << " , " << totalSbSb[i] << " , " << totalSbTe[i] << " , "
        <<  totalTeFe[i] << " , " << totalTeGe[i] << " , " << totalTeSb[i] << " , " << totalTeTe[i] << " , "
        <<  fractionofwrongbonds[i] << "\n";
        
    }
}


void dynamicanalysis::VMDviewer(Fileread& file, int atomtocentre){
    
    vector < vector<Atomcoord> > newframes;
    vector<Atomcoord> newframe (file.modelsize, Atomcoord(0.0, 0.0, 0.0));
    
    //calculate the new frame 1 placing Fe at the centre
    //  ---** This only works for cubic unit cell!!!----**
    for (int j = 0; j < file.frames.size(); ++j){
        
        for (int i = 0; i < file.modelsize; ++i) {
            
            
            double xmove, ymove, zmove;
            
            xmove = file.cellvectors[0].x/2 - file.frames[j][atomtocentre].x;
            ymove = file.cellvectors[1].y/2 - file.frames[j][atomtocentre].y;
            zmove = file.cellvectors[2].z/2 - file.frames[j][atomtocentre].z;
            
            
            if (file.frames[j][i].x +xmove > file.cellvectors[0].x)
                newframe[i].x = file.frames[j][i].x +xmove - file.cellvectors[0].x;
            else
                newframe[i].x = file.frames[j][i].x +xmove;
            
            
            if (file.frames[j][i].y +ymove > file.cellvectors[1].y)
                newframe[i].y = file.frames[j][i].y +ymove - file.cellvectors[0].y;
            else
                newframe[i].y = file.frames[j][i].y +ymove;
            
            
            if (file.frames[j][i].z +zmove > file.cellvectors[2].z)
                newframe[i].z = file.frames[j][i].z +zmove - file.cellvectors[0].z;
            else
                newframe[i].z = file.frames[j][i].z +zmove;
            
        }
        
        newframes.push_back(newframe);
        
        newframe.clear();
        newframe.resize(file.modelsize, Atomcoord(0.0, 0.0, 0.0));
    }
    
    
    
    string name = "/Users/joekirk/Desktop/Part III Project - Initial Results/" +file.name + "/VMDviewer/AtomcentredVMD_" + file.name;
    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    ost << "AMNISTEPS " << newframes.size() << endl;
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
            ost << " " << file.atomtypes[j] << "\t" << newframes[i][j].x << "\t" <<newframes[i][j].y <<"\t" << newframes[i][j].z << endl;
    }
}
    
   

void dynamicanalysis::maxfourieranalysis(Fileread& file)
{
    int Numberofcells = 32;
    
    fftw_complex signal[Numberofcells][Numberofcells][Numberofcells];
    fftw_complex result[Numberofcells][Numberofcells][Numberofcells];
    
    fftw_plan plan = fftw_plan_dft_3d(Numberofcells, Numberofcells, Numberofcells, signal[0][0], result[0][0], FFTW_FORWARD, FFTW_PATIENT);
    
    vector<double> maxcoefficient;
    double max = 0;
    
    double x = (file.cellvectors[0].x/32.0);
    double y = (file.cellvectors[1].y/32.0);
    double z = (file.cellvectors[2].z/32.0);
    
    for (int k = 0; k < file.frames.size(); ++k) {
        
        for (int l = 0; l < Numberofcells; ++l){
            for (int m = 0; m < Numberofcells; ++m){
                for (int n = 0; n < Numberofcells; ++n){
                    
                    signal[l][m][n][REAL] = 0;
                    signal[l][m][n][IMAG] = 0;
                    result[l][m][n][REAL] = 0;
                    result[l][m][n][IMAG] = 0;
                }
            }
        }
        
        for (int j = 0; j < file.modelsize; ++j){
            
            int indexx = floor((file.frames[k][j].x)/x);
            if (indexx >= 32)
                indexx -= 1;
            
            int indexy = floor((file.frames[k][j].y)/y);
            if (indexy >= 32)
                indexy -= 1;
            
            int indexz = floor((file.frames[k][j].z)/z);
            if (indexz >= 32)
                indexz -= 1;
            
            signal[indexx][indexy][indexz][REAL] += 1.0;
        }
        
        
        fftw_execute(plan);
        
        for (int p = 0; p < Numberofcells; ++p) {
            for (int q = 0; q < Numberofcells; ++q) {
                for (int r = 0; r < Numberofcells; ++r) {
                    if (p == 0 && q ==0 && r ==0){}
                    else{
                        double mag = sqrt((result[p][q][r][REAL] * result[p][q][r][REAL]) +
                                          (result[p][q][r][IMAG] * result[p][q][r][IMAG]));
                        
                        
                        
                        
                        if (mag > max)
                            max = mag;
                        
                    }
                }
            }
        }
        maxcoefficient.push_back(max);
        max = 0;
    }
    
    fftw_destroy_plan(plan);
    
    
    int numberofframes =  file.frames.size();
    double timestep = 0.005;
    double a = numberofframes;
    double b = file.framesinfile;
    
    int increment = floor (b/a);
    
    string name = "/Users/joekirk/Desktop/Part III Project - Initial Results/" +file.name + "/Fourieranalysis/Fouriermax" + file.name + ".csv";    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    ost << "Data extracted from file: " << file.name << "\n";
    
    ost <<  "i ,time / s, Max Fourier Coefficient\n";
    
    
    for (int i=0; i < maxcoefficient.size() ; ++i)
    {
        ost <<  i << " , " << (i*timestep*increment) + timestep << " , " << maxcoefficient[i] << "\n";
    }

        
}
    
    

void dynamicanalysis::nearestneighbourevolution(Fileread& file)
{
    
    int numberofneighbours = 10;
    int NFe =0, NGe = 0, NSb = 0, NTe = 0;
    
    for (int i = 0; i < file.modelsize; ++i)
    {
        if (file.atomtypes[i] == Fe)
            ++NFe;
        if (file.atomtypes[i] == Ge)
            ++NGe;
        if (file.atomtypes[i] == Sb)
            ++NSb;
        if (file.atomtypes[i] == Te)
            ++NTe;
    }
    
    //define a vector to be used over and over for sorting
    
    vector<double> tosort(file.modelsize);
    
    vector<double> Fenearest(numberofneighbours);
    
    vector<double> Genearest(numberofneighbours);
    
    vector<double> Sbnearest(numberofneighbours);
    
    vector<double> Tenearest(numberofneighbours);
    
    vector <vector<double>> Fedata, Gedata, Sbdata, Tedata;
    
    
    
    for (int k = 0; k < file.neighbourtables.size(); ++k) {
        
        for (int i = 0; i < file.modelsize; ++i)
        {
            
            if (file.atomtypes[i] == Fe)
            {
                for(int j = 0; j < file.modelsize; ++j){
                    tosort[j] = file.neighbourtables[k][i][j];}
                
                
                sort(tosort.begin(), tosort.end());
                Fenearest[0]    += tosort[1];
                Fenearest[1]    += tosort[2];
                Fenearest[2]    += tosort[3];
                Fenearest[3]    += tosort[4];
                Fenearest[4]    += tosort[5];
                Fenearest[5]    += tosort[6];
                Fenearest[6]    += tosort[7];
                Fenearest[7]    += tosort[8];
                Fenearest[8]    += tosort[9];
                Fenearest[9]    += tosort[10];
   
                
                tosort.clear();
                tosort.resize(file.modelsize);
                
            }
            
            
            if (file.atomtypes[i] == Ge)
            {
                for(int j = 0; j < file.modelsize; ++j){
                    tosort[j] = file.neighbourtables[k][i][j];}
                
                
                sort(tosort.begin(), tosort.end());
                Genearest[0]    += tosort[1];
                Genearest[1]    += tosort[2];
                Genearest[2]    += tosort[3];
                Genearest[3]    += tosort[4];
                Genearest[4]    += tosort[5];
                Genearest[5]    += tosort[6];
                Genearest[6]    += tosort[7];
                Genearest[7]    += tosort[8];
                Genearest[8]    += tosort[9];
                Genearest[9]    += tosort[10];
                
                tosort.clear();
                tosort.resize(file.modelsize);
            }
            
            
            
            if (file.atomtypes[i] == Sb)
            {
                for(int j = 0; j < file.modelsize; ++j){
                    tosort[j] = file.neighbourtables[k][i][j];}
                
                
                sort(tosort.begin(), tosort.end());
                Sbnearest[0]    += tosort[1];
                Sbnearest[1]    += tosort[2];
                Sbnearest[2]    += tosort[3];
                Sbnearest[3]    += tosort[4];
                Sbnearest[4]    += tosort[5];
                Sbnearest[5]    += tosort[6];
                Sbnearest[6]    += tosort[7];
                Sbnearest[7]    += tosort[8];
                Sbnearest[8]    += tosort[9];
                Sbnearest[9]    += tosort[10];
                
                tosort.clear();
                tosort.resize(file.modelsize);
            }
            
            
            if (file.atomtypes[i] == Te)
            {
                for(int j = 0; j < file.modelsize; ++j){
                    tosort[j] = file.neighbourtables[k][i][j];}
                
                
                sort(tosort.begin(), tosort.end());
                Tenearest[0]    += tosort[1];
                Tenearest[1]    += tosort[2];
                Tenearest[2]    += tosort[3];
                Tenearest[3]    += tosort[4];
                Tenearest[4]    += tosort[5];
                Tenearest[5]    += tosort[6];
                Tenearest[6]    += tosort[7];
                Tenearest[7]    += tosort[8];
                Tenearest[8]    += tosort[9];
                Tenearest[9]    += tosort[10];
                
                tosort.clear();
                tosort.resize(file.modelsize);
            }
            
            
        }
        
        for (int i =0; i < 10; ++i)
        {
            Fenearest[i] /= NFe;
            Genearest[i] /= NGe;
            Sbnearest[i] /= NSb;
            Tenearest[i] /= NTe;
        }
        
        Fedata.push_back(Fenearest);
        Gedata.push_back(Genearest);
        Sbdata.push_back(Sbnearest);
        Tedata.push_back(Tenearest);
        
        Fenearest.clear();
        Fenearest.resize(10);
        
        Genearest.clear();
        Genearest.resize(10);
        
        Sbnearest.clear();
        Sbnearest.resize(10);
        
        Tenearest.clear();
        Tenearest.resize(10);
        
    }
    
    int numberofframes =  file.frames.size();
    double timestep = 0.005;
    
    double a = numberofframes;
    double b = file.framesinfile;
    
    int increment = floor (b/a);
    
    
    string name = "/Users/joekirk/Desktop/Part III Project - Initial Results/" +file.name + "/Nearest_Neighbour_Evolution/FeEvolution" + file.name + ".csv";    
    ofstream ost(name.c_str());
    if (!ost) cout << "cannot open file " << name;
    
    ost << "Data extracted from file: " << file.name << "\n";
    ost <<  "i ,time / s, Fe1, Fe2, Fe3, Fe4, Fe5, Fe6, Fe7, Fe8, Fe9, Fe10 \n";
    
    
    for (int i=0; i < file.frames.size() ; ++i)
        
    {
        ost <<  i << " , " << (i*timestep*increment) + timestep << " , "
        <<  Fedata[i][0] << " , " << Fedata[i][1] << " , " << Fedata[i][2] << " , " << Fedata[i][3] << " , "
        <<  Fedata[i][4] << " , " << Fedata[i][5] << " , " << Fedata[i][6] << " , " << Fedata[i][7] << " , "
        <<  Fedata[i][8] << " , " << Fedata[i][9] <<  "\n";
    }
    
    ost.close();
    
    string name1 = "/Users/joekirk/Desktop/Part III Project - Initial Results/" +file.name + "/Nearest_Neighbour_Evolution/GeEvolution" + file.name + ".csv";
    ofstream ost1(name1.c_str());
    if (!ost1) cout << "cannot open file " << name1;
    
    ost1 << "Data extracted from file: " << file.name << "\n";
    ost1 <<  "i ,time / s, Ge1, Ge2, Ge3, Ge4, Ge5, Ge6, Ge7, Ge8, Ge9, Ge10 \n";
    
    
    for (int i=0; i < file.frames.size() ; ++i)
        
    {
        ost1 <<  i << " , " << (i*timestep*increment) + timestep << " , "
        <<  Gedata[i][0] << " , " << Gedata[i][1] << " , " << Gedata[i][2] << " , " << Gedata[i][3] << " , "
        <<  Gedata[i][4] << " , " << Gedata[i][5] << " , " << Gedata[i][6] << " , " << Gedata[i][7] << " , "
        <<  Gedata[i][8] << " , " << Gedata[i][9] <<  "\n";
    }
    
    ost1.close();
    
    string name2 = "/Users/joekirk/Desktop/Part III Project - Initial Results/" +file.name + "/Nearest_Neighbour_Evolution/SbEvolution" + file.name + ".csv";
    ofstream ost2(name2.c_str());
    if (!ost2) cout << "cannot open file " << name2;
    
    ost2 << "Data extracted from file: " << file.name << "\n";
    ost2 <<  "i ,time / s, Sb1, Sb2, Sb3, Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10 \n";
    
    
    for (int i=0; i < file.frames.size() ; ++i)
        
    {
        ost2 <<  i << " , " << (i*timestep*increment) + timestep << " , "
        <<  Sbdata[i][0] << " , " << Sbdata[i][1] << " , " << Sbdata[i][2] << " , " << Sbdata[i][3] << " , "
        <<  Sbdata[i][4] << " , " << Sbdata[i][5] << " , " << Sbdata[i][6] << " , " << Sbdata[i][7] << " , "
        <<  Sbdata[i][8] << " , " << Sbdata[i][9] <<  "\n";
    }
    
    ost2.close();
    
    string name3 = "/Users/joekirk/Desktop/Part III Project - Initial Results/" +file.name + "/Nearest_Neighbour_Evolution/TeEvolution" + file.name + ".csv";
    ofstream ost3(name3.c_str());
    if (!ost3) cout << "cannot open file " << name3;
    
    ost3 << "Data extracted from file: " << file.name << "\n";
    ost3 <<  "i ,time / s, Te1, Te2, Te3, Te4, Te5, Te6, Te7, Te8, Te9, Te10 \n";
    
    
    for (int i=0; i < file.frames.size() ; ++i)
        
    {
        ost3 <<  i << " , " << (i*timestep*increment) + timestep << " , "
        <<  Tedata[i][0] << " , " << Tedata[i][1] << " , " << Tedata[i][2] << " , " << Tedata[i][3] << " , "
        <<  Tedata[i][4] << " , " << Tedata[i][5] << " , " << Tedata[i][6] << " , " << Tedata[i][7] << " , "
        <<  Tedata[i][8] << " , " << Tedata[i][9] <<  "\n";
    }
    
    ost3.close();
    
    
}




 
