//AUTHOR Chris McGinn
#ifndef PYTHIATRACE_H
#define PYTHIATRACE_H

//cpp
#include <fstream>
#include <string>
#include <map>
#include <vector>

//ROOT
#include "TDatime.h"
#include "TMath.h"

//PYTHIA
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"

class pythiaTrace{
 public:
  //functions, external access
  pythiaTrace(Pythia8::Event inPytEvent, Pythia8::Info inPytInfo, int inSeed, int inNEvt);
  ~pythiaTrace(){};

  void processPrintShort();
  void printHistoryFromPos(unsigned long long inVal);
  void printHistoryFullHardScatt();
  void printLineFromVectDepth(std::vector<unsigned long long> inVect, std::vector<unsigned long long> prevLinePos, std::vector<std::string> prevLineName, unsigned long long maxDepth, unsigned long long currDepth, int curr, int total);
  void printHardScatt3(){printHistoryFromPos(hardScatt3); return;}
  void printHardScatt4(){printHistoryFromPos(hardScatt4); return;}

  void print1and2to3and4();
  void recursivePrint(int pos, int depth);
  void recursiveToOrigin(int pos, std::vector<unsigned long long>* inVect);
  
  void printOriginPoints();
  
 protected:
  //Objects, internal
  double ptDeltaZero = 0.1;
  double beamerWidthCM = 12.8;
  double beamerHeightCM = 9.6;

  const static int nNodesX = 28;
  const static int nNodesY = 22;

  std::vector<double> nodesX, nodesY;
  
  Pythia8::Event pytEvent;

  std::ofstream outFile;
  std::string dateStr;  
  
  int seed;
  int nEvt;
  
  double pthat;

  unsigned long long hardScatt1 = 0;
  unsigned long long hardScatt2 = 0;

  unsigned long long hardScatt3 = 0;
  unsigned long long hardScatt4 = 0;

  std::map<unsigned long long, std::vector<unsigned long long> > motherToDaughterMap;
  std::map<unsigned long long, std::vector<unsigned long long> > motherToDaughterMapInv;
  std::map<unsigned long long, std::vector<unsigned long long> > motherToDaughterMapInv2;
  std::map<unsigned long long, std::vector<unsigned long long> > daughterToMotherMap;
  std::vector<unsigned long long> originPoints;
  std::map<unsigned long long, std::vector<unsigned long long> > posToOriginPoints;
  std::map<unsigned long long, std::vector<unsigned long long> > cleanedMotherToDaughterMap;
  
  //Functions, internal
  std::vector<unsigned long long> firstTransverseDaughter(unsigned long long particle, unsigned long long* counter);
  unsigned long long getMaxDepth(unsigned long long particle);
  
  template<typename T>
    std::vector<T> removeDuplicates(std::vector<T> inVect);  
};

#endif
