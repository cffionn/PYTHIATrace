//cpp
#include <iostream>

//PYTHIA
#include "Pythia8/Pythia.h"

//Local
#include "include/doGlobalDebug.h"
#include "include/pythiaTrace.h"

int main()
{
  if(doGlobalDebug) std::cout << "FILE, LiNE: " <<__FILE__ << ", " << __LINE__ << std::endl;

  Pythia8::Pythia pythia;
  pythia.readString("Beams:eCM = 5020.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 1563583");
  //  pythia.readString("Random:seed = 0");

  if(doGlobalDebug) std::cout << "FILE, LiNE: " <<__FILE__ << ", " << __LINE__ << std::endl;

  pythia.readString("PhaseSpace:pTHatMin = 15.");
  
  pythia.init();

  int totEvt = 0;
  const int nTot = 1;
  int globalCounter = 0;

  if(doGlobalDebug) std::cout << "FILE, LiNE: " <<__FILE__ << ", " << __LINE__ << std::endl;

  while(totEvt < nTot){
    ++globalCounter;
    if(!pythia.next()) continue;

    if(doGlobalDebug) std::cout << "FILE, LiNE: " <<__FILE__ << ", " << __LINE__ << std::endl;
    pythiaTrace binaryTree(pythia.event, pythia.info, pythia.getSeed(), globalCounter);
    if(doGlobalDebug) std::cout << "FILE, LiNE: " <<__FILE__ << ", " << __LINE__ << std::endl;
    binaryTree.print1and2to3and4();
    binaryTree.printOriginPoints();
    
    ++totEvt;
  }
  
  return 0;
}
