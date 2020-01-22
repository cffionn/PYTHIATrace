//cpp
#include <algorithm>

//Include header of class
#include "include/pythiaTrace.h"

//Local
#include "include/getLinBins.h"
#include "include/doGlobalDebug.h"
#include "include/pdgToPrettyString.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

pythiaTrace::pythiaTrace(Pythia8::Event inPytEvent, Pythia8::Info inPytInfo, int inSeed = -1, int inNEvt = -1)
{
  //Setup our grid
  Double_t tempX[nNodesX+1];
  Double_t tempY[nNodesY+1];
  getLinBins(-(beamerWidthCM/2.), beamerWidthCM/2., nNodesX, tempX);
  getLinBins(-(beamerHeightCM/2.), beamerHeightCM/2., nNodesY, tempY);

  for(Int_t xI = 0; xI < nNodesX; ++xI){
    nodesX.push_back((tempX[xI] + tempX[xI+1])/2.);
  }
  for(Int_t yI = 0; yI < nNodesY; ++yI){
    nodesY.push_back((tempY[yI] + tempY[yI+1])/2.);
  }

  std::cout << "CHECK NODES X: " << std::endl;
  for(auto const & x : nodesX){
    std::cout << x << ", ";
  }
  std::cout << std::endl;

  std::cout << "CHECK NODES Y: " << std::endl;
  for(auto const & y : nodesY){
    std::cout << y << ", ";
  }
  std::cout << std::endl;

  //Lets document the seed and nEvt
  seed = inSeed;
  nEvt = inNEvt;

  TDatime* date = new TDatime();
  dateStr = std::to_string(date->GetDate());
  delete date;
  
  //Keep an internal copy of this event
  pytEvent = inPytEvent;
  
  //Build out the map from the input pythia event
  for(int pI = 0; pI < pytEvent.size(); ++pI){
    unsigned long long d1 = (unsigned long long)pytEvent[pI].daughter1();
    unsigned long long d2 = (unsigned long long)pytEvent[pI].daughter2();

    unsigned long long m1 = (unsigned long long)pytEvent[pI].mother1();
    unsigned long long m2 = (unsigned long long)pytEvent[pI].mother2();

    std::vector<unsigned long long> dVect;
    if(d1 != 0) dVect.push_back(d1);
    if(d2 != 0 && d2 != d1) dVect.push_back(d2);

    std::vector<unsigned long long> mVect;
    if(m1 != 0) mVect.push_back(m1);
    if(m2 != 0 && m2 != m1) mVect.push_back(m2);

    motherToDaughterMap[pI] = dVect;
    daughterToMotherMap[pI] = mVect;

    if(m1 == 1 || m1 == 2 || m2 == 1 || m2 == 2) originPoints.push_back(pI);
    
    //Dump final particle info for debug
    if(pytEvent[pI].isFinal()){
      //      std::cout << " FINAL pos, id, pt, eta: " << pI << ", " << pytEvent[pI].id() << ", " << TMath::Sqrt(pytEvent[pI].py()*pytEvent[pI].py() + pytEvent[pI].px()*pytEvent[pI].px()) << ", " << pytEvent[pI].eta() << std::endl;      
    }
  }

  //check that 3, 4 has trivial singular origin in proton, and substitute hard scattering as pseudo-origin point
  std::vector<unsigned long long> tempVect, tempVect2;
  recursiveToOrigin(5, &tempVect);
  if(tempVect.size() != 2){
    std::cout << "WARNING: position 5 has more than 2 origin point! return" << std::endl;
    for(unsigned int tI = 0; tI < tempVect.size(); ++tI){
      std::cout << tempVect[tI] << std::endl;
    }
    
    return;
  }
  else std::sort(std::begin(tempVect), std::end(tempVect));    

  recursiveToOrigin(6, &tempVect2);
  if(tempVect2.size() != 2){
    std::cout << "WARNING: position 6 has more than 2 origin point! return" << std::endl;
    for(unsigned int tI = 0; tI < tempVect2.size(); ++tI){
      std::cout << tempVect2[tI] << std::endl;
    }
    return;
  }
  else{
    std::sort(std::begin(tempVect2), std::end(tempVect2));

    for(unsigned int tI = 0; tI < tempVect.size(); ++tI){
      if(tempVect[tI] != tempVect2[tI]){
	std::cout << "WARNING: position 6 has different origin points than 5! return" << std::endl;
	return;
      }
    }
    
    for(unsigned int oI = 0; oI < originPoints.size(); ++oI){
      if(originPoints[oI] == tempVect[0]) originPoints[oI] = 5;
      else if(originPoints[oI] == tempVect[1]) originPoints[oI] = 6;
    }
  }
  
  std::sort(std::begin(originPoints), std::end(originPoints));
  for(auto const & iter : motherToDaughterMap){
    std::vector<unsigned long long> tempLocal;
    recursiveToOrigin(iter.first, &tempLocal);

    posToOriginPoints[iter.first] = tempLocal;
  }

  //build the inverse map
  for(auto const & iter : motherToDaughterMap){
    std::vector<unsigned long long> dVect = iter.second;

    for(unsigned int dI = 0; dI < dVect.size(); ++dI){
      if(motherToDaughterMapInv.count(dVect[dI]) == 0) motherToDaughterMapInv[dVect[dI]] = {};

      motherToDaughterMapInv[dVect[dI]].push_back(iter.first);
    }
  }
  
  //Lets dump every instance of more than 2 parents
  for(auto const & iter : motherToDaughterMapInv){
    if(iter.second.size() > 2){
      std::cout << "Particle at position \'" << iter.first << "\' has more than 2 mothers: " << std::endl;
      for(unsigned int iter2 = 0; iter2 < iter.second.size(); ++iter2){
	std::cout << " " << iter2 << "/" << iter.second.size() << ": " << iter.second[iter2] << std::endl;
      }
    }
  }
  
  //We have a mother map, an inverse mother map, and a daughter map - clean mother based on daughter
  for(auto const & iter : motherToDaughterMapInv){
    if(iter.second.size() > 2){
      std::cout << "CLEANING " << iter.first << std::endl;
      
      for(unsigned int iter2 = 0; iter2 < iter.second.size(); ++iter2){
	bool isGoodMatch = false;

	std::cout << " SEARCHING FOR: " << iter.second[iter2] << std::endl;
	for(unsigned int pI = 0; pI < daughterToMotherMap[iter.first].size(); ++pI){
	  std::cout << "  " << daughterToMotherMap[iter.first][pI] << std::endl;
	  if(iter.second[iter2] == daughterToMotherMap[iter.first][pI]){	    
	    isGoodMatch = true;
	    break;
	  }
	}

	std::cout << "MATCH FOUND: " << isGoodMatch << std::endl;
	if(!isGoodMatch){
	  int daughtPos = -1;
	  
	  for(unsigned int iter3 = 0; iter3 < motherToDaughterMap[iter.second[iter2]].size(); ++iter3){
	    if(motherToDaughterMap[iter.second[iter2]][iter3] == iter.first){
	      daughtPos = iter3;
	      break;
	    }
	  }

	  if(daughtPos < 0){
	    std::cout << "Despite redundancy, cant find daughtPos" << std::endl;
	  }
	  else{
	    std::cout << " Removing daughter " << motherToDaughterMap[iter.second[iter2]][daughtPos] << " from mother " << iter.second[iter2] << std::endl;
	    motherToDaughterMap[iter.second[iter2]].erase(motherToDaughterMap[iter.second[iter2]].begin()+daughtPos);
	    if(cleanedMotherToDaughterMap.count(iter.second[iter2]) == 0) cleanedMotherToDaughterMap[iter.second[iter2]] = {};

	    cleanedMotherToDaughterMap[iter.second[iter2]].push_back(iter.first);
	  }
	}
      }      
    }
  }

  //build the inverse map again
  for(auto const & iter : motherToDaughterMap){
    std::vector<unsigned long long> dVect = iter.second;

    for(unsigned int dI = 0; dI < dVect.size(); ++dI){
      if(motherToDaughterMapInv2.count(dVect[dI]) == 0) motherToDaughterMapInv2[dVect[dI]] = {};

      motherToDaughterMapInv2[dVect[dI]].push_back(iter.first);
    }
  }

  //Lets dump every instance of more than 2 parents again
  for(auto const & iter : motherToDaughterMapInv2){
    if(iter.second.size() > 2){
      std::cout << "Particle at position \'" << iter.first << "\' has more than 2 mothers (AGAIN): " << std::endl;
      for(unsigned int iter2 = 0; iter2 < iter.second.size(); ++iter2){
	std::cout << " " << iter2 << "/" << iter.second.size() << ": " << iter.second[iter2] << std::endl;
      }
    }
  }

  
  /*
  for(auto const & m : motherToDaughterMap){
    for(auto const & m2 : m.second){
      std::cout << " " << m2 << std::endl;
    }
  }  
  */

  for(int pI = 0; pI < pytEvent.size(); ++pI){
    if(!pytEvent[pI].isFinal()) continue;
    std::cout << " FINAL pos, id, pt, eta: " << pI << ", " << pytEvent[pI].id() << ", " << TMath::Sqrt(pytEvent[pI].py()*pytEvent[pI].py() + pytEvent[pI].px()*pytEvent[pI].px()) << ", " << pytEvent[pI].eta() << std::endl;

    for(unsigned int pI2 = 0; pI2 < posToOriginPoints[pI].size(); ++pI2){
      std::cout << "  Origin: " << posToOriginPoints[pI][pI2] << std::endl;
    }
  }
  
  //Lets fill out the hardscattering
  pthat = inPytInfo.pTHat();

  std::vector<unsigned long long> tempDVect = motherToDaughterMap[1];
  if(tempDVect.size() != 1) std::cout << "BIG WARNING: MORE THAN 1 BEAM DAUGHTER" << std::endl;
  std::vector<unsigned long long> hardScattV1 = firstTransverseDaughter(tempDVect[0], NULL);

  tempDVect = motherToDaughterMap[2];
  if(tempDVect.size() != 1) std::cout << "BIG WARNING: MORE THAN 1 BEAM DAUGHTER" << std::endl;
  std::vector<unsigned long long> hardScattV2 = firstTransverseDaughter(tempDVect[0], NULL);

  hardScattV1.insert(hardScattV1.end(), hardScattV2.begin(), hardScattV2.end());  
  hardScattV1 = removeDuplicates(hardScattV1);

  if(hardScattV1.size() != 2) std::cout << "BIG WARNING: MORE THAN 2 HARD SCATTERED IN FINAL" << std::endl;
  hardScatt3 = hardScattV1[0];
  hardScatt4 = hardScattV1[1];

  hardScatt1 = pytEvent[hardScatt3].mother1();
  hardScatt2 = pytEvent[hardScatt3].mother2();

  processPrintShort();
  printHistoryFullHardScatt();
  
  return;
}

void pythiaTrace::processPrintShort()
{
  std::cout << "This event is hard scattering of " << pdgToPrettyString(pytEvent[hardScatt1].id()) << " + " << pdgToPrettyString(pytEvent[hardScatt2].id()) << " --> " << pdgToPrettyString(pytEvent[hardScatt3].id()) << " + " << pdgToPrettyString(pytEvent[hardScatt4].id()) << ", pTHat = " << prettyString(pthat, 2, false) << "." << std::endl;
  
  return;
}

void pythiaTrace::printHistoryFromPos(unsigned long long inVal)
{
  unsigned long long maxDepth = getMaxDepth(inVal)+1;//Misses the very first instance - should be 
  std::cout << " Depth of " << maxDepth << std::endl;

  //Lets create a unique output file name 
  std::string outFileName = "particle" + std::to_string(inVal) + "_Seed";
  if(seed < 0) outFileName = outFileName + "NA_NEvt";
  else outFileName = outFileName + std::to_string(seed) + "_NEvt";
  if(nEvt < 0) outFileName = outFileName + "NA_";
  else outFileName = outFileName + std::to_string(nEvt) + "_";
  outFileName = outFileName + dateStr + ".tex";

  outFile.open(outFileName.c_str());

  outFile << "\\begin{frame}" << std::endl;
  outFile << "\\begin{tikzpicture}[remember picture, overlay]" << std::endl;
  outFile << "\\fontsize{4}{4}\\selectfont" << std::endl;
  outFile << std::endl;
  
  std::vector<unsigned long long> processVect = {inVal};
  printLineFromVectDepth(processVect, {}, {}, maxDepth, 0, 1, 0);

  outFile << "\\end{tikzpicture}" << std::endl;
  outFile << "\\end{frame}" << std::endl;
  outFile.close();
  
  return;
}


void pythiaTrace::printHistoryFullHardScatt()
{
  unsigned long long maxDepth = TMath::Max(getMaxDepth(5), getMaxDepth(6))+1;//Misses the very first instance - should be 
  std::cout << " Depth of " << maxDepth << std::endl;

  //Lets create a unique output file name 
  std::string outFileName = "hardScatt_Seed";
  if(seed < 0) outFileName = outFileName + "NA_NEvt";
  else outFileName = outFileName + std::to_string(seed) + "_NEvt";
  if(nEvt < 0) outFileName = outFileName + "NA_";
  else outFileName = outFileName + std::to_string(nEvt) + "_";
  outFileName = outFileName + dateStr + ".tex";

  outFile.open(outFileName.c_str());

  outFile << "\\begin{frame}" << std::endl;
  outFile << "\\begin{tikzpicture}[remember picture, overlay]" << std::endl;
  outFile << "\\fontsize{4}{4}\\selectfont" << std::endl;
  outFile << std::endl;
  
  std::vector<unsigned long long> processVect = {5};
  std::vector<unsigned long long> processVect2 = {6};
  printLineFromVectDepth(processVect, {}, {}, maxDepth, 0, 0, 2);
  printLineFromVectDepth(processVect2, {}, {}, maxDepth, 0, 1, 2);

  for(auto const & iter : cleanedMotherToDaughterMap){
    bool goodOrigin = false;
    for(unsigned int iter2 = 0; iter2 < posToOriginPoints[iter.first].size(); ++iter2){
      if(posToOriginPoints[iter.first][iter2] == 5 || posToOriginPoints[iter.first][iter2] == 6){
	goodOrigin = true;
	break;
      }
    }

    if(!goodOrigin) continue;

    for(unsigned int iter2 = 0; iter2 < iter.second.size(); ++iter2){      

      goodOrigin = false;
      for(unsigned int iter3 = 0; iter3 < posToOriginPoints[iter.second[iter2]].size(); ++iter3){
	if(posToOriginPoints[iter.second[iter2]][iter3] == 5 || posToOriginPoints[iter.second[iter2]][iter3] == 6){
	  goodOrigin = true;
	  break;
	}
      }

      
      if(goodOrigin) outFile << "\\path [line] (Part" << iter.first << ") -- (Part" << iter.second[iter2] << ");" << std::endl;
    }
  }
    
  outFile << "\\end{tikzpicture}" << std::endl;
  outFile << "\\end{frame}" << std::endl;
  outFile.close();
  
  return;
}

void pythiaTrace::printLineFromVectDepth(std::vector<unsigned long long> inVect, std::vector<unsigned long long> prevLinePos, std::vector<std::string> prevLineName, unsigned long long maxDepth, unsigned long long currDepth, int curr = 0, int total = 1)
{
  if(total > 2){
    std::cout << "PYTHIATRACE::printLineFromVectDepth ERROR: YOU HAVENT TAUGHT ME HOW TO DO THREE PARTONS. return" << std::endl;
    return;
  }
  
  std::vector<int> takenXNodes;
  Double_t yVal = nodesY[nNodesY-1-currDepth];
  std::string yValStr = "0.0";
  if(TMath::Abs(yVal) > 0.000001) yValStr = prettyString(yVal, 5, false);
  
  std::vector<unsigned long long> prevPos, prevPos2;
  std::vector<std::string> prevName, prevName2;

  double centerOfMass = 0;
  if(prevLinePos.size() > 0){
    for(auto const & val : prevLinePos){
      centerOfMass += val;
    }
    centerOfMass /= (double)prevLinePos.size();
    centerOfMass -= nNodesX/2;
  }
  
  std::cout << "COM: " << centerOfMass << std::endl;
  
  for(unsigned int pI = 0; pI < inVect.size(); ++pI){
    std::string nameStr = "Part" + std::to_string(inVect[pI]);
    if(currDepth == 0){
      int nodeToTake = (nNodesX/2)+1;
      if(total == 2){
	if(curr == 0) nodeToTake = 0;
	if(curr == 1) nodeToTake = nNodesX-1;
      }
      double xVal = nodesX[nodeToTake];
      takenXNodes.push_back(nodeToTake);

      std::string xValStr = "0.0";
      if(TMath::Abs(xVal) > 0.0000001) xValStr = prettyString(xVal, 5, false);

      std::string partStr = pdgToPrettyStringTex(pytEvent[inVect[pI]].id());
      if(partStr.find("\\") != std::string::npos) partStr = "$" + partStr + "$";
      else if(partStr.find("^") != std::string::npos) partStr = "$" + partStr + "$";
      else if(partStr.find("_") != std::string::npos) partStr = "$" + partStr + "$";

      outFile << "\\node[yshift=" << yValStr << "cm,xshift=" << xValStr << "] at (current page.center) [block] (" << nameStr << ") {" << partStr << ", "<< inVect[pI] <<"};" << std::endl;
    }
    else{
      int nodeToTake = 1000;
      int adj = 0;
      while(true){
	int nodeToTake1 = prevLinePos[pI];
	int nodeToTake2 = prevLinePos[pI];

	if(centerOfMass > 0){
	  nodeToTake1 -= adj;
	  nodeToTake2 += adj;	  	 
	}
	else{
	  nodeToTake1 += adj;
	  nodeToTake2 -= adj;	  	 
	}
	

	if(nodeToTake1 < (int)nodesX.size() && nodeToTake1 >= 0){
	  bool goodNode = true;
	  for(unsigned int tI = 0; tI < takenXNodes.size(); ++tI){
	    if(nodeToTake1 == takenXNodes[tI]){
	      goodNode = false;
	      break;
	    }
	  }
	  if(goodNode){
	    nodeToTake = nodeToTake1;
	    break;
	  }
	}

	if(nodeToTake2 >= 0 && nodeToTake2 < (int)nodesX.size()){
	  bool goodNode = true;
	  for(unsigned int tI = 0; tI < takenXNodes.size(); ++tI){
	    if(nodeToTake2 == takenXNodes[tI]){
	      goodNode = false;
	      break;
	    }
	  }
	  if(goodNode){
	    nodeToTake = nodeToTake2;
	    break;
	  }
	}
	if((nodeToTake1 < 0 || nodeToTake1 >= (int)nodesX.size()) && (nodeToTake2 < 0 || nodeToTake2 >= (int)nodesX.size())){
	  std::cout << "MORE NODES THAN SPACE SHIT" << std::endl;
	}
	
	++adj;
      }
      
      double xVal = nodesX[nodeToTake];
      takenXNodes.push_back(nodeToTake);

      std::string xValStr = "0.0";
      if(TMath::Abs(xVal) > 0.0000001) xValStr = prettyString(xVal, 5, false);

      std::string partStr = pdgToPrettyStringTex(pytEvent[inVect[pI]].id());
      if(partStr.find("\\") != std::string::npos) partStr = "$" + partStr + "$";
      else if(partStr.find("^") != std::string::npos) partStr = "$" + partStr + "$";
      else if(partStr.find("_") != std::string::npos) partStr = "$" + partStr + "$";

      outFile << "\\node[yshift=" << yValStr << "cm,xshift=" << xValStr << "cm] at (current page.center) [block] (" << nameStr << ") {"<< partStr << ", " << inVect[pI] << "};" << std::endl;

      outFile << "\\path [line] (" << prevLineName[pI] << ") -- (" << nameStr << ");" << std::endl;      
    }

  
    outFile << std::endl;
    prevName.push_back(nameStr);
    prevPos.push_back(takenXNodes[takenXNodes.size()-1]);    
  }

  std::vector<unsigned long long> nextVect;

  for(unsigned int pI = 0; pI < inVect.size(); ++pI){
    std::vector<unsigned long long> daught = motherToDaughterMap[inVect[pI]];

    if(inVect[pI] == 908){
      std::cout << "WE HIT 908 at depth " << currDepth << ": " << pytEvent[inVect[pI]].id() << std::endl;
    }
    
    for(auto const & dI : daught){         
      if(dI == 908){
	std::cout << "WE HIT 908 DAUGHTER AT " << currDepth << ", PARENT " << inVect[pI] << ", " << pytEvent[inVect[pI]].id() << std::endl;
      }

      nextVect.push_back(dI);
      prevPos2.push_back(prevPos[pI]);
      prevName2.push_back(prevName[pI]);
    }
  }
  
  if(nextVect.size() > 0) printLineFromVectDepth(nextVect, prevPos2, prevName2, maxDepth, currDepth + 1);
  
  return;
}


//Functions, internal
std::vector<unsigned long long> pythiaTrace::firstTransverseDaughter(unsigned long long particle, unsigned long long* counter)
{
  if(counter == NULL) counter = new unsigned long long(0);
  Pythia8::Particle input = pytEvent[particle];

  std::vector<unsigned long long> particles;
  
  //Check that this particle isn't already transverse
  if(TMath::Abs(input.pT() - pthat) < ptDeltaZero) particles.push_back(particle);
  else{
    ++(*counter);
    
    //Check for transverse daughters
    std::vector<unsigned long long> tempDVect = motherToDaughterMap[particle];  
    
    //Set our minCounter arbitrary large, minDaughter ID of 0
    unsigned long long minCounter = 100000;
    
    for(auto const & daught : tempDVect){
      unsigned long long tempCounter = 0;
      std::vector<unsigned long long> tempDaughter = firstTransverseDaughter(daught, &tempCounter);
      
      if(tempCounter < minCounter){
	particles = tempDaughter;
	minCounter = tempCounter;
      }
      else if(tempCounter == minCounter){
	particles.insert(particles.end(), tempDaughter.begin(), tempDaughter.end());
      }
    }

    (*counter) += minCounter;
  }
  
  return particles;
}

unsigned long long pythiaTrace::getMaxDepth(unsigned long long particle)
{
  std::vector<unsigned long long> tempVect = motherToDaughterMap[particle];
  unsigned long long counter = 0;
  for(auto const & daught : tempVect){
    unsigned long long depth = getMaxDepth(daught);
    if(depth > counter) counter = depth;
  }    
  if(tempVect.size() > 0) ++counter;
  return counter;
}


template<typename T>
std::vector<T> pythiaTrace::removeDuplicates(std::vector<T> inVect)
{
  unsigned int pos = 0;
  while(pos < inVect.size()){
    unsigned int pos2 = pos+1;
    while(pos2 < inVect.size()){
      if(inVect[pos] == inVect[pos2]) inVect.erase(inVect.begin()+pos2);
      else ++pos2;
    }

    ++pos;
  }
  
  return inVect;
}

void pythiaTrace::print1and2to3and4()
{
  std::cout << "TRACE 3: " << std::endl;
  recursivePrint(3, 0);

  std::cout << "TRACE 4: " << std::endl;
  recursivePrint(4, 0);

  return;
}


void pythiaTrace::recursivePrint(int pos, int depth)
{
  std::string depthStr = "";
  for(int dI = 0; dI < depth; ++dI){
    depthStr = depthStr + " ";
  }
  
  std::vector<unsigned long long> temp = daughterToMotherMap[pos];
  for(unsigned int tI = 0; tI < temp.size(); ++tI){
    std::cout << depthStr << "Mother of pos \'" << pos << "\' is " << temp[tI] << "\'." << std::endl;
    if(temp[tI] == 1 || temp[tI] == 2) break;
    else recursivePrint(temp[tI], depth+1);
  }
 
  return;
}


void pythiaTrace::printOriginPoints()
{
  std::cout << "Printing origin points: " << std::endl;
  for(unsigned int pI = 0; pI < originPoints.size(); ++pI){
    std::cout << " " << pI << "/" << originPoints.size() << ": " << originPoints[pI] << std::endl;
  }
  return;
}

void pythiaTrace::recursiveToOrigin(int pos, std::vector<unsigned long long>* inVect)
{
  std::cout << "RECURSIVE CALL ON POS " << pos << std::endl;
  
  for(unsigned int mI = 0; mI < daughterToMotherMap[pos].size(); ++mI){
    unsigned long long m1 = daughterToMotherMap[pos][mI];
    std::cout << " CHECKING " << m1 << std::endl;

    bool isOriginPoint = false;
    for(unsigned int oI = 0; oI < originPoints.size(); ++oI){
      if(originPoints[oI] == m1){isOriginPoint = true; break;}
    }    

    if(isOriginPoint){
      bool isInVect = false;
      for(unsigned int pI = 0; pI < inVect->size(); ++pI){
	if(inVect->at(pI) == m1){
	  isInVect = true;
	  break;
	}
      }
      if(!isInVect) inVect->push_back(m1);
    }
    else recursiveToOrigin(m1, inVect);
  }
  
  return;
}
