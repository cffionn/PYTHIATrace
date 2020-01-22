#ifndef PDGTOPRETTYSTRING_H
#define PDGTOPRETTYSTRING_H

#include <string>

std::string pdgToPrettyString(int inVal)
{
  std::string outStr = std::to_string(inVal);

  if(inVal == 1) outStr = "d";
  else if(inVal == -1) outStr = "dbar";
  else if(inVal == 2) outStr = "u";
  else if(inVal == -2) outStr = "ubar";
  else if(inVal == 3) outStr = "s";
  else if(inVal == -3) outStr = "sbar";
  else if(inVal == 4) outStr = "c";
  else if(inVal == -4) outStr = "cbar";
  else if(inVal == 5) outStr = "b";
  else if(inVal == -5) outStr = "bbar";
  else if(inVal == 6) outStr = "t";
  else if(inVal == -6) outStr = "tbar";
  else if(inVal == 22) outStr = "gamma";
  else if(inVal == 23) outStr = "Z";
  else if(inVal == 24) outStr = "W^{+}";
  else if(inVal == -24) outStr = "W^{-}";
  else if(inVal == 11) outStr = "e^{-}";
  else if(inVal == -11) outStr = "e^{+}";
  else if(inVal == 12) outStr = "nu_{e}";
  else if(inVal == -12) outStr = "nu_{e}bar";
  else if(inVal == 13) outStr = "mu^{-}";
  else if(inVal == -13) outStr = "mu^{+}";
  else if(inVal == 14) outStr = "nu_{mu}";
  else if(inVal == -14) outStr = "nu_{mu}bar";
  else if(inVal == 15) outStr = "tau^{-}";
  else if(inVal == -15) outStr = "tau^{+}";
  else if(inVal == 16) outStr = "nu_{tau}";
  else if(inVal == -16) outStr = "nu_{tau}bar";
  else if(inVal == 21) outStr = "g";
  else if(inVal == 22) outStr = "gamma";
  else if(inVal == 23) outStr = "Z";
  else if(inVal == 24) outStr = "W^{+}";
  else if(inVal == -24) outStr = "W^{-}";
  else if(inVal == 111) outStr = "pi^{0}";
  else if(inVal == -111) outStr = "pi^{0}bar";
  else if(inVal == 113) outStr = "rho^{0}";
  else if(inVal == 130) outStr = "K^{0}_{L}";
  else if(inVal == 211) outStr = "pi^{+}";
  else if(inVal == -211) outStr = "pi^{-}";
  else if(inVal == 213) outStr = "rho^{+}";
  else if(inVal == -213) outStr = "rho^{-}";
  else if(inVal == 221) outStr = "eta";
  else if(inVal == 223) outStr = "omega";
  else if(inVal == -223) outStr = "omegabar";
  else if(inVal == 310) outStr = "K^{0}_{S}";
  else if(inVal == 311) outStr = "K^{0}";
  else if(inVal == -311) outStr = "K^{0}bar";
  else if(inVal == 313) outStr = "K^{*0}";
  else if(inVal == 321) outStr = "K^{+}";
  else if(inVal == -321) outStr = "K^{-}";
  else if(inVal == 323) outStr = "K^{*+}";
  else if(inVal == -323) outStr = "K^{*-}";
  else if(inVal == 411) outStr = "D^{+}";
  else if(inVal == -411) outStr = "D^{-}";
  else if(inVal == 421) outStr = "D^{0}";
  else if(inVal == 423) outStr = "D^{*0}";
  else if(inVal == 511) outStr = "B^{0}";
  else if(inVal == 521) outStr = "B^{+}";
  else if(inVal == -521) outStr = "B^{-}";
  else if(inVal == 2112) outStr = "n";
  else if(inVal == -2112) outStr = "nbar";
  else if(inVal == 2212) outStr = "p";
  else if(inVal == -2212) outStr = "pbar";
  else if(inVal == 2114) outStr = "Delta^{0}";
  else if(inVal == 2214) outStr = "Delta^{+}";
  else if(inVal == -2214) outStr = "Delta^{+}bar";
  else if(inVal == 3222) outStr = "Sigma^{+}";
  else if(inVal == -3222) outStr = "Sigma^{+}bar";

  return outStr;
}

std::string pdgToPrettyStringTex(int inVal)
{
  std::string outStr = pdgToPrettyString(inVal);
  if(outStr.find("gamma") != std::string::npos) outStr.replace(outStr.find("gamma"), 5, "\\gamma");
  else if(outStr.find("pi") != std::string::npos) outStr = "\\" + outStr;
  else if(outStr.find("omega") != std::string::npos) outStr = "\\" + outStr;
  else if(outStr.find("eta") != std::string::npos) outStr = "\\" + outStr;
  else if(outStr.find("rho") != std::string::npos) outStr = "\\" + outStr;
  else if(outStr.find("Delta") != std::string::npos) outStr = "\\" + outStr;
  else if(outStr.find("Sigma") != std::string::npos) outStr = "\\" + outStr;

  if(outStr.find("bar") != std::string::npos) outStr = "\\overline{" + outStr.substr(0, outStr.find("bar")) + "}";  
  
  return outStr;
}

#endif
