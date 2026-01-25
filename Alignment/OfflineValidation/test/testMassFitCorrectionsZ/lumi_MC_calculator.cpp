#include <ROOT/RDataFrame.hxx>
#include "TFile.h"

using namespace std;
using namespace ROOT;
using ROOT::RDF::RNode;

void lumi_MC_calculator(std::string pathToMCFiles) {
    ROOT::EnableImplicitMT();

    // Zmumu cross section in fb from Run 2
    // TODO confirm we can use same xsec as run2  xsec_ZmmPostVFP = 2001.9
    double xsec = 2001.9e+03;  

    //TFile* fout = TFile::Open("./inoutfiles/results/lumi_MC.root", "RECREATE");
    ROOT::RDataFrame d( "preCutsTree", pathToMCFiles );
    
    auto dlast = std::make_unique<RNode>(d);
    std::cout <<"Total initial entries count is " << *(dlast->Count()) << std::endl;

    // Define MC weight
    dlast = std::make_unique<RNode>(dlast->Define("weight", [](float weight) -> float
	{
	  return std::copysign(1.0, weight);
    }, {"genweight"} ));   

    // Book weight histo
    auto myHist = dlast->Histo1D({"h_genweight", "h_genweight", 2, -1., 1.1}, "weight");
    std::cout << "Equivalent integrated lumi of MC in fb^-1 = number of weighted events / xsec " << ( myHist->GetBinContent(2) - myHist->GetBinContent(1) ) / xsec << std::endl;
    // myHist->Write();

}