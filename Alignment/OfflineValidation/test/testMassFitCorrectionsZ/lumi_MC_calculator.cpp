#include <ROOT/RDataFrame.hxx>
#include "TFile.h"

using namespace std;
using namespace ROOT;
using ROOT::RDF::RNode;

int main() {
    ROOT::EnableImplicitMT();

    // Zmumu cross section in fb
    double xsec = 2001.9e+03;  

    //TFile* fout = TFile::Open("lumi_MC_out.root", "RECREATE");
    ROOT::RDataFrame d( "preCutsTree", "/afs/cern.ch/user/c/calexe/CMSSW_14_0_18/src/RecoTracker/TrackProducer/test/globalcor_0_preCutsTree.root" );
    
    auto dlast = std::make_unique<RNode>(d);
    std::cout <<"Total initial entries count is " << *(dlast->Count()) << std::endl;

    // Define MC weight
    dlast = std::make_unique<RNode>(dlast->Define("weight", [](float weight) -> float
	{
	  return std::copysign(1.0, weight);
    }, {"genweight"} ));   

    // Book weight histo
    auto myHist = dlast->Histo1D({"h_genweight", "h_genweight", 2, -1., 1.1}, "weight");
    std::cout << ( myHist->GetBinContent(2) - myHist->GetBinContent(1) ) / xsec ; //TODO confirm with someone the weights are right and xsec is same as run2  xsec_ZmmPostVFP = 2001.9

}