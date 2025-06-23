// Takes a TTree with one filename per track and reshapes it to a TTree with one filename per event
// Needed before running massscales_data on the outputs of the CVH plugin
// Author: Cristina-Andreea Alexe

#include "TFile.h"
#include "TString.h"
#include "TVector.h"
#include "TVectorT.h"

#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TROOT.h"

#include <iostream>
#include <filesystem>
#include <chrono>

using namespace std;
using namespace ROOT;
namespace fs = std::filesystem;

typedef ROOT::VecOps::RVec<double> RVecD;
typedef ROOT::VecOps::RVec<unsigned int> RVecUI;
typedef ROOT::VecOps::RVec<int> RVecI;
typedef ROOT::VecOps::RVec<float> RVecF;
typedef ROOT::VecOps::RVec<bool> RVecB;
using ROOT::RDF::RNode;



int main() {
    ROOT::EnableImplicitMT();
    auto start = std::chrono::high_resolution_clock::now();

    // Read input and output directory names
    const std::string input_directory = "/gpfs/ddn/srm/cms/store/user/calexe/Muon/Run2022F-22Sep2023-v2-with-CVH/250606_161622/";
    // "/gpfs/ddn/srm/cms/store/user/calexe/DYto2Mu_MLL-50to120_keepMDSHits_TuneCP5_13p6TeV_powheg-pythia8/CVH_refit_MC/250507_094615/";
    const std::string output_directory = "/home/users/alexe/workingarea/CMSSW_15_0_0_pre1/src/Alignment/OfflineValidation/test/testMassFitCorrectionsZ/inoutfiles/Run2022F-22Sep2023-v2-with-CVH/";
    // "/home/users/alexe/workingarea/CMSSW_15_0_0_pre1/src/Alignment/OfflineValidation/test/testMassFitCorrectionsZ/inoutfiles/DYto2Mu_MLL-50to120_keepMDSHits_TuneCP5_13p6TeV_powheg-pythia8_with_CVH/";
    const fs::path output_dir_path{output_directory};

    // Create output directory if needed
    if (!fs::exists(output_dir_path)) {
        if (fs::create_directory(output_dir_path)) {
            std::cout << "Directory created: " << output_directory << std::endl;
        } else {
            std::cerr << "Failed to create directory: " << output_directory << std::endl;
        }
    }
    
    // For each file in the input directory
    for (const auto& entry : fs::recursive_directory_iterator(input_directory)) {
        if (entry.is_regular_file() && entry.path().extension() == ".root") {

            // Open input file
            TFile* f_in = TFile::Open(entry.path().string().c_str(), "READ");
            // Read input TTree
            auto tree_in = f_in->Get<TTree>("tree");
            
            ULong64_t event_in;
            int muon_trigger_in;
            float pt_in, eta_in, phi_in, charge_in, gen_pt_in, gen_eta_in, gen_phi_in, gen_weight_in;
            bool muon_loose_in, muon_is_global_in, track_purity_in, muon_medium_in;

            tree_in->SetBranchAddress("event", &event_in);
            tree_in->SetBranchAddress("UpdPt", &pt_in);
            tree_in->SetBranchAddress("UpdEta", &eta_in);
            tree_in->SetBranchAddress("UpdPhi", &phi_in);
            tree_in->SetBranchAddress("genPt", &gen_pt_in);
            tree_in->SetBranchAddress("genEta", &gen_eta_in);
            tree_in->SetBranchAddress("genPhi", &gen_phi_in);
            tree_in->SetBranchAddress("genweight", &gen_weight_in);
            tree_in->SetBranchAddress("trackCharge", &charge_in);
            tree_in->SetBranchAddress("muonLoose", &muon_loose_in);
            tree_in->SetBranchAddress("muonIsGlobal", &muon_is_global_in);
            tree_in->SetBranchAddress("trackHighPurity", &track_purity_in);
            tree_in->SetBranchAddress("muonMedium", &muon_medium_in);
            tree_in->SetBranchAddress("HLT_IsoMu24", &muon_trigger_in);

            // Open output file
            std::string filename = entry.path().filename().string();
            filename.replace(filename.size() - 5, 5, "_reshaped.root");
            TFile* f_out = TFile::Open((output_directory+filename).c_str(), "RECREATE");

            // Define output TTree
            auto tree_out = std::make_unique<TTree>("tree", "tree_reshaped");
            const int basketSize = 4*1024*1024;

            // define RVecs
            ULong64_t event_out;
            RVecI muon_trigger_out;
            RVecF pt_out, eta_out, phi_out, charge_out, gen_pt_out, gen_eta_out, gen_phi_out;
            float gen_weight_out;
            RVecB muon_loose_out, muon_is_global_out, track_purity_out, muon_medium_out;

            // define output tree and branches
            tree_out->Branch("event", &event_out, basketSize);
            tree_out->Branch("UpdPt", &pt_out, basketSize);
            tree_out->Branch("UpdEta", &eta_out, basketSize);
            tree_out->Branch("UpdPhi", &phi_out, basketSize);
            tree_out->Branch("genPt", &gen_pt_out, basketSize);
            tree_out->Branch("genEta", &gen_eta_out, basketSize);
            tree_out->Branch("genPhi", &gen_phi_out, basketSize);  
            tree_out->Branch("genweight", &gen_weight_out, basketSize);
            tree_out->Branch("trackCharge", &charge_out, basketSize);
            tree_out->Branch("muonLoose",&muon_loose_out, basketSize);
            tree_out->Branch("muonIsGlobal",&muon_is_global_out, basketSize);
            tree_out->Branch("trackHighPurity",&track_purity_out, basketSize);
            tree_out->Branch("muonMedium",&muon_medium_out, basketSize);
            tree_out->Branch("HLT_IsoMu24", &muon_trigger_out, basketSize);

            // Reshape tree
            tree_in->GetEntry(0);
            event_out = event_in;
            gen_weight_out = gen_weight_in; 

            int nEntries = tree_in->GetEntries();

            for (int iEntry = 0; iEntry < nEntries; ++iEntry) {
                tree_in->GetEntry(iEntry);
                if (event_in != event_out) { // if you are reading a new event
                    tree_out->Fill(); // fill tree out with the last full event
                    // resize vectors
                    pt_out.resize(0);
                    eta_out.resize(0);
                    phi_out.resize(0);
                    gen_pt_out.resize(0);
                    gen_eta_out.resize(0);
                    gen_phi_out.resize(0);
                    charge_out.resize(0);
                    muon_loose_out.resize(0);
                    muon_is_global_out.resize(0);
                    track_purity_out.resize(0);
                    muon_medium_out.resize(0);
                    muon_trigger_out.resize(0);

                    // read the 1st track of the new event
                    pt_out.emplace_back(pt_in);
                    eta_out.emplace_back(eta_in);
                    phi_out.emplace_back(phi_in);
                    gen_pt_out.emplace_back(gen_pt_in);
                    gen_eta_out.emplace_back(gen_eta_in);
                    gen_phi_out.emplace_back(gen_phi_in);
                    charge_out.emplace_back(charge_in);           
                    muon_loose_out.emplace_back(muon_loose_in);
                    muon_is_global_out.emplace_back(muon_is_global_in);
                    track_purity_out.emplace_back(track_purity_in);
                    muon_medium_out.emplace_back(muon_medium_in);
                    muon_trigger_out.emplace_back(muon_trigger_in);

                    gen_weight_out = gen_weight_in;
                    event_out = event_in;
                } else { // continue reading the current event
                    pt_out.emplace_back(pt_in);
                    eta_out.emplace_back(eta_in);
                    phi_out.emplace_back(phi_in);
                    gen_pt_out.emplace_back(gen_pt_in);
                    gen_eta_out.emplace_back(gen_eta_in);
                    gen_phi_out.emplace_back(gen_phi_in);
                    charge_out.emplace_back(charge_in);           
                    muon_loose_out.emplace_back(muon_loose_in);
                    muon_is_global_out.emplace_back(muon_is_global_in);
                    track_purity_out.emplace_back(track_purity_in);  
                    muon_medium_out.emplace_back(muon_medium_in);
                    muon_trigger_out.emplace_back(muon_trigger_in);
                }
            }

            tree_out->Fill(); // Fill the last event after the loop

            // Write and close files
            tree_out->Write();
            f_out->Write();
            f_in->Close();
        break; // do only 1 file for debugging
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    return 0;
}
