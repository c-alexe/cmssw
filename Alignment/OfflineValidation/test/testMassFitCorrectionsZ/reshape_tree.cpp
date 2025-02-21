// Takes a TTree with one entry per track and reshapes it to a TTree with one entry per event
// Needed before running massscales_data on the outputs of the CVH plugin
// Author: Cristina-Andreea Alexe

#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TVector.h"
#include "TVectorT.h"

#include <iostream>

using namespace std;
using namespace ROOT;

typedef ROOT::VecOps::RVec<double> RVecD;
typedef ROOT::VecOps::RVec<unsigned int> RVecUI;
typedef ROOT::VecOps::RVec<int> RVecI;
typedef ROOT::VecOps::RVec<float> RVecF;
typedef ROOT::VecOps::RVec<bool> RVecB;
using ROOT::RDF::RNode;

int reshape_tree() {

    auto tree_out = std::make_unique<TTree>("tree", "tree_reshaped");
    const int basketSize = 4*1024*1024;

   // define RVecs
   ULong64_t event_out;
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

   // read initial tree
   TFile* fin = TFile::Open("/afs/cern.ch/user/c/calexe/CMSSW_14_0_18/src/RecoTracker/TrackProducer/test/globalcor_0_data.root", "READ");
   auto tree_in = fin->Get<TTree>("tree");

   ULong64_t event_in;
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

   tree_in->GetEntry(0);
   event_out = event_in;
   gen_weight_out = gen_weight_in; 

   for (int iEntry = 0; tree_in->LoadTree(iEntry) >= 0; ++iEntry) {

    tree_in->GetEntry(iEntry);
        if (event_in != event_out) {
            tree_out->Fill();

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

            gen_weight_out = gen_weight_in;
            event_out = event_in;
        }
        else { 
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
        }

   }

   TFile* fout = TFile::Open("./globalcor_0_data_reshaped.root", "RECREATE");
   tree_out->Write();

   return 0;

}