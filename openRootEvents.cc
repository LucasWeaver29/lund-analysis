
#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"


// ROOT - files
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <iostream>

// ROOT - histogram
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TGraph.h"

using namespace Pythia8;

int main() {

    TString root_file_name = "events5000x9.root";
    //TString root_file_name = "Lund, bg, area sub - 45000.root";

    //TString output_name = "jet energy plots 45k, match, resized";
    //TString output_name = "Lund, bg, area sub - 45000, resized";

    // Open the event ROOT file
    TFile rootFile(root_file_name, "READ");
    
    if (rootFile.IsZombie()) {
        std::cout << "Error: Could not open ROOT file" << std::endl;
        return 1;
    }

    
    TTree* eventTree = (TTree*)rootFile.Get("events");
    if (!eventTree) {
        std::cout << "Error: Could not find tree in ROOT file" << std::endl;
        return 1;
    }

    // Define variables to store event data
    double ptHatMin, ptHatMax, crossSection, weight;
    int numParticles;
    std::vector<double>* px = nullptr;
    std::vector<double>* py = nullptr;
    std::vector<double>* pz = nullptr;
    std::vector<double>* energy = nullptr;
    std::vector<double>* eta = nullptr;
    std::vector<int>* pdgId = nullptr;

    eventTree->SetBranchAddress("ptHatMin", &ptHatMin);
    eventTree->SetBranchAddress("ptHatMax", &ptHatMax);
    eventTree->SetBranchAddress("crossSection", &crossSection);
    eventTree->SetBranchAddress("weight", &weight);
    eventTree->SetBranchAddress("numParticles", &numParticles);
    eventTree->SetBranchAddress("px", &px);
    eventTree->SetBranchAddress("py", &py);
    eventTree->SetBranchAddress("pz", &pz);
    eventTree->SetBranchAddress("energy", &energy);
    eventTree->SetBranchAddress("eta", &eta);
    eventTree->SetBranchAddress("pdgId", &pdgId);

    
    for (int iEvent = 0; iEvent < 20; ++iEvent) {

        eventTree->GetEntry(iEvent);

        vector<fastjet::PseudoJet> event_particles;

        cout << numParticles <<" particles in this event";

        for (int iPart = 0; iPart < numParticles; ++iPart) {
            if ((*eta)[iPart] > .9) continue;
            event_particles.push_back(fastjet::PseudoJet((*px)[iPart], (*py)[iPart], (*pz)[iPart], (*energy)[iPart]));
        }

        cout << event_particles.size() << " particles in |eta| < .9" << endl;

    }


    /*
    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);


    c1->cd();
    c1->SetLogz(1);
    hist->Draw("COLZ");
    c1->Print(output_name+".pdf(","pdf");
    c1->Clear();

    
    c1->cd();
    jes->Draw();
    c1->Print(output_name+".pdf","pdf");
    c1->Clear();

    c1->cd();
    jer->Draw();
    c1->Print(output_name+".pdf)","pdf");
    c1->Clear();
    */
    


    return 0;
}
