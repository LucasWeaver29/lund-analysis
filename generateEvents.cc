// Generates events and stores them in ROOT

#include <Pythia8/Pythia.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>

#include "TGraph.h"
#include "TCanvas.h"

using namespace Pythia8;

int main() {
    
    // Create a ROOT output file
    TString filename = "event_scaled3Large.root";
    TFile file(filename, "RECREATE");

    // number of events per bin
    // do 5000

    double etaMax = 5;
    
    // Give bounds for each bin
    
    // Scaled3 Lite
    std::vector<vector<int>> pt_hat_bins = {
        {20, 30, 160000},
        {30, 50, 32500},
        {50, 70, 13500},
        {70, 90, 9500},
        {90, 110, 7750},
        {110, 130, 7000},
        {130, 150, 6000},
        {150, 250, 6000}
    };
    
    
    
    /*
    // Scaled3 
    std::vector<vector<int>> pt_hat_bins = {
        {20, 30, 64000},
        {30, 50, 12000},
        {50, 70, 5100},
        {70, 90, 3800},
        {90, 110, 2800},
        {110, 130, 2700},
        {130, 150, 2500},
        {150, 250, 2500}
    };
    
    */


    /*
    vector<vector<int>> pt_hat_bins = {
        {10, 20, 2000},
        {20, 30, 2000},
        {30, 50, 2000},
        {50, 70, 2000},
        {70, 90, 2000},
        {90, 110, 2000},
        {110, 130, 2000},
        {130, 150, 2000},
        {150, 250, 2000}
    };
    */
    
    




    //vector<double> ptHatMins = {10, 20, 30, 50, 70, 90, 110, 130, 150};
    //vector<double> ptHatMaxs = {20, 30, 50, 70, 90, 110, 130, 150, 250};
    //vector<double> ptHatMins = {10, 50, 70, 90, 110, 130, 150};
    //vector<double> ptHatMaxs = {50, 70, 90, 110, 130, 150, 250};
    
    


    // Create a TTree
    TTree tree("events", "Pythia Events Tree");
    
    // Define variables to store event data
    int ptHatMin, ptHatMax; 
    int event_iBin;
    double pTHat;
    int numParticles;
    std::vector<double> px;
    std::vector<double> py;
    std::vector<double> pz;
    std::vector<double> energy;
    std::vector<double> eta;
    std::vector<int> id;
    
    // Create branches for the tree
    tree.Branch("ptHatMin", &ptHatMin, "ptHatMin/I");
    tree.Branch("ptHatMax", &ptHatMax, "ptHatMax/I");
    tree.Branch("event_iBin", &event_iBin, "event_iBin/I");
    //tree.Branch("crossSection", &crossSection, "crossSection/D");
    tree.Branch("numParticles", &numParticles, "numParticles/I");
    tree.Branch("px", &px);
    tree.Branch("py", &py);
    tree.Branch("pz", &pz);
    tree.Branch("energy", &energy);
    tree.Branch("eta", &eta);
    tree.Branch("pdgId", &id);
    tree.Branch("pTHat", &pTHat);


    vector<int> pt_hat_bin;
    TTree ptHatBins_tree("ptHatBins", "Pythia pT Hat Bins Tree");
    ptHatBins_tree.Branch("pt_hat_bin", &pt_hat_bin);
    for (int iBin = 0; iBin < pt_hat_bins.size(); iBin++) {
        pt_hat_bin = pt_hat_bins[iBin];
        ptHatBins_tree.Fill();
    }
    ptHatBins_tree.Write();

    // For testing how pythia.info.sigmaGen() evolves over time. 
    TGraph* xsecs = new TGraph();
    xsecs->SetTitle("Pythia event cross sections over events");
    
    // To store bin weights such that they can be added back to the tree entries later
    vector<double> bin_weights(pt_hat_bins.size());

    // =================================================================

    // Initialize Pythia
    Pythia pythia;
    pythia.readString("Beams:eCM = 5360");
    pythia.readString("HardQCD:all = on");
    pythia.readString("Next:numbershowEvent = 0");

    int tot_events = 0

    // Bin loop
    for (int iBin = 0; iBin < pt_hat_bins.size(); ++iBin) {

        ptHatMin = pt_hat_bins[iBin][0];
        ptHatMax = pt_hat_bins[iBin][1];
        int nEvents = pt_hat_bins[iBin][2];
        tot_events++;

        pythia.settings.readString("PhaseSpace:pTHatMin = " + to_string(ptHatMin));
        pythia.settings.readString("PhaseSpace:pTHatMax = " + to_string(ptHatMax));

        // Proton - proton collision is default
        pythia.init();

        event_iBin = iBin;

        // Event loop
        for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
            
            if (!pythia.next()) continue; // Skip failed events

            
            // Fill event data
                    
            //crossSection = pythia.info.sigmaGen();
            //weight = pythia.info.sigmaGen() / nEvents;
            xsecs->AddPoint(tot_events, pythia.info.sigmaGen());

            // Clear vectors for new event
            px.clear();
            py.clear();
            pz.clear();
            energy.clear();
            eta.clear();
            id.clear();

            // Fill particle data
            for (int iPart = 0; iPart < pythia.event.size(); ++iPart) {
                
                // Only take final, charged particles in correct eta bin
                if (not pythia.event[iPart].isFinal()) continue;

                if (abs(pythia.event[iPart].eta()) > etaMax) continue;

                if (not pythia.event[iPart].isCharged()) continue;
                

                const Pythia8::Particle& particle = pythia.event[iPart];
                //std::cout << "    Particle " << iPart << ":" << std::endl;
                //std::cout << "      Pythia says px is " << pythia.event[iPart].px() << std::endl;
                px.push_back(particle.px());
                //std:cout << "      px vector says its "<< px.back() << endl;
                py.push_back(particle.py());
                pz.push_back(particle.pz());
                energy.push_back(particle.e());
                eta.push_back(particle.eta());
                id.push_back(particle.id());
                
                // Doesn't work bc iPart doesn't match px idexes
                /*
                std::cout << "      px: " << px[iPart] << std::endl;
                std::cout << "      py: " << py[iPart] << std::endl;
                std::cout << "      pz: " << pz[iPart] << std::endl;
                std::cout << "      energy: " << energy[iPart] << std::endl;
                std::cout << "      pdgId: " << id[iPart] << std::endl;
                */
            }

            numParticles = px.size();

            // Fill the tree
            tree.Fill();
        } // end event loop
        
        bin_weights[iBin] = pythia.info.sigmaGen() / nEvents;
        
    }


    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    c1->cd();
    c1->SetLogy(1);
    xsecs->Draw();
    c1->Print("Pythia Cross Sections Over Time.pdf","pdf");
    c1->Clear();


    // Write the tree to the file
    tree.Write();
    //weight_tree.Write();
    file.Close();




    // Reopen file to add bin_weight
    TFile *reopened_file = new TFile(filename, "UPDATE");
    if (reopened_file->IsZombie()) {
        cout << "Error: Could not reopen file " << filename << endl;
        return 1;
    }

    // Reaccess the tree
    TTree* reopened_tree = (TTree*)reopened_file->Get("events");
    if (!reopened_tree) {
        std::cout << "Error: Could not find 'events' tree in ROOT file" << std::endl;
    }

    // Add new branch for bin weights
    double bin_weight = 0.0;
    TBranch *weight_branch = reopened_tree->Branch("bin_weight", &bin_weight, "bin_weight/D");

    // Set branch address for iBin
    int reopened_iBin;
    reopened_tree->SetBranchAddress("event_iBin", &reopened_iBin);

    int nEntries = reopened_tree->GetEntries();
    for (int iEntry = 0; iEntry < nEntries; iEntry++) {
        reopened_tree->GetEntry(iEntry);
        bin_weight = bin_weights[reopened_iBin];
        weight_branch->Fill();
    }

    // Save with overwrite flag
    reopened_tree->Write("", TObject::kOverwrite);
    reopened_file->Close();


    return 0;
}
