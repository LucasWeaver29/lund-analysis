#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "vector"

#include "TCanvas.h"

using namespace Pythia8;

/*
This script uses the Pythia8 event simulator to generate proton-proton hard scattering events, 
and store all charged, final-state particles in a ROOT file. 

Primary parameters:
vector<vector<int>> pt_hat_bins - the pT hat bins used for event generation. The suggested binning (provided 
by Zoltan Varga) generates 2,700,000 events, which has been found to adequately populate the phase
space of the primary Lund plane.
*/


int main() {
    
    // Create a ROOT output file
    TString output_file_name = "event_test.root";
    TFile output_file(output_file_name, "RECREATE");

    // Maximum particle eta
    double eta_max = 5;
    
    
    // Declare pT hat bins in the form {pT min, pT max, # of events}
    
    // For a smaller test run
    /*
    std::vector<vector<int>> pt_hat_bins = {
        {20, 30, 10000},
        {30, 50, 8000},
        {50, 70, 4000},
        {70, 100, 2000},
        {100, 250, 4000},
    };
    */

    // For a complete analysis. These bins are courtesy of Zoltan Varga.
    std::vector<vector<int>> pt_hat_bins = {
        {20, 30, 1000000},
        {30, 50, 800000},
        {50, 70, 400000},
        {70, 100, 200000},
        {100, 130, 100000},
        {130, 180, 100000},
        {180, 250, 100000}
    };
    

    // ======================================================================

    // Create a TTree to store data on each particle for every event.
    TTree tree("events", "Pythia Events Tree");
    
    // Define variables to store event data
    int ptHatMin, ptHatMax, event_iBin, numParticles;
    double pTHat;
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
    tree.Branch("numParticles", &numParticles, "numParticles/I");
    tree.Branch("px", &px);
    tree.Branch("py", &py);
    tree.Branch("pz", &pz);
    tree.Branch("energy", &energy);
    tree.Branch("eta", &eta);
    tree.Branch("pdgId", &id);
    tree.Branch("pTHat", &pTHat);


    // Create a TTree to store the binning information in the same ROOT file.
    TTree ptHatBins_tree("ptHatBins", "Pythia pT Hat Bins Tree");

    vector<int> pt_hat_bin;

    ptHatBins_tree.Branch("pt_hat_bin", &pt_hat_bin);

    for (int iBin = 0; iBin < pt_hat_bins.size(); iBin++) {
        pt_hat_bin = pt_hat_bins[iBin];
        ptHatBins_tree.Fill();
    }
    ptHatBins_tree.Write();

    // ======================================================================
    
    /*
    The weighting (derived from event cross section) for each pT hat bin is not finalized until 
    all events in that bin are run. bin_weights stores these weights so they can be added back to
    the corresponding "Pythia Events Tree" entries later.
    */
    vector<double> bin_weights(pt_hat_bins.size());

    // =================================================================

    // Initialize Pythia
    Pythia pythia;
    pythia.readString("Beams:eCM = 5360");
    pythia.readString("HardQCD:all = on");
    pythia.readString("Next:numbershowEvent = 0");

    int tot_events = 0;

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
            
            // Clear data vectors for new event
            px.clear();
            py.clear();
            pz.clear();
            energy.clear();
            eta.clear();
            id.clear();

            // Fill particle data
            for (int iPart = 0; iPart < pythia.event.size(); ++iPart) {
                
                // Only take final, charged particles with eta < eta_max
                if (not pythia.event[iPart].isFinal()) continue;

                if (not pythia.event[iPart].isCharged()) continue;

                if (abs(pythia.event[iPart].eta()) > eta_max) continue;
                
                // record particle data
                const Pythia8::Particle& particle = pythia.event[iPart];
                
                px.push_back(particle.px());
                py.push_back(particle.py());
                pz.push_back(particle.pz());
                energy.push_back(particle.e());
                eta.push_back(particle.eta());
                id.push_back(particle.id());
            }

            numParticles = px.size();

            // Fill the tree
            tree.Fill();

        } // end event loop
        
        // Bin weight can be calculated once all events in the bin have been generated
        bin_weights[iBin] = pythia.info.sigmaGen() / nEvents;
        
    }

    
    // Write the tree to the file and close the file
    tree.Write();
    output_file.Close();


    // ======================================================================

    // Now the file must be reopened so the correct bin weights can be added

    // Reopen file to add bin_weight
    TFile *reopened_file = new TFile(output_file_name, "UPDATE");
    if (reopened_file->IsZombie()) {
        cout << "Error: Could not reopen file " << output_file_name << endl;
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
    // Loop through all entries in Pythia Events Tree to add the correct bin weight
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
