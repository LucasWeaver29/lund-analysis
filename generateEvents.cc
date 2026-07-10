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
    TFile file("events1000x9, bin weight.root", "RECREATE");

    // number of events per bin
    // do 5000
    int nEvents = 1000;

    double etaMax = 5;
    
    // Give bounds for each bin
    vector<double> ptHatMins = {10, 20, 30, 50, 70, 90, 110, 130, 150};
    vector<double> ptHatMaxs = {20, 30, 50, 70, 90, 110, 130, 150, 250};
    //vector<double> ptHatMins = {10, 50, 70, 90, 110, 130, 150};
    //vector<double> ptHatMaxs = {50, 70, 90, 110, 130, 150, 250};
    
    //vector<double> ptHatMins = {5};
    //vector<double> ptHatMaxs = {300};

    //vector<double> ptHatMins = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
    //vector<double> ptHatMaxs = {20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250};


    // Create a TTree
    TTree tree("events", "Pythia Events Tree");
    
    // Define variables to store event data
    double ptHatMin, ptHatMax; 
    double crossSection, weight;
    double iBin;
    double pTHat;
    int numParticles;
    std::vector<double> px;
    std::vector<double> py;
    std::vector<double> pz;
    std::vector<double> energy;
    std::vector<double> eta;
    std::vector<int> id;
    
    // Create branches for the tree
    tree.Branch("ptHatMin", &ptHatMin, "ptHatMin/D");
    tree.Branch("ptHatMax", &ptHatMax, "ptHatMax/D");
    tree.Branch("iBin", &iBin, "iBin/D");
    //tree.Branch("crossSection", &crossSection, "crossSection/D");
    //tree.Branch("weight", &weight, "weight/D");
    tree.Branch("numParticles", &numParticles, "numParticles/I");
    tree.Branch("px", &px);
    tree.Branch("py", &py);
    tree.Branch("pz", &pz);
    tree.Branch("energy", &energy);
    tree.Branch("eta", &eta);
    tree.Branch("pdgId", &id);
    tree.Branch("pTHat", &pTHat);


    
    TTree weight_tree("weights", "Pythia Bin Weights Tree");
    double bin_weight;
    weight_tree.Branch("bin_weight", &bin_weight, "bin_weight/D");
    

    // For testing how pythia.info.sigmaGen() evolves over time. 
    TGraph* xsecs = new TGraph();
    xsecs->SetTitle("Pythia event cross sections over events");
    

    // =================================================================

    // Initialize Pythia
    Pythia pythia;
    pythia.readString("Beams:eCM = 5360");
    pythia.readString("HardQCD:all = on");
    pythia.readString("Next:numbershowEvent = 0");

    // Bin loop
    for (int iBin = 0; iBin < ptHatMins.size(); ++iBin) {

        pythia.settings.readString("PhaseSpace:pTHatMin = " + to_string(ptHatMins[iBin]));
        pythia.settings.readString("PhaseSpace:pTHatMax = " + to_string(ptHatMaxs[iBin]));

        // Proton - proton collision is default
        pythia.init();

        // Event loop
        for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
            
            if (!pythia.next()) continue; // Skip failed events

            
            // Fill event data
            ptHatMin = ptHatMins[iBin];
            ptHatMax = ptHatMaxs[iBin];

            pTHat = pythia.info.pTHat();
            
            //crossSection = pythia.info.sigmaGen();
            //weight = pythia.info.sigmaGen() / nEvents;
            xsecs->AddPoint(iEvent + iBin * nEvents, pythia.info.sigmaGen());

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
        bin_weight = pythia.info.sigmaGen() / nEvents;
        weight_tree.Fill();
    }


    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    c1->cd();
    c1->SetLogy(1);
    xsecs->Draw();
    c1->Print("Pythia Cross Sections Over Time.pdf","pdf");
    c1->Clear();


    // Write the tree to the file
    tree.Write();
    weight_tree.Write();
    file.Close();


    return 0;
}
