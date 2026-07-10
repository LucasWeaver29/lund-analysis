#ifndef EVENTDATA_H
#define EVENTDATA_H

#include "fastjet/PseudoJet.hh"

// ROOT - files
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include <vector>



// This file will read from Root event files, for easy access to pythia and background data


class my_event_tree {

public:
    // Variables to access event data from TTrees
    double ptHatMin, ptHatMax, crossSection, weight;
    double pTHat;
    double iBin; // bin_weight
    int numParticles;
    std::vector<double>* px = nullptr;
    std::vector<double>* py = nullptr;
    std::vector<double>* pz = nullptr;
    std::vector<double>* energy = nullptr;
    std::vector<double>* eta = nullptr;
    std::vector<int>* pdgId = nullptr;
    
    TFile event_file;
    TTree* event_tree;

    TTree* weight_tree; // bin_weight
    double bin_weight; // bin_weight

    // constructor
    my_event_tree(TString event_file_name):
    event_file(event_file_name, "READ")
    {
        // Check the event ROOT file
        if (event_file.IsZombie()) {
            std::cout << "Error: Could not open ROOT event file" << std::endl;
        }
        // Get event TTree
        event_tree = (TTree*)event_file.Get("events");
        if (!event_tree) {
            std::cout << "Error: Could not find 'events' tree in ROOT file" << std::endl;
        }
        // Get weight TTree // bin_weight
        weight_tree = (TTree*)event_file.Get("weights");
        if (!event_tree) {
            std::cout << "Error: Could not find 'weights' tree in ROOT file" << std::endl;
        }


        // set branches
        event_tree->SetBranchAddress("ptHatMin", &ptHatMin);
        event_tree->SetBranchAddress("ptHatMax", &ptHatMax);
        //event_tree->SetBranchAddress("crossSection", &crossSection);
        //event_tree->SetBranchAddress("weight", &weight);
        event_tree->SetBranchAddress("numParticles", &numParticles);
        event_tree->SetBranchAddress("px", &px);
        event_tree->SetBranchAddress("py", &py);
        event_tree->SetBranchAddress("pz", &pz);
        event_tree->SetBranchAddress("energy", &energy);
        event_tree->SetBranchAddress("eta", &eta);
        event_tree->SetBranchAddress("pdgId", &pdgId);
        if (event_file_name == "events1000x9_ptHat.root") event_tree->SetBranchAddress("pTHat", &pTHat);
        weight_tree->SetBranchAddress("bin_weight", &bin_weight); // bin_weight

    } // end constructor

    int GetEntries() {
        return event_tree->GetEntries();
    }
    void GetEntry(int iEvent) {
        event_tree->GetEntry(iEvent);
        weight_tree->GetEntry(iBin);
    }
    
    
    std::vector<fastjet::PseudoJet> get_particles(int iEvent, double part_eta_max = .9) { // Returns the final, charged, |etc| < max particles saved by storeWithRoot
        
        std::vector<fastjet::PseudoJet> event_particles;
        GetEntry(iEvent);
        //event_tree->GetEntry(iEvent);
        //weight_tree->GetEntry(iBin);
        for (int iPart=0; iPart<numParticles; ++iPart) {

            if(abs((*eta)[iPart]) > part_eta_max) continue;
            
            event_particles.push_back(fastjet::PseudoJet((*px)[iPart], (*py)[iPart],(*pz)[iPart],(*energy)[iPart]));
        }
        return event_particles;
    }


    void close_file() {
        event_file.Close();
    }



};


class my_background_tree {

public:

    int numBGparts;
    std::vector<double>* bgPts = nullptr;
    std::vector<double>* bgPhis = nullptr;
    std::vector<double>* bgEtas = nullptr;
    
    TFile bg_file;
    TTree* bg_tree;

    // constructor
    my_background_tree(TString background_file_name):
    bg_file(background_file_name, "READ")
    {
        if (bg_file.IsZombie()) { // check the ROOT file
            std::cout << "Error: Could not open ROOT background file" << std::endl;
        }

        bg_tree = (TTree*)bg_file.Get("backgrounds");
        if (!bg_tree) {
            std::cout << "Error: Could not find 'backgrounds' tree in ROOT file (probably bc it was called events before or vv)" << std::endl;
        }

        // set branches
        bg_tree->SetBranchAddress("numParticles", &numBGparts);
        bg_tree->SetBranchAddress("pts", &bgPts);
        bg_tree->SetBranchAddress("phis", &bgPhis);
        bg_tree->SetBranchAddress("etas", &bgEtas);

    } // end constructor

    int GetEntries() {
        return bg_tree->GetEntries();
    }
    void GetEntry(int iEvent) {
        bg_tree->GetEntry(iEvent);
    }

    std::vector<fastjet::PseudoJet> get_particles(int iEvent) {
        // Particle loop. These are the final, charged, |etc| < max particles saved by storeWithRoot
        std::vector<fastjet::PseudoJet> bg_particles;
        bg_tree->GetEntry(iEvent);
        for (int iPart=0; iPart<numBGparts; ++iPart) {
            fastjet::PseudoJet bg_prtcl;
            bg_prtcl.reset_PtYPhiM((*bgPts)[iPart],(*bgEtas)[iPart],(*bgPhis)[iPart], 0);
            bg_particles.push_back(bg_prtcl);        
        }
        return bg_particles;
    }

    void close_file() {
        bg_file.Close();
    }


};

#endif



    
    
    


