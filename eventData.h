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
    int ptHatMin, ptHatMax;
    double bin_weight;
    //double pTHat;
    int event_iBin; 
    int numParticles;
    std::vector<double>* px = nullptr;
    std::vector<double>* py = nullptr;
    std::vector<double>* pz = nullptr;
    std::vector<double>* energy = nullptr;
    std::vector<double>* eta = nullptr;
    std::vector<int>* pdgId = nullptr;
    
    TFile event_file;
    TTree* event_tree;

    std::vector<int>* pt_hat_bin = nullptr;
    TTree* bins_tree;

    bool debug = false;

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
        
        // set branches
        event_tree->SetBranchAddress("ptHatMin", &ptHatMin);
        event_tree->SetBranchAddress("ptHatMax", &ptHatMax);
        //event_tree->SetBranchAddress("crossSection", &crossSection);
        event_tree->SetBranchAddress("bin_weight", &bin_weight);
        event_tree->SetBranchAddress("numParticles", &numParticles);
        event_tree->SetBranchAddress("px", &px);
        event_tree->SetBranchAddress("py", &py);
        event_tree->SetBranchAddress("pz", &pz);
        event_tree->SetBranchAddress("energy", &energy);
        event_tree->SetBranchAddress("eta", &eta);
        event_tree->SetBranchAddress("pdgId", &pdgId);
        event_tree->SetBranchAddress("event_iBin", &event_iBin);
        //if (event_file_name == "events1000x9_ptHat.root") event_tree->SetBranchAddress("pTHat", &pTHat);

        // Get pT Hat Bins TTree
        bins_tree = (TTree*)event_file.Get("ptHatBins");
        if (!bins_tree) {
            std::cout << "ERROR: Could not find 'ptHatBins' tree in ROOT file" << std::endl;
        }
        bins_tree->SetBranchAddress("pt_hat_bin", &pt_hat_bin);



    } // end constructor

    int GetEntries() {
        return event_tree->GetEntries();
    }
    void GetEntry(int iEvent) {
        if(debug) std::cout << "Calling GetEntry" << std::endl;
        event_tree->GetEntry(iEvent);
    }

    int GetNumBins() {
        return bins_tree->GetEntries();
    }
    void GetBinEntry(int iBin) {
        bins_tree->GetEntry(iBin);
    }
    
    
    std::vector<fastjet::PseudoJet> get_particles(int iEvent, double part_eta_max, double part_pt_min = 0) { // Returns the final, charged, |etc| < max particles saved by storeWithRoot
        
        std::vector<fastjet::PseudoJet> event_particles;
        GetEntry(iEvent);
        //event_tree->GetEntry(iEvent);
        //weight_tree->GetEntry(event_iBin);
        //if(debug) std::cout << "in get_particles, looping through event particles" << std::endl;
        //std::cout << "numparticles: " << numParticles << std::endl;
        //std::cout << "Array size" << (*px).size() << std::endl;
        for (int iPart=0; iPart < numParticles; iPart++) {
            if (debug) std::cout << "iPart = " << iPart << std::endl;

            if(std::fabs((*eta)[iPart]) > part_eta_max) continue;
            //if (debug) std::cout << "Got eta" << std::endl;
            if((pow((*px)[iPart],2) + pow((*py)[iPart],2)) < pow(part_pt_min,2)) continue;
            //if (debug) std::cout << "iPart = " << iPart << ", b" << std::endl;
            //if (debug) std::cout << "   pushing back to event_particles" << std::endl;
            event_particles.push_back(fastjet::PseudoJet((*px)[iPart], (*py)[iPart], (*pz)[iPart], (*energy)[iPart]));
            //std::cout << "pushed" << std::endl;
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

    std::vector<fastjet::PseudoJet> get_particles(int iEvent, double part_pt_min = 0) {
        // Particle loop. These are the final, charged, |etc| < max particles saved by storeWithRoot
        std::vector<fastjet::PseudoJet> bg_particles;
        bg_tree->GetEntry(iEvent);
        for (int iPart=0; iPart<numBGparts; ++iPart) {
            if ((*bgPts)[iPart] < part_pt_min) continue;
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


//=============================================
// For many files
//=============================================

class my_background_trees {

    public:

    int numBGparts;
    std::vector<double>* bgPts = nullptr;
    std::vector<double>* bgPhis = nullptr;
    std::vector<double>* bgEtas = nullptr;
    
    TFile* bg_file = nullptr;
    TTree* bg_tree = nullptr;

    // For keeping track of which file to use
    int file_index = -1;
    std::vector<int> entries_per_file;
    std::vector<TString> bg_file_names;

    bool debug = false;

    // constructor
    my_background_trees(std::vector<TString> bg_file_names_in):
    entries_per_file(bg_file_names_in.size()),
    bg_file_names(bg_file_names_in.size())
    {
        bg_file_names = bg_file_names_in;

        for (int iFile = 0; iFile < bg_file_names.size(); iFile++) {

            if(debug) std::cout << "Opening " << bg_file_names[iFile] << std::endl;

            bg_file = TFile::Open(bg_file_names[iFile], "READ");
            check_file(iFile);

            bg_tree = (TTree*)bg_file->Get("backgrounds");
            check_tree(iFile);

            entries_per_file[iFile] = bg_tree->GetEntries();
            if(debug) std::cout << entries_per_file[iFile] << " entries in " << bg_file_names[iFile] << std::endl;


            bg_file->Close();


        }
        
        
    } // end constructor

    int GetEntries() {
        int all_entries = 0;
        for (int entries : entries_per_file) {
            all_entries += entries;
        }
        return all_entries;
    }

    void GetEntry(int iEvent) {
        for (int i = 0; i < entries_per_file.size(); i++) {
            if (iEvent > entries_per_file[i]) {
                iEvent -= entries_per_file[i];
            }
            else {
                if (i > file_index) open_next_file();
                break;
            }
        }
        
        bg_tree->GetEntry(iEvent);
    }

    void check_file(int iFile) { // takes iFile so cout can print correct name
        if (!bg_file || bg_file->IsZombie()) std::cout << "Error: Could not open ROOT background file " << bg_file_names[iFile] << std::endl;
    }

    void check_tree(int iFile) {
        if (!bg_tree) std::cout << "Error: Could not find 'backgrounds' tree in " << bg_file_names[iFile] << "(probably bc it was called events before or vv)" << std::endl;
    }

    void open_next_file() {
        
        file_index++;
        if (file_index >= entries_per_file.size()) std::cout << "ERROR: background file_index has exceeded number of files" << std::endl;

        if (bg_file && bg_file->IsOpen()) {
            bg_file->Close();
        }
        
        bg_file = TFile::Open(bg_file_names[file_index], "READ");
        check_file(file_index);

        bg_tree = (TTree*)bg_file->Get("backgrounds");
        check_tree(file_index);

        bg_tree->SetBranchAddress("numParticles", &numBGparts);
        bg_tree->SetBranchAddress("pts", &bgPts);
        bg_tree->SetBranchAddress("phis", &bgPhis);
        bg_tree->SetBranchAddress("etas", &bgEtas);
    }

    std::vector<fastjet::PseudoJet> get_particles(int iEvent, double part_pt_min = 0) {
        // Particle loop. These are the final, charged, |etc| < max particles saved by storeWithRoot
        if (debug) std::cout << "get_particles in my bg trees called" << std::endl;
        std::vector<fastjet::PseudoJet> bg_particles;
        GetEntry(iEvent);
        for (int iPart=0; iPart<numBGparts; ++iPart) {
            if ((*bgPts)[iPart] < part_pt_min) continue;
            fastjet::PseudoJet bg_prtcl;
            bg_prtcl.reset_PtYPhiM((*bgPts)[iPart],(*bgEtas)[iPart],(*bgPhis)[iPart], 0);
            bg_particles.push_back(bg_prtcl);        
        }
        return bg_particles;
    }

    void close_file() {
        bg_file->Close();
    }


};





#endif



    
    
    

