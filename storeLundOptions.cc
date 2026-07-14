#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"


// ROOT - files
//#include <TFile.h>
//#include <TTree.h>
//#include <TBranch.h>

#include <vector>
#include <iostream>
#include <filesystem>

// ROOT - histogram
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLatex.h"

#include "eventData.h"
#include "storeLund.h"

using namespace Pythia8;


int main() {
    
    //TString output_file_name = "9k, groomed";
    TString output_folder = "LundPlanes/LundRootFiles";

    cout << "Make sure you have created the needed subfolder!" << endl;

    TString notes = "";
    //"z_cut = .2, beta = 3";

    
    vector<lund_sub_options> pythia_lund_options = 
        {
            lund_sub_options{.name = "Pythia, 2.3mil"},
            //lund_sub_options{.name = "Pythia, SD (mine, all)", .groom_option =  "SD_mine_all"},
        };

    vector<lund_sub_options> bg_lund_options = 
        {
            //lund_sub_options{.name = "Background, SD (mine, all)", .groom_option = "SD_mine_all"},
            //lund_sub_options{.name = "Constituent subtraction, SD (mine, all), .event_sub = "ConSub", .groom_option = "SD_mine_all"}
        };


    TString eventFileName = "event_2.3mil.root";
    TString backgroundFileName = "backgrounds50k.root";
    //"thermalBackgrounds9000,etaMax=2.root";

    double part_eta_max = .9;
    double jet_eta_max = .5;
    double part_pt_min = .15;

    double Rparam = .4;
    bool leading_jet_only = false;
    if (part_eta_max - Rparam != jet_eta_max) cout << "WARNING: Rparam, part_eta_max, and jet_eta_max are not in agreement" << endl;

    bool debug = false;

    TString space = " ";

    //==================================================================
    // Preparing to read from ROOT files
    my_event_tree met(eventFileName);
    my_background_tree mbt(backgroundFileName);
    // Checking if output subfolder exists
    //std::filesystem::path current_path = std::filesystem::current_path();
    //std::filesystem::path full_path = current_path / output_folder / subfolder

    //=======================================================================
    // Declare TTrees for storing Lund points, and total jet pt
    vector<TTree*> pythia_lund_coords_trees;

    // lund_coords store the variables that the TTrees reference for filling branches. Each lund_coords in the pythia_lund_coords_holders vector is for a different option
    vector<lund_coords> pythia_lund_coords_holders;

    vector<TH1F*> pythia_weighted_jet_pt_hists;



    if (debug) cout << "Abt to assemble pythia_kinvars_trees and pythia_jetpt_trees" << endl;
    for (int iOption = 0; iOption < pythia_lund_options.size(); iOption++) {
        
        TTree* kinvars_tree = new TTree(pythia_lund_options[iOption].name + space + "kinvars", pythia_lund_options[iOption].name + space + "Kinematic Variables");
        
        pythia_lund_coords_holders.push_back(lund_coords{});

        kinvars_tree->Branch("log_inv_R", &pythia_lund_coords_holders[iOption].log_inv_R, "log_inv_R/D");
        kinvars_tree->Branch("log_kt", &pythia_lund_coords_holders[iOption].log_kt, "log_kt/D");
        kinvars_tree->Branch("jet_pt", &pythia_lund_coords_holders[iOption].jet_pt, "jet_pt/D");
        kinvars_tree->Branch("bin_weight", &pythia_lund_coords_holders[iOption].bin_weight, "bin_weight/D");
        kinvars_tree->Branch("iBin", &pythia_lund_coords_holders[iOption].iBin, "iBin/I");


        pythia_lund_coords_trees.push_back(kinvars_tree);


        TH1F* weighted_jet_pt_hist = new TH1F(pythia_lund_options[iOption].name + space + "weighted_jet_pt_hist", "Weighted Jet pT Hist", 300, 0, 300);
        pythia_weighted_jet_pt_hists.push_back(weighted_jet_pt_hist);
        pythia_lund_options[iOption].weighted_jet_pt_hist = pythia_weighted_jet_pt_hists[iOption];

    }
    

    // Set up fastjet analysis
    LundGroomer lund_groomer(0, 400, leading_jet_only, Rparam, jet_eta_max, part_eta_max);
    
    // =======================================================================

    // Reading events from ROOT
    int numEvents = met.GetEntries();
    //std::cout << "Number of events: " << numEvents << std::endl;

    
    int numBackgrounds;
    if (bg_lund_options.size() != 0) {
        numBackgrounds = mbt.GetEntries();
        if (numBackgrounds < numEvents) {
            cout << "Warning: There are " << numEvents << " provided events but only" << numBackgrounds << " provided backgrounds." << endl;
            cout << "Lund plane creation will stop once backgrounds run out" << endl;
            numEvents = numBackgrounds;
        }
    
    }
    
        

    int counter = 1000;

    // Loop over events in ROOT TTree
    for (int iEvent = 0; iEvent < numEvents; ++iEvent) {
        

        if (debug) cout << "iEvent: " << iEvent << endl;
        if (iEvent > counter) {
            cout << counter << " events analyzed from ROOT file" << endl;
            counter += 1000;
        }

        if (debug) cout << "Calling met.get_particles" << endl;
        vector<fastjet::PseudoJet> event_particles = met.get_particles(iEvent, part_eta_max, part_pt_min); // iEvent, part_eta_max
        
        
        if (debug) cout << "Looping through pythia_lund_options" << endl;
        
        
        for (int iOption = 0; iOption < pythia_lund_options.size(); ++iOption) {
            if (debug) cout << "iOption " << iOption << endl;
            
            pythia_lund_options[iOption].bin_weight = met.bin_weight; // Need to give bin weight to lund_groomer so weighted jet pt hist can be filled with correct bin weight

            vector<lund_kin_vars> all_kin_vars = lund_groomer.get_kin_vars(event_particles, pythia_lund_options[iOption]);
            
            for (lund_kin_vars vars : all_kin_vars) {
                pythia_lund_coords_holders[iOption].log_inv_R = log(Rparam/vars.delta);
                pythia_lund_coords_holders[iOption].log_kt = log(vars.kt);
                pythia_lund_coords_holders[iOption].jet_pt = vars.jet_pt;
                pythia_lund_coords_holders[iOption].bin_weight = met.bin_weight;
                pythia_lund_coords_holders[iOption].iBin = met.event_iBin;


                pythia_lund_coords_trees[iOption]->Fill();

            }
        }
        

        // Embedding in the background
        /*
        if (bg_lund_options.size() != 0) {
            vector<fastjet::PseudoJet> background_prtcls = mbt.get_particles(iEvent, part_pt_min);
            move(background_prtcls.begin(), background_prtcls.end(), back_inserter(event_particles));

            if (debug) cout << "Looping through background_lund_options" << endl;
            for (int iOption = 0; iOption < bg_lund_options.size(); ++iOption) {
                if (debug) cout << "iOption " << iOption << endl;
                bg_lund_options[iOption].sub_options.bin_weight = met.bin_weight;
                //bg_lund_options[iOption].sub_options.jet_pt_hist
                vector<lund_kin_vars> all_kin_vars = lund_groomer.get_kin_vars(event_particles, bg_lund_options[iOption].sub_options);
                for (lund_kin_vars vars : all_kin_vars) {
                    double logR0_R = log(Rparam/vars.delta);
                    double logkt = log(vars.kt);
                    bg_lund_cuts[iOption][0]->Fill(logR0_R,logkt, met.bin_weight); // The inclusive pt lund
                    bg_lund_cuts[iOption][vars.num_cut + 1]->Fill(logR0_R,logkt, met.bin_weight); // the jet pt cut lund (which is at iCut + 1 in pythia_lund_cuts);
                }
            }
        }
        */

    }  // End reading events from ROOT
    

    // =========================================================
    // Save these plots in a root file, so they can be modified/resized without having to run the whole thing again
    

    for (int iOption = 0; iOption < pythia_lund_options.size(); iOption ++) {

        TFile output_file(output_folder + "/" + pythia_lund_options[iOption].name + ".root", "RECREATE");

        pythia_lund_coords_trees[iOption]->Write();

        TTree jet_pt_tree(pythia_lund_options[iOption].name + space + "jet_pt_tree", "Jet pT Tree");
        jet_pt_tree.Branch("jet_pt_hist", &pythia_weighted_jet_pt_hists[iOption]);
        jet_pt_tree.Fill();
        jet_pt_tree.Write();
        output_file.Close();
    
    }


    met.close_file();
    mbt.close_file();

    return 0;
}



