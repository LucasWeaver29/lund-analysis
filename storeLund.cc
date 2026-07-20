#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"


// ROOT - files
#include <vector>
#include <iostream>

// ROOT - histogram 
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"

#include "eventData.h"
#include "storeLund.h"

using namespace std;


int main() {
    
    //TString output_file_name = "9k, groomed";
    //TString output_folder = "LundPlanes/LundRootFiles";

    cout << "Make sure you have created the needed subfolder!" << endl;

    TString notes = "";
    //"z_cut = .2, beta = 3";

    
    lund_sub_options lund_ops = lund_sub_options{.name = "Embedded, SD (mine, all), from event_Zoltans", .groom_option = "SD_mine_all"};
        
    bool embed_in_bg = true;

    TString eventFileName = "event_Zoltans.root";
    vector<TString> backgroundFileNames = {
        "backgrounds2.7m.root", 
        "backgrounds2.7m_1.root",
        "backgrounds2.7m_2.root",
        "backgrounds2.7m_3.root",
        "backgrounds2.7m_4.root",
    };
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
    if (debug) cout << "Creating my_event_tree" << endl;
    my_event_tree met(eventFileName);
    if (debug) cout << "Creating my_background_trees" << endl;
    my_background_trees mbt(backgroundFileNames);
    // Checking if output subfolder exists
    //std::filesystem::path current_path = std::filesystem::current_path();
    //std::filesystem::path full_path = current_path / output_folder / subfolder

    //=======================================================================
    // 
    TFile output_file(lund_ops.name + ".root", "RECREATE");
    
    // Declare TTrees for storing Lund points, and total jet pt
    TTree* lund_coords_tree = new TTree("lund_coords_tree", "Lund Plane Points");
        
    // Variables for tree branches
    double log_inv_R, log_kt, jet_pt, bin_weight;
    int iBin;

    lund_coords_tree->Branch("log_inv_R", &log_inv_R, "log_inv_R/D");
    lund_coords_tree->Branch("log_kt", &log_kt, "log_kt/D");
    lund_coords_tree->Branch("jet_pt", &jet_pt, "jet_pt/D");
    lund_coords_tree->Branch("bin_weight", &bin_weight, "bin_weight/D");
    lund_coords_tree->Branch("iBin", &iBin, "iBin/I");


    TH1F* weighted_jet_pt_hist = new TH1F("weighted_jet_pt_hist", "Weighted Jet p_{T}; jet p_{T}; Weighted counts", 300, 0, 300);
    lund_ops.weighted_jet_pt_hist = weighted_jet_pt_hist;


    // Set up fastjet analysis
    LundGroomer lund_groomer(0, 400, leading_jet_only, Rparam, jet_eta_max, part_eta_max); // 0, 400 for no jet pt cut
    
    // =======================================================================

    // Reading events from ROOT
    int numEvents = met.GetEntries();
    //std::cout << "Number of events: " << numEvents << std::endl;

    
    int numBackgrounds;
    if (embed_in_bg) {
        numBackgrounds = mbt.GetEntries();
        cout << "There are " << numBackgrounds << " available backgrounds" << endl;
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
        vector<fastjet::PseudoJet> event_particles = met.get_particles(iEvent, part_eta_max, part_pt_min); // get_particles also calls getEntry(iEvent)
        
        if (embed_in_bg) {
            vector<fastjet::PseudoJet> background_prtcls = mbt.get_particles(iEvent, part_pt_min);
            move(background_prtcls.begin(), background_prtcls.end(), back_inserter(event_particles));
        }

        lund_ops.bin_weight = met.bin_weight; // Need to give bin weight to lund_groomer so weighted jet pt hist can be filled with correct bin weight

        if (debug) cout << "Calling lund_groomed.get_kin_vars()" << endl;
        vector<lund_kin_vars> all_kin_vars = lund_groomer.get_kin_vars(event_particles, lund_ops);
        
        for (lund_kin_vars vars : all_kin_vars) {
            log_inv_R = log(Rparam/vars.delta);
            log_kt = log(vars.kt);
            jet_pt = vars.jet_pt;
            bin_weight = met.bin_weight;
            iBin = met.event_iBin;

            if (debug) cout << "Filling lund_coords_tree" << endl;
            lund_coords_tree->Fill();

        }
        
    }  // End reading events from ROOT
    

    // =========================================================
    // Save these plots in a root file, so they can be modified/resized without having to run the whole thing again
    

    lund_coords_tree->Write();

    TTree jet_pt_tree("jet_pt_tree", "Jet pT Tree");
    jet_pt_tree.Branch("jet_pt_hist", &weighted_jet_pt_hist);
    jet_pt_tree.Fill();
    jet_pt_tree.Write();

    output_file.Close();

    met.close_file();
    mbt.close_file();



    // Make a card with specs about what's in the lund plane
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    // Adding a text information sheet with stats about this generation
    TString title = "Lund plane points: " + lund_ops.name;
    TString line1 = "Number of events: " + to_string(numEvents);
    TString line2 = "Events from: " + eventFileName;
    TString line3 = embed_in_bg? "Backgrounds from: " + print_vec(backgroundFileNames) : "No background.";
    TString line4 = "Leading Jet Only: " + bool2Str(leading_jet_only);

    TString line5 = "Jet radius: " + to_string_round(Rparam);
    TString line6 = "Particle eta max: " + to_string_round(part_eta_max);
    TString line7 = "Jet eta max: " + to_string_round(jet_eta_max);

    vector<TString> lines = {title, notes, line1, line2, line3, line4, line5, line6, line7};

    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }
    
    c1->Print(lund_ops.name + "-Data_card.pdf","pdf");



    return 0;
}



