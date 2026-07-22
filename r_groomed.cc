
#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
//#include "fastjet/tools/Filter.hh"
//#include "fastjet/tools/Pruner.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh" // jet by jet constituent-based subtraction
#include "fastjet/contrib/SoftKiller.hh"

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
#include "TLegend.h"

// For the background
#include "backSubTools.h"
#include "eventData.h"

using namespace Pythia8;

struct r_groomed_vars {
    double delta, jet_pt;
}

class delta_getter {

public:
    
    // General
    double part_eta_max;
    double jet_eta_max;
    double Rparam;
    bool leading_jet_only;

    fastjet::JetDefinition jetDef_akt; // for jet identification
    fastjet::JetDefinition jetDef_ca_recluster; // for reclustering the antikt-identified jets to make Lund plane
    //fastjet::JetDefinition jetDef_kt_recluster(fastjet::kt_algorithm, Rparam, recombScheme, strat); //  ???

    fastjet::Selector jet_eta_selector;


    // Softdrop
    double z_cut = .2;
    double beta = 0; // beta = 0 does just a normal z cut

    // my pc
    my_pc_subtractor pc_subtractor;

    // Filtering:
    /*
    double Rfilt = .3;
    double nfilt = 3;
    fastjet::Filter filter;

    // Pruning
    double prune_zcut = .2;
    double Rcut_factor = .2;
    fastjet::Pruner pruner;
    */

    // Jet by jet constituent subtraction
    fastjet::GridMedianBackgroundEstimator bge; // GridMedianBackgroundEstimator is faster than JetMedianBackgroundEstimator and "performs equally well in nearly all cases"
    fastjet::contrib::ConstituentSubtractor subtractor;

    // SoftKiller
    fastjet::contrib::SoftKiller soft_killer;

    
    // RecursiveSoftDrop - fastjet contrib
    double z_cut_rsd = .2;
    double beta_rsd = 0;
    fastjet::contrib::RecursiveSoftDrop rsd;
    

    // constructor
    delta_getter(bool leading_jet_only_in, double Rparam_in, double jet_eta_max_in, double part_eta_max_in):
    jetDef_akt(
        fastjet::antikt_algorithm, 
        Rparam, 
        fastjet::E_scheme, 
        fastjet::Best),
    jetDef_ca_recluster(
        fastjet::cambridge_algorithm, 
        1, 
        fastjet::E_scheme, 
        fastjet::Best),
    pc_subtractor(Rparam_in),
    /*
    pruner(
        fastjet::cambridge_algorithm, 
        prune_zcut, 
        Rcut_factor),
    filter(Rfilt, fastjet::SelectorNHardest(nfilt)),
    */
    bge(part_eta_max_in, .5), // particle eta max, grid spacing
    subtractor(&bge),
    rsd(beta_rsd, z_cut_rsd, Rparam_in),
    soft_killer(part_eta_max_in, .4) // rapidity max, grid_size
        {
        Rparam = Rparam_in;
        jet_eta_max = jet_eta_max_in;
        leading_jet_only = leading_jet_only_in;
        

        jet_eta_selector = fastjet::SelectorAbsRapMax(jet_eta_max_in);

        // Set up subtractor for jet by jet constituent subtraction
        subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
        subtractor.set_max_distance(0.3);  // R_max for ghost-particle pairing
        subtractor.set_max_eta(part_eta_max_in); // These two are for whole-event mode
        subtractor.initialize(); //
        }

    // Takes event particles. Finds jets with akt. Reclusters those jets using CA.
    // returns a vector of the delta of the first splitting for all the jets.
    vector<r_groomed_vars> get_deltas(vector<fastjet::PseudoJet> particles, bg_sub_options ops) {

        vector<r_groomed_vars> deltas;
        
        if (ops.event_sub == "ConSub") {
            bge.set_particles(particles);
            particles = subtractor.subtract_event(particles);
        }
        else if (ops.event_sub == "SoftKill") {
            vector<fastjet::PseudoJet> soft_killed_event;
            double pt_thresh = 0; // pt threshold for killed particles
            soft_killer.apply(particles, soft_killed_event, pt_thresh);
            particles = soft_killed_event;
        }
        else if (ops.event_sub != "null") cout << "Warning: Unknown event subtraction request" << endl;




        fastjet::ClusterSequence clust_seq(particles, jetDef_akt);
        vector<fastjet::PseudoJet> jets = sorted_by_pt(jet_eta_selector(clust_seq.inclusive_jets(ptMin)));


        for (const fastjet::PseudoJet& jet : jets) { 
           
            vector<fastjet::PseudoJet> constituents;
            
           if (ops.jet_sub == "MyPCkT") {
                constituents = pc_subtractor.kt_subtract_constit(jet, particles);
                if(debug && constituents.empty()) cout << "No remaining constituents after my pc subtraction" << endl;
            }
            else if (ops.jet_sub == "MyPCGM") {
                constituents = pc_subtractor.geometric_subtract_constit(jet, particles); 
                if(debug && constituents.empty()) cout << "No remaining constituents after my pc geometric match subtraction" << endl;
            }
            else if (ops.jet_sub =="PC_parts") {
                constituents = pc_subtractor.get_pc_particles(jet, particles);
                if(debug && constituents.empty()) cout << "No remaining constituents in perpendicular cone" << endl;
            }
            else if (ops.jet_sub == "RSD_mine") {
                fastjet::ClusterSequence reclusterSeq_temp(jet.constituents(), jetDef_ca_recluster);
                fastjet::PseudoJet reclustered_jet = sorted_by_pt(reclusterSeq_temp.inclusive_jets())[0];
                constituents = recursive_soft_drop_constit(reclustered_jet, Rparam);
            }
            else if (ops.jet_sub != "RSD_contrib") {
                cout << "Warning: Unknown jet subtraction request. No subtraction will be performed" << endl;
                constituents = jet.constituents();
            }
            else {
                if (debug) cout << "Calling jet.constituents()" << endl;
                constituents = jet.constituents();
            }

            if (constituents.empty()) {
                break; // Breaks are used because jets are in pt order, so once one has no remaining particles after subtraction, we've reacehed the jets that are all/mostly backgroun
            }
            
            if (debug) cout << "reclusterSeq with new jet constituents" << endl;
            fastjet::ClusterSequence reclusterSeq(constituents, jetDef_ca_recluster);
            jet = sorted_by_pt(reclusterSeq.inclusive_jets())[0];

            // The subtractions that work on an already CA reclustered jet
            if (ops.jet_sub = "RSD_contrib") {
                jet = rsd(jet);
            }

            fastjet::PseudoJet parent1, parent2;
            // For each jet, iteratively compare branchings
            while (reclustered_jet.has_parents(parent1, parent2)) { 
                // In each case, identify the higher pt parent 
                fastjet::PseudoJet harder  = (parent1.pt() > parent2.pt()) ? parent1 : parent2;
                fastjet::PseudoJet softer  = (parent1.pt() > parent2.pt()) ? parent2 : parent1;
                
                // rap-phi distance, delta = delta_ab
                double delta = harder.delta_R(softer);
                
                double z = softer.pt() / (harder.pt() + softer.pt());

                // Softdrop as described in https://arxiv.org/pdf/1402.2657
                if (z < z_cut * pow(delta/Rparam, beta)) { // does not meet soft drop condition
                    reclustered_jet = harder;
                    continue;
                }
                
                deltas.push_back(r_groomed_vars{.delta = delta, .jet_pt = jet.pt()});
                break;

                
            } // end while loop

            if (leading_jet_only) break;

        } // end jet loop 
        return deltas;
    } // end get_deltas() def
};


TString bool2Str(bool b) {return b? "True" : "False";}

int main() {
    
    TString output_folder_name = "";
    TString output_file_name = "r_groomed";
    //TString root_output_name = output_name;

    vector<back_sub_options> all_sub_ops = {
        back_sub_options{.name = "Pythia", .embed_in_bg = false},
        back_sub_options{.name = "Embedded", .embed_in_bg = true},
        back_sub_options{.name = "ConSub", .event_sub = "ConSub", .embed_in_bg = true},
        back_sub_options{.name = "SoftKiller", .event_sub = "SoftKill", .embed_in_bg = true},
        back_sub_options{.name = "My PC", .jet_sub = "MyPCGM", .embed_in_bg = true}
    };

    
    TString event_file_name = "event_test.root";
    vector<TString> background_file_names = {"backgrounds28k.root"};


    // For jet pt cuts of the lund plane
    vector<vector<double>> jet_pt_cuts = {
        {20, 120}, // first cuts is the inclusive plane
        {20, 40},
        {40, 60},
        {60, 80},
        {80, 100},
        {100, 120}
    };


    TString notes = ""; // leading jets only
    double part_eta_max = .9; // Particle eta max;
    double jet_eta_max = .5; // jet eta max;
    double Rparam = 0.4;
    if (part_eta_max - Rparam != jet_eta_max) cout << "WARNING: Rparam, part_eta_max, and jet_eta_max are not in agreement" << endl;
    double part_pt_min = .15;

    bool leading_jet_only = false;

    bool normalize = true; // normalize histograms
    
    TString space = " ";
    
    my_event_tree met(event_file_name);
    my_background_trees mbt(background_file_names);


    //=======================================================================
    // Declare 2D ROOT histograms
    double nBins = 20;
    double Rg_min = 0;
    double Rg_max = .5;
    double Rg_plot_length = (Rg_max - Rg_min) / nBins;

    vector<vector<TH1F*>> all_Rg_hists;

    for (int iOp = 0; iOp < all_sub_ops.size(); iOp++) {
        for (int iCut = 0; iCut < jet_pt_cuts.size(); iCut++) {
            all_Rg_hists[iOp].push_back(new TH1F(sub_ops.name + space + iCut, "R_{groomed}: " + jet_pt_cuts[iCut][0] " < Jet p_{T} < " jet_pt_cuts[iCut][1] +"; delta R; weighted counts", nBins, Rg_min, Rg_max));
        }       
    }

    for (back_sub_options ops : all_sub_ops) {
        ops.weighted_jet_pt_hist = new TH1F("weighted_jet_pt_hist" + space + ops.name, "Weighted Jet p_{T} for" + space + ops.name + "; jet p_{T}; Weighted counts", 400, 0, 200);
    }

    /*
    TH1F *true_hist = new TH1F("true_hist", "pythia; delta R; weighted counts",nBins, 0, .5);
    TH1F *bg_hist = new TH1F("bg_hist", "bg; delta R; weighted counts",nBins, 0, .5);
    TH1F *pc_hist = new TH1F("pc_hist", "my pc; delta R; weighted counts",nBins, 0, .5);
    // Filtering clusters a jet's constituents with a smaller-than-original radius R_filt, and keepts the n_filt hardest of the resulting subjets
    TH1F *filter_hist = new TH1F("filter_hist", "filtered; delta R; weighted counts",nBins, 0, .5);
    // Pruning reclusters a jet's constituents with some algorithm, vetoeing soft and large-angle recombinations 
    TH1F *prune_hist = new TH1F("prune_hist", "pruned; delta R; weighted counts",nBins, 0, .5);
    TH1F *cs_hist = new TH1F("cs_hist", "Constituent Subtractor; delta R; weighted counts",nBins, 0, .5);
    TH1F *softKill_hist = new TH1F("softKill_hist", "SoftKiller; delta R; weighted counts",nBins, 0, .5);
    */
    //TH1F *recursive_sd_hist = new TH1F("recursive_sd_hist", "Recursive SoftDrop; delta R; weighted counts",nBins, 0, .5);



    //TH1F *true_deltas_hist = new TH1F("true_deltas_hist", "z_pass plot, " + output_file_name + "; z_pass; weighted counts", 100, 0, .5);
    //TH1F *bg_deltas_hist = new TH1F("bg_deltas_hist", "z_pass plot, " + output_file_name + "; z_pass; weighted counts", 100, 0, .5);
    //TH1F *true_deltas_hist = new TH1F("true_deltas_hist", "jet constituent multiplicity, " + output_file_name + "; jet constituent multiplicity; weighted counts", 100, 0, 300);
    //TH1F *bg_deltas_hist = new TH1F("bg_deltas_hist", "jet constituent multiplicity, " + output_file_name + "; jet constituent multiplicity; weighted counts", 100, 0, 300);


    // ROOT hists for pt cuts
    /*
    vector<double> ptCutMins = {10, 20, 30, 50, 70, 100, 130, 160, 190};
    vector<double> ptCutMaxs = {20, 30, 50, 70, 100, 130, 160, 190, 250};
    vector<TH1F*> ptCuts;
    
    TString space = " ";
    for (int iCut = 0; iCut < ptCutMins.size(); ++iCut) {
        ptCuts.push_back(new TH1F("", "jet pT" + space + ptCutMins[iCut] + " to " + ptCutMaxs[iCut] + "; delta of first splitting; weighted counts",
            100, 0, .5));
    }
    */

    //=========================================================================
    // Set up fastjet analysis
    my_pc_subtractor pc_subtractor(Rparam);
    delta_getter dg(leading_jet_only, Rparam, jet_eta_max);
    rho_estimator rho_finder(Rparam, jet_eta_max);
    
    //=========================================================================
    // Event loop

    int numEvents = met.GetEntries();
    int numBackgrounds = mbt.GetEntries();
    
    int counter = 1000;

    for (int iEvent = 0; iEvent < numEvents; ++iEvent) {
        
        if (iEvent > counter) {
            cout << counter << " events analyzed from ROOT file" << endl;
            counter += 1000;
        }

        //met.GetEntry(iEvent); get_particles makes this same call

        vector<fastjet::PseudoJet> event_particles = met.get_particles(iEvent, part_eta_max);
        
        vector<fastjet::PseudoJet> all_particles = event_particles;

        // ============================================================
        // find first declusters for options that don't require a background
        
        for (int iOp = 0; iOp < all_sub_ops.size(); iOp++) {
            if (all_sub_ops[iOp].embed_in_bg) continue; 
            
            vector<r_groomed_vars> deltas = dg.get_deltas(event_particles, all_sub_ops[iOp]);
            for (r_groomed_vars vars : deltas) {
                
                for (int iCut = 0; iCut < jet_pt_cuts.size(); iCut++) { // Find which cuts to add this to
                    
                    if ((vars.jet_pt > jet_pt_cuts[iCut][0]) && (vars.jet_pt < jet_pt_cuts[iCut][1])) {
                        all_Rg_hists[iOp][iCut]->Fill(vars.delta, met.bin_weight);
                    }
                }
            }
        }

        // =====================================================================================================
        // Embedding in the background
         
        //mbt.GetEntry(iEvent); get_particles also makes GetEntry(iEvent) call
        vector<fastjet::PseudoJet> background_prtcls = mbt.get_particles(iEvent);
        
        move(background_prtcls.begin(), background_prtcls.end(), back_inserter(all_particles));

        // =======================================================================================================================
        // declusterings for options that do require a background
        
        for (int iOp = 0; iOp < all_sub_ops.size(); iOp++) {
            if (!all_sub_ops[iOp].embed_in_bg) continue; 
            
            vector<r_groomed_vars> deltas = dg.get_deltas(event_particles, all_sub_ops[iOp]);
            for (r_groomed_vars vars : deltas) {
                
                for (int iCut = 0; iCut < jet_pt_cuts.size(); iCut++) { // Find which cuts to add this to
                    
                    if ((vars.jet_pt > jet_pt_cuts[iCut][0]) && (vars.jet_pt < jet_pt_cuts[iCut][1])) {
                        all_Rg_hists[iOp][iCut]->Fill(vars.delta, met.bin_weight);
                    }
                }
            }
        }
        
    

    }  // End reading events from ROOT
    //==========================================================
    

    // =========================================================
    // Save these plots in a root file, so they can be modified/resized without having to run the whole thing again
    /*
    TFile output_root_file(root_output_name + ".root", "RECREATE");

    TH2F* lund_pointer = lund_inclusive;

    TTree lund_tree("lund_tree", "Lund Tree");
    lund_tree.Branch("lund", &lund_pointer);
 
    // First entry is inclusive lund
    lund_tree.Fill();
    // subsequent entries are cuts
    for (int iCut = 0; iCut < ptCutMins.size(); ++iCut) {
        lund_pointer = lundCuts[iCut];
        lund_tree.Fill();
    }

    lund_tree.Write();
    output_root_file.Close();
    */
    //=======================================================
    
    // Normalizing
    if (normalize) {
        for (int iOp = 0; iOp < all_sub_ops.size(); iOp++) {
            for (int iCut = 0; iCut < jet_pt_cuts.size(); iCut++) {
                double nJets_weighted = 0;
                for (int iBin  = all_sub_ops[iOp].weighted_jet_pt_hist->FindBin(jet_pt_cuts[iCut][0]);
                        iBin <= all_sub_ops[iOp].weighted_jet_pt_hist->FindBin(jet_pt_cuts[iCut][1]);
                        iBin ++) { // Add up all weighted jets in this bin
                    nJets_weighted += all_sub_ops[iOp].weighted_jet_pt_hist->GetBinContent(iBin);
                }
                all_Rg_hists[iOp][iCut]->Scale(1 / (Rg_plot_length * nJets_weighted))
                
                
            }
        }
    }


    //========================================================
    // Printing histograms

    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    
    // Adding a text information sheet with stats about this generation
    TString title = "Splitting Scale Analysis: " + output_file_name;
    TString line0 = notes;
    TString line1 = "Number of events: " + to_string(numEvents);
    TString line2 = "Events from: " + event_file_name;
    TString line3 = "Backgrounds from: " + print_vec(background_file_names);
    TString line4 = "Leading jet only: " + bool2Str(leading_jet_only);
    TString line5 = "Normalized: " + bool2str(normalize);

    TString line6 = "Jet radius: " + to_string(Rparam);
    TString line7 = "Particle eta max: " + to_string(part_eta_max);
    


    vector<TString> lines = {title, line0, line1, line2, line3, line4, line5, line6, line7};

    c1->cd();
    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }
    
    c1->Print(output_folder_name + output_file_name + ".pdf(","pdf");
    c1->Clear();
    

    // put hist in vector
    //vector<TH1F*> hists = {bg_hist, sd_hist, pc_hist, prune_hist, filter_hist, cs_hist, softKill_hist, true_hist};
    vector<TH1F*> hists = {bg_hist, pc_hist, cs_hist, softKill_hist, true_hist};
    
    double display_yMax = 0;

    for (int iCut = 0; iCut < jet_pt_cuts.size(); iCut++) {
        TLegend* legend = new TLegend(.1, .7, .28, .9);

        for (int iOp = 0; iOp < all_sub_ops.size(); iOp++) {

            display_yMax = max(display_yMax, all_Rg_hists[iOp][iCut]->GetMaximum());

            if (iOp == 0) all_Rg_hists[iOp][iCut]->Draw("HIST PLC");
            else all_Rg_hists[iOp][iCut]->Draw("HIST PLC SAME");

            legend->AddEntry(all_Rg_hists[iOp][iCut], all_sub_ops[iOp].name, "l");
        }

        all_Rg_hists[0][iCut]->SetStats(false); // no stat box
        all_Rg_hists[0][iCut]->GetYaxis()->SetRangeUser(0, display_yMax * 1.1);
        c1->Update();

        legend->Draw();
        if (iCut == jet_pt_cuts.size() - 1) c1->Print(output_folder_name + output_file_name + ".pdf)", "pdf");
        else c1->Print(output_folder_name + output_file_name + ".pdf", "pdf");


    }

    
    // Since ROOT histograms point to the data they store, once the eventRootFile
    // closes they no longer have access to their data
    met.close_file();
    mbt.close_file();

    return 0;
}


