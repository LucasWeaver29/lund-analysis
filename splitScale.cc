
#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
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

class delta_getter {

public:
    
    // General
    double ptMin = 60;
    double ptMax = 80;
    double jet_eta_max = .5;
    double Rparam = .4;
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
    double Rfilt = .3;
    double nfilt = 3;
    fastjet::Filter filter;

    // Pruning
    double prune_zcut = .2;
    double Rcut_factor = .2;
    fastjet::Pruner pruner;

    // Jet by jet constituent subtraction
    fastjet::GridMedianBackgroundEstimator bge; // GridMedianBackgroundEstimator is faster than JetMedianBackgroundEstimator and "performs equally well in nearly all cases"
    fastjet::contrib::ConstituentSubtractor subtractor;

    // SoftKiller
    fastjet::contrib::SoftKiller soft_killer;

    
    /*  // RecursiveSoftDrop - fastjet contrib
    double z_cut_rsd = .2;
    double beta_rsd = 0;
    double Rparam_rsd = .4;
    fastjet::contrib::RecursiveSoftDrop rsd(beta_rsd, z_cut_rsd, Rparam_rsd);
    */

    // constructor
    delta_getter(double ptMin_in, double ptMax_in, bool leading_jet_only_in, double Rparam_in = .4, double jet_eta_max_in = .5, double part_eta_max_in = .9):
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
    pruner(
        fastjet::cambridge_algorithm, 
        prune_zcut, 
        Rcut_factor),
    filter(Rfilt, fastjet::SelectorNHardest(nfilt)),
    bge(part_eta_max_in, .5), // particle eta max, grid spacing
    subtractor(&bge),
    soft_killer(part_eta_max_in, .4) // rapidity max, grid_size
        {
        Rparam = Rparam_in;
        jet_eta_max = jet_eta_max_in;
        leading_jet_only = leading_jet_only_in;
        ptMin = ptMin_in;
        ptMax = ptMax_in;

        jet_eta_selector = fastjet::SelectorAbsRapMax(jet_eta_max_in);

        // Set up subtractor for jet by jet constituent subtraction
        subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
        subtractor.set_max_distance(0.3);  // R_max for ghost-particle pairing
        subtractor.set_max_eta(part_eta_max_in); // These two are for whole-event mode
        subtractor.initialize(); //
        }

    // Takes event particles. Finds jets with akt. Reclusters those jets using CA.
    // returns a vector of the delta of the first splitting for all the jets.
    vector<double> get_deltas(vector<fastjet::PseudoJet> particles, bg_sub_options subs = bg_sub_options{}) {

        vector<double> deltas;
        
        if (subs.jet_by_jet_constit) {
            bge.set_particles(particles);
            particles = subtractor.subtract_event(particles);
        }
        if (subs.soft_kill) {
            vector<fastjet::PseudoJet> soft_killed_event;
            double pt_thresh = 0; // pt threshold for killed particles
            soft_killer.apply(particles, soft_killed_event, pt_thresh);
            particles = soft_killed_event;
        }

        fastjet::ClusterSequence clust_seq(particles, jetDef_akt);

        vector<fastjet::PseudoJet> jets = sorted_by_pt(jet_eta_selector(clust_seq.inclusive_jets(ptMin)));

        for (const fastjet::PseudoJet& jet : jets) { 
           
            vector<fastjet::PseudoJet> constituents;
            
            if (subs.my_pc) {
                constituents = pc_subtractor.subtract_constit(jet, particles);
                if (constituents.empty()) continue;
            }
            else if (subs.filter) {
                constituents = filter(jet).constituents();
            }
            else if (subs.prune) {
                constituents = pruner(jet).constituents();
            }
            else if (subs.recursive_soft_drop) {
                constituents = recursive_soft_drop_constit(jet);
            }
            else constituents = jet.constituents();

            fastjet::ClusterSequence reclusterSeq(constituents, jetDef_ca_recluster);
            fastjet::PseudoJet reclustered_jet = sorted_by_pt(reclusterSeq.inclusive_jets())[0];
            //fastjet::PseudoJet reclustered_jet = jet;

            if ((reclustered_jet.pt() < ptMin) || (reclustered_jet.pt() > ptMax)) continue; // jets are in pt order


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
                
                deltas.push_back(delta);
                break;

                
            } // end while loop

            if (leading_jet_only) break;

        } // end jet loop 
        return deltas;
    } // end get_deltas() def
};


TString bool2Str(bool b) {return b? "True" : "False";}

int main() {
    
    TString output_folder_name = "splitScale9000/";
    TString output_file_name = "compare";
    //TString root_output_name = output_name;

    TString notes = ""; // leading jets only
    double ptMin = 60; // Minimum jet pT
    double ptMax = 80;
    double part_eta_max = .9; // Particle eta max;
    double jet_eta_max = .5; // jet eta max;
    
    double Rparam = 0.4;

    bool leading_jet_only = false;

    bool normalize = true; // normalize histograms
    
    TString eventFileName = "events1000x9.root";
    TString backgroundFileName = "thermalBackgrounds9000.root";

    my_event_tree met(eventFileName);
    my_background_tree mbt(backgroundFileName);


    //=======================================================================
    // Declare 2D ROOT histograms
    double nBins = 20;


    TH1F *true_hist = new TH1F("true_hist", "pythia; delta R; weighted counts",nBins, 0, .5);
    TH1F *bg_hist = new TH1F("bg_hist", "bg; delta R; weighted counts",nBins, 0, .5);
    TH1F *pc_hist = new TH1F("pc_hist", "my pc; delta R; weighted counts",nBins, 0, .5);
    // Filtering clusters a jet's constituents with a smaller-than-original radius R_filt, and keepts the n_filt hardest of the resulting subjets
    TH1F *filter_hist = new TH1F("filter_hist", "filtered; delta R; weighted counts",nBins, 0, .5);
    // Pruning reclusters a jet's constituents with some algorithm, vetoeing soft and large-angle recombinations 
    TH1F *prune_hist = new TH1F("prune_hist", "pruned; delta R; weighted counts",nBins, 0, .5);
    TH1F *cs_hist = new TH1F("cs_hist", "Constituent Subtractor; delta R; weighted counts",nBins, 0, .5);
    TH1F *softKill_hist = new TH1F("softKill_hist", "SoftKiller; delta R; weighted counts",nBins, 0, .5);
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
    delta_getter dg(ptMin, ptMax, leading_jet_only, Rparam, jet_eta_max);
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
        // find first decluster of true jets
        
        vector<double> true_deltas = dg.get_deltas(event_particles);

        for (double delta : true_deltas) {
            true_hist->Fill(delta, met.weight);
        }

        // =====================================================================================================
        // Embedding in the background
         
        //mbt.GetEntry(iEvent); get_particles also makes GetEntry(iEvent) call
        vector<fastjet::PseudoJet> background_prtcls = mbt.get_particles(iEvent);
        
        move(background_prtcls.begin(), background_prtcls.end(), back_inserter(all_particles));

        // =======================================================================================================================
        // Decluster embedded jets. No softdrop, pc sub, etc
        vector<double> bg_deltas = dg.get_deltas(all_particles);
        for (double bg_delta : bg_deltas) {
            bg_hist->Fill(bg_delta, met.weight);
        }

        // Now with my perpendicular cone method // Currenlt with ca recluster at r = .4
        vector<double> bg_pc_deltas = dg.get_deltas(all_particles, bg_sub_options{.my_pc = true}); 
        for (double bg_pc_delta : bg_pc_deltas) {
            pc_hist->Fill(bg_pc_delta, met.weight);
        }

        /*
        // Filtering
        vector<double> filter_deltas = dg.get_deltas(all_particles, bg_sub_options{.filter=true});
        for (double delta : filter_deltas) {
            filter_hist->Fill(delta, met.weight);
        }

        // Pruning
        vector<double> prune_deltas = dg.get_deltas(all_particles, bg_sub_options{.prune=true});
        for (double delta : prune_deltas) {
            prune_hist->Fill(delta, met.weight);
        }
        */

        // jet by jet constituent subtraction
        vector<double> cs_deltas = dg.get_deltas(all_particles, bg_sub_options{.jet_by_jet_constit=true});
        for (double delta : cs_deltas) {
            cs_hist->Fill(delta, met.weight);
        }

        vector<double> softKill_deltas = dg.get_deltas(all_particles, bg_sub_options{.soft_kill=true});
        for (double delta : softKill_deltas) {
            softKill_hist->Fill(delta, met.weight);
        }

        /*
        vector<double> rsd_deltas = dg.get_deltas(all_particles, bg_sub_options{.recursive_soft_drop=true});
        for (double delta : rsd_deltas) {
            recursive_sd_hist->Fill(delta, met.weight);
        }
        */

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
    //========================================================
    // Printing histograms

    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    
    // Adding a text information sheet with stats about this generation
    TString title = "Splitting Scale Analysis: " + output_file_name;
    TString line0 = notes;
    TString line1 = "Number of events: " + to_string(numEvents);
    TString line2 = "Events from: " + eventFileName;
    TString line3 = "Backgrounds from: " + backgroundFileName;
    TString line4 = "Leading jet only: " + bool2Str(leading_jet_only);
    TString line5 = "Normalized: " + bool2str(normalize);

    TString line6 = "Jet pT min: " + to_string(ptMin);
    TString line7 = "Jet pT max: " + to_string(ptMax);
    TString line8 = "Jet radius: " + to_string(Rparam);
    TString line9 = "Particle eta max: " + to_string(part_eta_max);
    


    vector<TString> lines = {title, line0, line1, line2, line3, line4, line5, line6, line7, line8, line9};

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
    TLegend* legend = new TLegend(.1, .7, .28, .9);
    
    c1->SetLogy(0);

    for (int iHist = 0; iHist < hists.size(); ++iHist) {
        cout << hists[iHist]->GetEntries() << endl;
        if (normalize) hists[iHist]->Scale(1.0/hists[iHist]->Integral());
        c1->cd();
        if (iHist == 0) hists[iHist]->Draw("HIST PLC");
        else hists[iHist]->Draw("HIST SAME PLC");
        legend->AddEntry(hists[iHist], hists[iHist]->GetTitle(), "l");
    }

    legend->Draw();

    c1->Print(output_folder_name + output_file_name + ".pdf)", "pdf");

    
    // Since ROOT histograms point to the data they store, once the eventRootFile
    // closes they no longer have access to their data
    met.close_file();
    mbt.close_file();

    return 0;
}



