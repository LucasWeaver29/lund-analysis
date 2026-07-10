// Make jet spectra (pt, y, phi) histograms
// Use events from root file and backgrounds from root file
// Attempt first background subtraction techniques


#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"

// ROOT - files
//#include <TFile.h>
//#include <TTree.h>
//#include <TBranch.h>
#include <vector>
//#include <iostream>

// ROOT - histogram
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

// ROOT - interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT - saving file.
#include "TFile.h"

#include "backSubTools.h"
#include "eventData.h"

#include <algorithm>

using namespace Pythia8;

TString bool2Str(bool b) {return b? "True" : "False";}

/*
double median(vector<double> v) {
    sort(v.begin(), v.end());
    if (v.size()%2 == 1) {
        return v[v.size()/2];
    }
    else {
        return .5 * (v[v.size()/2] + v[(v.size()/2)-1]);
    }
}
*/

int main(int argc, char* argv[]) {

    TString output_name = "jetSpectraBg - 9000, jet_eta_max, test eventData w Part";
    TString root_output_name = output_name;
    
    TString eventFileName = "events1000x9.root";
    TString backgroundFileName = "thermalBackgrounds9000.root";
    
    TString notes = "Added jet_eta_max = .5 ";
    // search for '// leading jet only' to change this

    // Unconditionally does background, maybe globalBgSub
    bool globalBgSub = false;
    bool my_perpCone = true;
    bool lilys_perpCone = true;
    bool leading_jet_only = true;
    
    double ptMin = 0;
    double part_eta_max = .9;
    double Rparam = .4;
    double jet_eta_max = part_eta_max - Rparam;


    //int toSkip = 1; // Only analyze 1/toSkip of the events in root file

    //Get bin specifications here from storeWithRoot
    // Make sure recreation of iBin later down also matches
    //vector<double> ptHatMins = {10, 20, 30, 50, 70, 90, 110, 130, 150};
    //vector<double> ptHatMaxs = {20, 30, 50, 70, 90, 110, 130, 150, 250};



    //==================================================================
    // Preparing to read from ROOT event file

    my_event_tree met(eventFileName);
    my_background_tree mbt(backgroundFileName);

    // ======================================================
    // Prepare ROOT analysis
    
    // Declare ROOT histograms
    TH1F *ptHist = new TH1F("ptHist", "Jet pT (anti-kt) - Weights; pT (GeV); weighted counts", 400, -70, 200);
    //TH1F *pt_noWeights = new TH1F("pt", "Jet pT (anti-kt) - No weights; pT (GeV); counts", 400, 0, 200);
    TH1F *etaHist = new TH1F("etaHist", "Jet pseudorapidity (anti-kt); pseudorapidity ; weighted counts", 100, -5.5, 5.5);
    TH1F *phiHist = new TH1F("phiHist", "Jet azimuthal angle (anti-kt); azimuthal angle; weighted counts", 100, 0, 2*M_PI);

    TH1F *pt_bgHist = new TH1F ("pt_bgHist", "Jet pT, Background; pT (GeV); weighted counts", 400, -70, 200);
    
    TH1F *pt_bg_subHist = new TH1F ("pt_bg_subHist", "Jet pT, Background, global background subtraction; pT (GeV); weighted counts", 400, -70, 200);
    
    // Make sure these are over the same range, with the same # of bins
    double pt_track_min = 0;
    double pt_track_max = 10;
    int pt_track_num_bins = 50;
    
    // For Lily's pc method
    TH1F *leadingJet_pt_tracks = new TH1F ("leadingJet_pt_tracksHist", "pT track spectrum of leading jet; pT track; counts", pt_track_num_bins, pt_track_min, pt_track_max);
    TH1F *pc_pt_tracks = new TH1F ("pc_pt_tracks", "Perpendicular cone pT tracks spectrum; pT track; counts", pt_track_num_bins, pt_track_min , pt_track_max);
    TH1F *jet_minus_pc_tracks = new TH1F ("jet_minus_pc_track", "jet track pt - pc track pt; pT track; counts", pt_track_num_bins, pt_track_min , pt_track_max);
    TH1F *jet_minus_Cpc_tracks = new TH1F ("jet_minus_Cpc_track", "jet track pt - pc track pt, with C correction; pT track; counts", pt_track_num_bins, pt_track_min , pt_track_max);
    

    TH1F *pt_pc_sub = new TH1F ("pt_pc_sub", "Jet pT, Background, perpendicular cone subtraction; pT (GeV); weighted counts", 400, -70, 200);
    
    // ROOT histogram for each bin
    /*
    vector<TH1F*> binHists = {};
    for (int i = 0; i < ptHatMins.size(); ++i) {
        TString num = i;
        TString name = "pT bin range: " + to_string(ptHatMins[i]) + "to" + to_string(ptHatMaxs[i]);
        binHists.push_back(new TH1F(num, name,400, 0, 200));
    }
    
    vector<TH1F*> binHists_noW = {};
    for (int i = 0; i < ptHatMins.size(); ++i) {
        TString num = i;
        TString name = "pT bin range: " + to_string(ptHatMins[i]) + "to" + to_string(ptHatMaxs[i]);
        binHists_noW.push_back(new TH1F(num, name,400, 0, 200));
    }
    */
    // ======================================================


    // Setup fastjet analysis
    fastjet::Strategy strat = fastjet::Best;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;

    fastjet::JetDefinition jetDef_akt(fastjet::antikt_algorithm, Rparam, recombScheme,strat);
    // For global background estimation
    fastjet::JetDefinition jetDef_kt(fastjet::kt_algorithm, Rparam, recombScheme,strat);

    // Area definition
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(part_eta_max));
     
    fastjet::Selector jet_eta_selector = fastjet::SelectorEtaRange(-1 * jet_eta_max, jet_eta_max);

    my_pc_subtractor pc_subtractor(Rparam, fastjet::cambridge_algorithm);
    RhoEstimator rho_finder(Rparam, .5, .9); // rparam, jet_eta_max, part_eta_max

    //=======================================
    // Reading events from ROOT

    int numEvents = met.GetEntries();
    
    int numBackgrounds = mbt.GetEntries();
    
    if (numBackgrounds < numEvents) {
        cout << "Warning: There are " << numEvents << " provided events but only" << numBackgrounds << " provided backgrounds." << endl;
        cout << "Jet Spectra creation will stop once backgrounds run out" << endl;
        numEvents = numBackgrounds;
    }


    //int iBin = 0;
    //double currentBin = ptHatMins[0];
    int counter = 1000;

    //cout << "Remember that we set iEvent = 15000 and had event loop stop after iEvent = 20000" << endl;
    //numEvents = 2000;
    // Event loop for events stored in ROOT file
    for (int iEvent = 0; iEvent < numEvents; ++iEvent) {

        if (iEvent > counter) {
            cout << counter  << " events analyzed from ROOT file" << endl;
            counter += 1000;
        }

        // met.GetEntry(iEvent); because get_particles also does getEntry

        /*
        if (ptHatMin>currentBin) {
            ++iBin;
            currentBin = ptHatMin;
        }
        */

        vector<fastjet::PseudoJet> pythia_particles = met.get_particles(iEvent, part_eta_max);

        //============================
        // stats for no background

        // Fastjet analysis
        fastjet::ClusterSequence pythia_clustSeq(pythia_particles, jetDef_akt);
        
        vector<fastjet::PseudoJet> pythia_jets = sorted_by_pt(jet_eta_selector(pythia_clustSeq.inclusive_jets(ptMin)));
        
        // Add each jet's stats to histograms
        for (fastjet::PseudoJet jet : pythia_jets) {            
            
            // ROOT histogram
            ptHist->Fill(jet.pt(), met.bin_weight);
            //pt_noWeights->Fill(jet.pt());
            etaHist->Fill(jet.eta(), met.bin_weight);
            phiHist->Fill(jet.phi(), met.bin_weight);
            
            //binHists[iBin]->Fill(jet.pt(), weight);
            //binHists_noW[iBin]->Fill(jet.pt());
            if (leading_jet_only) break; // leading jet only
        }
        
        // ==============================================
        // Embedding in the background

        //mbt.GetEntry(iEvent); // because get_particles also does GetEntry

        // Create a vector to store all new particles
        vector<fastjet::PseudoJet> background_prtcls = mbt.get_particles(iEvent);
        
        vector<fastjet::PseudoJet> all_particles = pythia_particles;

        // Append background particle vector to fjInputs
        move(background_prtcls.begin(), background_prtcls.end(), back_inserter(all_particles));

        // end embedding in background
        //=============================================
        // global background density estimation
        
        double rho;
        
        if (globalBgSub) {
            
            rho = rho_finder.rho(all_particles);
        }

        //=====================================================

        // Stats with background

        // Fastjet analysis
        fastjet::ClusterSequenceArea embedded_clustSeq(all_particles, jetDef_akt, areaDef);
        
        vector<fastjet::PseudoJet> embedded_jets = sorted_by_pt(jet_eta_selector(embedded_clustSeq.inclusive_jets(ptMin)));
        
        // Add each jet's stats to background histograms
        for (fastjet::PseudoJet jet : embedded_jets) {            

            // ROOT histogram
            pt_bgHist->Fill(jet.pt(), met.bin_weight);
            
            if (globalBgSub) {
                pt_bg_subHist->Fill(jet.pt() - rho*jet.area(), met.bin_weight);
            }

            if (leading_jet_only) break; // leading jet only

        } // End stats with background

        //================================================
        // (my) Perpendicular cone background estimation and constituent-based subtraction algorithm
        
        if (lilys_perpCone) {
            // Find spectrum of constituents in leading jet:
            fastjet::PseudoJet leadingJet = embedded_jets[0];
            for (fastjet::PseudoJet constituent : leadingJet.constituents()) {
                leadingJet_pt_tracks->Fill(constituent.pt(), met.bin_weight);
            }

            // Find pt spectrum of pc particles
            fastjet::PseudoJet pc_axis;
            pc_axis.reset_PtYPhiM(0, leadingJet.eta(), fmod(leadingJet.phi() + M_PI/2, 2*M_PI));

            double pc_pt_sum = 0; // Tracking total scalar momentum in pc

            fastjet::Selector pc_selector = fastjet::SelectorCircle(Rparam);
            pc_selector.set_reference(pc_axis);
            vector<fastjet::PseudoJet> pc_tracks = pc_selector(all_particles);
            for (fastjet::PseudoJet track :: pc_tracks) {
                pc_pt_tracks->Fill(track.pt(), met.bin_weight);
                pc_pt_sum += track.pt();
            }
        }

        // Try to use perpendicular cone to measure multiplicities of particles with certain pt, then match those particles by pt to constituents of the real jet
        if (my_perpCone) {   
            
            for (fastjet::PseudoJet jet : inclusiveJets_bg) {
                fastjet::PseudoJet subtracted_jet = pc_subtractor.subtract(jet, fjInputs);
                if (subtracted_jet.pt() != 0) pt_pc_subHist->Fill(subtracted_jet.pt(), met.bin_weight);
                if (leading_jet_only) break;
            }
            
        }

    } // end event loop

    
    
    
    // Applying bin by bin correction factor for Lily's pc method
    if (lilys_perpCone) {
        
        for (int iBin = 1; iBin < pt_track_num_bins + 1; iBin) {


        }
        
        
        
        TH1F *corrected_ptHist = new TH1F ("corrected_ptHist", "Perpendicular cone pT tracks spectrum", pt_track_num_bins, pt_track_min , pt_track_max);
        if(perpCone) {
            // bin 0 is underflow bin
            for (int iBin = 1; iBin < pt_track_num_bins + 1; ++iBin) {
                double N_jet = leadingJet_pt_tracks->GetBinContent(iBin);
                double N_pc = pc_pt_tracks->GetBinContent(iBin);
                double correction_factor = N_jet / N_pc;
                jet_minus_pc_tracks->Fill()
                //corrected_ptHist->Fill();
            }
        }
    }

    /*
    // Store these plots in a root file
    TFile output_root_file(root_output_name + ".root", "RECREATE"); 

    TTree spectra_tree("spectra_tree", "spectra_tree");
    spectra_tree.Branch("ptHist", &*ptHist);
    spectra_tree.Branch("pt_bgHist", &*pt_bgHist);
    spectra_tree.Branch("pt_bg_subHist", &*pt_bg_subHist);

    spectra_tree.Fill();
    spectra_tree.Write();
    output_root_file.Close();
    */

    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    

    // Adding a text information sheet with stats about this generation
    // Use Form() here, and just /n for new line
    TString title = "Jet Spectra from events in ROOT folder: " + output_name;
    TString line0 = notes;
    TString line1 = "Number of events: " + to_string((int)(numEvents));
    TString line2 = "Events from: " + eventFileName;
    
    TString line3 = "Embedded in background: True";
    TString line4 = "Backgrounds from: " + backgroundFileName;
    TString line5 = "Global Background Estimation + Subtraction: " + bool2Str(globalBgSub);
    TString line6 = "Perpendicular cone method: " + bool2Str(perpCone);
    TString line7 = "Leading jet only: " + bool2Str(leading_jet_only);

    TString line8 = "Jet pT min: " + to_string(ptMin);
    TString line9 = "Jet radius: " + to_string(Rparam);
    TString line10 = "Particle eta max: " + to_string(part_eta_max);
    TString line11 = "Jet eta max: " + to_string(jet_eta_max);

    vector<TString> lines = {title, line0, line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11};

    c1->cd();
    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }
    
    c1->Print(output_name+".pdf(","pdf");
    c1->Clear();


    TLegend* legend = new TLegend(.7, .7, .9, .9);
    
    
    c1->cd();
    
    c1->SetLogy(1);
    
    // pt, with background
    pt_bgHist->Draw("HIST" "PLC");
    legend->AddEntry(pt_bgHist, "Jet pt, (background)", "l");

    // pt, no background
    ptHist->Draw("HIST" "SAME" "PLC");
    legend->AddEntry(ptHist, "Jet pt, no background","l");
    
    // pt, with background after subtraction
    pt_bg_subHist->Draw("HIST" "SAME" "PLC");
    legend->AddEntry(pt_bg_subHist, Form("Jet pt (area-based sub)"), "l");


    pt_pc_subHist->Draw("HIST SAME PLC");
    legend->AddEntry(pt_pc_subHist, Form("Jet pt (pc sub)"), "l");


    legend->Draw();
    c1->Print(output_name + ".pdf","pdf");
    c1->Clear();
    
    
    /*
// Looping through histograms for each weighted bin
    TLegend* legend = new TLegend(.1, .7, .28, .9);
    legend->AddEntry(pt, "Total", "l");
    for (int i = 0; i < binHists.size(); ++i) {
        binHists[i]->Draw("HIST" "SAME" "PLC");
        TString caption = "" + to_string((int)ptHatMins[i]) + " to " + to_string((int)ptHatMaxs[i]);
        legend->AddEntry(binHists[i], caption, "l");
    }
    legend->Draw();
    c1->Print(output_name + "(","pdf");
    c1->Clear();
    

    // Unweighted pt
    c1->SetLogy(1);
    c1->cd();
    pt_noWeights->Draw("HIST" "PLC");

    // Looping through histograms for each unweighted bin
    TLegend* legend2 = new TLegend(.1, .7, .28, .9);
    legend2->AddEntry(pt, "Total", "l");
    for (int i = 0; i < binHists_noW.size(); ++i) {
        binHists_noW[i]->Draw("HIST" "SAME" "PLC");
        TString caption = "" + to_string((int)ptHatMins[i]) + " to " + to_string((int)ptHatMaxs[i]);
        legend2->AddEntry(binHists[i], caption, "l");
    }
    legend2->Draw();
    c1->Print(output_name, "pdf");
    c1->Clear();
    */
//

    c1->SetLogy(0);
    c1->cd();
    etaHist->Draw("HIST");
    c1->Print(output_name + ".pdf","pdf");
    c1->Clear();

    c1->cd();
    phiHist->Draw("HIST");
    c1->Print(output_name+".pdf)","pdf");
    c1->Clear();
    
    met.close_file();
    mbt.close_file();
    
    return 0;
}