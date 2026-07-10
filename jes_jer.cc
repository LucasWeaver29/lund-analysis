#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/Selector.hh"

// ROOT - files
#include <vector>
#include <iostream>

// ROOT - histogram
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"

#include "backSubTools.h"
#include "eventData.h"

#include <chrono>

using namespace Pythia8;


int main() {

    TString output_file_name = "45k, jet_eta_max";
    TString output_folder_name = "jetEnergyPlots/";

    TString event_file_name = "events5000x9.root";
    TString background_file_name = "thermalBackgrounds45000.root";

    TString notes = "Added jet_eta_max = .5";

    //cout << "Remember that this is pt sub - pt true, NOT divided by pt true" << endl;

    // Create a root output file
    TString root_output_name = output_file_name;
    
    double ptMin = 0; // jet pt min
    double partEtaMax = .9;
    double jet_eta_max = .5; // = partEtaMax - Rparam

    bool leading_jet_only = true;

    // Graphs 'n stuff
    double nBinsX = 400;
    double xMin = 0;
    double xMax = 200;
    TH2F* sub_vs_true_hist = new TH2F("sub_vs_true_hist", "Area-based Subtraction: pT sub vs pT True; pt True; (pt_sub - pt_true) / pt_true", nBinsX, xMin, xMax, 400, -20, 20);
    TH2F* sub_minus_true_hist = new TH2F("sub_minus_true_hist", "Area-based Subtraction: pT sub minus pT True; pt True; pt_sub - pt_true", nBinsX, xMin, xMax, 400, -100, 100);

    TH2F* fastjet_sub_vs_true_hist = new TH2F("fastjet_sub_vs_true_hist", "Area-based Subtraction, with fastjet's methods: pT sub vs pT True; pt True; (pt_sub - pt_true)/pt_true", nBinsX, xMin, xMax, 400, -20, 20);
    TH2F* fastjet_sub_minus_true_hist = new TH2F("fastjet_sub_minus_true_hist", "Area-based Subtraction, with fastjet's methods: pT sub minus pT True; pt True; pt_sub - pt_true", nBinsX, xMin, xMax, 400, -100, 100);

    /*
    TGraph* jes = new TGraph();
    jes->SetTitle("Jet Energy Scale vs pt True; pt True; Jet energy scale");
    TGraph* jer = new TGraph();
    jer->SetTitle("Jet Energy Resolution vs pt True; pt True; Jet energy resolution");
    */

    TH1F* delta_pt = new TH1F("delta_pt", "Delta pT; pT_raw - pT_sub; counts", 50, -10, 150);
    TH1F* my_rho_hist = new TH1F("my_rho_hist", "Rho from my method; rho; counts", 100, -1, 250);
    TH1F* fastjet_rho_hist = new TH1F("fastjet_rho_hist", "Rho from fastjet; rho; counts", 100, -50, 250);
    TH1F* akt_area_hist = new TH1F("kt_area_hist", "Area from akt jets; area' counts", 100, -.1, 1);
   

    //==================================================================
    // Preparing to read from ROOT files

    my_event_tree met(event_file_name);
    my_background_tree mbt(background_file_name);

    //=======================================================================

    // Setup fastjet analysis
    fastjet::Strategy strat = fastjet::Best;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
    double Rparam = .4;

    fastjet::JetDefinition jetDef_akt(fastjet::antikt_algorithm, Rparam, recombScheme,strat);
    // For global background estimation
    fastjet::JetDefinition jetDef_kt(fastjet::kt_algorithm, Rparam, recombScheme,strat);

    // Area definition
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(partEtaMax + Rparam));
      
    fastjet::Selector jet_eta_selector = fastjet::SelectorAbsRapMax(.5);
    fastjet::Selector part_eta_selector = fastjet::SelectorAbsRapMax(.9);


    fastjet::JetMedianBackgroundEstimator bge(jet_eta_selector, jetDef_kt, areaDef);
    fastjet::Subtractor subtractor(&bge);

    RhoEstimator rho_finder(Rparam, .5, .9); // rparam, jet_eta_max, part_eta_max


    // For diferentiating background particles
    //MyInfo* bg_true = new MyInfo(true); // Attatcht this user info to bakcground particles
    //MyInfo* bg_false = new MyInfo(false); // to not background particles
    //=====================================================================
    // Reading events from ROOT

    int numEvents = met.GetEntries();

    int numBackgrounds = mbt.GetEntries();

    if (numBackgrounds < numEvents) {
        cout << "Warning: There are " << numEvents << " provided events but only" << numBackgrounds << " provided backgrounds." << endl;
        cout << "Jet Spectra creation will stop once backgrounds run out" << endl;
        numEvents = numBackgrounds;
    }

    int counter = 1000;

    //std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now(); // for timing first iteraction of event loop
    //std::chrono::time_point<std::chrono::high_resolution_clock> end;
    //cout << "Remember that iEvent set to 10000, numEvents t0 15000 to speed things up" << endl;
    //numEvents = 27000;
    // Event loop for events from ROOT file
    for (int iEvent = 0; iEvent < numEvents; ++iEvent) {


        if (iEvent > counter) {
            cout << counter  << " events analyzed from ROOT file" << endl;
            /* For timing 
            end = std::chrono::high_resolution_clock::now();
            auto duration  = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            cout << counter  << " events analyzed from ROOT file in " << duration.count() << " ms" << endl;
            //start = std::chrono::high_resolution_clock::now();
            */
            counter += 1000;
        }

        // met.GetEntry(iEvent); // get_particles() does this too.

        vector<fastjet::PseudoJet> event_particles = part_eta_selector(met.get_particles(iEvent));


        // Finding truth level jets
        fastjet::ClusterSequence clustSeq_true(event_particles, jetDef_akt);
        // What should pt min of truth level jets?
        vector<fastjet::PseudoJet> jets_true = sorted_by_pt(jet_eta_selector(clustSeq_true.inclusive_jets(ptMin)));

        // make a copy
        vector<fastjet::PseudoJet> all_particles = event_particles;

        // ==============================================
        // Embedding in the background

        //mbt.GetEntry(iEvent);

        // Create a vector to store all new particles
        vector<fastjet::PseudoJet> background_prtcls = mbt.get_particles(iEvent); // no eta selection required because bg particles are only up to .9 eta
        
        // Append background particle vector to fjInputs_bg
        move(background_prtcls.begin(), background_prtcls.end(), back_inserter(all_particles));

        // end embedding in background
        
        //=============================================
        // global background density estimation
        // my method
        /*
        fastjet::ClusterSequenceArea clustSeq_kt(all_particles, jetDef_kt, areaDef);

        vector<fastjet::PseudoJet> kt_jets_bg = sorted_by_pt(jet_eta_selector(clustSeq_kt.inclusive_jets(0)));

        // remove the two leading jets
        kt_jets_bg.erase(kt_jets_bg.begin());
        kt_jets_bg.erase(kt_jets_bg.begin());

        // My method for rho
        vector<double> rhos;

        for (fastjet::PseudoJet jet : kt_jets_bg) {
            rhos.push_back(jet.pt()/jet.area());
        }

        double my_rho = median(rhos);
        */

        double my_rho = rho_finder.rho(all_particles);
        my_rho_hist->Fill(my_rho);

        // fastjet methods for global background density estimation
        bge.set_particles(all_particles);
        double fastjet_rho = bge.rho();
        fastjet_rho_hist->Fill(fastjet_rho);

        //=====================================================
        // Stats with background

        // Fastjet analysis
        fastjet::ClusterSequenceArea clustSeq_bg(all_particles, jetDef_akt, areaDef);

        vector<fastjet::PseudoJet> jets_background = jet_eta_selector(clustSeq_bg.inclusive_jets(0));

        // Need to match jets geometrically between akt_jets_bg, jets_true. Done simply by rap phi distance. 
        fastjet::Selector sel = fastjet::SelectorCircle(.4);// The maximum R two jets can be from each other to be matched

        if (jets_true.empty()) continue;
        //fastjet::PseudoJet jet_true = jets_true[0];
        for (fastjet::PseudoJet jet_true : jets_true) {
            
            if (abs(jet_true.eta()) > jet_eta_max) continue;

            sel.set_reference(jet_true);
            vector<fastjet::PseudoJet> potential_matches = sel(jets_background);

            double delta = .4;
            double temp_delta;
            fastjet::PseudoJet closest_match;
            for (fastjet::PseudoJet potential_match : potential_matches) {

                temp_delta = potential_match.delta_R(jet_true);
                if (temp_delta < delta) {
                    delta = temp_delta;
                    closest_match = potential_match;
                }
                
            }

            if (closest_match.pz() != 0) { // if a match was found ...
                //cout << "Found valid match for leading jet at delta = " << delta << endl;

                double pt_true = jet_true.pt();

                double pt_sub = closest_match.pt() - my_rho*closest_match.area();

                akt_area_hist->Fill(closest_match.area());

                sub_vs_true_hist->Fill(pt_true, (pt_sub - pt_true)/pt_true, met.weight);
                //sub_vs_true_hist->Fill(pt_true, (pt_sub - pt_true), weight);

                sub_minus_true_hist->Fill(pt_true, pt_sub - pt_true, met.weight);
                //sub_minus_true_unweighted_hist->Fill(pt_sub - pt_true);

                delta_pt->Fill(closest_match.pt() - pt_sub);


                // fastjet method
                fastjet::PseudoJet subtracted_jet = subtractor(closest_match);
                fastjet_sub_vs_true_hist->Fill(pt_true, (subtracted_jet.pt() - pt_true) / pt_true, met.weight);
                fastjet_sub_minus_true_hist->Fill(pt_true, subtracted_jet.pt() - pt_true, met.weight);


            }
            
            if (leading_jet_only) break; // for leading jet only

        } // end jets_true loop
        
//       */
    } // end event loop


    // Loop through 2D hist to extract scale and resolution
    /*
    double binWidth = (xMax - xMin) / nBinsX;
    for (int iBin = 1; iBin <= sub_vs_true_hist->GetNbinsX(); ++iBin) {
        // Use weight for jes, jer vs pt true? -> Don't think so, because weight was accounted for in creation of 2D sub_vs_true
        TH1D* slice = sub_vs_true_hist->ProjectionY("slice", iBin, iBin);
        //if (slice->IsEmpty()) continue;
        // plots point in the middle of the bin
        jes->AddPoint(iBin*binWidth + xMin + binWidth*.5, slice->GetMean());
        jer->AddPoint(iBin*binWidth + xMin + binWidth*.5, slice->GetStdDev());
    }
    */

    // Store these plots in a root file
    TFile output_root_file(output_file_name + root_output_name + ".root", "RECREATE");

    TTree jet_energy_tree("jet_energy_tree", "Output Tree");
    jet_energy_tree.Branch("2dhist", &*sub_vs_true_hist);
    //jet_energy_tree.Branch("jes", &*jes);
    //jet_energy_tree.Branch("jer", &*jer);

    jet_energy_tree.Fill();
    jet_energy_tree.Write();

    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);


    // Adding a text information sheet with stats about this generation
    TString title = "JES, JER plots: " + output_file_name;
    TString line0 = notes;
    TString line1 = "Number of events: " + to_string(numEvents);
    TString line2 = "Events from: " + event_file_name;
    TString line3 = "Backgrounds from: " + background_file_name;


    TString line4 = "Jet pT min: " + to_string(ptMin);
    TString line5 = "Jet radius: " + to_string(Rparam);
    TString line6 = "Particle eta max: " + to_string(partEtaMax);
    TString line7 = "Jet eta max: " + to_string(jet_eta_max);
    TString line8 = "Leading jet only: " + to_string(leading_jet_only);

    vector<TString> lines = {title, notes, line1, line2, line3, line4, line5, line6, line7, line8};

    c1->cd();
    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }
    
    c1->Print(output_folder_name + output_file_name+".pdf(","pdf");
    c1->Clear();



    vector<TH2F*> hists = {sub_vs_true_hist, sub_minus_true_hist, fastjet_sub_vs_true_hist, fastjet_sub_minus_true_hist};
    
    TString space = " ";

    double binWidth = (xMax - xMin) / nBinsX;
    for (TH2F* hist : hists) {
        
        TGraph* jes = new TGraph();
        jes->SetTitle("Jet Energy Scale vs pt True:" + space + hist->GetTitle() + "; pt True; Jet energy scale");
        TGraph* jer = new TGraph();
        jer->SetTitle("Jet Energy Resolution vs pt True" + space + hist->GetTitle() + "; pt True; Jet energy resolution");
        
        for (int iBin = 1; iBin <= hist->GetNbinsX(); ++iBin) {
            // Use weight for jes, jer vs pt true? -> Don't think so, because weight was accounted for in creation of 2D sub_vs_true
            TH1D* slice = hist->ProjectionY("slice", iBin, iBin);
            //if (slice->IsEmpty()) continue;
            // plots point in the middle of the bin
            jes->AddPoint(iBin*binWidth + xMin + binWidth*.5, slice->GetMean());
            jer->AddPoint(iBin*binWidth + xMin + binWidth*.5, slice->GetStdDev());
        }

        c1->cd();
        hist->Draw("COLZ");
        c1->Print(output_folder_name + output_file_name+".pdf","pdf");
        c1->Clear();

        c1->cd();
        jes->Draw();
        c1->Print(output_folder_name + output_file_name+".pdf","pdf");
        c1->Clear();

        c1->cd();
        jer->Draw();
        c1->Print(output_folder_name + output_file_name+".pdf","pdf");
        c1->Clear();

    }


    // Loop through 2D hist to extract scale and resolution
    /*
    double binWidth = (xMax - xMin) / nBinsX;
    for (int iBin = 1; iBin <= sub_vs_true_hist->GetNbinsX(); ++iBin) {
        // Use weight for jes, jer vs pt true? -> Don't think so, because weight was accounted for in creation of 2D sub_vs_true
        TH1D* slice = sub_vs_true_hist->ProjectionY("slice", iBin, iBin);
        //if (slice->IsEmpty()) continue;
        // plots point in the middle of the bin
        jes->AddPoint(iBin*binWidth + xMin + binWidth*.5, slice->GetMean());
        jer->AddPoint(iBin*binWidth + xMin + binWidth*.5, slice->GetStdDev());
    }
    */


    /*
    c1->cd();
    sub_vs_true_hist->Draw("COLZ");
    c1->Print(output_file_name+".pdf","pdf");
    c1->Clear();

    c1->cd();
    jes->Draw();
    c1->Print(output_file_name+".pdf","pdf");
    c1->Clear();

    c1->cd();
    jer->Draw();
    c1->Print(output_file_name+".pdf","pdf");
    c1->Clear();
    */

    c1->cd();
    delta_pt->Draw();
    c1->Print(output_folder_name + output_file_name+".pdf", "pdf");
    c1->Clear();

    c1->cd();
    my_rho_hist->Draw();
    c1->Print(output_folder_name + output_file_name+".pdf", "pdf");
    c1->Clear();

    c1->cd();
    fastjet_rho_hist->Draw();
    c1->Print(output_folder_name + output_file_name+".pdf", "pdf");
    c1->Clear();

    c1->cd();
    akt_area_hist->Draw();
    c1->Print(output_folder_name + output_file_name+".pdf)", "pdf");
    c1->Clear();

    

    met.close_file();
    mbt.close_file();
    output_root_file.Close();

    return 0;
}