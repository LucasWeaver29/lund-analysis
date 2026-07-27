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

#include "unifiedSubtractors.h"
#include "eventData.h"

#include <chrono>

using namespace std;


int main() {

    TString output_file_name = "jes_jer";
    TString output_folder_name = "";

    //TString event_file_name = "event_test.root";
    //vector<TString> background_file_names = {"backgrounds14ka.root", "backgrounds14kb.root"};

    TString event_file_name = "event_Zoltans.root";

    vector<TString> background_file_names = {        
        "backgrounds2.7m.a.root", 
        "backgrounds2.7m.b.root",
        "backgrounds2.7m.c.root",
        "backgrounds2.7m.d.root",
        "backgrounds2.7m.e.root",
        "backgrounds2.7m.f.root" 
    };

    vector<unified_sub_options> all_sub_ops = {
        //unified_sub_options{.name = "Pythia", .embed_in_bg = false},
        unified_sub_options{.name = "Embedded, area sub", .jet_sub = jet_subtraction::Area_pT, .embed_in_bg = true},
        //unified_sub_options{.name = "Embedded", .embed_in_bg = true},
        unified_sub_options{.name = "ConSub", .event_sub = event_subtraction::ConSub, .embed_in_bg = true},
        unified_sub_options{.name = "SoftKiller", .event_sub = event_subtraction::SoftKill, .embed_in_bg = true},
        //unified_sub_options{.name = "My PC", .jet_sub = jet_subtraction::MyPCGM, .embed_in_bg = true},
        unified_sub_options{.name = "RSD (contrib)", .jet_sub = jet_subtraction::RSD_contrib, .embed_in_bg = true}
    
    };


    TString notes = "";


    double Rparam = .5;
    double part_eta_max = .9;
    double jet_eta_max = .5; // = part_eta_max - Rparam
    double part_pt_min = .15;

    bool leading_jet_only = false;

    bool debug = false;

    // Graphs 'n stuff
    int jes_nBinsX = 100;
    double jes_xMin = 0;
    double jes_xMax = 200;
    double jes_xbinWidth = (jes_xMax - jes_xMin) / jes_nBinsX;

    int jes_nBinsY = 50;
    double jes_yMin = -10;
    double jes_yMax = 10;

    TString space = " ";

    vector<TH2F*> jes_plots;

    for (int iOp = 0; iOp < all_sub_ops.size(); iOp++) {
        jes_plots.push_back(new TH2F("jes plot" + space + iOp, "JES for " + all_sub_ops[iOp].name + "; pT_{true}; pT_{sub} - pT_{true} / pT_{true}", jes_nBinsX, jes_xMin, jes_xMax, jes_nBinsY, jes_yMin, jes_yMax));
        jes_plots[iOp]->SetStats(false);

    }

    /*
    TGraph* jes = new TGraph();
    jes->SetTitle("Jet Energy Scale vs pt True; pt True; Jet energy scale");
    TGraph* jer = new TGraph();
    jer->SetTitle("Jet Energy Resolution vs pt True; pt True; Jet energy resolution");
    */

    /*
    TH1F* delta_pt = new TH1F("delta_pt", "Delta pT; pT_raw - pT_sub; counts", 50, -10, 150);
    TH1F* my_rho_hist = new TH1F("my_rho_hist", "Rho from my method; rho; counts", 100, -1, 250);
    TH1F* fastjet_rho_hist = new TH1F("fastjet_rho_hist", "Rho from fastjet; rho; counts", 100, -50, 250);
    TH1F* akt_area_hist = new TH1F("kt_area_hist", "Area from akt jets; area' counts", 100, -.1, 1);
    */

    //==================================================================
    // Preparing to read from ROOT files

    my_event_tree met(event_file_name);
    my_background_trees mbt(background_file_names);

    //=======================================================================
    // Setup fastjet analysis
    double rho, jet_area;

    EventSubtractor event_subtractor(part_eta_max);
    JetSubtractor jet_subtractor(Rparam, part_eta_max);

    fastjet::GridMedianBackgroundEstimator bge(part_eta_max, .5);

    std::unique_ptr<fastjet::ClusterSequence> cs_ptr;
    std::unique_ptr<fastjet::ClusterSequenceArea> csa_ptr;

    fastjet::Selector jet_eta_selector = fastjet::SelectorAbsRapMax(jet_eta_max);

    // Need to match jets geometrically between akt_jets_bg, jets_true. Done simply by rap phi distance. 
    fastjet::Selector sel = fastjet::SelectorCircle(.4);// The maximum R two jets can be from each other to be matched

    fastjet::JetDefinition jet_def_akt(fastjet::antikt_algorithm, Rparam, fastjet::E_scheme,fastjet::Best);
    fastjet::AreaDefinition area_def(fastjet::active_area, fastjet::GhostedAreaSpec(part_eta_max));


    /*

    fastjet::Strategy strat = fastjet::Best;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
    double Rparam = .4;

    fastjet::JetDefinition jet_def_akt(fastjet::antikt_algorithm, Rparam, recombScheme,strat);
    // For global background estimation
    fastjet::JetDefinition jetDef_kt(fastjet::kt_algorithm, Rparam, recombScheme,strat);

    // Area definition
      
    fastjet::Selector jet_eta_selector = fastjet::SelectorAbsRapMax(.5);
    fastjet::Selector part_eta_selector = fastjet::SelectorAbsRapMax(.9);


    fastjet::JetMedianBackgroundEstimator bge(jet_eta_selector, jetDef_kt, areaDef);
    fastjet::Subtractor subtractor(&bge);

    RhoEstimator rho_finder(Rparam, .5, .9); // rparam, jet_eta_max, part_eta_max
    */

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

    // Event loop for events from ROOT file
    for (int iEvent = 0; iEvent < numEvents; ++iEvent) {

        if (iEvent > counter) {
            cout << counter  << " events analyzed from ROOT file" << endl;
            counter += 1000;
        }

        vector<fastjet::PseudoJet> event_particles = met.get_particles(iEvent, part_eta_max, part_pt_min);


        // Finding truth level jets
        fastjet::ClusterSequence clustSeq_pythia(event_particles, jet_def_akt);
        vector<fastjet::PseudoJet> jets_pythia = fastjet::sorted_by_pt(jet_eta_selector(clustSeq_pythia.inclusive_jets()));

        // make a copy
        vector<fastjet::PseudoJet> all_particles = event_particles;

        // Embed the background
        vector<fastjet::PseudoJet> background_prtcls = mbt.get_particles(iEvent, part_pt_min); // no eta selection required because bg particles are only up to .9 eta
        move(background_prtcls.begin(), background_prtcls.end(), back_inserter(all_particles));
    

        // ==============================================
        // Loop through option
        if (debug) cout << "Looping through options" << endl;
        for (int iOp = 0; iOp < all_sub_ops.size(); iOp++) {

            unified_sub_options sub_ops = all_sub_ops[iOp];

            // Need to make a copy because event_subtractor edits these directly. Different form storeLundNext and r_groomed where a function call ensures a copy is being used
            vector<fastjet::PseudoJet> particles;
            sub_ops.embed_in_bg? particles = all_particles : particles = event_particles;

          

            // Apply event subtractions
            if (sub_ops.event_sub != event_subtraction::null) {
                event_subtractor.subtract(&particles, all_sub_ops[iOp].event_sub);
            }

            if (particles.empty()) continue;

            if(debug) cout << "Clustering jets" << endl;
            vector<fastjet::PseudoJet> jets_embedded;
            if (sub_ops.jet_sub == jet_subtraction::Area_pT) { // Use cluster sequence area, and record rho
                bge.set_particles(particles);
                rho = bge.rho();
                csa_ptr = std::make_unique<fastjet::ClusterSequenceArea>(particles, jet_def_akt, area_def);
                jets_embedded = jet_eta_selector(csa_ptr->inclusive_jets());
            }
            else { // Use regular cluster sequence
                cs_ptr = std::make_unique<fastjet::ClusterSequence>(particles, jet_def_akt);
                jets_embedded = jet_eta_selector(cs_ptr->inclusive_jets());
            }

            // perform jet subtractions
            if(debug) cout << "Performing jet subtraction" << endl;
            if (sub_ops.jet_sub != jet_subtraction::Area_pT) {// Don't need to subtract anything and don't want to recluster if we're doing area-based pT subtraction
                for (fastjet::PseudoJet& jet_embedded : jets_embedded) {
                    jet_subtractor.subtract_recluster(&jet_embedded, sub_ops.jet_sub, particles);
                }
            }

            // match pythia jets with subtracted jets
            if (debug) cout << "Matching pythia jets with subtracted jet" << endl;
            for (fastjet::PseudoJet jet_pythia : jets_pythia) {

                sel.set_reference(jet_pythia);
                vector<fastjet::PseudoJet> potential_matches = sel(jets_embedded);

                double best_delta = .4;
                double temp_delta;
                fastjet::PseudoJet* closest_match = nullptr;
                for (fastjet::PseudoJet& potential_match : potential_matches) {
                    temp_delta = jet_pythia.delta_R(potential_match);
                    if (temp_delta < best_delta) {
                        best_delta = temp_delta;
                        closest_match = &potential_match;
                    }
                }

                if (closest_match == nullptr) break; // If not match was found break because that means there's not enough pythia jets remaining

                if (debug) cout << "Getting jet pT" << endl;
                double pt_true = jet_pythia.pt();
                double pt_sub = closest_match->pt();
                if (sub_ops.jet_sub == jet_subtraction::Area_pT) pt_sub -= closest_match->area() * rho;
                
                // find cut
                jes_plots[iOp]->Fill(pt_true, (pt_sub - pt_true)/pt_true, met.bin_weight);


                if (leading_jet_only) break; // for leading jet only

            } // end jet loop

        } // end option loop  

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

    // Usually I do 800, 600. Triying 900, 600 so z-axis label isn't cut off
    TCanvas *c1 = new TCanvas("c1", "Canvas", 900, 600);


    // Adding a text information sheet with stats about this generation
    TString title = "JES, JER plots: " + output_file_name;
    TString line0 = notes;
    TString line1 = "Number of events: " + to_string(numEvents);
    TString line2 = "Events from: " + event_file_name;
    TString line3 = "Backgrounds from: " + print_vec(background_file_names);

    TString line4 = "Jet Radius: " + to_string(Rparam);
    TString line5 = "Particle eta max: " + to_string(part_eta_max);
    TString line6 = "Jet eta max: " + to_string(jet_eta_max);
    TString line7 = "Leading jet only: " + to_string(leading_jet_only);

    vector<TString> lines = {title, notes, line1, line2, line3, line4, line5, line6, line7};

    c1->cd();
    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }
    
    c1->Print(output_folder_name + output_file_name+".pdf(","pdf");
    c1->Clear();

    
    for (int iOp = 0; iOp < all_sub_ops.size(); iOp++) {

        
        TGraph* jes = new TGraph();
        jes->SetTitle("Jet Energy Scale vs pt True:" + space + jes_plots[iOp]->GetTitle() + "; pt True; Jet energy scale");
        TGraph* jer = new TGraph();
        jer->SetTitle("Jet Energy Resolution vs pt True" + space + jes_plots[iOp]->GetTitle() + "; pt True; Jet energy resolution");
        
        for (int iBin = 1; iBin <= jes_plots[iOp]->GetNbinsX(); ++iBin) {
            // Use weight for jes, jer vs pt true? -> Don't think so, because weight was accounted for in creation of 2D sub_vs_true
            TH1D* slice = jes_plots[iOp]->ProjectionY("slice", iBin, iBin);
            //if (slice->IsEmpty()) continue;
            // plots point in the middle of the bin
            jes->AddPoint(iBin*jes_xbinWidth + jes_xMin + jes_xbinWidth*.5, slice->GetMean());
            jer->AddPoint(iBin*jes_xbinWidth + jes_xMin + jes_xbinWidth*.5, slice->GetStdDev());
        }
        

        c1->cd();
        jes_plots[iOp]->Draw("COLZ");
        /*
        if ((iOp == all_sub_ops.size() - 1)) {
            c1->Print(output_folder_name + output_file_name+".pdf)","pdf");
        }
        */
        c1->Print(output_folder_name + output_file_name+".pdf","pdf");
                    c1->Clear();

        
        c1->cd();
        jes->Draw();
        c1->Print(output_folder_name + output_file_name+".pdf","pdf");
        c1->Clear();

        c1->cd();
        jer->Draw();
        if (iOp == all_sub_ops.size() - 1) {
            c1->Print(output_folder_name + output_file_name+".pdf)","pdf");
        }
        else c1->Print(output_folder_name + output_file_name+".pdf","pdf");
        c1->Clear();
        

        






    }



    /*
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
    */

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

    /*
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
    */
    

    met.close_file();
    mbt.close_file();
    //output_root_file.Close();

    return 0;
}