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
#include "backSubTools.h"

using namespace Pythia8;


struct lund_option {
    TString name;
    bg_sub_options sub_options;
};

int main() {
    
    //TString output_file_name = "9k, groomed";
    TString output_folder = "LundPlanes";
    TString subfolder = "9k, groomed";

    cout << "Make sure you have created the needed subfolder!" << endl;

    TString notes = "";
    //"z_cut = .2, beta = 3";

    vector<lund_option> pythia_lund_options = 
        {
            //lund_option{.name = "Pythia", .sub_options = bg_sub_options{}},
            //lund_option{"Pythia with Contrib's SoftDrop", bg_sub_options{.contrib_softDrop = true}}
            lund_option{"Pythia with my SoftDrop", bg_sub_options{.my_softDrop = true}}

        };

    vector<lund_option> bg_lund_options = 
        {
            //lund_option{.name = "Constituent subtraction, PC particles", bg_sub_options{.pc_particles = true, .constit_sub = true}},
            //lund_option{.name = "Contrib's Recursive SoftDrop", .sub_options = bg_sub_options{.contrib_rsd = true}},
            //lund_option{.name = "Embedded in bg, log", bg_sub_options{}},
            //lund_option{.name = "Constituent subtraction, RSD", bg_sub_options{.constit_sub = true, .contrib_rsd = true}},
            //lund_option{.name = "Constituent subtraction, SD", bg_sub_options{.constit_sub = true, .contrib_softDrop = true}},
            //lund_option{.name = "softKiller", bg_sub_options{.soft_kill = true}},
            //lund_option{.name = "Recursive SoftDrop", bg_sub_options{.my_rsd = true}},
            //lund_option{.name = "My Perpendicular Cone", bg_sub_options{.my_pc= true}},
            //lund_option{.name = "My Perpendicular Cone, with SD", bg_sub_options{.my_pc= true, .contrib_softDrop = true}}
            //lund_option{.name = "My Perpendicular Cone, geometric match, log", bg_sub_options{.my_pc_geometric_match= true}}
            //lund_option{.name = "Constituent subtraction", bg_sub_options{.constit_sub=true}}
        };


    TString eventFileName = "events1000x9, bin weight.root";
    TString backgroundFileName = "thermalBackgrounds9000.root";
    //"thermalBackgrounds9000,etaMax=2.root";

    double ptMin = 20; // Minimum jet pT. 20 for Alice
    double ptMax = 120; // Maximum jet pt. 120 for Alice
    double part_eta_max = .9;
    double jet_eta_max = .5;
    double Rparam = .4;
    bool leading_jet_only = false;
    if (part_eta_max - Rparam != jet_eta_max) cout << "WARNING: Rparam, part_eta_max, and jet_eta_max are not in agreement" << endl;

    bool debug = false;

    //==================================================================
    // Preparing to read from ROOT files
    my_event_tree met(eventFileName);
    my_background_tree mbt(backgroundFileName);
    // Checking if output subfolder exists
    //std::filesystem::path current_path = std::filesystem::current_path();
    //std::filesystem::path full_path = current_path / output_folder / subfolder

    //=======================================================================
    // Graphs 'n stuff
    // Declare 2D ROOT histograms
    TString xAxisName = "ln((R0/delta)";
    TString yAxisName = "ln(kt)";
    TString space = " ";
    double lund_xMin = 0;
    double lund_xMax = 5;
    int lund_nBinsX = 50;
    double lund_yMin = -4;
    double lund_yMax = 3.5;
    int lund_nBinsY = 50;

    int pt_nBins = 150;

    vector<double> ptCutMins = {10, 20, 30, 50, 70, 100, 130, 160, 190};
    vector<double> ptCutMaxs = {20, 30, 50, 70, 100, 130, 160, 190, 250};


    vector<vector<TH2F*>> pythia_lund_cuts;
    vector<TH1F*> pt_pythia;
    for (int iOption = 0; iOption < pythia_lund_options.size(); ++iOption) {

        vector<TH2F*> cuts;

        // Push back inclusive lund plane to index 0 of cuts
        cuts.push_back(new TH2F(pythia_lund_options[iOption].name, pythia_lund_options[iOption].name + "; " + xAxisName + "; " + yAxisName, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax));

        // Add a histogram for tracking weighted jet pt counts to pt_bg
        pt_pythia.push_back(new TH1F("pt_pythia" + space + to_string(iOption), pythia_lund_options[iOption].name + ": jet pt spectrum; jet pt; weighted counts", pt_nBins, -.5, 299.5));
        pythia_lund_options[iOption].sub_options.jet_pt_hist = pt_pythia[iOption];

        // Add lund plane jet pt cuts
        for (int iCut = 0; iCut < ptCutMins.size(); ++iCut) {
            cuts.push_back(new TH2F("pythia_lund_cuts" + pythia_lund_options[iOption].name + to_string(iCut), "jet pT" + space + ptCutMins[iCut] + " to " + ptCutMaxs[iCut] + "; " + xAxisName + "; " + yAxisName, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax));
        }
        pythia_lund_cuts.push_back(cuts);

    }


    /*
    TH1F *pt_pythia = new TH1F("pt_pythia", "Pythia; Jet pt spectrum; weighted counts", pt_nBins, -.5, 299.5);

    if (do_pythia) {
        // Push inclusive jet pt pythia lund to index 0 of pythia_lund_cts
        pythia_lund_cuts.push_back(new TH2F("lund_pythia", "Pythia; " + xAxisName + "; " + yAxisName, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax)); 
        
        for (int iCut = 0; iCut < ptCutMins.size(); ++iCut) {
            pythia_lund_cuts.push_back(new TH2F("plc" + space + to_string(iCut), "jet pT" + space + ptCutMins[iCut] + " to " + ptCutMaxs[iCut] + "; " + xAxisName + "; " + yAxisName, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax));
        }
    }
    */

    //vector<TString> bg_file_names;
    vector<vector<TH2F*>> bg_lund_cuts; // index 0 is inclusive
    vector<TH1F*> pt_bg;
    for (int iOption = 0; iOption < bg_lund_options.size(); ++iOption) {

        vector<TH2F*> cuts;

        // Push back inclusive lund plane to index 0 of cuts
        cuts.push_back(new TH2F(bg_lund_options[iOption].name, bg_lund_options[iOption].name + "; " + xAxisName + "; " + yAxisName, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax));

        // Add a histogram for tracking weighted jet pt counts to pt_bg
        pt_bg.push_back(new TH1F("pt_bg" + space + to_string(iOption), bg_lund_options[iOption].name + ": jet pt spectrum; jet pt; weighted counts", pt_nBins, -.5, 299.5));
        bg_lund_options[iOption].sub_options.jet_pt_hist = pt_bg[iOption];

        // Add lund plane jet pt cuts
        for (int iCut = 0; iCut < ptCutMins.size(); ++iCut) {
            cuts.push_back(new TH2F("blc" + bg_lund_options[iOption].name + to_string(iCut), "jet pT" + space + ptCutMins[iCut] + " to " + ptCutMaxs[iCut] + "; " + xAxisName + "; " + yAxisName, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax));
        }
        bg_lund_cuts.push_back(cuts);

    }
    
    // Set up fastjet analysis
    LundGroomer lund_groomer(ptMin, ptMax, leading_jet_only, Rparam, jet_eta_max, part_eta_max);
    lund_groomer.set_cuts(ptCutMins, ptCutMaxs);
    
    // =======================================================================

    // Reading events from ROOT
    int numEvents = met.GetEntries();
    //std::cout << "Number of events: " << numEvents << std::endl;

    int numBackgrounds = mbt.GetEntries();
    if (numBackgrounds < numEvents) {
        cout << "Warning: There are " << numEvents << " provided events but only" << numBackgrounds << " provided backgrounds." << endl;
        cout << "Lund plane creation will stop once backgrounds run out" << endl;
        numEvents = numBackgrounds;
    }

    int counter = 1000;

    // Loop over events in ROOT TTree
    for (int iEvent = 0; iEvent < numEvents; ++iEvent) {
        
        if (debug) cout << "iEvent: " << iEvent << endl;

        if (iEvent > counter) {
            cout << counter << " events analyzed from ROOT file" << endl;
            counter += 1000;
        }

        vector<fastjet::PseudoJet> event_particles = met.get_particles(iEvent, part_eta_max); // iEvent, part_eta_max

        for (int iOption = 0; iOption < pythia_lund_options.size(); ++iOption) {
            pythia_lund_options[iOption].sub_options.bin_weight = met.bin_weight;
            vector<vector<double>> kin_vars = lund_groomer.get_kin_vars(event_particles, pythia_lund_options[iOption].sub_options);
            for (vector<double> vec : kin_vars) {
                double logR0_R = log(Rparam/vec[0]);
                double logkt = log(vec[1]);
                pythia_lund_cuts[iOption][0]->Fill(logR0_R,logkt, met.bin_weight); // The inclusive pt lund
                pythia_lund_cuts[iOption][vec[2] + 1]->Fill(logR0_R,logkt, met.bin_weight); // the jet pt cut lund (which is at iCut + 1 in pythia_lund_cuts);
            }
        }
        

        // Embedding in the background
        vector<fastjet::PseudoJet> background_prtcls = mbt.get_particles(iEvent);
        move(background_prtcls.begin(), background_prtcls.end(), back_inserter(event_particles));

        for (int iOption = 0; iOption < bg_lund_options.size(); ++iOption) {
            bg_lund_options[iOption].sub_options.bin_weight = met.bin_weight;
            //bg_lund_options[iOption].sub_options.jet_pt_hist
            vector<vector<double>> kin_vars = lund_groomer.get_kin_vars(event_particles, bg_lund_options[iOption].sub_options);
            for (vector<double> vec : kin_vars) {
                double logR0_R = log(Rparam/vec[0]);
                double logkt = log(vec[1]);
                bg_lund_cuts[iOption][0]->Fill(logR0_R,logkt, met.bin_weight); // The inclusive pt lund
                bg_lund_cuts[iOption][vec[2] + 1]->Fill(logR0_R,logkt, met.bin_weight); // the jet pt cut lund (which is at iCut + 1 in pythia_lund_cuts);
            }
        }

    }  // End reading events from ROOT
    

// =========================================================
    // Save these plots in a root file, so they can be modified/resized without having to run the whole thing again
    /*
    TFile output_root_file(output_folder_name + root_output_name + ".root", "RECREATE");

    TH2F* lund_pointer = lund_pythia;

    TTree lund_tree("lund_tree", "Lund Tree");
    lund_tree.Branch("lund", &lund_pointer);
 
    // First entry is inclusive lund
    lund_tree.Fill();
    // subsequent entries are cuts
    for (int iLund = 0; iLund < lunds.size(); ++iLund) {
        lund_pointer = lunds[iLund];
        lund_tree.Fill();
    }

    lund_tree.Write();
    output_root_file.Close();
    */
//========================================================
    // Output files

    double area = ((lund_xMax - lund_xMin)/lund_nBinsX) * ((lund_yMax-lund_yMin)/lund_nBinsY); // For normalization


    // Pythia ================================================
    for (int iOption = 0; iOption < pythia_lund_options.size(); ++iOption) {
        
        TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

        // Adding a text information sheet with stats about this generation
        TString title = "Groomed Lund planes:; " + subfolder + "; " + pythia_lund_options[iOption].name;
        TString line1 = "Number of events: " + to_string(numEvents);
        TString line2 = "Events from: " + eventFileName;
        TString line3 = "Backgrounds from: " + backgroundFileName;
        TString line4 = "Leading Jet Only: " + bool2Str(leading_jet_only);

        TString line5 = "Jet pT min: " + to_string((int)ptMin);
        TString line6 = "Jet pT max: " + to_string((int)ptMax);

        TString line7 = "Jet radius: " + to_string_round(Rparam);
        TString line8 = "Particle eta max: " + to_string_round(part_eta_max);
        TString line9 = "Jet eta max: " + to_string_round(jet_eta_max);

        vector<TString> lines = {title, notes, line1, line2, line3, line4, line5, line6, line7, line8, line9};

        c1->cd();
        for (int iLine = 0; iLine < lines.size(); ++iLine) {
            TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
            text->SetTextSize(.04);
            text->Draw();
        }
        
        c1->Print(output_folder + "/" + subfolder + "/"  + pythia_lund_options[iOption].name + ".pdf(","pdf");
        c1->Clear();

        // Normalizing: Divide each bin by bin area and number of jets
        //cout << "nJets: " << nJets_bg[iFile] << endl;
        double nbgJets_weighted = 0;
        for(int iBin = pt_pythia[iOption]->FindBin(ptMin); iBin < pt_pythia[iOption]->FindBin(ptMax); ++iBin) {
            nbgJets_weighted += pt_pythia[iOption]->GetBinContent(iBin);
        }
        // Normalizing: Divide each bin by bin area and number of jets
        //area = ((lund_xMax - lund_xMin)/lund_nBinsX) * ((lund_yMax-lund_yMin)/lund_nBinsY);
        double factor = 1/(area * nbgJets_weighted);
        for (TH2F* lund : pythia_lund_cuts[iOption]) {
            lund->Scale(factor);
        }
         
        for (int iLund = 0; iLund < pythia_lund_cuts[iOption].size(); ++iLund) {
            
            c1->cd();
            pythia_lund_cuts[iOption][iLund]->Draw("COLZ");
            c1->Print(output_folder + "/" + subfolder + "/"  + pythia_lund_options[iOption].name + ".pdf","pdf");
            c1->Clear();
        }

        c1->cd();
        c1->SetLogy(1);
        pt_pythia[iOption]->Draw("HIST");
        c1->Print(output_folder + "/" + subfolder + "/"  + pythia_lund_options[iOption].name + ".pdf)","pdf");
        c1->Clear();
        c1->SetLogy(0);

    } // End pythia file loop
    


    // bg Lunds ======================================================
    for (int iOption = 0; iOption < bg_lund_options.size(); ++iOption) {
        
        TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);


        // Adding a text information sheet with stats about this generation
        TString title = "Groomed Lund planes:; " + subfolder + "; " + bg_lund_options[iOption].name;
        TString line1 = "Number of events: " + to_string(numEvents);
        TString line2 = "Events from: " + eventFileName;
        TString line3 = "Backgrounds from: " + backgroundFileName;
        TString line4 = "Leading Jet Only: " + to_string(leading_jet_only);

        TString line5 = "Jet pT min: " + to_string(ptMin);
        TString line6 = "Jet pT max: " + to_string(ptMax);

        TString line7 = "Jet radius: " + to_string(Rparam);
        TString line8 = "Particle eta max: " + to_string(part_eta_max);
        TString line9 = "Jet eta max: " + to_string(jet_eta_max);

        vector<TString> lines = {title, notes, line1, line2, line3, line4, line5, line6, line7, line8, line9};

        c1->cd();
        for (int iLine = 0; iLine < lines.size(); ++iLine) {
            TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
            text->SetTextSize(.04);
            text->Draw();
        }
        
        c1->Print(output_folder + "/" + subfolder + "/"  + bg_lund_options[iOption].name + ".pdf(","pdf");
        c1->Clear();

        // Normalizing: Divide each bin by bin area and number of jets
        //cout << "nJets: " << nJets_bg[iFile] << endl;
        double nbgJets_weighted = 0;
        for(int iBin = pt_bg[iOption]->FindBin(ptMin); iBin < pt_bg[iOption]->FindBin(ptMax); ++iBin) {
            nbgJets_weighted += pt_bg[iOption]->GetBinContent(iBin);
        }
        // Normalizing: Divide each bin by bin area and number of jets
        //area = ((lund_xMax - lund_xMin)/lund_nBinsX) * ((lund_yMax-lund_yMin)/lund_nBinsY);
        double factor = 1/(area * nbgJets_weighted);
        for (TH2F* lund : bg_lund_cuts[iOption]) {
            lund->Scale(factor);
        }
         
        for (int iLund = 0; iLund < bg_lund_cuts[iOption].size(); ++iLund) {
            
            c1->cd();
            bg_lund_cuts[iOption][iLund]->Draw("COLZ");
            c1->Print(output_folder + "/" + subfolder + "/"  + bg_lund_options[iOption].name + ".pdf","pdf");
            c1->Clear();
        }

        c1->cd();
        c1->SetLogy(1);
        pt_bg[iOption]->Draw("HIST");
        c1->Print(output_folder + "/" + subfolder + "/"  + bg_lund_options[iOption].name + ".pdf)","pdf");
        c1->Clear();
        c1->SetLogy(0);

    } // End file loop


    // =======================================
    // Messing with it manually to directly subtract two lund planes
    /*
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    TString manual_file_name = "Constituent Subtraction, Jet - PC";

    // Adding a text information sheet with stats about this generation
    TString title = "Groomed Lund planes: " + subfolder + "; " + manual_file_name;
    TString line1 = "Number of events: " + to_string(numEvents);
    TString line2 = "Events from: " + eventFileName;
    TString line3 = "Backgrounds from: " + backgroundFileName;
    TString line4 = "Leading Jet Only: " + to_string(leading_jet_only);

    TString line5 = "Jet pT min: " + to_string(ptMin);
    TString line6 = "Jet pT max: " + to_string(ptMax);

    TString line7 = "Jet radius: " + to_string(Rparam);
    TString line8 = "Particle eta max: " + to_string(part_eta_max);
    TString line9 = "Jet eta max: " + to_string(jet_eta_max);

    vector<TString> lines = {title, notes, line1, line2, line3, line4, line5, line6, line7, line8, line9};

    c1->cd();
    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }

    c1->Print(output_folder + "/" + subfolder + "/"  + manual_file_name + ".pdf(","pdf");
    c1->Clear();


    for (int iLund = 0; iLund < bg_lund_cuts[0].size(); ++iLund) {
        c1->cd();
        TH2F* sub = (TH2F*)bg_lund_cuts[1][iLund]->Clone("sub");
        sub->Add(bg_lund_cuts[0][iLund], -1);
        sub->Draw("COLZ");
        if (iLund == bg_lund_cuts[0].size()-1) c1->Print(output_folder + "/" + subfolder + "/"  + manual_file_name + ".pdf)","pdf");
        else c1->Print(output_folder + "/" + subfolder + "/"  + manual_file_name + ".pdf","pdf");
        c1->Clear();
    }
    */
    // end manually messing


    // Since ROOT histograms point to the data they store, once the eventRootFile
    // closes they no longer have access to their data
    met.close_file();
    mbt.close_file();

    return 0;
}



