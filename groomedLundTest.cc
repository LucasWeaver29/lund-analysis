#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/LundGenerator.hh"


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





// Integrate on jet pt spectrum for n jets
// 1d hist of event weihts
/*
vector<vector<double>> Get_kin_vars(vector<fastjet::PseudoJet> particles, double ptMin, double ptMax, double Rparam = .4, TH1F* pt_weighted = nullptr, TH1F* pt_unweighted = nullptr, double weight = 0, double ptHatMax = 0) {

    bool fastjets_lund = true;

    vector<vector<double>> kin_vars;

    fastjet::contrib::LundGenerator lund; // for fastjet::contrib lunds
    fastjet:: JetDefinition jetDef_akt(fastjet::antikt_algorithm, Rparam, fastjet::E_scheme, fastjet::Best);
    fastjet:: JetDefinition jetDef_ca_recluster(fastjet::cambridge_algorithm, 1, fastjet::E_scheme, fastjet::Best);
    fastjet::ClusterSequence cs(particles, jetDef_akt);
    vector<fastjet::PseudoJet> jets = cs.inclusive_jets();
    for (fastjet::PseudoJet jet : jets) {

        fastjet::ClusterSequence recluster_sequence(jet.constituents(), jetDef_ca_recluster);
        fastjet::PseudoJet reclustered_jet = recluster_sequence.inclusive_jets()[0];

        if ((reclustered_jet.pt() < ptMin) || (reclustered_jet.pt() > ptMax)) continue;

        // Slope/smoothness cut
        if (jet.pt() > 1.5*ptHatMax) {
            cout << "Smoothness cut applied" << endl;
            continue;
        }


        if (pt_weighted != nullptr) {
            pt_weighted->Fill(jet.pt(), weight);
            pt_unweighted->Fill(jet.pt());
        }


        if (fastjets_lund) {
            vector<fastjet::contrib::LundDeclustering> declusterings = lund(reclustered_jet);
            for (const auto& d : declusterings) {
                vector<double> vars = {d.Delta(), d.kt()};
                kin_vars.push_back(vars);
            }

        }
        else {

            fastjet::PseudoJet j1, j2;

            while (reclustered_jet.has_parents(j1, j2)) {

                fastjet::PseudoJet harder = (j1.pt() > j2.pt())? j1 : j2;
                fastjet::PseudoJet softer = (j1.pt() > j2.pt())? j2 : j1;
                
                if (harder.pt() > 100) break;

                double delta = j1.delta_R(j2);
                double kt = sin(delta) * softer.pt();
                vector<double> vars = {delta, kt};
                kin_vars.push_back(vars);
                reclustered_jet = harder;

            }
        } // end my lund 
    } // end jet loop

    return kin_vars;
    
}

*/


bool check4outliers(double x, double y, vector<vector<double>> outlier_ranges_x, vector<vector<double>> outlier_ranges_y) {

    for (int iOutlier = 0; iOutlier < outlier_ranges_x.size(); ++iOutlier) {
        if ((x > outlier_ranges_x[iOutlier][0]) && 
            (x < outlier_ranges_x[iOutlier][1]) &&
            (y > outlier_ranges_y[iOutlier][0]) &&
            (y < outlier_ranges_y[iOutlier][1])) {
                return true;
            }
    }
    return false;

}


int main() {
    
    TString output_file_name = "9k, debug, bin weights";
    TString output_folder = "LundPlanes/";
    
    //TString root_output_name = output_file_name;

    TString eventFileName = "events1000x9, bin weight.root";
    TString backgroundFileName = "thermalBackgrounds45000.root";

    
    double ptMin = 20; // Minimum jet pT. 20 for Alice
    double ptMax = 120; // Maximum jet pt. 120 for Alice
    double part_eta_max = .9;
    double jet_eta_max = .5;
    double Rparam = 0.4;

    bool leading_jet_only = false;

    bool find_outliers = false;
    // Found that the first 5 spikes for 20-120 come from event 191

    // Get bin specifications here from storeWithRoot
    // Make sure recreation of iBin later down also matches
    /*
    vector<double> ptHatMins = {10, 20, 30, 50, 70, 90, 110, 130, 150};
    vector<double> ptHatMaxs = {20, 30, 50, 70, 90, 110, 130, 150, 250};
    int iBin = 0;
    double currentBin = ptHatMins[0];
    */

    //==================================================================
    // Preparing to read from ROOT files
    my_event_tree met(eventFileName);
    //my_background_tree mbt(backgroundFileName);

    // Checking if output subfolder exists
    //std::filesystem::path current_path = std::filesystem::current_path();
    //std::filesystem::path full_path = current_path / output_folder / subfolder

    //=======================================================================
    // Graphs 'n stuff
    // Declare 2D ROOT histograms
    TString xAxisName = "ln((R=.4/delta)";
    TString yAxisName = "ln(kt)";
    TString space = " ";
    double xMin = 0;
    double xMax = 5;
    int nBinsX = 50;
    double yMin = -4;
    double yMax = 3.5;
    int nBinsY = 50;

    TH2F *lund = new TH2F("lund", "Pythia; " + xAxisName + "; " + yAxisName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    TH1F *pt_weighted = new TH1F("pt_weighted", "Pythia: Weighted Jet pt spectrum; jet pt; weighted counts", 125 , -.5, 249.5);

    // For debugging
    TH1F *pt_unweighted = new TH1F("pt_pythia_unweighted", "Pythia: Unweighted Jet pt spectrum; jet pt; counts", 125 , -.5, 249.5);
    TH2F *lund_eventWeighted_unnormalized = new TH2F("lund_eventWeighted_unnormalized", "Event weights, unnormalized; " + xAxisName + "; " + yAxisName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    TH2F *lund_eventWeighted_normalized = new TH2F("lund_eventWeighted_normalized", "Event weights, normalized; " + xAxisName + "; " + yAxisName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    TH2F *lund_unWeighted_unnormalized = new TH2F("lund_unWeighted_unnormalized", "Unweighted, unnormalized; " + xAxisName + "; " + yAxisName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    TH2F *lund_unWeighted_normalized = new TH2F("lund_unWeighted_normalized", "Unweighted, normalized; " + xAxisName + "; " + yAxisName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    TH1F *event_weights = new TH1F("event_weights", "Event Weights; log(event weight); Counts", 100, -22, -4);
    

    // Implementing jet pt cuts
    vector<double> ptCutMins = {10, 20, 30, 50, 70, 100, 130, 160, 190};
    vector<double> ptCutMaxs = {20, 30, 50, 70, 100, 130, 160, 190, 250};
    
    vector<TH2F*> pythia_lund_cuts;
    for (int iCut = 0; iCut < ptCutMins.size(); ++iCut) {
            pythia_lund_cuts.push_back(new TH2F("plc" + space + to_string(iCut), "jet pT" + space + ptCutMins[iCut] + " to " + ptCutMaxs[iCut] + "; " + xAxisName + "; " + yAxisName, nBinsX, xMin, xMax, nBinsY, yMin, yMax));
    }

    /*
    TH2F *lund_bg = new TH2F("lund_bg", "Background; " + xAxisName + "; " + yAxisName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    TH2F *lund_constit_sub = new TH2F("lund_constit_sub", "Constituent Subtraction; " + xAxisName + "; " + yAxisName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    TH2F *lund_softKill = new TH2F("lund_softKill", "SoftKiller; " + xAxisName + "; " + yAxisName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    TH2F *lund_my_rsd = new TH2F("lund_my_rsd", "Recursive Soft Drop - mine; " + xAxisName + "; " + yAxisName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    //TH2F *lund_contrib_rsd = new TH2F("lund_contrib_rsd", "Recursive Soft Drop - fj contrib; " + xAxisName + "; " + yAxisName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    TH2F *lund_my_pc = new TH2F("lund_my_pc", "My pc; " + xAxisName + "; " + yAxisName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);

    vector<TH2F*> bg_lunds =  {lund_bg, lund_constit_sub, lund_softKill, lund_my_rsd, lund_my_pc};
    // For keeping track of number of jets, weighted by pt, for average lund density.
    vector<TH1F*> pt_bg; 
    for (int i = 0; i < bg_lunds.size(); ++i) {
        pt_bg.push_back(new TH1F("pt_bg" + to_string(i), space + bg_lunds[i]->GetTitle() + ": jet pt spectrum; jet pt; weighted counts", 125, 0, 250));
    }


    vector<bg_sub_options> options =   {bg_sub_options{}, 
                                        bg_sub_options{.jet_by_jet_constit = true},
                                        bg_sub_options{.soft_kill = true},
                                        bg_sub_options{.my_rsd = true},
                                        bg_sub_options{.my_pc= true}};
    

    */

    // Implementing jet pt cuts
    
    // Figuring out why there's spikes
    // Outliers for 9k, debug, 20 - 120
    vector<vector<double>> outlier_ranges_x = {{.2, .3}, {.6, .7}, {1.2, 1.3}, {1.4, 1.5}, {2, 2.1}};
    vector<vector<double>> outlier_ranges_y = {{-.25, -.1}, {-.4, -.25}, {-1.15, -1}, {-1, -.85}, {-2.65, -2.5}};



    // Set up fastjet analysis
    LundGroomer lund_groomer(ptMin, ptMax, leading_jet_only, Rparam, jet_eta_max, part_eta_max);
    int iCut = -1;
    int* iCut_ptr = &iCut;
    lund_groomer.set_cuts(iCut_ptr, ptCutMins, ptCutMaxs);

    // =======================================================================

    // Reading events from ROOT
    int numEvents = met.GetEntries();
    //std::cout << "Number of events: " << numEvents << std::endl;

    /*
    int numBackgrounds = mbt.GetEntries();
    if (numBackgrounds < numEvents) {
        cout << "Warning: There are " << numEvents << " provided events but only" << numBackgrounds << " provided backgrounds." << endl;
        cout << "Lund plane creation will stop once backgrounds run out" << endl;
        numEvents = numBackgrounds;
    }
    */

    int counter = 1000;

    // Loop over events in ROOT TTree
    for (int iEvent = 0; iEvent < numEvents; ++iEvent) {
         
        //if (iEvent == 191) continue;

        
        if (iEvent > counter) {
            cout << counter << " events analyzed from ROOT file" << endl;
            counter += 1000;
        }
        

        vector<fastjet::PseudoJet> event_particles = met.get_particles(iEvent, .9); // iEvent, part_eta_max

        event_weights->Fill(log(met.bin_weight));

        vector<vector<double>> pythia_kin_vars = lund_groomer.get_kin_vars(event_particles, bg_sub_options{.jet_pt_hist = pt_weighted, .jet_pt_hist_unweighted = pt_unweighted, .event_weight = met.bin_weight});
        //vector<vector<double>> pythia_kin_vars = lund_groomer.get_kin_vars(event_particles, bg_sub_options{.jet_pt_hist = pt_weighted, .jet});
        
        
        for (vector<double> vec : pythia_kin_vars) {
            double log_inv_R0 = log(Rparam/vec[0]);
            double log_pt = log(vec[1]);
            lund_eventWeighted_unnormalized->Fill(log_inv_R0,log_pt, met.bin_weight);
            lund_eventWeighted_normalized->Fill(log_inv_R0, log_pt, met.bin_weight);
            lund_unWeighted_unnormalized->Fill(log_inv_R0, log_pt);
            lund_unWeighted_normalized->Fill(log_inv_R0, log_pt);
                
            if (iCut != -1) pythia_lund_cuts[iCut]->Fill(log_inv_R0, log_pt, met.bin_weight);
            
            if (find_outliers) {
                if(check4outliers(log_inv_R0, log_pt, outlier_ranges_x, outlier_ranges_y) && iCut == 1) {
                    cout << "Outlier from iEvent " << iEvent << ". Weight: " << met.bin_weight << endl;
                }
            }
            //lund_pythia->Fill(log(Rparam/vec[0]),log(vec[1]));
        }

        /*
        // Embedding in the background
        vector<fastjet::PseudoJet> background_prtcls = mbt.get_particles(iEvent);
        move(background_prtcls.begin(), background_prtcls.end(), back_inserter(event_particles));

        
        for (int iOption = 0; iOption < options.size(); ++iOption) {
            vector<vector<double>> kin_vars = lund_groomer.get_kin_vars(event_particles, options[iOption]);
            for (vector<double> vec : kin_vars) {
                bg_lunds[iOption]->Fill(log(Rparam/vec[0]),log(vec[1]), met.weight);
                //bg_lunds[iOption]->Fill(log(Rparam/vec[0]),log(vec[1]));
            }
        }
        */
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

    // Pythia
    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    // Adding a text information sheet with stats about this generation
    TString title = "Groomed Lund planes: " + output_file_name;
    TString line1 = "Number of events: " + to_string(numEvents);
    TString line2 = "Events from: " + eventFileName;
    TString line3 = "Backgrounds from: " + backgroundFileName;
    TString line4 = "Leading jet only: " + bool2Str(leading_jet_only);
    TString line5 = "Jet pT min: " + to_string(ptMin);
    TString line6 = "Jet pT max: " + to_string(ptMax);

    TString line7 = "Jet radius: " + to_string(Rparam);
    TString line8 = "Particle eta max: " + to_string(part_eta_max);
    TString line9 = "Jet eta max: " + to_string(jet_eta_max);

    vector<TString> lines = {title, line1, line2, line3, line4, line5, line6, line7, line8, line9};

    c1->cd();
    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }
    
    c1->Print(output_folder + output_file_name +".pdf(","pdf");
    c1->Clear();


    // Finding weighted njets
    double nJets_pythia = pt_weighted->Integral();
    //cout << "FindBin(ptMin): " << pt_pythia->FindBin(ptMin) << endl;
    //cout << "FindBin(ptMax): " << pt_pythia->FindBin(ptMax) << endl;
    /*
    for(int iBin = pt_weighted->FindBin(ptMin); iBin < pt_weighted->FindBin(ptMax); ++iBin) {
        nJets_pythia += pt_weighted->GetBinContent(iBin);
    }
    */

    // Normalizing: Divide each bin by bin area and number of jets
    double area = ((xMax - xMin)/nBinsX) * ((yMax-yMin)/nBinsY);
    double factor = 1/(area * nJets_pythia);
    //lund_pythia->Scale(factor);
    lund_eventWeighted_normalized->Scale(factor);
    lund_unWeighted_normalized->Scale(factor);
                    
    //cout << "Factor:" << factor << endl;
    /*
    for (TH2F* cut : pythia_lund_cuts) {
        cut->Scale(factor);
    }
    */
    
    // Find the cause of the 5 spots
    if (find_outliers) {
        cout << "Bin indices of outliers" << endl;
        for (int x = 1; x <= lund_eventWeighted_normalized->GetNbinsX(); ++x) {
            for (int y = 1; y <= lund_eventWeighted_normalized->GetNbinsY(); ++y) {
                if (lund_eventWeighted_normalized->GetBinContent(x,y) > 3) {
                    cout << "Indices: (" << x << "," << y << ")" << endl;
                    double xlow = lund_eventWeighted_normalized->GetXaxis()->GetBinLowEdge(x);
                    double xhigh = lund_eventWeighted_normalized->GetXaxis()->GetBinUpEdge(x);
                    double ylow = lund_eventWeighted_normalized->GetYaxis()->GetBinLowEdge(y);
                    double yhigh = lund_eventWeighted_normalized->GetYaxis()->GetBinUpEdge(y);

                    cout << "Edges: X: " << xlow << ", " << xhigh <<". Y: " << ylow << ", " << yhigh << endl;
                }
            }
        }
    }
    /*
    Bin indices of outliers
    (3,26)
    (7,25)
    (13,20)
    (15,21)
    (21,10)
    */


    c1->cd();
    lund_unWeighted_unnormalized->Draw("COLZ");
    c1->Print(output_folder + output_file_name +".pdf","pdf");
    c1->Clear();

    //c1->SetLogz(1);
    c1->cd();
    lund_eventWeighted_unnormalized->Draw("COLZ");
    c1->Print(output_folder + output_file_name +".pdf","pdf");
    c1->Clear();

    c1->cd();
    lund_eventWeighted_normalized->Draw("COLZ");
    c1->Print(output_folder + output_file_name +".pdf","pdf");
    c1->Clear();

    //c1->SetLogz(0);

    c1->cd();
    lund_unWeighted_normalized->Draw("COLZ");
    c1->Print(output_folder + output_file_name +".pdf","pdf");
    c1->Clear();

    c1->cd();
    c1->SetLogy(1);
    event_weights->Draw("HIST");
    c1->Print(output_folder + output_file_name + ".pdf", "pdf");
    c1->Clear();

    c1->cd();
    pt_weighted->Draw("HIST");
    c1->Print(output_folder + output_file_name +".pdf","pdf");
    c1->Clear();

    c1->cd();
    pt_unweighted->Draw("HIST");
    c1->Print(output_folder + output_file_name +".pdf","pdf");
    c1->Clear();


    c1->SetLogy(0);
    for (int iLund = 0; iLund < pythia_lund_cuts.size(); ++iLund) {
        
        c1->cd();
        pythia_lund_cuts[iLund]->Draw("COLZ");
        if (iLund + 1 == pythia_lund_cuts.size()) {
            c1->Print(output_folder + output_file_name +".pdf)","pdf");
        } 
        else {
            c1->Print(output_folder + output_file_name +".pdf","pdf");
        }
        c1->Clear();
    }
    
    /*
    // bg Lunds
    for (int iFile = 0; iFile < bg_file_names.size(); ++iFile) {
        
        // Create a ROOT canvas
        TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

        // Adding a text information sheet with stats about this generation
        TString title = "Groomed Lund planes:; " + subfolder + "; " + bg_file_names[iFile];
        TString line1 = "Number of events: " + to_string(numEvents);
        TString line2 = "Events from: " + eventFileName;
        TString line3 = "Backgrounds from: " + backgroundFileName;

        TString line4 = "Jet pT min: " + to_string(ptMin);
        TString line5 = "Jet pT max: " + to_string(ptMax);

        TString line6 = "Jet radius: " + to_string(Rparam);
        TString line7 = "Particle eta max: " + to_string(part_eta_max);
        TString line8 = "Jet eta max: " + to_string(jet_eta_max);

        vector<TString> lines = {title, line1, line2, line3, line4, line5, line6, line7, line8};

        c1->cd();
        for (int iLine = 0; iLine < lines.size(); ++iLine) {
            TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
            text->SetTextSize(.04);
            text->Draw();
        }
        
        c1->Print(output_folder + subfolder + bg_file_names[iFile] + ".pdf(","pdf");
        c1->Clear();

        // Normalizing: Divide each bin by bin area and number of jets
        
        cout << "nJets: " << nJets_bg[iFile] << endl;
        double factor = 1/(area * nJets_bg[iFile]);
        bg_lunds[iFile]->Scale(factor);
        for (TH2F* lund : bg_lund_cuts[iFile]) {
            lund->Scale(factor);
        }
        *

        c1->cd();
        bg_lunds[iFile]->Draw("COLZ");
        c1->Print(output_folder + subfolder + bg_file_names[iFile] + ".pdf","pdf");
        c1->Clear();

        for (int iLund = 0; iLund < bg_lund_cuts[iFile].size(); ++iLund) {
            
            c1->cd();
            bg_lund_cuts[iFile][iLund]->Draw("COLZ");
            if (iLund + 1 == bg_lund_cuts[iFile].size()) {
                c1->Print(output_folder + subfolder + bg_file_names[iFile] + ".pdf)","pdf");
            } 
            else {
                c1->Print(output_folder + subfolder + bg_file_names[iFile] + ".pdf","pdf");
            }
            c1->Clear();
        }
    } // End file loop

    */

    // Since ROOT histograms point to the data they store, once the eventRootFile
    // closes they no longer have access to their data
    met.close_file();
    //mbt.close_file();

    return 0;
}



