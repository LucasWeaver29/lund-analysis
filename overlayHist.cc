// pythia, fastjet
#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"


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
#include "TGraph.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace Pythia8;

int main() {

    TString file1_name = "LundPlanes/45k, pythia.root";
    TString file2_name = "LundPlanes/45k, bg.root";
    TString file3_name = "LundPlanes/45k, area sub.root";

    vector<TString> file_names = {file2_name, file3_name, file1_name};
    vector<TString> file_names_to_print = {"bg embedded", "area sub", "pythia"};

    TString output_file_name = "45k, overlay"; //  a pdf
    TString output_folder_name = "LundPlanes/";
    //TString output_name = "Lund, bg, area sub - 45000, resized";

    // cuts of cuts
    vector<double> kt_cut_params = {-2, -1, 0, 2, 2.5}; // Lines of ln(kt) = x along which to take cuts
    vector<double> delta_cut_params = {.1, .2, .3, .5, 1};

    // pt lund cuts, for labeling. From makeLund. These are lund 9000 series
    //vector<double> ptCutMins = {10, 20, 30, 50, 70, 100, 130, 160, 190};
    //vector<double> ptCutMaxs = {20, 30, 50, 70, 100, 130, 160, 190, 250};
    
        
    // Use bin sizes, max, min in makeLund to convert from kt, delta values to bin numbers

    // This was the TH2F decleration for the Lund 9000 series: TH2F *lund_inclusive = new TH2F("lund_inclusive", "Inclusive pt; " + xAxisName + "; " + yAxisName, 400, 0, 5, 400, -4, 3.5);
    // For 45k series:  TH2F *lund_inclusive = new TH2F("lund_inclusive", "Inclusive pt; " + xAxisName + "; " + yAxisName, 50, 0, 5, 50, -4, 3.5);

    vector<double> delta_iBins;
    for (int iDelta = 0; iDelta < delta_cut_params.size(); ++iDelta) {
        delta_iBins.push_back((int) ((delta_cut_params[iDelta] / 5) * 50));
    }
    vector<double> kt_iBins;
    for (int i_kt = 0; i_kt < kt_cut_params.size(); ++i_kt) {
        kt_iBins.push_back((int) (((kt_cut_params[i_kt] + 4) / 7.5) * 50));
    }


    // Opening root files, trees
    vector<TFile*> files;
    vector<TTree*> trees;
    vector<TH2F*> hists(file_names.size(), nullptr);

    /*
    for (int iHist = 0; iHist < hists.size(); ++iHist) {
        // Again using known TH2F configuration from Lund 9000 series
        hists[iHist] = new TH2F(Form("hist_%d", iHist), "Lund Histogram", 400, 0, 5, 400, -4, 3.5);
    }
    */


    for (int iFile = 0; iFile < file_names.size(); ++iFile) {
        files.push_back(new TFile(file_names[iFile], "READ"));
        if(files[iFile]->IsZombie()) {
            cout << "Error: Coult not open Root file " << file_names[iFile] << endl;
            return 1;
        }

        trees.push_back((TTree*)files[iFile]->Get("lund_tree"));
        if (!trees[iFile]) {
            cout << "Error: Could not find tree in ROOT file" << endl;
            return 1;
        }

        trees[iFile]->SetBranchAddress("lund", &hists[iFile]);
    }
    


    // =====================================================================
    // Set up root canvas
    
    
    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    // Adding a text information sheet with stats about this generation
    TString title = "kt, delta cuts of primary lund plane: " + output_file_name;
    TString line1 = "File 1: " + file1_name;
    TString line2 = "File 2: " + file2_name;
    TString line3 = "File 3: " + file3_name;
    

    vector<TString> lines = {title, line1, line2, line3};

    c1->cd();
    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }
    
    c1->Print(output_folder_name+output_file_name+".pdf(","pdf");
    c1->Clear();


    // =====================================================================
    // try to overlay the th2d directly on top of each other
    
    int numCuts = trees[0]->GetEntries(); // Number of lund cuts

    for (int iCut = 0; iCut < numCuts; ++iCut) {

        c1->cd();

        for (int iTree = 0; iTree < trees.size(); ++iTree) { // To inclusive lund planes from TTrees
            trees[iTree]->GetEntry(iCut);
            //hists[iTree]->Rebin2D(8,8);
            //hists[iTree]->SetFillStyle(1001);
            hists[iTree]->SetLineColor(kBlack);

        }

        hists[0]->SetFillColorAlpha(kRed, .5);
        hists[0]->Draw("BOX");

        hists[1]->SetFillColorAlpha(kGreen, .5);
        hists[1]->Draw("BOX SAME");

        hists[2]->SetFillColorAlpha(kBlue, .5);
        hists[2]->Draw("BOX SAME");

        c1->SetTheta(45);  // Elevation angle
        c1->SetPhi(30);    // Azimuthal angle

        TLegend *legend = new TLegend(.1, .7, .28, .9);
        for (int iHist = 0; iHist < hists.size(); ++iHist) {
            legend->AddEntry(hists[iHist], file_names_to_print[iHist], "f");
        }
        legend->Draw();

        if(iCut == numCuts - 1) c1->Print(output_folder_name+output_file_name + ".pdf)", "pdf");
        else c1->Print(output_folder_name+output_file_name + ".pdf", "pdf");
        c1->Clear();

    }
    

    // =====================================================================
    // For cuts along kt in the inclusive lund plane
    
    TString space = " ";
    
    int numCuts = trees[0]->GetEntries(); // Number of lund cuts

    for (int iCut = 0; iCut < numCuts; ++iCut) {

        for (int iTree = 0; iTree < trees.size(); ++iTree) { // To inclusive lund planes from TTrees
            trees[iTree]->GetEntry(iCut);
        }

        
        // kt cuts
        for (int i_kt = 0; i_kt < kt_iBins.size(); ++i_kt) { 

            TLegend* legend = new TLegend(.1,.7,.28, .9);

            c1->cd();

            for (int iTree = 0; iTree < trees.size(); ++iTree) {
                
                
                TH1D* kt_cut = hists[iTree]->ProjectionX(Form("kt_cut_%d", iTree), kt_iBins[i_kt] - 3, kt_iBins[i_kt] + 3, "e");
                legend->AddEntry(kt_cut, file_names[iTree], "l");
                
                if (iTree==0) {
                    kt_cut->SetTitle(space + hists[iTree]->GetTitle() + ", ln(kt) =" + space + kt_cut_params[i_kt] + "; ln(R=.4/delta); weighted counts");
                    kt_cut->Draw("HIST" "PLC");
                }
                else kt_cut->Draw("HIST" "SAME" "PLC");
            
            }

            legend->Draw();
            c1->Print(output_folder_name+output_file_name + ".pdf", "pdf");
            c1->Clear();

        }

        // delta cuts
        for (int iDelta = 0; iDelta < delta_iBins.size(); ++iDelta) { 

            TLegend* legend = new TLegend(.1,.7,.28,.9);

            c1->cd();
            
            for (int iTree = 0; iTree < trees.size(); ++iTree) {

                TH1D* delta_cut = hists[iTree]->ProjectionY(Form("delta_cut_%d",iTree), delta_iBins[iDelta] - 3, delta_iBins[iDelta] + 3, "e");
                legend->AddEntry(delta_cut, file_names[iTree], "l");

                if (iTree == 0) {
                    delta_cut->SetTitle(space + hists[iTree]->GetTitle() + ", ln(.4/delta) =" + space + delta_cut_params[iDelta] + "; ln(kt); weighted counts");
                    delta_cut->Draw("HIST" "PLC");
                }
                else delta_cut->Draw("HIST" "SAME" "PLC");
            
            }

            legend->Draw();
            if ((iDelta == delta_iBins.size() - 1) && (iCut == numCuts - 1)) c1->Print(output_folder_name+output_file_name + ".pdf)", "pdf");
            else c1->Print(output_folder_name + output_file_name + ".pdf", "pdf");
            c1->Clear();

        }

    }
    

    
    return 0;
}
