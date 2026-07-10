// Read from a root file, make lund plane


#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"


// ROOT - files
//#include <TFile.h>
//#include <TTree.h>
//#include <TBranch.h>
#include <vector>
//#include <iostream>

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

using namespace Pythia8;

// Added eta option for storeWithRoot, hopefully this still works

// Takes a CA clustered jet, returns a vector of kinematic 
// variables (delta, kt) vectors for each primary splitting
vector<vector<double>> decluster(fastjet::PseudoJet jet, bool softDrop = false) {
    
    double z_cut = .05; // z cut > .02
    double beta = 2;
    double Rparam = .4;

    vector<vector<double>> kinematic_vars;

    fastjet::PseudoJet parent1;
    fastjet::PseudoJet parent2;
    // Pseudojets a and b 
    fastjet::PseudoJet* pa = &parent1;
    fastjet::PseudoJet* pb = &parent2;

    // For each jet, iteratively compare branchings
    while (jet.has_parents(parent1, parent2)) { 
        // In each case, identify the higher pt parent 
        if(parent1.pt() > parent2.pt()) {
            pa = &parent1;
            pb = &parent2;
        }
        else {
            pa = &parent2;
            pb = &parent1;
        }

        // Classify pb as the emission, and pa + pb as the emitter
        // and define kinematic variables
        
        // rap-phi distance, delta = delta_ab
        double delta = pa->delta_R(*pb);


        // from https://arxiv.org/pdf/1402.2657
        // In this paper, softDrop stops once the condition is met.
        // Here, I continue applying it to every branching
        if (softDrop) {
            // If softdrop condition is not met, discard lower pt subjet. Otherwise return kinematic variables and continue
            if ((pb->pt() / (pa->pt() + pb->pt())) > z_cut * pow(delta/Rparam, beta)) {// soft drop condition
                jet = *pa;
                continue;
            }
        }


        double kt = sin(delta) * pb->pt();

        vector<double> vars = {delta, kt};

        kinematic_vars.push_back(vars);
    
        jet = *pa;
    }
    return kinematic_vars;
}

double median(vector<double> v) {
    sort(v.begin(), v.end());
    if (v.size()%2 == 1) {
        return v[v.size()/2];
    }
    else {
        return .5 * (v[v.size()/2] + v[(v.size()/2)-1]);
    }
}

TString bool2Str(bool b) {return b? "True" : "False";}

int main() {
    
    TString output_file_name = "9k, makeLund";
    TString output_folder_name = "LundPlanes/";
    //cout << "Decluster() has been modified to only return the first branching!" << endl;

    TString root_output_name = output_file_name;

    TString eventFileName = "events1000x9.root";
    TString backgroundFileName = "thermalBackgrounds9000.root";


    // Minimum jet pT
    double ptMin = 0;
    // Particle eta max;
    double etaMax = .9;
    // jet eta max;
    double jet_eta_max = .5;
    
    double Rparam = 0.4;

    //int toSkip = 1000;
    // For not using all the events in the ROOT file. Will only use 1 of every toSkip events
   
    bool background = false;

    // Area based background subtraction with rho
    bool area_sub = false;
    if (!area_sub) cout << "Remember to not use clustersequenceArea for faster performance" << endl;
    if (area_sub && !background) {
        cout << "Can't do global background estimation and subtraction with no background! Setting area_sub = false";
        area_sub = false;
    } 

    bool soft_drop = false;
    if (soft_drop && !background) {
        cout << "Can't do soft drop with no background! Setting soft_drop = false";
        soft_drop = false;
    } 

    // Note: An outlier found in pt spectra in jetSpectraF from "pythia_events_20000x9.root"
    // is manually removed below. This required bin specification
    // Get bin specifications here from storeWithRoot
    // Make sure recreation of iBin later down also matches
    vector<double> ptHatMins = {10, 20, 30, 50, 70, 90, 110, 130, 150};
    vector<double> ptHatMaxs = {20, 30, 50, 70, 90, 110, 130, 150, 250};
    int iBin = 0;
    double currentBin = ptHatMins[0];

    //==================================================================
    // Preparing to read from ROOT files

    my_event_tree met(eventFileName);
    my_background_tree mbt(backgroundFileName);

    //=======================================================================
    // Graphs 'n stuff
    // Declare 2D ROOT histograms
    TString xAxisName = "ln((R=.4/delta)";
    TString yAxisName = "ln(kt)";

    
    TH2F *lund_inclusive = new TH2F("lund_inclusive", "Inclusive pt; " + xAxisName + "; " + yAxisName, 50, 0, 5, 50, -4, 3.5);
   
    //TH3F *lund_3d = new TH3F("lund", "3D Lund plane; " + xAxisName + "; " + yAxisName + "; jet pT", 50, 0, 5, 50, -4, 3.5, 50, 0, 200);

    // ROOT hists for pt cuts
    vector<double> ptCutMins = {10, 20, 30, 50, 70, 100, 130, 160, 190};
    vector<double> ptCutMaxs = {20, 30, 50, 70, 100, 130, 160, 190, 250};
    vector<TH2F*> lundCuts;
    
    TString space = " ";
    for (int iCut = 0; iCut < ptCutMins.size(); ++iCut) {
        lundCuts.push_back(new TH2F("", "jet pT" + space + ptCutMins[iCut] + " to " + ptCutMaxs[iCut] + "; " + xAxisName + "; " + yAxisName,
            50, 0, 5, 50, -4, 3.5));
    }

   

    //=========================================================================
    // Set up fastjet analysis
    fastjet::Strategy strat = fastjet::Best;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;

    fastjet::JetDefinition jetDef_akt(fastjet::antikt_algorithm, Rparam, recombScheme, strat); // for jet identification
    fastjet::JetDefinition jetDef_kt(fastjet::kt_algorithm, Rparam, recombScheme, strat); // for global background and area estimation
    fastjet::JetDefinition jetDef_ca_recluster(fastjet::cambridge_algorithm, 3, recombScheme, strat); // for reclustering the antikt-identified jets to make Lund plane

    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(etaMax + Rparam)); // for jet area estimation
            
    //=========================================================================
    // Reading events from ROOT


    int numEvents = met.GetEntries();
    //std::cout << "Number of events: " << numEvents << std::endl;

    int numBackgrounds;
    if(background) {
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
        
        if (iEvent > counter) {
            cout << counter << " events analyzed from ROOT file" << endl;
            counter += 1000;
        }

        // Get the current event. This updates all branch addresses
        //met.GetEntry(iEvent);

        vector<fastjet::PseudoJet> event_particles = met.get_particles(iEvent, .9); // iEvent, part_eta_max
        
        // =====================================================================================================
        // Embedding in the background
        if (background) {
            
            //mbt.GetEntry(iEvent);

            // Create a vector to store all new particles
            vector<fastjet::PseudoJet> background_prtcls = mbt.get_particles(iEvent);

            // Append background particle vector to fjInputs
            move(background_prtcls.begin(), background_prtcls.end(), back_inserter(event_particles));

        
        } // End if (background)

        // =====================================================
        // area based subtraction - global background density estimation and subtraction. Only works if background was done
        double rho;

        if (area_sub) {

            fastjet::ClusterSequenceArea clustSeq_area(event_particles, jetDef_kt, areaDef);
        
            vector<fastjet::PseudoJet> kt_jets_bg = sorted_by_pt(clustSeq_area.inclusive_jets());

            // remove the two leading jets
            kt_jets_bg.erase(kt_jets_bg.begin());
            kt_jets_bg.erase(kt_jets_bg.begin());

            vector<double> rhos;

            for (fastjet::PseudoJet jet : kt_jets_bg) {
                rhos.push_back(jet.pt()/jet.area());
            }

            rho = median(rhos);
        }

        // =======================================================================================================================
        // Making primary Lund plane from 
        // https://hal.science/hal-01851158/document
        
        // Only use area cluster sequence if doing area subtraction
        fastjet::ClusterSequenceArea clustSeq(event_particles, jetDef_akt, areaDef);
        //fastjet::ClusterSequence clustSeq(event_particles, jetDef_akt);
        vector<fastjet::PseudoJet> event_jets = clustSeq.inclusive_jets(ptMin);
        
        // Now recluster the jets identified by anti-kt according to CA, add points to lund plane
        bool cut = false;
        int iCut = 0;

        double jet_pt;

        //vector<fastjet::PseudoJet> reclusteredJets;
        for (fastjet::PseudoJet jet : event_jets) {
            
            if (jet.eta() > jet_eta_max) continue;

            // Try to take out that one outlier
            //if (iBin == 0 && jet.pt() > 100) continue;

            cut = false;

            if (area_sub) jet_pt = jet.pt() - rho*jet.area();
            else jet_pt = jet.pt(); 

            // Find which cut this jet belongs to
            for (int i = 0; i < ptCutMins.size(); ++i) {
                                
                if (jet_pt > ptCutMins[i] && jet_pt < ptCutMaxs[i]) {
                    cut = true;
                    iCut = i;
                    break;
                }
            }
            
            fastjet::ClusterSequence reclustered_jet(clustSeq.constituents(jet), jetDef_ca_recluster);
            if (reclustered_jet.inclusive_jets().size() != 1) {
                cout << "Warning, an anti-kt jet reclustered into " << reclustered_jet.inclusive_jets().size() << " CA jets." << endl;
            }
            
            vector<vector<double>> kin_vars = decluster(reclustered_jet.inclusive_jets()[0], soft_drop);
            for (vector<double> vec : kin_vars) {
                lund_inclusive->Fill(log(Rparam/vec[0]),log(vec[1]), met.weight);
                //lund_3d->Fill(log(Rparam/vec[0]),log(vec[1]), jet_pt, weight);
                if (cut) lundCuts[iCut]->Fill(log(Rparam/vec[0]),log(vec[1]), met.weight);
            }


        } // End Lund construction
        // =======================================================================================================================
     
    }  // End reading events from ROOT
    //==========================================================
    

    // =========================================================
    // Save these plots in a root file, so they can be modified/resized without having to run the whole thing again
    TFile output_root_file(output_folder_name + root_output_name + ".root", "RECREATE");

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

    //========================================================
    // Printing histograms

    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    // Adding a text information sheet with stats about this generation
    TString title = "Lund plane with pT cuts: " + output_file_name;
    TString line1 = "Number of events: " + to_string(numEvents);
    TString line2 = "Events from: " + eventFileName;
    
    
    TString line3 = "Embedded in background: " + bool2Str(background);
    TString line4 = "Backgrounds from: " + backgroundFileName;
    TString line5 = "Global Background Estimation + Area Based Subtraction: " + bool2Str(area_sub);
    TString line6 = "Soft Drop: " + bool2Str(soft_drop);

    TString line7 = "Jet pT min: " + to_string(ptMin);
    TString line8 = "Jet radius: " + to_string(Rparam);
    TString line9 = "Particle eta max: " + to_string(etaMax);



    vector<TString> lines = {title, line1, line2, line3, line4, line5, line6, line7, line8, line9};

    c1->cd();
    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }
    
    c1->Print(output_folder_name + output_file_name +".pdf(","pdf");
    c1->Clear();


    c1->cd();
    lund_inclusive->Draw("COLZ");
    c1->Print(output_folder_name + output_file_name +".pdf","pdf");
    c1->Clear();

    for (int iCut = 0; iCut < ptCutMins.size(); ++iCut) {
        
        c1->cd();
        lundCuts[iCut]->Draw("COLZ");
        if (iCut + 1 == ptCutMins.size()) {
            c1->Print(output_folder_name + output_file_name + ".pdf)","pdf");
        } 
        else {
            c1->Print(output_folder_name + output_file_name + ".pdf","pdf");
        }
        c1->Clear();
    }

    // Since ROOT histograms point to the data they store, once the eventRootFile
    // closes they no longer have access to their data
    met.close_file();
    mbt.close_file();

    return 0;
}



