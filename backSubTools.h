#ifndef BACKSUBTOOLS_H
#define BACKSUBTOOLS_H

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"

#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "fastjet/contrib/ConstituentSubtractor.hh" // jet by jet constituent-based subtraction
#include "fastjet/contrib/SoftKiller.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "algorithm"

#include "simpleTools.h"

using namespace std;


struct bg_sub_options {
    TString name;
    TString event_sub = "null";
    TString jet_sub = "null";
    TString groom_option = "null";
    TH1F* weighted_jet_pt_hist = nullptr;
    double bin_weight;
    bool embed_in_bg = false;
};


//====================================================================
// Background subtraction objects and methods

// Takes a CA clustered jet, returns constituents after applying recursive softdrop. 
// from https://arxiv.org/pdf/1804.03657
vector<fastjet::PseudoJet> recursive_soft_drop_constit(const fastjet::PseudoJet& jet, double Rparam, int N = 5, double z_cut = .2, double beta = 0){

    vector<fastjet::PseudoJet> branches;
    double n = 0;
    branches.push_back(jet);

    while (true) {

        fastjet::PseudoJet j1,j2;

        // "2. Take the remaining branch whose two parent subjets have the widest separation in delta R"
        // Find greatest delta:
        double widest_delta = Rparam;
        int iWidest = -1;
        for (int iBranch = 0; iBranch < branches.size(); ++iBranch) {
            if (branches[iBranch].has_parents(j1, j2)) {
                double current_delta = j1.delta_R(j2);
                if (current_delta > widest_delta) {
                    widest_delta = current_delta;
                    iWidest = iBranch;
                }
            }
            else {
                continue;
            }
        }

        if (iWidest == -1 || (not branches[iWidest].has_parents(j1, j2))) break; // "4. Iterate this process until ... the C/A tree has been fully recursed throug"
        // Make another call to .has_parents to correctly update j1 and j2


        // "Remove that branch from the list of branches"
        branches.erase(branches.begin() + iWidest);

        // "3. If the two subjets pass the SD condition in Eq. (2.1), keep both subjets as new branches; 
        // otherwise, remove the softer of the two subjets and keep the hardest as a new branch."
        double z_12 = min(j1.pt(), j2.pt())/(j1.pt() + j2.pt());
        if (z_12 > z_cut * pow((widest_delta/Rparam), beta)) {
            branches.push_back(j1);
            branches.push_back(j2);
            ++n;
            if (n >= N) break; // "4. Iterate until the SD condition has been met N times"
        }
        else {
            branches.push_back((j1.pt() > j2.pt())? j1 : j2);
        }
    } // end while loop

    // "4. The groomed jet is made of all the remaining branches"
    vector<fastjet::PseudoJet> constituents;
    for (fastjet::PseudoJet branch : branches) {
        vector<fastjet::PseudoJet> branch_constituents = branch.constituents();
        constituents.insert(constituents.end(), branch_constituents.begin(), branch_constituents.end());
    }

    return constituents;
    
}



class RhoEstimator {
    public:

    fastjet::JetDefinition jet_def;
    fastjet::AreaDefinition area_def;
    fastjet::Selector jet_eta_selector;
    //double jet_eta_max;
    //double part_eta_max;

    // constructor
    RhoEstimator(double Rparam, double jet_eta_max, double part_eta_max):
    area_def(fastjet::active_area, fastjet::GhostedAreaSpec(part_eta_max)),
    jet_def(fastjet::kt_algorithm, Rparam, fastjet::E_scheme, fastjet::Best)
    {
    jet_eta_selector = fastjet::SelectorEtaRange(-1 * jet_eta_max, jet_eta_max);
    //jet_eta_max = jet_eta_max_in;
    //part_eta_max = part_eta_max_in;
    }

    double rho(vector<fastjet::PseudoJet> particles) {

        fastjet::ClusterSequenceArea clust_seq(particles, jet_def, area_def);

        vector<fastjet::PseudoJet> kt_jets = sorted_by_pt(jet_eta_selector(clust_seq.inclusive_jets()));

        // remove the two leading jets
        kt_jets.erase(kt_jets.begin());
        kt_jets.erase(kt_jets.begin());

        vector<double> rhos;

        for (fastjet::PseudoJet jet : kt_jets) {
            rhos.push_back(jet.pt()/jet.area());
        }

        return median(rhos);
    }
};






// Try to use perpendicular cone to measure multiplicities of particles with certain pt, then match those particles by pt to constituents of the real jet
class my_pc_subtractor {
  
    public:

    // constructor - Rparam so pc radius matches jet radius.
    my_pc_subtractor(double Rparam_in, double part_eta_max_in, fastjet::JetAlgorithm algo = fastjet::cambridge_algorithm):
    r_selector(fastjet::SelectorCircle(Rparam_in)),
    jet_recluster_def(algo, Rparam_in + .3, fastjet::E_scheme, fastjet::Best),
    area_def(fastjet::active_area, fastjet::GhostedAreaSpec(part_eta_max))
    {
        max_pt_diff = 100; // We're not using this currently
        Rparam = Rparam_in;
        cone_area = M_PI * pow(Rparam_in, 2);
        part_eta_max = part_eta_max_in;
        //r_recluster = r_recluster_in;

        set_pc_cfs();
    }

    fastjet::Selector r_selector;
    fastjet::JetDefinition jet_recluster_def;
    double max_pt_diff;
    double Rparam;
    double cone_area;
    double part_eta_max;
    //double r_recluster;

    fastjet::AreaDefinition area_def;

    // For using correction factor
    //TString cf_file_name;
    TFile* cf_file = nullptr;
    TTree* cf_tree = nullptr;
    TH1F* cf_hist = nullptr;

    void set_pc_cfs(TString cf_file_name = "PC_Correction_Factors.root") {
        
        if (cf_file != nullptr && cf_file->IsOpen()) cf_file->Close();

        cf_file = TFile::Open(cf_file_name, "READ");
        if (cf_file->IsZombie()) {
            cout << "Error: Could not find PC CF Root file " << cf_file_name << endl;
            return;
        }
        cf_tree = (TTree*)cf_file->Get("pc_cf_tree");
        if (!cf_tree) {
            cout << "Error: Could not find pc_cf_tree in ROOT file " << cf_file_name << endl;
            return;
        }
        cf_tree->SetBranchAddress("cf_hist", &cf_hist);
        cf_tree->GetEntry(0); // Only the 1 entry;

    }

    void close_cf_file() {
        cf_file->Close();
    }


    // Takes in a jet to apply perpendicular cone subtraction to, and all particles in the event
    // Returns that jet's constituents (vector<PseudoJet>) after applying perpendicular cone subtraction.
    vector<fastjet::PseudoJet> kt_subtract_constit(const fastjet::PseudoJet& jet, vector<fastjet::PseudoJet> particles) {

        fastjet::PseudoJet pc_axis;

        pc_axis.reset_PtYPhiM(0,jet.eta(), fmod(jet.phi() + M_PI/2 , 2*M_PI));
            
        // selector to find all particles with Rparam of leading jet
        //fastjet::Selector r_selector = fastjet::SelectorCircle(Rparam);
        r_selector.set_reference(pc_axis);
            
        vector<fastjet::PseudoJet> pc_tracks = r_selector(particles);
        vector<fastjet::PseudoJet> constituents = jet.constituents();

        for (fastjet::PseudoJet track : pc_tracks) {
            
            double track_pt = track.pt();

            // Subtract track in leading jet with closest pt. 
            int i_match = -1;
            double pt_diff = max_pt_diff;
            for (int iCon = 0; iCon < constituents.size(); ++iCon) {
                
                double temp_diff = abs(constituents[iCon].pt() - track_pt);

                if (temp_diff < pt_diff) {
                    pt_diff = temp_diff;
                    i_match = iCon;
                }
            }

            // if we found a match, remove from jet
            if (i_match != -1) {
                constituents.erase(constituents.begin() + i_match);
            }

        }

        return constituents;

    }

    // Subtracts as kt_subtract_constit(), but then reclusters into single pseduojet
    fastjet::PseudoJet kt_subtract(const fastjet::PseudoJet& jet, vector<fastjet::PseudoJet> particles) {
        vector<fastjet::PseudoJet> constituents = kt_subtract_constit(jet, particles);
        vector<fastjet::PseudoJet> pc_reclust_jets = fastjet::ClusterSequence(constituents, jet_recluster_def).inclusive_jets(0);
        if (pc_reclust_jets.size() == 0) {
            //cout << "Warning: pc-subtracted jet reclustered into 0 jets" << endl; //Some jets recluster to 0 jets
            return fastjet::PseudoJet(0,0,0,0);
        }
        else if (pc_reclust_jets.size() > 1) {
            //cout << "Warning: pc subtracted jet reclustered into " << pc_reclust_jets.size() << "jets." << endl; // Some jets recluster into more than 1 if R_recluster = Rparam
        }
        
        return pc_reclust_jets[0];

    }

    vector<fastjet::PseudoJet> get_pc_particles(const fastjet::PseudoJet& jet, vector<fastjet::PseudoJet> particles) {
        fastjet::PseudoJet pc_axis;
        pc_axis.reset_PtYPhiM(0,jet.eta(), fmod(jet.phi() + M_PI/2 , 2*M_PI));
        r_selector.set_reference(pc_axis);    
        return r_selector(particles);
    }

    // This one matches pc particles geometrically with event particles 
    // Takes in a jet to apply perpendicular cone subtraction to, and all particles in the event
    // Returns that jet's constituents (vector<PseudoJet>) after applying perpendicular cone subtraction.
    vector<fastjet::PseudoJet> geometric_subtract_constit(const fastjet::PseudoJet& jet, vector<fastjet::PseudoJet> particles) {        

        fastjet::PseudoJet pc_axis;

        pc_axis.reset_PtYPhiM(0,jet.eta(), fmod(jet.phi() + M_PI/2 , 2*M_PI));
            
        // selector to find all particles with Rparam of leading jet
        //fastjet::Selector r_selector = fastjet::SelectorCircle(Rparam);
        r_selector.set_reference(pc_axis);
            
        vector<fastjet::PseudoJet> pc_tracks = r_selector(particles);
        for (fastjet::PseudoJet& track : pc_tracks) { // make a set of "ghost particles" by shifting pc tracks into original jet cone
            track.reset_PtYPhiM(track.pt(), track.eta(), fmod(track.phi() - M_PI/2 , 2*M_PI), track.m());
        }

        vector<fastjet::PseudoJet> constituents = jet.constituents();

        for (fastjet::PseudoJet track : pc_tracks) {
            while (true) {
                double smallest_delta = Rparam;
                
                // Find jet constituent closest to ghost particles 
                int i_match = -1;
                for (int iCon = 0; iCon < constituents.size(); ++iCon) {
                    
                    double temp_delta = track.delta_R(constituents[iCon]);

                    if (temp_delta < smallest_delta) {
                        smallest_delta = temp_delta;
                        i_match = iCon;
                    }
                }

                if (i_match == -1) break; // move to next track
                
                if (constituents[i_match].pt() - track.pt() >= 0) {
                    constituents[i_match] -= track;
                    break;
                }
                else {
                    track -= constituents[i_match];
                    constituents.erase(constituents.begin() + i_match);
                    continue; // Match this subtracted track to another jet constituent
                }

            }

        }

        return constituents;
    }


    // jet_subtraction code is PC_with_cf
    // Same as above, but takes particles from both pcs. Halves the kt of each ghost to account for the second cone. Applies correction factors
    vector<fastjet::PseudoJet> sub_with_cf_constit(const fastjet::PseudoJet& jet, vector<fastjet::PseudoJet> particles) {        

        //cout << "Jet area from sub_with_cf: " << jet.area() << endl;

        fastjet::PseudoJet pc_axis1, pc_axis2;

        pc_axis1.reset_PtYPhiM(0,jet.eta(), fmod(jet.phi() + M_PI/2 , 2*M_PI));
        pc_axis2.reset_PtYPhiM(0,jet.eta(), fmod(jet.phi() - M_PI/2 , 2*M_PI));

            
        // selector to find all particles with Rparam of pc axis
        //fastjet::Selector r_selector = fastjet::SelectorCircle(Rparam);
        vector<fastjet::PseudoJet> axes = {pc_axis1, pc_axis2};
        vector<fastjet::PseudoJet> all_pc_tracks;
        for (int i_axis = 0; i_axis < axes.size(); i_axis++) {
            r_selector.set_reference(axes[i_axis]);
            vector<fastjet::PseudoJet> pc_tracks = r_selector(particles);
            for (fastjet::PseudoJet& track : pc_tracks) { // make a set of "ghost particles" by shifting pc tracks into original jet cone
                double ghost_pt = track.pt() * (.5) * (jet.area() / (cone_area)) * cf_hist->GetBinContent(cf_hist->FindBin(track.pt()));
                //double ghost_pt = track.pt() * (.5) * cf_hist->GetBinContent(cf_hist->FindBin(track.pt())); // Not using area correction here
                track.reset_PtYPhiM(ghost_pt, track.eta(), fmod(((i_axis == 0)? track.phi() - M_PI/2 : track.phi() + M_PI/2), 2*M_PI), track.m()); 
            }
            move(pc_tracks.begin(), pc_tracks.end(), back_inserter(all_pc_tracks));
        }
        /*
        r_selector.set_reference(pc_axis1);
        vector<fastjet::PseudoJet> pc1_tracks = r_selector(particles);
        for (fastjet::PseudoJet& track : pc1_tracks) { // make a set of "ghost particles" by shifting pc tracks into original jet cone
            //double ghost_pt = track.pt() * (.5) * (jet.area() / (cone_area)) * cf_hist->GetBinContent(cf_hist->FindBin(track.pt()));
            double ghost_pt = track.pt() * (.5) * cf_hist->GetBinContent(cf_hist->FindBin(track.pt()));
            track.reset_PtYPhiM(ghost_pt, track.eta(), fmod(track.phi() - M_PI/2 , 2*M_PI), track.m()); 
        }

        r_selector.set_reference(pc_axis2);
        vector<fastjet::PseudoJet> pc2_tracks = r_selector(particles);
        for (fastjet::PseudoJet& track : pc2_tracks) { // make a set of "ghost particles" by shifting pc tracks into original jet cone
            double ghost_pt = track.pt() * (.5) * (jet.area() / (cone_area)) * cf_hist->GetBinContent(cf_hist->FindBin(track.pt()));
            // double ghost_pt = track.pt() * (.5) * (jet.area() / (cone_area)) * cf_hist->GetBinContent(cf_hist->FindBin(track.pt()));

            track.reset_PtYPhiM(ghost_pt / 2, track.eta(), fmod(track.phi() + M_PI/2 , 2*M_PI), track.m());
        }
        

        vector<fastjet::PseudoJet> all_pc_tracks;
        move(pc1_tracks.begin(), pc1_tracks.end(), back_inserter(all_pc_tracks));
        move(pc2_tracks.begin(), pc2_tracks.end(), back_inserter(all_pc_tracks));
        */


        vector<fastjet::PseudoJet> constituents = jet.constituents();

        for (fastjet::PseudoJet track : all_pc_tracks) {
            while (true) {
                double smallest_delta = Rparam;
                
                // Find jet constituent closest to ghost particles 
                int i_match = -1;
                for (int iCon = 0; iCon < constituents.size(); ++iCon) {
                    
                    double temp_delta = track.delta_R(constituents[iCon]);

                    if (temp_delta < smallest_delta) {
                        smallest_delta = temp_delta;
                        i_match = iCon;
                    }
                }

                if (i_match == -1) break; // move to next track
                
                if (constituents[i_match].pt() >= track.pt()) {
                    constituents[i_match] -= track;
                    break;
                }
                else {
                    track -= constituents[i_match];
                    constituents.erase(constituents.begin() + i_match);
                    continue; // Match this subtracted track to another jet constituent
                }

            }

        }

        return constituents;
    }

};



//pc_correction_factor() {

//}





//========================================================




#endif