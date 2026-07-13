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

using namespace Pythia8;


// ================================================================
// Tools
double median(vector<double> v) {
    sort(v.begin(), v.end());
    if (v.size()%2 == 1) {
        return v[v.size()/2];
    }
    else {
        return .5 * (v[v.size()/2] + v[(v.size()/2)-1]);
    }
}

bool contains_TString(vector<TString> vec, TString x) {
    return (std::find(vec.begin(), vec.end(), x) != vec.end());
}

// return a TString representation of a double, rounded to nearest .1
TString to_string_round(double num) {
    if (num < 1) {
        return "." + to_string((int)(10*num));
    }
    else {
        return to_string((int)num) + "." + to_string((int)(10*fmod(num, 1.0)));
    }
}

void print_bin (std::vector<int> bin) {
    cout << "pt of " << bin[0] << " to " << bin[1] << ", " << bin[2] << " events." << endl;
}

void print_bins (std::vector<vector<int>> bins) {
    for (int iBin = 0; iBin < bins.size(); iBin++) {
        cout << "Bin " << iBin << ": ";
        print_bin(bins[iBin]);
    }
}

void print_bins2code (std::vector<vector<int>> bins) {
    cout << "vector<vector<int>> bins = {" << endl;
    for (int iBin = 0; iBin < bins.size(); iBin++) {
        cout << "{" << bins[iBin][0] << ", " << bins[iBin][1] << ", " << bins[iBin][2] << (iBin == bins.size()-1? "}" : "},") << endl;
    }
    cout << "};" << endl;

}


TString bool2Str(bool b) {return b? "True" : "False";}

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
    my_pc_subtractor(double Rparam_in, fastjet::JetAlgorithm algo = fastjet::cambridge_algorithm):
    r_selector(fastjet::SelectorCircle(Rparam_in)),
    jet_recluster_def(algo, Rparam_in + .3, fastjet::E_scheme, fastjet::Best)
    {
        max_pt_diff = 100; // We're not using this currently
        Rparam = Rparam_in;
        //r_recluster = r_recluster_in;
    }

    fastjet::Selector r_selector;
    fastjet::JetDefinition jet_recluster_def;
    double max_pt_diff;
    double Rparam;
    //double r_recluster;


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


};


/*
Event wide subtraction options
"ConSub" - Contrib's Constituent Subtraction
"ConSubRho" - Contrib's Constituent Subtraction, with manual rho
"SoftKill" - Contrib's SoftKill

Jet subtraction options
"MyPCkT" - My perpendindicular Cone, with only kt match
"MyPCGM" - My perpendindicular Cone, geometric match
"PCParts" - Use PC particles instead of jet particles
"Filter" - filter
"Prune" - Prune
"RSD_mine" - my recursive softdrop
"RSD_contrib" - contrib's softdrop
"Constituents" - return just the jet's constituents
"Recluster_CA" - Recluster according to CA

Recluster options
"SD_mine_first"
"SD_mine_all"
"SD_contrib_first"
*/


struct bg_sub_options {
    TString event_sub = "null";
    TString jet_sub = "null";
    vector<TString> groom_options = {};
    
    
    /*
    bool my_pc = false;
    bool my_pc_geometric_match = false;
    bool pc_particles = false;
    bool filter = false;
    bool prune = false;
    bool constit_sub = false; // Uses rho to generate background ghost particles, who's pt are subtracted from jet particles
    bool constit_sub_manual = false; // Use your own rho for fastjet's constituent subtraction
    bool soft_kill = false; // Removes particles until rho of event = 0
    bool my_rsd = false; // Recursive Soft Drop: Applies softdrop to branches until soft drop condition has been met N times. Makes no difference for R_g because it only matters that the soft drop condition is being applied to the first clustering. 
    bool contrib_rsd = false;
    bool my_softDrop = false; // Apply softDrop cut while declustering CA jet
    bool contrib_softDrop = false;
    */

    TH1F* jet_pt_hist = nullptr;
    TH1F* jet_pt_hist_unweighted = nullptr;
    double bin_weight = 0;
    //double event_pTHat = 0;
    //double nJets_weighted = -1;
};

struct lund_kin_vars {
    double delta, kt;
    int num_cut, num_jet;
};



// Making primary Lund plane from 
// https://hal.science/hal-01851158/document
class LundGroomer {

public:
    
    // General
    double ptMin;
    double ptMax;
    double jet_eta_max;
    double Rparam;
    bool leading_jet_only;
    bool debug = false;

    fastjet::JetDefinition jetDef_akt; // for jet identification
    fastjet::JetDefinition jetDef_ca_recluster; // for reclustering the antikt-identified jets to make Lund plane
    //fastjet::JetDefinition jetDef_kt_recluster(fastjet::kt_algorithm, Rparam, recombScheme, strat); //  ???
    fastjet::Selector jet_eta_selector;

    //int* iCut;

    RhoEstimator rho_estimator;

    // Recording cuts
    bool cuts = false;
    vector<double> ptCutMins;
    vector<double> ptCutMaxs;

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

    // RecursiveSoftDrop - fastjet contrib
    fastjet::contrib::RecursiveSoftDrop rsd;

    fastjet::contrib::SoftDrop sd;
    

    // constructor
    LundGroomer(double ptMin_in, double ptMax_in, bool leading_jet_only_in, double Rparam_in, double jet_eta_max_in, double part_eta_max_in):
    jetDef_akt(
        fastjet::antikt_algorithm, 
        Rparam_in, 
        fastjet::E_scheme, 
        fastjet::Best),
    jetDef_ca_recluster(
        fastjet::cambridge_algorithm, 
        Rparam_in + .4, 
        fastjet::E_scheme, 
        fastjet::Best),
    rho_estimator(Rparam_in, jet_eta_max_in, part_eta_max_in),
    pc_subtractor(Rparam_in),
    pruner(
        fastjet::cambridge_algorithm, 
        prune_zcut, 
        Rcut_factor),
    filter(Rfilt, fastjet::SelectorNHardest(nfilt)),
    bge(part_eta_max_in, .5), // particle eta max, grid spacing
    subtractor(&bge),
    soft_killer(part_eta_max_in, .4), // rapidity max, grid_size
    rsd(beta, z_cut, 5), // beta, z_cut, N (# of iterations)
    sd(beta, z_cut)
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

    
    //void set_cuts(int* iCut_in, vector<double>& ptCutMins_in, vector<double>& ptCutMaxs_in) {
    void set_cuts(vector<double>& ptCutMins_in, vector<double>& ptCutMaxs_in) {

        cuts = true;
        //iCut = iCut_in;
        ptCutMins = ptCutMins_in;
        ptCutMaxs = ptCutMaxs_in;
    }
    /*
    void set_nJets(int& nJets_in) {
        nJets = nJets_in;
        nJets = 0;
    }
    */
    
        // Takes event particles, int to count number of jets, int to record iCut, bg_sub_options.
    // Finds jets with akt. Reclusters those jets using CA.
    // returns a vector of the delta of the kinematic variables of those jets' splittings for the lund plane.
    vector<lund_kin_vars> get_kin_vars(vector<fastjet::PseudoJet> particles, bg_sub_options ops = bg_sub_options{}) {

        vector<lund_kin_vars> all_kinematic_vars;

        if (debug) cout << "Beginning event wide subtractions" << endl;
        if (ops.event_sub == "ConSub") {
            bge.set_particles(particles);
            particles = subtractor.subtract_event(particles);
            if (particles.empty()) cout << "No remaining particles after contrib constituent subtraction" << endl;
        }
        else if (ops.event_sub == "ConSubRho") {
            bge.set_particles(particles);
            //double rho = rho_estimator.rho();
            fastjet::contrib::ConstituentSubtractor manual_subtractor(bge.rho());
            particles = manual_subtractor.subtract_event(particles);
            if (particles.empty()) {
                cout << "No remaining particles after bigger rho constituent subtraction" << endl;
                return all_kinematic_vars;
            }
        }
        else if (ops.event_sub == "SoftKill") {
            vector<fastjet::PseudoJet> soft_killed_event;
            double pt_thresh = 0; // pt threshold for killed particles
            soft_killer.apply(particles, soft_killed_event, pt_thresh);
            particles = soft_killed_event;
        }
        else if (ops.event_sub == "null") {;}
        else cout << "Warning: Unknown event subtraction request" << endl;

        
        

        
        fastjet::ClusterSequence clust_seq(particles, jetDef_akt);
        vector<fastjet::PseudoJet> jets = sorted_by_pt(jet_eta_selector(clust_seq.inclusive_jets()));

        int jet_counter = -1;

        if (debug) cout << "Looping through akt jets" << endl;
        for (fastjet::PseudoJet& jet : jets) {  

            jet_counter ++;
            if (debug) {
                cout << "Jet num " << jet_counter << endl;
            }

            vector<fastjet::PseudoJet> constituents;
            
            // First, things that act based on the akt jet, and return constituents
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
            else {
                constituents = jet.constituents();
            }

            if (constituents.empty()) {
                break; // Breaks are used because jets are in pt order, so once one has no remaining particles after subtraction, we've reacehed the jets that are all/mostly backgroun
            }

            fastjet::ClusterSequence reclusterSeq(constituents, jetDef_ca_recluster);
            jet = sorted_by_pt(reclusterSeq.inclusive_jets())[0];
            
            // Next, the things that work only on an already CA reclustered jet, and return a new jet
            if  (ops.jet_sub == "Filter") {
                jet = filter(jet);
            }
            else if (ops.jet_sub == "Prune") {
                jet = pruner(jet);
            }
            else if (ops.jet_sub == "RSD_contrib") {
                jet = rsd(jet);
            }
            
            
            
            /*
            // First, things that act based on the akt jet, and return constituents
            if (ops.jet_sub == "MyPCkT") {
                vector<fastjet::PseudoJet> constituents = pc_subtractor.kt_subtract_constit(jet, particles);
                if (debug && constituents.empty()) cout << "No remaining constituents after my pc subtraction" << endl;
                if (constituents.empty()) break;
                fastjet::ClusterSequence reclusterSeq(constituents, jetDef_ca_recluster);
                jet = sorted_by_pt(reclusterSeq.inclusive_jets())[0];
            }
            else if (ops.jet_sub == "MyPCGM") {
                vector<fastjet::PseudoJet> constituents = pc_subtractor.geometric_subtract_constit(jet, particles); 
                if(debug && constituents.empty()) cout << "No remaining constituents after my pc geometric match subtraction" << endl;
                if (constituents.empty()) break;
                fastjet::ClusterSequence reclusterSeq(constituents, jetDef_ca_recluster);
                jet = sorted_by_pt(reclusterSeq.inclusive_jets())[0];
            }
            else if (ops.jet_sub =="PC_parts") {
                vector<fastjet::PseudoJet> constituents = pc_subtractor.get_pc_particles(jet, particles);
                if(debug && constituents.empty()) cout << "No remaining constituents in perpendicular cone" << endl;
                if (constituents.empty()) break;
                fastjet::ClusterSequence reclusterSeq(constituents, jetDef_ca_recluster);
                jet = sorted_by_pt(reclusterSeq.inclusive_jets())[0];
            }
            // These next ones need to be reclusterd with CA before the jet subtraction technique is applied
            else if  (ops.jet_sub == "Filter") {
                //vector<fastjet::PseudoJet> constituents = jet.constituents();
                fastjet::ClusterSequence reclusterSeq(jet.constituents(), jetDef_ca_recluster);
                jet = sorted_by_pt(reclusterSeq.inclusive_jets())[0];
                jet = filter(jet);
            }
            else if (ops.jet_sub == "Prune") {
                fastjet::ClusterSequence reclusterSeq(jet.constituents(), jetDef_ca_recluster);
                jet = sorted_by_pt(reclusterSeq.inclusive_jets())[0];
                jet = pruner(jet);
            }
            else if (ops.jet_sub == "RSD_contrib") {
                fastjet::ClusterSequence reclusterSeq(jet.constituents(), jetDef_ca_recluster);
                jet = sorted_by_pt(reclusterSeq.inclusive_jets())[0];
                jet = rsd(jet);
            }
            // My recursive softdrop needs to be clustered both before and after subtraction
            else if (ops.jet_sub == "RSD_mine") {
                fastjet::ClusterSequence reclusterSeq1(jet.constituents(), jetDef_ca_recluster);
                jet = sorted_by_pt(reclusterSeq1.inclusive_jets())[0];
                vector<fastjet::PseudoJet> constituents = recursive_soft_drop_constit(jet, Rparam);
                if (constituents.empty()) break;
                fastjet::ClusterSequence reclusterSeq2(constituents, jetDef_ca_recluster);
                jet = sorted_by_pt(reclusterSeq2.inclusive_jets())[0];
            }
            // Default is to just recluster by CA
            else if (ops.jet_sub == "null") {
                if (debug) cout << "ops.jet_sub == null: Reclustering jet with CA" << endl;
                fastjet::ClusterSequence reclusterSeq(jet.constituents(), jetDef_ca_recluster);
                jet = sorted_by_pt(reclusterSeq.inclusive_jets())[0];
            }
            else cout << "Warning: Unknown jet subtraction request" << endl;
            */
            


            if (debug) cout << "Reclustered jet pt" << endl;
            if ((jet.pt() < ptMin) || (jet.pt() > ptMax)) continue; // jets are in pt order
            // already selected for eta by jet_eta_selector

            int iCut;
            if (cuts) { // find which cut this jet is in
                //int iCut_temp = -1;
                for (int i = 0; i < ptCutMins.size(); ++i) {
                    if ((jet.pt() > ptCutMins[i]) && (jet.pt() < ptCutMaxs[i])) {
                        //iCut_temp = i;
                        iCut = i;
                        break;
                    }
                }
                //if (iCut_temp != -1) *iCut = iCut_temp;
            }

            // Tracking #jets with weighted pt
            if (ops.jet_pt_hist != nullptr && ops.bin_weight != 0) {
                ops.jet_pt_hist->Fill(jet.pt(), ops.bin_weight);
                //ops.jet_pt_hist_unweighted->Fill(jet.pt());
            }
            
            if (debug) cout << "Beginning grooming + declustering" << endl;
            if (contains_TString(ops.groom_options, "SD_contrib_first")) {
                jet = sd(jet);
            }

            fastjet::PseudoJet parent1, parent2;
            bool first = true;
            // For each jet, iteratively compare branchings
            while (jet.has_parents(parent1, parent2)) { 
                
                // In each case, identify the higher pt parent 
                fastjet::PseudoJet harder  = (parent1.pt() > parent2.pt()) ? parent1 : parent2;
                fastjet::PseudoJet softer  = (parent1.pt() > parent2.pt()) ? parent2 : parent1;
                
                //if (harder.pt() > 100) break; // Event 3992 in Events9x1000 seems to be the problem

                double delta = harder.delta_R(softer);
                double kt = sin(delta) * softer.pt();
                
                if (first && contains_TString(ops.groom_options, "SD_mine_first")) {
                    first = false;
                    // If softdrop condition is not met, discard lower pt subjet. Otherwise return kinematic variables and continue
                    if ((softer.pt() / (harder.pt() + softer.pt())) < z_cut * pow(delta/Rparam, beta)) {// soft drop condition
                        jet = harder;
                        continue;
                    }
                }


                if (contains_TString(ops.groom_options, "SD_mine_all")) {
                    first = false;
                    // If softdrop condition is not met, discard lower pt subjet. Otherwise return kinematic variables and continue
                    if ((softer.pt() / (harder.pt() + softer.pt())) < z_cut * pow(delta/Rparam, beta)) {// soft drop condition
                        jet = harder;
                        continue;
                    }
                }

                lund_kin_vars vars = {.delta = delta, .kt = kt};
                vars.num_jet = jet_counter;

                //vector<double> vars = {delta, kt};
                //if (cuts) vars.push_back(iCut);
                if (cuts) vars.num_cut = iCut;

                all_kinematic_vars.push_back(vars);

                jet = harder;
                first = false;
            } // end while loop

            if (leading_jet_only) break;

        } // end jet loop 
        return all_kinematic_vars;
    } // end decluster() def
};





//pc_correction_factor() {

//}





//========================================================




#endif