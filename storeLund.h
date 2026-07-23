#ifndef STORELUND_H
#define STORELUND_H

#include "backSubTools.h"


struct lund_kin_vars {
    double delta, kt, jet_pt;
};

struct lund_coords {
    double log_inv_R, log_kt, jet_pt, bin_weight;
    int iBin;
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

    // For area based jet pT subtraction
    fastjet::AreaDefinition area_def;

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
    area_def(fastjet::active_area, fastjet::GhostedAreaSpec(part_eta_max_in)),
    pc_subtractor(Rparam_in, part_eta_max_in),
    pruner(fastjet::cambridge_algorithm, prune_zcut, Rcut_factor),
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

    
    
    // Takes event particles, int to count number of jets, int to record iCut, bg_sub_options.
    // Finds jets with akt. Reclusters those jets using CA.
    // returns a vector of the delta of the kinematic variables of those jets' splittings for the lund plane.
    vector<lund_kin_vars> get_kin_vars(vector<fastjet::PseudoJet> particles, bg_sub_options ops) {

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

        double rho;
        if (ops.jet_sub == "pT_rho_sub") rho = rho_estimator.rho(particles);
        
        
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
            else if ((ops.jet_sub != "pT_rho_sub") && (ops.jet_sub != "Filter") && (ops.jet_sub != "Prune") && (ops.jet_sub != "RSD_contrib")) {
                cout << "Warning: Unknown jet subtraction request. No subtraction will be performed" << endl;
                constituents = jet.constituents();
            }
            else {
                if (debug) cout << "Calling jet.constituents()" << endl;
                constituents = jet.constituents();
            }

            if (constituents.empty()) {
                break; // Breaks are used because jets are in pt order, so once one has no remaining particles after subtraction, we've reacehed the jets that are all/mostly backgroun
            }

            if (debug) cout << "reclusterSeq with new jet constituents" << endl;
            //fastjet::ClusterSequenceArea reclusterSeq(constituents, jetDef_ca_recluster, area_def);
            fastjet::ClusterSequence reclusterSeq(constituents, jetDef_ca_recluster);
            jet = fastjet::sorted_by_pt(reclusterSeq.inclusive_jets())[0];
            
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

            double jet_pt;
            if (ops.jet_sub == "pT_rho_sub") {
                jet_pt = jet.pt() - rho*jet.area();
            }
            else jet_pt =  jet.pt();
            

            if (ops.weighted_jet_pt_hist != nullptr) {
                if(debug) cout << "Filling weighted jet pt hist: Jet pt =  " << jet.pt() << ", bin weight = " << ops.bin_weight << endl;
                ops.weighted_jet_pt_hist->Fill(jet_pt, ops.bin_weight);
            }

            // jets already selected for eta by eta selector

            
            if (debug) cout << "Beginning grooming + declustering" << endl;
            if (ops.groom_option == "SD_contrib_first") {
                jet = sd(jet);
            }

            fastjet::PseudoJet parent1, parent2;
            bool first = true;
            // For each jet, iteratively compare branchings
            while (jet.has_parents(parent1, parent2)) { 
                if (debug) cout << "New declustering branch" << endl;

                // In each case, identify the higher pt parent 
                fastjet::PseudoJet harder  = (parent1.pt() > parent2.pt()) ? parent1 : parent2;
                fastjet::PseudoJet softer  = (parent1.pt() > parent2.pt()) ? parent2 : parent1;
                
                double delta = harder.delta_R(softer);
                double kt = sin(delta) * softer.pt();
                
                if (first && ops.groom_option == "SD_mine_first") {
                    first = false;
                    // If softdrop condition is not met, discard lower pt subjet. Otherwise return kinematic variables and continue
                    if ((softer.pt() / (harder.pt() + softer.pt())) < z_cut * pow(delta/Rparam, beta)) {// soft drop condition
                        jet = harder;
                        continue;
                    }
                }


                if (ops.groom_option == "SD_mine_all") {
                    first = false;
                    // If softdrop condition is not met, discard lower pt subjet. Otherwise return kinematic variables and continue
                    if ((softer.pt() / (harder.pt() + softer.pt())) < z_cut * pow(delta/Rparam, beta)) {// soft drop condition
                        jet = harder;
                        continue;
                    }
                }

                lund_kin_vars vars = {.delta = delta, .kt = kt, .jet_pt = jet_pt};

                all_kinematic_vars.push_back(vars);

                jet = harder;
                first = false;
            } // end while loop

            if (leading_jet_only) break;

        } // end jet loop 
        return all_kinematic_vars;
    } // end decluster() def
};


#endif