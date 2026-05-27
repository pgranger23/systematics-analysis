#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Conventions/Constants.h" //for calculating event kinematics

Int_t fIsCC, fPdg;
Float_t fLep_Px, fLep_Py, fLep_Pz, fLep_E;
Float_t fNu_Px, fNu_Py, fNu_Pz, fNu_E, fR_nucleus;
Float_t fHitNuc_Px, fHitNuc_Py, fHitNuc_Pz, fHitNuc_E, fRemovalE, Q2, y, v, q0, q3;
std::string fProc;

std::vector<int> fPreFSI_Pdg;
std::vector<int> fPreFSI_Index, fPreFSI_Mother1, fPreFSI_Mother2, fPreFSI_Daughter1, fPreFSI_Daughter2;
std::vector<float> fPreFSI_Px, fPreFSI_Py, fPreFSI_Pz, fPreFSI_E;

std::vector<int> fPostFSI_Pdg;
std::vector<int> fPostFSI_Index, fPostFSI_Mother1, fPostFSI_Mother2, fPostFSI_Daughter1, fPostFSI_Daughter2;
std::vector<float> fPostFSI_Px, fPostFSI_Py, fPostFSI_Pz, fPostFSI_E;

Int_t fNPiPlus_Pre, fNPiMinus_Pre, fNPi0_Pre;
Int_t fNPiPlus_Post, fNPiMinus_Post, fNPi0_Post;
Int_t fNCex, fNAbs, fNInel, fNPiProd;

TTree* setup_otree(){
    TTree *otree = new TTree("genie_dump", "GENIE event tree");
    otree->Branch("isCC", &fIsCC, "isCC/I");
    otree->Branch("pdg", &fPdg, "pdg/I");
    otree->Branch("lep_Px", &fLep_Px, "lep_Px/F");
    otree->Branch("lep_Py", &fLep_Py, "lep_Py/F");
    otree->Branch("lep_Pz", &fLep_Pz, "lep_Pz/F");
    otree->Branch("lep_E", &fLep_E, "lep_E/F");
    otree->Branch("nu_Px", &fNu_Px, "nu_Px/F");
    otree->Branch("nu_Py", &fNu_Py, "nu_Py/F");
    otree->Branch("nu_Pz", &fNu_Pz, "nu_Pz/F");
    otree->Branch("nu_E", &fNu_E, "nu_E/F");
    otree->Branch("hitnuc_Px", &fHitNuc_Px, "hitnuc_Px/F");
    otree->Branch("hitnuc_Py", &fHitNuc_Py, "hitnuc_Py/F");
    otree->Branch("hitnuc_Pz", &fHitNuc_Pz, "hitnuc_Pz/F");
    otree->Branch("hitnuc_E", &fHitNuc_E, "hitnuc_E/F");
    otree->Branch("removalE", &fRemovalE, "removalE/F");
    otree->Branch("R_nucleus", &fR_nucleus, "R_nucleus/F");
    otree->Branch("Q2", &Q2, "Q2/F");
    otree->Branch("y", &y, "y/F");
    otree->Branch("v", &v, "v/F");
    otree->Branch("q0", &q0, "q0/F");
    otree->Branch("q3", &q3, "q3/F");

    otree->Branch("preFSI_Pdg", &fPreFSI_Pdg);
    otree->Branch("preFSI_Index", &fPreFSI_Index);
    otree->Branch("preFSI_Mother1", &fPreFSI_Mother1);
    otree->Branch("preFSI_Mother2", &fPreFSI_Mother2);
    otree->Branch("preFSI_Daughter1", &fPreFSI_Daughter1);
    otree->Branch("preFSI_Daughter2", &fPreFSI_Daughter2);
    otree->Branch("preFSI_Px", &fPreFSI_Px);
    otree->Branch("preFSI_Py", &fPreFSI_Py);
    otree->Branch("preFSI_Pz", &fPreFSI_Pz);
    otree->Branch("preFSI_E", &fPreFSI_E);
    otree->Branch("postFSI_Pdg", &fPostFSI_Pdg);
    otree->Branch("postFSI_Index", &fPostFSI_Index);
    otree->Branch("postFSI_Mother1", &fPostFSI_Mother1);
    otree->Branch("postFSI_Mother2", &fPostFSI_Mother2);
    otree->Branch("postFSI_Daughter1", &fPostFSI_Daughter1);
    otree->Branch("postFSI_Daughter2", &fPostFSI_Daughter2);
    otree->Branch("postFSI_Px", &fPostFSI_Px);
    otree->Branch("postFSI_Py", &fPostFSI_Py);
    otree->Branch("postFSI_Pz", &fPostFSI_Pz);
    otree->Branch("postFSI_E", &fPostFSI_E);

    otree->Branch("npiplus_pre", &fNPiPlus_Pre, "npiplus_pre/I");
    otree->Branch("npiminus_pre", &fNPiMinus_Pre, "npiminus_pre/I");
    otree->Branch("npi0_pre", &fNPi0_Pre, "npi0_pre/I");
    otree->Branch("npiplus_post", &fNPiPlus_Post, "npiplus_post/I");
    otree->Branch("npiminus_post", &fNPiMinus_Post, "npiminus_post/I");
    otree->Branch("npi0_post", &fNPi0_Post, "npi0_post/I");
    otree->Branch("ncex", &fNCex, "ncex/I");
    otree->Branch("nabs", &fNAbs, "nabs/I");
    otree->Branch("ninel", &fNInel, "ninel/I");
    otree->Branch("npiprod", &fNPiProd, "npiprod/I");

    otree->Branch("proc", &fProc);
    return otree;
}

void extract_genie(std::string ifilename, std::string ofilename){
    TFile *ifile = new TFile(ifilename.c_str(), "READ");
    TFile *ofile = new TFile(ofilename.c_str(), "RECREATE");

    TTree *otree = setup_otree();
    
    TTree * genie_tree;
    ifile->GetObject("genieEvt", genie_tree);
    if(genie_tree == nullptr){
        std::cerr << "Error: Could not find the genieEvt tree in the input file." << std::endl;
        return;
    }

    genie::NtpMCEventRecord *fEventRecord = nullptr;
    genie_tree->SetBranchAddress("genie_record", &fEventRecord);

    uint nentries = genie_tree->GetEntries();
    uint pct = nentries/100;

    for (int i = 0; i < nentries; i++){
        if(i % pct == 0){
            std::cout << "Processing event " << i << " of " << nentries << std::endl;
        }
        genie_tree->GetEntry(i);
        const genie::Interaction *inter = fEventRecord->event->Summary();
        // get the different components making up the interaction
        const genie::InitialState &initState  = inter->InitState();
        const genie::ProcessInfo  &procInfo   = inter->ProcInfo();

        fIsCC = !procInfo.IsWeakNC();
        fProc = procInfo.AsString();

        const TLorentzVector v4_null;
        genie::GHepParticle* probe = fEventRecord->event->Probe();
        genie::GHepParticle* finallepton = fEventRecord->event->FinalStatePrimaryLepton();
        genie::GHepParticle * hitnucl = fEventRecord->event->HitNucleon();
        const TLorentzVector & k1 = ( probe ? *(probe->P4()) : v4_null );
        const TLorentzVector & k2 = ( finallepton ? *(finallepton->P4()) : v4_null );
        const TLorentzVector & k3 = ( hitnucl ? *(hitnucl->P4()) : v4_null );

        fPdg = probe->Pdg();

        fLep_Px = k2.Px();
        fLep_Py = k2.Py();
        fLep_Pz = k2.Pz();
        fLep_E = k2.E();
        fNu_Px = k1.Px();
        fNu_Py = k1.Py();
        fNu_Pz = k1.Pz();
        fNu_E = k1.E();
        fHitNuc_Px = k3.Px();
        fHitNuc_Py = k3.Py();
        fHitNuc_Pz = k3.Pz();
        fHitNuc_E = k3.E();
        if(hitnucl)
            fRemovalE = hitnucl->RemovalEnergy();
        else
            fRemovalE = 0;

        TLorentzVector q  = k1-k2;

        Q2 = -q.M2();
        v = q.Energy();
        y = v/k1.E();
        q0 = q.Energy();
        q3 = q.P();

        // Clear vectors and counters
        fPreFSI_Pdg.clear();
        fPreFSI_Index.clear();
        fPreFSI_Mother1.clear(); fPreFSI_Mother2.clear();
        fPreFSI_Daughter1.clear(); fPreFSI_Daughter2.clear();
        fPreFSI_Px.clear(); fPreFSI_Py.clear(); fPreFSI_Pz.clear(); fPreFSI_E.clear();
        
        fPostFSI_Pdg.clear();
        fPostFSI_Index.clear();
        fPostFSI_Mother1.clear(); fPostFSI_Mother2.clear();
        fPostFSI_Daughter1.clear(); fPostFSI_Daughter2.clear();
        fPostFSI_Px.clear(); fPostFSI_Py.clear(); fPostFSI_Pz.clear(); fPostFSI_E.clear();

        fNPiPlus_Pre = 0; fNPiMinus_Pre = 0; fNPi0_Pre = 0;
        fNPiPlus_Post = 0; fNPiMinus_Post = 0; fNPi0_Post = 0;
        fNCex = 0; fNAbs = 0; fNInel = 0; fNPiProd = 0;

        // Fill vectors with pre and post FSI particles
        for(int j = 0; j < fEventRecord->event->GetEntries(); ++j){
            genie::GHepParticle * p = fEventRecord->event->Particle(j);
            if(!p) continue;
            int status = p->Status();
            int pdg = p->Pdg();

            if(status == 12 || status == 13 || status == 14){
                fPreFSI_Pdg.push_back(pdg);
                fPreFSI_Index.push_back(j);
                fPreFSI_Mother1.push_back(p->FirstMother());
                fPreFSI_Mother2.push_back(p->LastMother());
                fPreFSI_Daughter1.push_back(p->FirstDaughter());
                fPreFSI_Daughter2.push_back(p->LastDaughter());
                fPreFSI_Px.push_back(p->Px());
                fPreFSI_Py.push_back(p->Py());
                fPreFSI_Pz.push_back(p->Pz());
                fPreFSI_E.push_back(p->E());

                if(pdg == 211) fNPiPlus_Pre++;
                else if(pdg == -211) fNPiMinus_Pre++;
                else if(pdg == 111) fNPi0_Pre++;

            } else if (status == 1){
                fPostFSI_Pdg.push_back(pdg);
                fPostFSI_Index.push_back(j);
                fPostFSI_Mother1.push_back(p->FirstMother());
                fPostFSI_Mother2.push_back(p->LastMother());
                fPostFSI_Daughter1.push_back(p->FirstDaughter());
                fPostFSI_Daughter2.push_back(p->LastDaughter());
                fPostFSI_Px.push_back(p->Px());
                fPostFSI_Py.push_back(p->Py());
                fPostFSI_Pz.push_back(p->Pz());
                fPostFSI_E.push_back(p->E());

                if(pdg == 211) fNPiPlus_Post++;
                else if(pdg == -211) fNPiMinus_Post++;
                else if(pdg == 111) fNPi0_Post++;
            }
        }

        // Simple FSI categorization for pions
        // Tracing: for each pre-FSI pion, look at its children
        for(size_t ipre = 0; ipre < fPreFSI_Pdg.size(); ++ipre){
            int pdg = fPreFSI_Pdg[ipre];
            if(abs(pdg) != 211 && pdg != 111) continue;

            int d1 = fPreFSI_Daughter1[ipre];
            int d2 = fPreFSI_Daughter2[ipre];
            
            if(d1 == -1){ 
                // No FSI happened to this particle (it just passed through or was already final)
                continue; 
            }

            bool found_pion = false;
            bool same_pion = false;
            int n_pions_desc = 0;

            for(int id = d1; id <= d2; ++id){
                genie::GHepParticle * desc = fEventRecord->event->Particle(id);
                if(!desc) continue;
                if(abs(desc->Pdg()) == 211 || desc->Pdg() == 111){
                    found_pion = true;
                    n_pions_desc++;
                    if(desc->Pdg() == pdg) same_pion = true;
                }
            }

            if(!found_pion) fNAbs++;
            else {
                if(!same_pion && n_pions_desc == 1) fNCex++;
                else if(same_pion && n_pions_desc == 1) fNInel++;
                else if(n_pions_desc > 1) fNPiProd++;
            }
        }

        otree->Fill();

        delete fEventRecord->event;
    }
    
    ifile->Close();
    ofile->cd();
    otree->Write();
    ofile->Close();
}