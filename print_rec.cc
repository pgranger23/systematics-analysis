#include <iostream>
#include <string>

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Conventions/Constants.h" //for calculating event kinematics

void print_rec(std::string ifilename, int eid=0){
    TFile *ifile = new TFile(ifilename.c_str(), "READ");
    TTree * genie_tree;
    ifile->GetObject("genieEvt", genie_tree);
    if(genie_tree == nullptr){
        std::cerr << "Error: Could not find the genieEvt tree in the input file." << std::endl;
        return;
    }

    genie::NtpMCEventRecord *fEventRecord = nullptr;
    genie_tree->SetBranchAddress("genie_record", &fEventRecord);

    uint nentries = genie_tree->GetEntries();
    genie_tree->GetEntry(eid);
    std::cout << *fEventRecord << std::endl;
    std::cout << fEventRecord->event->HitNucleon()->RemovalEnergy()*1000 << std::endl;
    
    ifile->Close();
}