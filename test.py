import uproot
import polars as pl
import awkward as ak
from tqdm import tqdm
import ROOT
ROOT.gROOT.LoadMacro('/home/pgranger/software/Generator/src/scripts/gcint/genie_setup.C')
genie_record = ROOT.MakeNullPointer(ROOT.genie.NtpMCEventRecord)
f = ROOT.TFile('out.root')
t = f.Get('genieEvt')
t.Print()
t.SetBranchAddress('genie_record', ROOT.AddressOf(genie_record))