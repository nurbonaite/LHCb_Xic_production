import ROOT
from ROOT import TChain
from ctypes import c_double
import os

PIDfolder = "/project/bfys/jdevries/cmtuser/LHCb_Xic_production/pidcalib/UraniaDev_v7r0/"
PIDfilename = "PerfHists_Pi_Turbo16_MagDown_BHH_Binning_P_ETA_nTracks_Brunel.root"
PIDfile = ROOT.TFile.Open(PIDfolder + PIDfilename)
PIDhistbeforecut_name ="TotalHist_Pi_MC15TuneV1_ProbNNpi > 0.5_All__Pi_P_Pi_Eta_nTracks_Brunel"
PIDhistbeforecut = PIDfile.Get(PIDhistbeforecut_name)
PIDhistaftercut_name = "PassedHist_Pi_MC15TuneV1_ProbNNpi > 0.5_All__Pi_P_Pi_Eta_nTracks_Brunel"
PIDhistaftercut = PIDfile.Get(PIDhistaftercut_name)
'''
################################ Piplus_P ###################################

X = PIDhistbeforecut.ProjectionX('X')

#########################################
filedir = "/dcache/bfys/jdevries/ntuples/LcAnalysis/ganga/102"
subjobs = next(os.walk(filedir))[1]
filename = "MC_Lc2pKpiTuple_26103090.root"
excludedjobs = []

tree = TChain("tuple_Lc2pKpi/DecayTree")

#############################################
filedir = "/dcache/bfys/jdevries/ntuples/LcAnalysis/ganga/16"
subjobs = next(os.walk(filedir))[1]
filename = "MC_Lc2pKpiTuple_26103092.root"
excludedjobs = []

tree = TChain("tuple_Lc2pKpi/DecayTree")
#######################################


for job in subjobs:
  if not job in excludedjobs :
    print("- Adding subjobs {0}".format(job))
    tree.Add("{0}/{1}/{2}".format(filedir,job,filename))

axesbins = X.GetXaxis().GetXbins()
myhistogram = ROOT.TH1F('myhistogram','P of Piplus_Lc', axesbins.GetSize()-1, axesbins.GetArray())
tree.Draw("piplus_P>>+myhistogram")

XA = PIDhistaftercut.ProjectionX('XA')
ratio = XA.Clone()
ratio.Divide(X)

################################ Piplus_ETA ###################################

Y = PIDhistbeforecut.ProjectionY('Y')

filedir = "/dcache/bfys/jdevries/ntuples/LcAnalysis/ganga/102"
subjobs = next(os.walk(filedir))[1]
filename = "MC_Lc2pKpiTuple_26103090.root"
excludedjobs = []

tree = TChain("tuple_Lc2pKpi/DecayTree")

#############################################
filedir = "/dcache/bfys/jdevries/ntuples/LcAnalysis/ganga/16"
subjobs = next(os.walk(filedir))[1]
filename = "MC_Lc2pKpiTuple_26103092.root"
excludedjobs = []

tree = TChain("tuple_Lc2pKpi/DecayTree")
#######################################

for job in subjobs :
  if not job in excludedjobs :
    print("- Adding subjobs {0}".format(job))
    tree.Add("{0}/{1}/{2}".format(filedir,job,filename))
axesbins = Y.GetXaxis().GetXbins()
myhistogram = ROOT.TH1F('myhistogram','Eta of Piplus_Lc', axesbins.GetSize()-1, axesbins.GetArray())
tree.Draw("piplus_ETA>>+myhistogram")

YA = PIDhistaftercut.ProjectionY('YA')
ratio = YA.Clone()
ratio.Divide(Y)
'''
################################ Piplus_nTracks ###################################

Z = PIDhistbeforecut.ProjectionZ('Z')
'''
filedir = "/dcache/bfys/jdevries/ntuples/LcAnalysis/ganga/102"
subjobs = next(os.walk(filedir))[1]
filename = "MC_Lc2pKpiTuple_26103090.root"
excludedjobs = []

tree = TChain("tuple_Lc2pKpi/DecayTree")
'''
#############################################
filedir = "/dcache/bfys/jdevries/ntuples/LcAnalysis/ganga/16"
subjobs = next(os.walk(filedir))[1]
filename = "MC_Lc2pKpiTuple_26103092.root"
excludedjobs = []

tree = TChain("tuple_Lc2pKpi/DecayTree")
#######################################

for job in subjobs :
  if not job in excludedjobs :
    print("- Adding subjobs {0}".format(job))
    tree.Add("{0}/{1}/{2}".format(filedir,job,filename))

axesbins = Z.GetXaxis().GetXbins()
myhistogram = ROOT.TH1F('myhistogram','nTracks of Piplus_Lc', axesbins.GetSize()-1, axesbins.GetArray())
tree.Draw("nTracks>>+myhistogram")

ZA = PIDhistaftercut.ProjectionZ('ZA')
ratio = ZA.Clone()
ratio.Divide(Z)


efftotal_error = c_double(0.0)
nbinranges = [ axesbins.GetSize() -1 ]
myhistogram.Sumw2()
ntotal = myhistogram.Integral(0, nbinranges[0])
myhistogram.Multiply(ratio)
print(ntotal)
efftotal = myhistogram.IntegralAndError(0, nbinranges[0], efftotal_error, "")
print(efftotal)
ntotal = c_double(ntotal)
efftotal = c_double(efftotal)

efftotal = efftotal.value / ntotal.value
efftotal_error = efftotal_error.value / ntotal.value
print("PID eff Pplus_nTracks  \t= {0:.4f} +- {1:.4f}".format(efftotal, efftotal_error))


