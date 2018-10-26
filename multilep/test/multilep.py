import sys, os
import FWCore.ParameterSet.Config as cms

#function to return JSON file
def getJSON(is2017):
    if is2017: return "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt"
    else: return "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"

# Default arguments
#inputFile       = '/store/mc/RunIISummer16MiniAODv2/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/00A9113F-15D6-E611-9142-047D7B881D3A.root'
#inputFile       = '/store/mc/RunIISummer16MiniAODv2/TTGamma_Dilept_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/90000/003658EE-77E6-E611-ACB1-7CD30ABD295A.root'
#inputFile        = '/store/mc/RunIISummer16MiniAODv2/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/AABE3103-4FD5-E611-91BA-02163E01314F.root'
#inputFile       = '/store/data/Run2016D/DoubleMuon/MINIAOD/03Feb2017-v1/100000/52779EE0-F4ED-E611-BF87-70106F49CD3C.root'
#inputFile       = "root://cmsxrootd.fnal.gov///store/data/Run2017C/MuonEG/MINIAOD/PromptReco-v3/000/300/780/00000/86494C82-EA7E-E711-ACCC-02163E01441B.root"
#inputFile       = 'file:///pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/prompt/HeavyNeutrino_trilepton_M-100_V-0.01_2l_NLO/heavyNeutrino_1.root'
#inputFile       = "root://xrootd-cms.infn.it///store/mc/RunIISummer16MiniAODv2/SMS-TChiWZ_ZToLL_mZMin-0p1_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSummer16Fast_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/18589842-DCBD-E611-B8BF-0025905A48D8.root"
#inputFile       = '/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/005DC030-D3F4-E711-889A-02163E01A62D.root'
#inputFile       = '/store/mc/RunIIFall17MiniAOD/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/001C74A0-B4D6-E711-BD4B-FA163EB4F61D.root'
#inputFile       = '/store/mc/RunIIFall17MiniAOD/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/000DE38C-A903-E811-92F5-B083FED13C9E.root'
#inputFile       = '/store/mc/RunIISummer16MiniAODv2/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v3/70000/06598411-2DC8-E611-AFCC-02163E019D8F.root'
#inputFile       = 'file:///user/ikhvastu/CMSSW_9_4_4/src/heavyNeutrino/multilep/test/pickevents_DY_1_283_583420.root'
#inputFile       = '/store/mc/RunIISummer16MiniAODv2/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/00A25ADE-DFD4-E611-8EAC-0025905A48B2.root'
#inputFile       = '/store/data/Run2017B/DoubleMuon/MINIAOD/17Nov2017-v1/30000/3800AFDD-20D6-E711-87D6-02163E019D17.root'
#inputFile       = '/store/mc/RunIIFall17MiniAOD/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/02A1FA28-F1F9-E711-AFB2-002590E3A0FC.root'
#inputFile       = '/store/data/Run2017F/DoubleMuon/MINIAOD/17Nov2017-v1/50000/009DC3A2-A7DE-E711-99F7-02163E013717.root'
#inputFile       = '/store/mc/RunIIFall17MiniAOD/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/50000/004C666D-C0E0-E711-AADB-0CC47A6C183A.root'
#inputFile       = '/store/mc/RunIIFall17MiniAOD/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/02A1FA28-F1F9-E711-AFB2-002590E3A0FC.root'
#inputFile       = '/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/0695D53A-EEF4-E711-99BD-02163E019C44.root'
#inputFile       = '/store/mc/RunIIFall17MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v2/60000/0857E231-22EF-E711-BEAC-FA163E2DB7A6.root'
#inputFile       = '/store/mc/RunIIFall17MiniAOD/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v2/60000/063DA70C-EADA-E711-8493-008CFAFBF7C8.root'
inputFile       = '/store/data/Run2017F/SingleMuon/MINIAOD/17Nov2017-v1/00000/3E7C07F9-E6F1-E711-841A-0CC47A4C8E46.root'
#inputFile       = '/store/data/Run2016D/DoubleMuon/MINIAOD/03Feb2017-v1/100000/52779EE0-F4ED-E611-BF87-70106F49CD3C.root'
#inputFile       = '/store/mc/RunIISummer17MiniAOD/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v1/90000/FE1118B0-A192-E711-A33C-90B11C064AD8.root'
#inputFile       = 'file:///user/ikhvastu/CMSSW_9_4_4/src/heavyNeutrino/multilep/test/pickevents_data_275309_124_152717842.root'
nEvents         = 1000
outputFile      = 'singlelep.root'   # trilep    --> skim three leptons (basic pt/eta criteria)
                                 # dilep     --> skim two leptons
                                 # singlelep --> skim one lepton
                                 # ttg       --> skim two leptons + one photon
                                 # fakerate  --> not implemented

def getVal(arg):
    return arg.split('=')[-1]

# Loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i]
    if "outputFile"  in sys.argv[i]: outputFile = getVal(sys.argv[i])
    elif "inputFile" in sys.argv[i]: inputFile  = getVal(sys.argv[i])
    elif "events"    in sys.argv[i]: nEvents    = int(getVal(sys.argv[i]))

isData = not ('SIM' in inputFile or 'HeavyNeutrino' in inputFile)
is2017 = "Run2017" in inputFile or "17MiniAOD" in inputFile
isSUSY = "SMS-T" in inputFile

if 'pickevent' in inputFile.lower():
    isData = True
    is2017 = False
    isSUSY = False

print isData, is2017, isSUSY

process = cms.Process("BlackJackAndHookers")

# initialize MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Unscheduled mode makes no difference for us
#process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.source       = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFile.split(",")))
process.options      = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents    = cms.untracked.PSet(input = cms.untracked.int32(nEvents))
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if   isData and is2017: process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v6'
elif is2017:            process.GlobalTag.globaltag = '94X_mc2017_realistic_v13'
elif isData:            process.GlobalTag.globaltag = '94X_dataRun2_v10'
else:                   process.GlobalTag.globaltag = '94X_mcRun2_asymptotic_v3'

#
# Vertex collection
#
process.load('CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi')
process.goodOfflinePrimaryVertices.src    = cms.InputTag('offlineSlimmedPrimaryVertices')
process.goodOfflinePrimaryVertices.filter = cms.bool(False)                          #Don't use any EDFilter when relying on hCounter!

#
# Import some objectsequences sequence (details in cff files)
#
from heavyNeutrino.multilep.jetSequence_cff import addJetSequence
from heavyNeutrino.multilep.egmSequence_cff import addElectronAndPhotonSequence
addJetSequence(process, isData, is2017)
addElectronAndPhotonSequence(process)

#
# Read additional MET filters not stored in miniAOD
#
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
for module in [process.BadPFMuonFilter, process.BadChargedCandidateFilter]:
  module.muons        = cms.InputTag("slimmedMuons")
  module.PFCandidates = cms.InputTag("packedPFCandidates")

process.BadPFMuonFilter.filter = cms.bool(False)
process.BadChargedCandidateFilter.filter = cms.bool(False)

#clean 2016 data met from spurious muons and ECAL slew rate
metCollection = "slimmedMETs"
#if (not is2017) and isData : metCollection = "slimmedMETsMuEGClean" #No longer needed for new rereco

# Main Process
process.blackJackAndHookers = cms.EDAnalyzer('multilep',
  vertices                      = cms.InputTag("goodOfflinePrimaryVertices"),
  genEventInfo                  = cms.InputTag("generator"),
  lheEventInfo                  = cms.InputTag("externalLHEProducer"),
  pileUpInfo                    = cms.InputTag("slimmedAddPileupInfo"),
  genParticles                  = cms.InputTag("prunedGenParticles"),
  muons                         = cms.InputTag("slimmedMuons"),
  muonsEffectiveAreas           = cms.FileInPath('heavyNeutrino/multilep/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
  muonsEffectiveAreasFall17     = cms.FileInPath('heavyNeutrino/multilep/data/effAreas_cone03_Muons_Fall17.txt'),
  electrons                     = cms.InputTag("slimmedElectrons"),
  electronsEffectiveAreas       = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt'), # WARNING this is spring 15, following SUSY-standard, i.e. not the most up-to-date values
  electronsEffectiveAreasFall17 = cms.FileInPath('heavyNeutrino/multilep/data/effAreas_cone03_Electrons_Fall17.txt'),

  electronsMva                  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
  electronsMvaHZZ               = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
  electronMvaFall17Iso          = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values"),
  electronMvaFall17NoIso        = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
  electronsCutBasedVeto         = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
  electronsCutBasedLoose        = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
  electronsCutBasedMedium       = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
  electronsCutBasedTight        = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
  leptonMvaWeightsMuSUSY16      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_SUSY16_BDTG.weights.xml"),
  leptonMvaWeightsEleSUSY16     = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_SUSY16_BDTG.weights.xml"),
  leptonMvaWeightsMuttH16       = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_ttH16_BDTG.weights.xml"),
  leptonMvaWeightsElettH16      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_ttH16_BDTG.weights.xml"),
  leptonMvaWeightsMuSUSY17      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_SUSY17_BDTG.weights.xml"),
  leptonMvaWeightsEleSUSY17     = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_SUSY17_BDTG.weights.xml"),
  leptonMvaWeightsMuttH17       = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_ttH17_BDTG.weights.xml"),
  leptonMvaWeightsElettH17      = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_ttH17_BDTG.weights.xml"),
  leptonMvaWeightsEletZqTTV16   = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_tZqTTV16_BDTG.weights.xml"),
  leptonMvaWeightsMutZqTTV16    = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_tZqTTV16_BDTG.weights.xml"),
  leptonMvaWeightsEletZqTTV17   = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/el_tZqTTV17_BDTG.weights.xml"),
  leptonMvaWeightsMutZqTTV17    = cms.FileInPath("heavyNeutrino/multilep/data/mvaWeights/mu_tZqTTV17_BDTG.weights.xml"),
  JECtxtPath                    = cms.FileInPath("heavyNeutrino/multilep/data/JEC/dummy.txt"),
  #photons                       = cms.InputTag("slimmedPhotons"),
  #photonsChargedEffectiveAreas  = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfChargedHadrons_90percentBased.txt'),
  #photonsNeutralEffectiveAreas  = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased.txt'),
  #photonsPhotonsEffectiveAreas  = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt'),
  #photonsCutBasedLoose          = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose"),
  #photonsCutBasedMedium         = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium"),
  #photonsCutBasedTight          = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight"),
  #photonsMva                    = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
  #photonsChargedIsolation       = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
  #photonsNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
  #photonsPhotonIsolation        = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
  #photonsFull5x5SigmaIEtaIPhi   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIPhi"),
  taus                          = cms.InputTag("slimmedTaus"),
  packedCandidates              = cms.InputTag("packedPFCandidates"),
  rho                           = cms.InputTag("fixedGridRhoFastjetAll"),
  met                           = cms.InputTag(metCollection),
  jets                          = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
  #jetsSmeared                   = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmeared"),
  #jetsSmearedUp                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedUp"),
  #jetsSmearedDown               = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC" if isData else "slimmedJetsCorrectedAndSmearedDown"),
  jecUncertaintyFile16          = cms.FileInPath("heavyNeutrino/multilep/data/JEC/Summer16_07Aug2017_V9_MC_Uncertainty_AK4PFchs.txt"),
  jecUncertaintyFile17          = cms.FileInPath("heavyNeutrino/multilep/data/JEC/Fall17_17Nov2017_V6_MC_Uncertainty_AK4PFchs.txt"),
  prescales                     = cms.InputTag("patTrigger"),
  triggers                      = cms.InputTag("TriggerResults::HLT"),
  recoResultsPrimary            = cms.InputTag("TriggerResults::PAT"),
  recoResultsSecondary          = cms.InputTag("TriggerResults::RECO"),
  badPFMuonFilter               = cms.InputTag("BadPFMuonFilter"),
  badChargedCandFilter          = cms.InputTag("BadChargedCandidateFilter"),
  skim                          = cms.untracked.string(outputFile.split('/')[-1].split('.')[0].split('_')[0]),
  isData                        = cms.untracked.bool(isData),
  is2017                        = cms.untracked.bool(is2017),
  isSUSY                        = cms.untracked.bool(isSUSY),
  jetPtThreshold                = cms.untracked.double(15.0),

)

if isData:
  import FWCore.PythonUtilities.LumiList as LumiList
  process.source.lumisToProcess = LumiList.LumiList(filename = "../data/JSON/" + getJSON(is2017)).getVLuminosityBlockRange()

process.p = cms.Path(process.goodOfflinePrimaryVertices *
                     process.BadPFMuonFilter *
                     process.BadChargedCandidateFilter *
                     process.egmSequence *
                     process.jetSequence *
                     #process.fullPatMetSequence *
                     process.basicJetsForMETModifiedEMthreshold *
                     process.fullPatMetSequenceModifiedEMthreshold *
                     process.blackJackAndHookers)
