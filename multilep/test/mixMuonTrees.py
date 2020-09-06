#! /usr/bin/env python

# Original author: Tom Cornelis

#
# Takes the muon POG tag and probe trees
# Mix-in the leptonMva values from the tuples
#   Requires the ttg repository for some helper functions
#

import ROOT
import os,glob,argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel', action='store', default='INFO', help='Log level for logging', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE'])
argParser.add_argument('--jobId', action='store', default=None)
argParser.add_argument('--isChild', action='store_true', default=False, help='mark as subjob, will never submit subjobs by itself')
argParser.add_argument('--dryRun', action='store_true', default=False, help='do not launch subjobs, only show them')
argParser.add_argument('--runLocal', action='store_true', default=False, help='use local resources instead of Cream02')
argParser.add_argument('--subJob', action='store', default=None, help='The xth subjob for a sample')
args = argParser.parse_args()

from ttg.tools.logger import getLogger
log = getLogger(args.logLevel)

#
# Samples to be mixed
#
samplesToMix = []

for i in glob.glob('/pnfs/iihe/cms/store/user/tomc/tnpTuples_muons/POG/MC_Moriond17_DY_tranch4Premix_part*.root'):
  samplesToMix.append((i, '/pnfs/iihe/cms/store/user/tomc/heavyNeutrino/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_MiniAOD2016v3_ext1-v2_muonTuples/200722*/*/*.root'))

samplesToMix.append(('/pnfs/iihe/cms/store/user/tomc/tnpTuples_muons/POG/TnPTree_94X_DYJetsToLL_M50_Madgraph.root', 
    '/pnfs/iihe/cms/store/user/tomc/heavyNeutrino/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1_muonTuples/200718*/*/*.root'))

samplesToMix.append(('/pnfs/iihe/cms/store/user/tomc/tnpTuples_muons/POG/TnPTreeZ_102XAutumn18_DYJetsToLL_M50_MadgraphMLM.root', 
    '/pnfs/iihe/cms/store/user/tomc/heavyNeutrino/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_MiniAOD2018-v1_muonTuples/200718*/*/*.root'))

for run in ['Bver2', 'C', 'D', 'E', 'F', 'G', 'H']:
  samplesToMix.append(('/pnfs/iihe/cms/store/user/tomc/tnpTuples_muons/POG/TnPTreeZ_LegacyRereco07Aug17_SingleMuon_Run2016%s_GoldenJSON.root' % run,
      '/pnfs/iihe/cms/store/user/tomc/heavyNeutrino/SingleMuon/crab_Run2016%s-17Jul2018*_muonTuples/200718*/*/*.root' % run.replace('ver2','')))

for run in ['B', 'C', 'D', 'E', 'F']:
  samplesToMix.append(('/pnfs/iihe/cms/store/user/tomc/tnpTuples_muons/POG/TnPTree_17Nov2017_SingleMuon_Run2017%sv1_Full_GoldenJSON.root' % run,
      '/pnfs/iihe/cms/store/user/tomc/heavyNeutrino/SingleMuon/crab_Run2017%s*_muonTuples/200718*/*/*.root' % run))

samplesToMix.append(('/pnfs/iihe/cms/store/user/tomc/tnpTuples_muons/POG/TnPTreeZ_17Sep2018_SingleMuon_Run2018Av2_GoldenJSON.root',
    '/pnfs/iihe/cms/store/user/tomc/heavyNeutrino/SingleMuon/crab_Run2018A-17Sep2018-v2_muonTuples/200718*/*/*.root'))

for run in ['Bv1', 'Cv1']:
  samplesToMix.append(('/pnfs/iihe/cms/store/user/tomc/tnpTuples_muons/POG/TnPTreeZ_17Sep2018_SingleMuon_Run2018%s_GoldenJSON.root' % run,
      '/pnfs/iihe/cms/store/user/tomc/heavyNeutrino/SingleMuon/crab_Run2018%s-17Sep2018-v1_muonTuples/200718*/*/*.root' % run.replace('v2','').replace('v1', '')))

samplesToMix.append(('/pnfs/iihe/cms/store/user/tomc/tnpTuples_muons/POG/TnPTreeZ_SingleMuon_Run2018Dv2_GoldenJSON.root',
    '/pnfs/iihe/cms/store/user/tomc/heavyNeutrino/SingleMuon/crab_Run2018D-22Jan2019-v2_muonTuples/200718*/*/*.root'))


if not args.isChild:
  from ttg.tools.jobSubmitter import submitJobs
  jobs = []
  for i in (range(len(samplesToMix)) if not args.jobId else [args.jobId]):
    totalJobs = (100 if 'Run2018D' in samplesToMix[i][0] else 10)
    for j in (range(totalJobs) if not args.subJob else [args.subJob]):
      jobs.append((str(i),str(j)))
  submitJobs(__file__, ('jobId', 'subJob'), jobs, argParser, jobLabel="TP", wallTime='168')
  exit(0)


totalJobs = (100 if 'Run2018D' in samplesToMix[int(args.jobId)][0] else 10)

#
# Mixing function
#
def mix(input, subJob):
  log.info('Mixing %s with %s' % input)
  inputPOG, inputGhent = input
  output = inputPOG.replace('POG', 'updated7')

  if os.path.exists(output.replace('.root', '_%s.root' % subJob)):
    log.info('Finished - output already exists')
    return


  treeGhent  = ROOT.TChain('blackJackAndHookers/blackJackAndHookersTree')
  for path in glob.glob(inputGhent):
    treeGhent.Add(path)

  # You need to build with both run and event number, as apparently sometimes you have the same event number in different runs?
  treeGhent.BuildIndex('_runNb', '_eventNb');

  inputFile  = ROOT.TFile(inputPOG)
  treePOG    = inputFile.Get('tpTree/fitter_tree')

  outputFile = ROOT.TFile(output.split('/')[-1].replace('.root', '_%s.root' % subJob),"RECREATE")
  outputFile.mkdir('tpTree')
  outputFile.cd('tpTree')

  keepBranches = ['pt', 'mt', 'eta', 'phi', 'pair_probeMultiplicity', 'pair_dz', 'pair_nJets30',
                  'tag_pt', 'tag_abseta', 'tag_IsoMu24', 'tag_combRelIsoPF04dBeta', 'tag_dxyPVdzmin', 'tag_dzPV', 'tag_nVertices',
                  'lumi', 'run', 'event',
                  'abseta', 'mass', 'mcMass', 'Medium', 'charge',
                  'combRelIsoPF03', 'combRelIsoPF03dBeta', 'combRelIsoPF04', 'combRelIsoPF04dBeta',
                  'PF', 'Glb', 'TM', 'dxyPVdzmin', 'dzPV',
                  ]
                                       

  treePOG.SetBranchStatus("*", 0)
  for i in keepBranches: treePOG.SetBranchStatus(i, 1)
  outputTree = treePOG.CloneTree(0)

  newBranches=['leptonMvaTop/F', 'dxy/F', 'dz/F', 'SIP3D/F', 'miniIso/F']
  from ttg.tools.makeBranches import makeBranches
  newVars = makeBranches(outputTree, newBranches)

  log.info('Start loop...')
  entries = treePOG.GetEntries()
  thresholds = [i*entries/totalJobs for i in range(totalJobs)]+[entries]
  for i in xrange(thresholds[int(subJob)], thresholds[int(subJob)+1]):
    treePOG.GetEntry(i)
    if treePOG.pt < 10: continue
    if abs(treePOG.eta) > 2.4: continue

# do not expect this will give you the correct event (but often is does, sometimes it doesn't), do not expect ROOT will tell you that it is doing garbage here
# thank you ROOT developers for wasting again a lot of my precious time
#   treeGhent.GetEntryWithIndex(treePOG.event) 
# but ROOT is reasonably close in event number, so let's try to find it in the next events if it does not match
    indexNumber = treeGhent.GetEntryNumberWithIndex(treePOG.run, treePOG.event)
    treeGhent.GetEntry(indexNumber)

    while treeGhent._eventNb != treePOG.event:
      indexNumber += 1
      treeGhent.GetEntry(indexNumber)

    for j in range(treeGhent._nMu):
      if(abs(treeGhent._lPt[j]-treePOG.pt) > 0.1): continue
      if(abs(treeGhent._lEta[j]-treePOG.eta) > 0.1): continue
      newVars.leptonMvaTop = treeGhent._leptonMvaTOP[j]
      newVars.dxy          = treeGhent._dxy[j]
      newVars.dz           = treeGhent._dz[j]
      newVars.SIP3D        = treeGhent._3dIPSig[j]
      newVars.miniIso      = treeGhent._miniIso[j]
      outputTree.Fill()
      break
    else:
      newVars.leptonMvaTop = 0 # these muons are not PF muons, so make sure to also include the PFMuon requirement for tag- and probe 
      newVars.dxy          = 0 
      newVars.dz           = 0 
      newVars.SIP3D        = 0 
      newVars.miniIso      = 0
      outputTree.Fill()
      # this should never happen, as we should have all PF muons in miniAOD
      # nevertheless, you will still a few rare cases where the match didn't work because the pt between AOD and miniAOD is just too different
      # I guess it's better to let these few cases out, they are to little to matter anyway
      # If you see, however, below message way too often, then maybe the run/event matching is broken
      if treePOG.PF:
        log.warning('Found a PF muon in POG trees which is not in Ghent trees')
        log.info('Run:   %s %s' % (str(treePOG.run), str(treeGhent._runNb)))
        log.info('Lumi:  %s %s' % (str(treePOG.lumi), str(treeGhent._lumiBlock)))
        log.info('Event: %s %s' % (str(treePOG.event), str(treeGhent._eventNb)))
        log.info('  need:  %s %s' % (str(treePOG.pt), str(treePOG.eta)))
        for j in range(treeGhent._nMu):
          log.info('  found: %s %s' % (str(treeGhent._lPt[j]), str(treeGhent._lEta[j])))

  outputTree.Write()
  outputFile.Close()
  inputFile.Close()
  log.info('Files are mixed, now copying...')

  # Copy
  os.system('/user/$USER/production/proxyExpect.sh;gfal-copy -f -vvv file://%s srm://maite.iihe.ac.be:8443%s' % (os.path.join(os.getcwd(), output.split('/')[-1].replace('.root', '_%s.root' % subJob)), output.replace('.root', '_%s.root' % subJob)))
  log.info('Finished')

mix(samplesToMix[int(args.jobId)], args.subJob)
