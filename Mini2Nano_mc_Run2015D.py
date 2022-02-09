import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("MiniAOD2NanoAOD")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("MiniAOD2NanoAOD")

process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    wantSummary=cms.untracked.bool(True)
)

# Set the maximum number of events to be processed (-1 processes all events)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(500))

process.GlobalTag.globaltag = '76X_dataRun2_16Dec2015_v0'

process.source = cms.Source(
    "PoolSource", fileNames=cms.untracked.vstring(
        'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v2/00000/0052E102-10C8-E511-8B8D-003048357A8C.root',
        )
)

# Number of events to be skipped (0 by default)
process.source.skipEvents = cms.untracked.uint32(0)

# Register fileservice for output file
process.miniaod2nanoaod = cms.EDAnalyzer(
    "MiniAOD2NanoAOD",
    isData = cms.bool(False)
)

process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("TTJets_SingleLeptFromT_miniAOD2nanoAOD.root")
)

process.p = cms.Path(process.miniaod2nanoaod)
