import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/0A38C4EE-8F40-E811-84BF-0CC47AC52BAE.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/0CFA3F59-8140-E811-A657-0025905A60E0.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/0EAA5FDB-DC40-E811-8F4C-001E672480BB.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/12E09093-B440-E811-9AEE-0025905A6092.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/16616B72-7640-E811-AF07-0CC47A4D7658.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/208212D5-8740-E811-B073-0CC47A71F7B8.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/2E2B276E-9640-E811-8670-001E67248A25.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/2EC84DE7-A440-E811-9CE4-0CC47AC17550.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/2EE64E65-4841-E811-80ED-0CC47A78A3EC.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/4020F2B6-5D40-E811-9477-0025905A6118.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/48B68C48-CF40-E811-9267-001E67248566.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/54FFBEA6-5D40-E811-9F3B-0CC47A7C3612.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/6484E4E4-AF40-E811-91BB-0025905B85A2.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/683E2FE1-8740-E811-A425-001E67247968.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/6E3386D8-1041-E811-98F7-0CC47A4D7632.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7201E363-1541-E811-AFA3-001E672489C1.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/72231EA2-E340-E811-AA6C-0CC47AC577DC.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/749098EC-8740-E811-A97B-001E6724866F.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/781A3BD9-E740-E811-8A4A-0025905A608E.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7C780317-9940-E811-B1A4-0025905A612C.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7EBA1D20-C840-E811-A8DF-0025905B860E.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/887B5CEF-8F40-E811-AE43-0CC47AC52D78.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/8A0B61C9-AB40-E811-86A9-001E6724804D.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/8A43BB7D-5640-E811-B7D7-0CC47A7C35C8.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/8C92A909-EA40-E811-ABCB-0CC47AC17672.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/9047107C-D540-E811-A8C4-0CC47AC52FC2.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/90536E35-9440-E811-B31F-0025905B855A.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/94650235-B340-E811-B954-001E672488A4.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/9C2C3ED8-8740-E811-B1B5-0CC47AC52BDE.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/9C6B5F69-C340-E811-AC51-0025905A607A.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/9C8B587A-D540-E811-9D9A-001E67248566.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/A842BE38-A340-E811-8229-0CC47A4C8E38.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/AE15F7D8-8740-E811-9A13-001E6724816F.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/BAE79FD4-8740-E811-B0CB-0CC47AC52BAE.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/C4B6472D-B340-E811-8C4E-0CC47AC57942.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/C6DF8AA9-D740-E811-8240-0025905A60B2.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/C8938B93-B440-E811-A04C-0025905A6092.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/CA435BD5-8740-E811-8344-0CC47AC30F14.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/D0D535D6-8740-E811-8C71-0CC47AC17550.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/D8645FFA-8940-E811-8BC9-003048F3517A.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/D8D445D3-8740-E811-B85F-0CC47AC52A94.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/DC4B9FD6-8740-E811-871B-0CC47AC17678.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E43D58F9-8940-E811-A327-0CC47AC52CFE.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E48833DA-DC40-E811-8214-0CC47AC174BC.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E4F269B4-5D40-E811-B550-0025905A60F4.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/EE0C698A-B440-E811-81D7-0CC47A4C8E98.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F2FE13AF-7B40-E811-9F48-0CC47A78A3F4.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F634CF43-F140-E811-ACF6-001E672480BB.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F65EECB5-AB40-E811-BAF3-0CC47AC1757C.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F8FC1239-A340-E811-B6D6-0CC47A4D7690.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/FAACAE36-D040-E811-B35D-0CC47A4D7658.root',
       '/store/mc/RunIIFall17MiniAODv2/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/FAD83C11-C840-E811-98E1-0CC47A7139C4.root',
] )
