import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/0419776B-0A7D-E511-9CA4-0025905B8592.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/1401E97A-097D-E511-B4E2-002590E3A004.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/18F02667-0A7D-E511-8B59-002618943984.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/1E9BD5E4-097D-E511-8325-00261894392F.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/24AE1F78-097D-E511-A835-0025905A60F4.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/24BE0473-0A7D-E511-9B47-0025905B859E.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/2E9E6FDD-097D-E511-A480-00261894385A.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/2EAD7E88-097D-E511-A790-0025901D1668.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/307A8472-097D-E511-908F-0026189438A0.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/364144E0-097D-E511-A135-0026189438B3.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/3C6FFADC-097D-E511-9CF0-0026189438F9.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/3CEE6566-097D-E511-8BD0-0025905B858A.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/42840B6F-0A7D-E511-80F1-0025905A610C.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/42F2AC6E-0A7D-E511-B90B-003048FFCB74.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/443EC625-0A7D-E511-BD87-0025905B85AA.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/5075B8E1-097D-E511-AD1E-0025905A6056.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/5658DFE0-097D-E511-ADF2-003048D15E36.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/58B6AB73-097D-E511-9BC7-00261894396E.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/5CAA14DB-097D-E511-949D-00261894396E.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/5CB0F86D-0A7D-E511-B0ED-0025905A497A.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/621E226C-0A7D-E511-944D-0025905A610A.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/6454126C-0A7D-E511-9176-0025905A48D6.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/68EE11E3-097D-E511-8C92-0025905A611C.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/7C215B6B-0A7D-E511-A9CB-0025905B8592.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/7E8FE8E6-097D-E511-B62D-003048FFCB74.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/82886883-097D-E511-AF97-003048344B08.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/8820A0DF-097D-E511-86B6-0025905A6094.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/88BF73E0-097D-E511-8442-0025905B8598.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/90346DDD-097D-E511-A2C1-00261894385A.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/92C325DF-097D-E511-8EA7-002618943910.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/9A41D6DF-097D-E511-856C-003048FFCB9E.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/9CA97472-0A7D-E511-91AC-003048FFD7A2.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/A06E6FE1-097D-E511-9ED9-0025905B858A.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/A6B47EC2-097D-E511-846B-0025905A609A.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/B0E4DFE0-097D-E511-83C4-0025905A60EE.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/B62ED7E4-097D-E511-B350-00261894392F.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/BA105662-097D-E511-96FE-00248C0BE005.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/C41B4B3E-0A7D-E511-A0FA-00261894384F.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/C4A534DF-097D-E511-AE03-0025905A6056.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/D0AD6AE0-097D-E511-80D5-0025905A6094.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/DEF8D36D-0A7D-E511-A205-003048FFCB96.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/E2188D68-0A7D-E511-8B10-0026189438B3.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/ECCB2DE3-097D-E511-A0AF-0025905A48C0.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/F097DE70-0A7D-E511-9414-003048FFCC18.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/F0A08372-0A7D-E511-B3C5-0025905B85A2.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/F4BF28E2-097D-E511-B151-0026189438E8.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/F69115BD-097D-E511-9177-0025905A6066.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/FA351E62-097D-E511-849C-00248C0BE005.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/FAF687E1-097D-E511-B348-0025905B85A2.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/FC25CB6C-0A7D-E511-B7D2-0025905A60E0.root',
       '/store/mc/RunIISpring15MiniAODv2/SMS-T1bbbb_mGlu-1025-1050_mLSP-800to975-400to800-1000to1025_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/50000/FC8F8FDF-097D-E511-8743-0025905A6068.root' ] );


secFiles.extend( [
               ] )
