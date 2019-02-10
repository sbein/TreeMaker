import FWCore.ParameterSet.Config as cms

# import functions to be assigned as class methods
from TreeMaker.TreeMaker.makeTreeFromMiniAOD_cff import makeTreeFromMiniAOD
from TreeMaker.TreeMaker.JetDepot import JetVariations
from TreeMaker.TreeMaker.makeJetVars import makeJetVars, makeGoodJets, makeJetVarsAK8, makeMHTVars
from TreeMaker.TreeMaker.doHadTauBkg import doHadTauBkg, makeJetVarsHadTau
from TreeMaker.TreeMaker.doLostLeptonBkg import doLostLeptonBkg
from TreeMaker.TreeMaker.doZinvBkg import doZinvBkg, reclusterZinv


import os
def getAODfromMiniAODPath(datasetPathMiniAOD):
    print "getAODfromMiniAODPath" , datasetPathMiniAOD
    filepathbits = datasetPathMiniAOD.split('/store')
    if len(filepathbits)==2:
        filepathbits[1] = '/store'+filepathbits[1]
    # get AOD file(s) from MiniAOD file path
    # fix SSL site configuration for DAS client:
    os.environ['SSL_CERT_DIR'] = '/etc/pki/tls/certs:/etc/grid-security/certificates'
    #parentfiles = check_output(["./data/das_client.py", '--query=parent file=%s' % datasetPathMiniAOD, '--limit=0']) .split()
    parentfiles = []
    os.system('dasgoclient --query="parent file='+filepathbits[1]+'" --limit=0 > '+'tmp.txt')
    ftmp = open('tmp.txt')
    lines = ftmp.readlines()
    ftmp.close()
    os.system('rm tmp.txt')
    for line in lines: parentfiles.append(line.strip())

    aodFiles = []
    for parentfile in parentfiles:
        if "/AOD" in parentfile:
            # parent file is already AOD file
            aodFiles.append(parentfile)
        else:
            # parent file is e.g. RAW, so check its children:
            # childFiles = check_output(["./data/das_client.py", "--query=child file=%s" % parentfile, '--limit=0']).split()
            childFiles = []
            os.system('dasgoclient --query="child file='+parentfile+'" --limit=0 > '+'tmp.txt')
            ftmp = open('tmp.txt')
            lines = ftmp.readlines()
            ftmp.close()
            os.system('rm tmp.txt')
            for line in lines: childFiles.append(line.strip())

            for childFile in childFiles:
                if "/AOD" in childFile:
                    aodFiles.append(childFile)

    # avoid duplicate entries:
    print 'obtained ',len(aodFiles),',file-long aod set:', aodFiles
    return list(set(aodFiles))


def quitIfFileAlreadyProcessed(output_filename):

    import commands

    print "Checking if output file has been already processed"
    status_sam, output = commands.getstatusoutput("gfal-ls srm://dcache-se-cms.desy.de/pnfs/desy.de/cms/tier2/store/user/sbein/NtupleHub/Production2016v2/" + output_filename + ".root")
    status_viktor, output = commands.getstatusoutput("gfal-ls srm://dcache-se-cms.desy.de/pnfs/desy.de/cms/tier2/store/user/vkutzner/NtupleHub/" + output_filename + ".root")
    
    if status_sam == 0 or status_viktor == 0:
        print "Output file already exists. OK, exiting"
        quit(77)
    else:
        print "OK, proceeding"


class maker:
    def __init__(self,parameters):
        self.parameters = parameters
    
        # auto configuration for different scenarios
        self.scenarioName=self.parameters.value("scenario","")
        from TreeMaker.Production.scenarios import Scenario
        self.scenario = Scenario(self.scenarioName)

        # to keep track of MET fix, currently applied to all 2017 data and MC
        self.doMETfix = ("Fall17" in self.scenarioName or "2017" in self.scenarioName)
        self.getParamDefault("verbose",True)
        self.getParamDefault("inputFilesConfig","")
        self.getParamDefault("dataset",[])
        self.getParamDefault("sidecar",[])#SB     
        self.getParamDefault("privateSample",False)#SB
        self.getParamDefault("nstart",0)
        self.getParamDefault("nfiles",-1)
        self.getParamDefault("numevents",-1)
        self.getParamDefault("outfile","test_run")
        # get base sample name, assuming format is name or name_#
        outfilesplit = self.outfile.split('_')
        if outfilesplit[-1][-1].isdigit(): outfilesplit = outfilesplit[:-1]
        self.sample = '_'.join(outfilesplit)
        outfilesuff=self.parameters.value("outfilesuff","_RA2AnalysisTree")
        self.outfile += outfilesuff
        self.getParamDefault("treename","PreSelection")
        
        # background estimations on by default
        self.getParamDefault("lostlepton", True)
        self.getParamDefault("hadtau", False)
        self.getParamDefault("hadtaurecluster", 0)
        self.getParamDefault("doZinv", True)

        # special signal stuff
        self.getParamDefault("systematics",True);
        self.getParamDefault("semivisible",True);
        self.getParamDefault("deepAK8",True);
        self.getParamDefault("deepDoubleB",True);
        
        # compute the PDF weights
        self.getParamDefault("doPDFs", True);
        
        # other options off by default
        self.getParamDefault("debugtracks", False)
        self.getParamDefault("applybaseline", False)
        
        # take command line input (w/ defaults from scenario if specified)
        self.getParamDefault("globaltag",self.scenario.globaltag)
        self.getParamDefault("tagname",self.scenario.tagname)
        self.getParamDefault("hlttagname",self.scenario.hlttagname)
        self.getParamDefault("geninfo",self.scenario.geninfo)
        self.getParamDefault("pmssm",self.scenario.pmssm)
        self.getParamDefault("fastsim",self.scenario.fastsim)
        self.getParamDefault("signal",self.scenario.signal)
        self.getParamDefault("jsonfile",self.scenario.jsonfile)
        self.getParamDefault("jecfile",self.scenario.jecfile)
        self.getParamDefault("residual",self.scenario.residual)
        self.getParamDefault("jerfile",self.scenario.jerfile)
        self.getParamDefault("pufile",self.scenario.pufile)
        self.getParamDefault("wrongpufile",self.scenario.wrongpufile)
        self.getParamDefault("era",self.scenario.era)
        self.getParamDefault("localera",self.scenario.localera)
        
        # temporary redirector fix
        self.getParamDefault("redir", "root://cmsxrootd.fnal.gov/")
        # handle site name usage
        if self.redir[0]=="T":
            self.redir = "root://cmsxrootd.fnal.gov//store/test/xrootd/"+self.redir
        
        # Load input files
        self.readFiles = cms.untracked.vstring()
        self.readFiles_sidecar = cms.untracked.vstring()

        if self.inputFilesConfig!="" :
            readFilesImport = getattr(__import__("TreeMaker.Production."+self.inputFilesConfig+"_cff",fromlist=["readFiles"]),"readFiles")
            if self.nfiles==-1:
                self.readFiles.extend( readFilesImport )
            else:
                self.readFiles.extend( readFilesImport[self.nstart:(self.nstart+self.nfiles)] )

        if self.dataset!=[] :    
            self.readFiles.extend( [self.dataset] )
        for irf, rf in enumerate(self.readFiles):
            quitIfFileAlreadyProcessed(self.outfile)
            if '/store/' in rf: self.readFiles_sidecar += getAODfromMiniAODPath(rf)
            else: 
            	shpingy = rf.replace('mini','').replace('step4','step3')
            	self.readFiles_sidecar.append(shpingy)
        self.readFiles = [(self.redir if val[0:6]=="/store" else "")+val for val in self.readFiles]
        self.readFiles_sidecar = [(self.redir if val[0:6]=="/store" else "")+val for val in self.readFiles_sidecar]
        print 'readFiles', self.readFiles
        print 'readFiles_sidecar', self.readFiles_sidecar
        
        # branches for treemaker
        self.VectorRecoCand             = cms.vstring()
        self.VarsDouble                 = cms.vstring()
        self.VarsInt                    = cms.vstring()
        self.VarsBool                   = cms.vstring()
        self.VectorTLorentzVector       = cms.vstring()
        self.VectorDouble               = cms.vstring()
        self.VectorString               = cms.vstring()
        self.VectorInt                  = cms.vstring()
        self.VectorBool                 = cms.vstring()
        self.VectorVectorBool           = cms.vstring()
        self.VectorVectorInt            = cms.vstring()
        self.VectorVectorDouble         = cms.vstring()
        self.VectorVectorString         = cms.vstring()
        self.VectorVectorTLorentzVector = cms.vstring()

    def getParamDefault(self,param,default):
        setattr(self,param,self.parameters.value(param,default))
        
    def printSetup(self):
        print " readFiles: "+str(self.readFiles)
        print " outfile: "+self.outfile
        print " treename: "+self.treename
        print " "
        print " storing lostlepton variables: "+str(self.lostlepton)
        print " storing hadtau variables: "+str(self.hadtau)+" w/ reclustering "+str(self.hadtaurecluster)
        print " storing Zinv variables: "+str(self.doZinv)
        print " storing semi-visible jet variables: "+str(self.semivisible)
        print " storing deepAK8 variables: "+str(self.deepAK8)
        print " storing deepDoubleB variables: "+str(self.deepDoubleB)
        print " "
        print " storing JEC/JER systematics: "+str(self.systematics)
        print " storing PDF weights: "+str(self.doPDFs)
        print " "
        print " storing track debugging variables: "+str(self.debugtracks)
        print " Applying baseline selection filter: "+str(self.applybaseline)
        print " "
        print " scenario: "+self.scenarioName
        print " global tag: "+self.globaltag
        print " Instance name of tag information: "+self.tagname
        print " Instance name of HLT tag information: "+self.hlttagname
        print " Including gen-level information: "+str(self.geninfo)
        print " Including pMSSM-related information: "+str(self.pmssm)
        print " Using fastsim settings: "+str(self.fastsim)
        print " Running signal uncertainties: "+str(self.signal)
        if len(self.jsonfile)>0: print " JSON file applied: "+self.jsonfile
        if len(self.jecfile)>0: print " JECs applied: "+self.jecfile+(" (residuals)" if self.residual else "")
        if len(self.jerfile)>0: print " JERs applied: "+self.jerfile
        if len(self.pufile)>0: print " PU weights stored: "+self.pufile
        print " era of this dataset: "+self.era
        print " local era of this dataset: "+self.localera

    # assign class methods from functions
    makeTreeFromMiniAOD = makeTreeFromMiniAOD
    JetVariations = JetVariations
    makeGoodJets = makeGoodJets
    makeJetVars = makeJetVars
    makeJetVarsAK8 = makeJetVarsAK8
    makeMHTVars = makeMHTVars
    doHadTauBkg = doHadTauBkg
    makeJetVarsHadTau = makeJetVarsHadTau
    doLostLeptonBkg = doLostLeptonBkg
    doZinvBkg = doZinvBkg
    reclusterZinv = reclusterZinv
