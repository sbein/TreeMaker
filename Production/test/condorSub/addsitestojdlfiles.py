import os, sys
from glob import glob 

try: jdlnamekey = sys.argv[1]
except: jdlnamekey = 'jobExecCondor_Run2016F-17Jul2018-v1.MET.jdl'

jdlnamelist = glob(jdlnamekey)
if not os.path.exists('unusedjdls'): os.system('mkdir unusedjdls')


for jdlname in jdlnamelist:

	if not len(jdlname.split('_'))>1: continue
	print 'editing', jdlname
	fjdl = open(jdlname)
	block = fjdl.read()
	fjdl.close()
	if 'DESIRED_Sites' in block: continue
	scenariokey = block.split(' -j ')[-1].split()[0]
	print 'key:', scenariokey
	partialfilefilename = scenariokey.replace('.','/')+'_cff.py'
	filefilename = os.environ['CMSSW_BASE']+'/src/TreeMaker/Production/python/'+partialfilefilename
	filefile = open(filefilename)
	block2 = filefile.read()
	filefile.close()
	print 'fname:', filefilename
	examplefile = block2.split(',')[-2].strip().replace("'","")
	if not '.root' in examplefile[-7:]: 
		print 'problem', jdlname, filefilename, examplefile
		exit(0)
	print 'example file', examplefile
	os.system('dasgoclient --query="dataset file='+examplefile+'" > '+scenariokey+'.txt')
	datasetnamefile = open(scenariokey+'.txt')
	datasetname = datasetnamefile.read().strip()
	datasetnamefile.close()
	print 'datasetname=', datasetname
	os.system('dasgoclient --query="parent dataset='+datasetname+'" > '+scenariokey+'.txt')
	datasetnamefile = open(scenariokey+'.txt')
	datasetname = datasetnamefile.read().strip()
	datasetnamefile.close()
	print 'parent dataset=', datasetname
	os.system('dasgoclient --query="site dataset='+datasetname+'" > '+scenariokey+'.txt')
	sitefile = open(scenariokey+'.txt')
	sitelist = sitefile.readlines()
	sitefile.close()
	sites = ','.join(sitelist).replace('\n','')
	line4jdl = '+DESIRED_Sites = "'+sites+'"'
	block = block.replace('''# Queue''', \
	'''+DESIRED_Sites = "T1_UK_RAL_Disk,T2_IN_TIFR"
# Queue''')
	print 'new block', block
	os.system('mv '+jdlname+' unusedjdls/')
	newjdl = open(jdlname,'w')
	newjdl.write(block)
	newjdl.close()
	print 'just overwrote but backed up', jdlname, 'in', 'unusedjdls/'
	
	
	
	
	