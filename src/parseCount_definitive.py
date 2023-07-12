#!/usr/bin/env python

import subprocess
from optparse import OptionParser
import sys, os, random
import time 
import numpy as np
import pandas as pd
import re
import copy
import tempfile 
from tqdm import tqdm       #It creates a progress bar
from math import exp, log
from scipy.stats import poisson
#import statsmodels.api as sm
from function import *
import math
#import statsmodels.discrete.count_model as cm
#from statsmodels.distributions.empirical_distribution import ECDF
#from statsmodels.base.model import GenericLikelihoodModel
from scipy.stats import nbinom

#from statsmodels.genmod.generalized_linear_model import NegativeBinomial


import warnings
warnings.filterwarnings("ignore")





pathofexec  = os.path.abspath(sys.argv[0]) 
pathofexecarray=pathofexec.split('/')[:-2] 

pathofconfig=  ("/".join(pathofexecarray))+"/config/homosapiens_hg19.cfg" 
if(not os.path.exists(pathofconfig)):
    sys.stderr.write("\nERROR: The configuration file "+pathofconfig+" does not exist, re-clone the repository\n") 
    sys.exit(1) 


gcContent=  ("/".join(pathofexecarray))+"/data/hs37d5.gc" 
if(not os.path.exists(gcContent)):
    sys.stderr.write("\nERROR: The GC file "+gcContent+" does not exist\n") 
    sys.exit(1) 



mismappingrate=0.01
    
#print(pathofconfig) 
parser = OptionParser(
"\n\n"+

" ##     ##    ###    ########  ######## ##     ## ########\n"+
" ###   ###   ## ##   ##     ##    ##    ##     ## ##      \n"+ 
" #### ####  ##   ##  ##     ##    ##    ##     ## ##      \n"+      
" ## ### ## ##     ## ########     ##    ######### ######  \n"+   
" ##     ## ######### ##   ##      ##    ##     ## ##      \n"+       
" ##     ## ##     ## ##    ##     ##    ##     ## ##      \n"+       
" ##     ## ##     ## ##     ##    ##    ##     ## ########\n"+
"\n"+

"Bayesian inference of KARYOtype for low coverage data\n"+

"marthe [options] [bam file 1] [bam file 2]...") 



parser.add_option("-r", "--reference",      dest="ref",          help="Reference genome used for alignment",                 default=None,    type="string") 
parser.add_option("--tmpfolder",            dest="tmpfolder",    help="Temporary folder",                                    default="/tmp/",    type="string") 
parser.add_option("-o","--outfile",         dest="resultso",     help="Output prefix to store the results",                  default="results",    type="string") 
parser.add_option("--samtools",             dest="samtools",     help="Use this version of samtools",                        default="samtools",    type="string") 
parser.add_option("--tabix",                dest="tabix",        help="Use this version of tabix",                           default="tabix",    type="string") 
parser.add_option("--map",                  dest="map_mask",     help="Use a mappability map in BED format (Recommended)",   default=None,    type="string") 
parser.add_option("--mapscore",             dest="map_track",     help="Use a mappability track",   default=None,    type="string") 
parser.add_option("--hpc",                  dest="hpc",          help="Use high-performance computing (for queueing systems ex: SGE)",          action="store_true") 
parser.add_option("--resume",               dest="resume",       help="Resume by providing the temp directory used",                              type="string") 
parser.add_option("--nice",                 dest="nice",         help="Nice the jobs",                                                          action="store_true") 
parser.add_option("-t"  , "--threads",      dest="threads",      help="Number of threads to use if the local machine is used",  default=1,   type="int") 
parser.add_option("--bamidx",               dest="bamidx",       help="Index of the BAM file to process",                              type="int")
parser.add_option("--mismap",               dest="mismappingrate", help="Mismapping rate , default is "+str(mismappingrate)+"",            default=mismappingrate, type="float") 
parser.add_option("--winsize",              dest="winsize",      help="Size of genomic windows, default is 1Mbp",             default=1000000, type="int") 
parser.add_option("--minlen",               dest="minlen",       help="Minimum length of DNA fragment",             default=35, type="int") 
parser.add_option("--gc",                   dest="gcfile",       help="File containing the % of GC per window size (see above) default: "+gcContent+"",                  default=gcContent,    type="string") 
parser.add_option("-c", "--conf",           dest="configfile",   help="Configuration for various conditions file to use for species/build default: "+str(pathofconfig)+"", default=pathofconfig, type="string") 
parser.add_option("-F",                     dest="file",         help="There is a list with names of BAM files",             action="store_true") 
parser.add_option("-p", "--parallel",        dest="parallel",     help="List with nodes if parallelization is chosen",        default=None) 

#This has to be removed later
parser.add_option("--subsample",           dest="subsample",   help="Subsample value for the simulations", default=1, type="int") 


#TODO ?
#parser.add_option(""  , "--branchl",     dest="branchlscale", help="Seq-gen branch scale, default is 0.00045",            default=0.00045, type="float") 
#parser.add_option(""  , "--chrlen",      dest="lengthchr",    help="Chromosome length, default is 10kb",              default=10000,   type="int") 
#parser.add_option("-c", "--numcont",     dest="numcont",      help="Number of present-day human contaminants, default is 2", default=2, type="int")




(options,args) = parser.parse_args()

if( options.resume == None):
    if( len(args)==0 ):
        sys.stderr.write("\nPlease use -h to see options.\n") 
        sys.exit(1) 

mismappingrate=options.mismappingrate 

#detecting programs
samtoolscmd     = re.sub(b'\s+', b'', which(options.samtools))
samtoolscmd = samtoolscmd.decode()

results_path = options.tmpfolder

isExist = os.path.exists(results_path)
if not isExist:
    os.makedirs(results_path)

#detecting programs

rmcmd     = "rm" 
#rmcmd     = re.sub('\s+','',which("rm")) 
#awkcmd    = re.sub('\s+','',which("awk")) 
#cutcmd    = re.sub('\s+','',which("cut")) 

#seqgencmd = re.sub('\s+','',which("seq-gen")) 





chrnames   = {} 
#chrlengths = [] 
chrranges  = [] 
autoranges = [] 
sexranges  = [] 
sexrangesUnaffectedMale   = [] 
sexrangesUnaffectedFemale = [] 


autochr = [] 
sexchr = [] 

#
unaffemalechrs = [] 
unafmalechrs  = [] 

unaffemaleProb = 0.5 
unafmaleProb = 0.5 



referencefai = open(options.ref+".fai", "r") 

#Read the reference file and store information
listchr = list(range(1,23))+["X", "Y", "x", "y"]
listchr= [str(option) for option in listchr]


for linefai in referencefai:

    faia=linefai.split("\t") 
    if faia[0] in listchr:
        #chrnames.append(   faia[0] ) 
        chrnames[ faia[0] ] = 0 

        #chrlengths.append( faia[1] ) 
        for c in range(1, (int(faia[1])-options.winsize), options.winsize):
            rangedict =	{ "chr": faia[0], "start": c,   "end": c+(options.winsize-1), "ecov":1 }     
            chrranges.append(rangedict) 
    
referencefai.close() 


######################
#READING GC FILE #
######################

gcranges  = [] 
sys.stderr.write( "\nReading the GC content file "+str(options.gcfile) +"\n" ) 
if(not os.path.exists(options.gcfile)):
    sys.stderr.write("\nError gc content file:"+options.gcfile+" does not exist. This file is needed to quantify the GC content\n") 
    sys.exit(1) 
else:
    gcContentfile = open(options.gcfile, "r") 
    while True:

        linegc = gcContentfile.readline() 
        if not linegc :
            break 

        linegc=linegc.strip() 
        gcf=linegc.split("\t") 
        region=gcf[0] 
        region2=region.split(":") 
        chrgc   = region2[0] 
        coordgc = region2[1] 
        region3=coordgc.split("-") 
        startgc   = region3[0] 

        if(startgc.endswith("0")):
            startgc   = int( startgc )+1 
        else:      
            startgc   = int( startgc ) 
        
        endgc     = int(region3[1]) 
        
        gccont=gcf[1] 
        rangedict =	{ "chr": chrgc, "start": startgc,   "end": endgc, "gc": float(gcf[1]) }      
        gcranges.append( rangedict ) 


######################
#READING MAPPABILITY FILE 
#Correcting performance and storing it
######################

sys.stderr.write( "\nReading the mappability file "+str(options.map_track) +"\n" ) 
if options.map_track != None:
    map_file = open(options.map_track, "r") 
    line = map_file.readline()
    map_file.close()

    if line.startswith("fixedStep"):
        sys.stderr.write("The format of the mappability file is not correct.\nCorrecting format:")
        new_mappability_file = ("/".join(pathofexecarray))+"/data/map_ucsc.wig"
        convert_map_file(mappability_file, new_mappability_file)  #It doesn't return anything, it simply writes content to the new mappability file.
        mappability_file = new_mappability_file     #It gives a new value to the variable
        #Save the mappability scores
        infile = open(mappability_file, "r")
        map_score_list = [line.split()[2] for line in infile]
        infile.close()


    else:
        elements=line.split()
        len_line = len(elements)
        
        if len_line!=3:
            sys.stderr.write("The format of the mappability file is not correct, it should have 3 columns.")
            sys.exit(1)
        else:
            if not isInt(elements[0]):
                sys.stderr.write("The first column should be an integer referring to the genomic region")
                sys.exit(1)
            map_file = open(options.map_track, "r") 
            map_score_list = [line.split()[2] for line in map_file]
            map_file.close()




else:   #If mappability file is not used, all the mappability scores are set to 1
    map_score_list = [1]*len(chrranges)


        
######################
#READING CONFIG FILE #
######################

conditions=[] 

sys.stderr.write("\nReading the configuration file "+str(options.configfile) +"\n" ) 

if(not os.path.exists(options.configfile)):
    sys.stderr.write("\nError config file:"+options.configfile+" does not exist. This file is needed to define the different karyotypes\n") 
    sys.exit(1) 
else:
    configfile = open(options.configfile, "r") 
    
    ###################################################
    #this line defines the autosomes 1..22, stores them in autochr
    ###################################################
    lineconfig, autosomes  = readlineconfig(configfile, 'autosomes:', options)       #auto
 
    for autoidx in range(0,len(autosomes)):
        if autosomes[autoidx] in chrnames:
            autochr.append(autosomes[autoidx]) 

        else: 
            sys.stderr.write("\nParse error with config file:{} The autosome:{} was not found\n".format(options.configfile, str(autosomes[autoidx]))) 
            sys.exit(1) 

    sys.stderr.write("\nFound the following autosomes: "+(",".join(autochr))) 
                     
    ###################################################
    #this line defines which chromosomes are sexual, stores them in sexchr
    ###################################################

    lineconfig, sexsomes = readlineconfig(configfile, "sex:", options)   #sex

    for sexidx in range(0, len(sexsomes)):
        if sexsomes[sexidx] in chrnames:
            sexchr.append(sexsomes[sexidx]) 

        else:
            sys.stderr.write("\nParse error with config file:"+options.configfile+" The sex:"+str(sexsomes[sexidx])+" was not found\n") 
            sys.exit(1) 

    sys.stderr.write("\nFound the following sex chromosomes: "+(",".join(sexchr))) 

    #Check if the sex chromosome is also defined as autosome
    for sexidx in range(0,len(sexchr)):
        if sexchr[sexidx] in autochr:
            sys.stderr.write("\nParse error with config file:{} The sex chromosome:{} was also defined as autosome".format(options.configfile, str(sexchr[sexidx]))) 
            sys.exit(1) 


    ###################################################
    #defines an unaffected female (no chromosome conditions), stores the chromosome names in unafemalechrs and the probability in unafemaleProb
    ###################################################

    lineconfig, sexsomes = readlineconfig(configfile, 'unaffectedfemale:', options)   #unaffectedfemale

    if( len(sexsomes) != 3):
        sys.stderr.write("\nParse error with config file:"+options.configfile+" The third line needs to define the karyotype of an unaffectedfemale as such:\nunaffectedfemale: [chr#1] [chr#2] [prob]\nex:unaffectedfemale: X X 0.5\nFound "+str(len(sexsomes))+" fields instead of 3\n") 
        sys.exit(1) 

    for sexidx in range(0,len(sexsomes)-1):
        if not sexsomes[sexidx] in sexchr:
            sys.stderr.write("\nError with config file:"+options.configfile+" The chromosome :"+str(sexsomes[sexidx])+" was not found in sex chromosomes") 
            sys.exit(1) 
        unaffemalechrs.append(sexsomes[sexidx]) 

    unaffemaleProb = sexsomes[len(sexsomes)-1] 
    if(not isFloat(unaffemaleProb)):
        sys.stderr.write("\nParse error with config file:"+options.configfile+", last field in line "+lineconfig+" is not a probability\nThe line needs to define a condition as such:\nname: [aMF] [chr:region][dD] [prob]\nex:Down: a 21D 0.001\n") 
        sys.exit(1) 
    unaffemaleProb = float(unaffemaleProb) 


    sys.stderr.write("\nFound the following definition for an unaffected male: {} p= {}".format(",".join(unaffemalechrs), str(unaffemaleProb)) ) 


    ###################################################
    #defines an unaffected male (no chromosome conditions), stores the chromosome names in unafmalechrs and the probability in unafmaleProb
    ###################################################

    lineconfig, sexsomes = readlineconfig(configfile, 'unaffectedmale:', options)   #unaffectedmale

    if( len(sexsomes) != 3):
        sys.stderr.write("\nParse error with config file:"+options.configfile+" The third line needs to define the karyotype of an unaffectedmale as such:\nunaffectedmale: [chr#1] [chr#2] [prob]\nex:unaffectedmale: X X 0.5\nFound "+str(len(sexsomes))+" fields instead of 3\n") 
        sys.exit(1) 

    for sexidx in range(0,len(sexsomes)-1):
        if not sexsomes[sexidx] in sexchr:
            sys.stderr.write("\nError with config file:"+options.configfile+" The chromosome :"+str(sexsomes[sexidx])+" was not found in sex chromosomes\n") 
            sys.exit(1) 
        unafmalechrs.append(sexsomes[sexidx]) 

    unafmaleProb = sexsomes[len(sexsomes)-1] 
    if(not isFloat(unafmaleProb)):
        sys.stderr.write("\nParse error with config file:"+options.configfile+", last field in line "+lineconfig+" is not a probability\nThe line needs to define a condition as such:\nname: [aMF] [chr:region][dD] [prob]\nex:Down: a 21D 0.001\n") 
        sys.exit(1) 
    unafmaleProb = float(unafmaleProb) 

    sys.stderr.write("\nFound the following definition for an unaffected male: {} p= {}".format(",".join(unafmalechrs), str(unafmaleProb)) ) 


    ##################################################·
    #conditions
    ###################################################

 
                     
    sys.stderr.write("\nReading conditions to detect:\n") 
    while lineconfig:

        lineconfig = configfile.readline() 

        if(lineconfig.find("#")==0): #It is a comment line
            continue 

        if(lineconfig.find("#")!=-1):   #Don´t find #
            lineconfig=lineconfig[0:lineconfig.find("#")] #Until the next comment line
            lineconfig = lineconfig.strip() 

        if(len(lineconfig)==0): #Emptly line, skip
            continue 

        if(lineconfig.find(":")==-1):
            sys.stderr.write("\nParse error with config file:"+options.configfile+", did not find a colon\nThe line needs to define a condition as such:\nname: [aMF] [chr:region][dD] [prob]\nex:Down: a 21D 0.001\n") 
            sys.exit(1) 

        namecond_=lineconfig[0:lineconfig.find(":")] 
        namecond="" 
        namecond = "".join([ c if c.isalnum() else "_" for c in namecond_ ])

        print ("Name: "+namecond)
        fields=lineconfig[(lineconfig.find(":")+1):len(lineconfig)].split( ) #fields=['a', '22:1-15519412D', '0.00001351351']
        type_chr=-1 #0 = auto, 1 affects females 2 affects males
    
        if(fields[0].lower() == "a"):
            type_chr="auto"
            print ("Type: autosomal")

        if(fields[0].upper() == "M"):
            type_chr = "male"
            print ("Type: sexual chromosomes (M)")

        if(fields[0].upper() == "F"):
            type_chr = "female"
            print ("Type: sexual chromosomes (F)")

        priorprob=fields[len(fields)-1] 

        if(type_chr == -1):
            sys.stderr.write("\nParse error with config file:"+options.configfile+", first field must be a M or F\nThe line needs to define a condition as such:\nname: [aMF] [chr:region][dD] [prob]\nex:Down: a 21D 0.001\n") 
            sys.exit(1) 

        if(not isFloat(priorprob)):
            sys.stderr.write("\nParse error with config file:"+options.configfile+", last field is not a probability\nThe line needs to define a condition as such:\nname: [aMF] [chr:region][dD] [prob]\nex:Down: a 21D 0.001\n") 
            sys.exit(1) 

        priorprob=float(priorprob) 
        if(priorprob<0 or priorprob>1):
            sys.stderr.write("\nParse error with config file:"+options.configfile+", last field is not a probability 0-1\nThe line needs to define a condition as such:\nname: [aMF] [chr:region][dD] [prob]\nex:Down: a 21D 0.001\n") 
            sys.exit(1) 
            
        print ("Occurrence: {} (1 in every {})".format(str(priorprob), str(int(1/priorprob))))

        fields=fields[1:len(fields)-1] #fields = ['22:1-15519412D']
        defChr=[] 
        
        for defin in fields:    #defin = '22:1-15519412D'
            typedD=0   #0 = nothing; 1 = Duplication; -1 = deletion
            if(defin.endswith("d")):
                typedD="deletion"

            if(defin.endswith("D")):
                typedD= "duplication"

            if(typedD == 0):
                sys.stderr.write("\nParse error with config file:"+options.configfile+", each definition must end with d or D for deletions and duplications respectively, found: "+defin+"\nThe line needs to define a condition as such:\nname: [aMF] [chr:region][dD] [prob]\nex:Down: a 21D 0.001\n") 
                sys.exit(1) 

            defin=defin[0:len(defin)-1]     #defin = '22:1-15519412'
            
            chrnamec = "" 
            coordstc = -1 
            coordenc = -1 
            
            if( (defin.find(":") != -1) and (defin.find("-") != -1) and (defin.find(":") < defin.find("-") )):
                chrnamec = defin[0:defin.find(":")]     #chrnamec = '22'
                coordstc = defin[(defin.find(":")+1):defin.find("-")]   #coordstc = '1'
                coordenc = defin[(defin.find("-")+1):len(defin)]    #coordend = '15519412'

                if( (not isInt(coordstc)) or (not isInt(coordenc) ) ):
                    sys.stderr.write("\nParse error with config file:"+options.configfile+", the region has to be either [chr] or [chr]:[start]-[end], found: "+defin+"\nThe line needs to define a condition as such:\nname: [aMF] [chr:region][dD] [prob]\nex:criduchat_longarm a 22:15519412-50818468D 0.001\n") 
                    sys.exit(1) 
                coordstc = int(coordstc) 
                coordenc = int(coordenc) 

                if( coordstc > coordenc  ):
                    sys.stderr.write("\nParse error with config file:"+options.configfile+", the region  [chr]:[start]-[end] needs an end coordinate greater than the start, found: "+defin+"\nThe line needs to define a condition as such:\nname: [aMF] [chr:region][dD] [prob]\nex:criduchat_longarm a 22:15519412-50818468D 0.001\n") 
                    sys.exit(1) 
            else:
                chrnamec = defin 
            
            if(type_chr == "auto"): #auto
                if not chrnamec in autochr:
                    sys.stderr.write("\nParse error with config file:"+options.configfile+" The chromosome:"+str(chrnamec)+" was defined as autosome but we cannot find it in the autosomes\n") 
                    sys.exit(1) 

            if(type_chr == "female" or type_chr == "male"): #sex chr
                if not chrnamec in sexchr:
                    sys.stderr.write("\nParse error with config file:"+options.configfile+" The chromosome:"+str(chrnamec)+" was defined as a sex chromosome but we cannot find it in the sex chromosome\n") 
                    sys.exit(1) 
            defChr.append( {"type": typedD ,
                            "chrname":chrnamec,
                            "coordst":coordstc,
                            "coorden":coordenc }  ) 

        conditions.append({ "namecond":namecond,
                            "type": type_chr, #affects autosome, male or female
                            "priorprob":priorprob,
                            "defChr":defChr}) 
        


            
        print ("------------")

#check if a chr was both defined as auto and sex
for autoidx in range(0,len(autochr)):
    if autochr[autoidx] in sexchr:  #line starts with "sex"
        sys.stderr.write("\nProblem with "+str(autochr[autoidx])+" which was also defined as sex chromosome.\n") 
        sys.exit(1) 

#split chrranges into autosomes and sex chromosomes 
chrrangesnew  = [] 
for r in chrranges: #chrranges contains { "chr", "start","end", "ecov": }
    if r["chr"] in autochr:
        autoranges.append(r)    #autoranges contains { "chr", "start","end", "ecov": }
        chrrangesnew.append(r) 
    if r["chr"] in sexchr:
        sexranges.append(r)     #sexranges contains { "chr", "start","end", "ecov": }
        chrrangesnew.append(r) 


chrranges = chrrangesnew    #chrranges contains { "chr", "start","end", "ecov": }

#unafmalechrs and unaffemalechr
for r in sexranges:
    rmale   = copy.deepcopy(r) 
    rfemale = copy.deepcopy(r) 
    rmale["ecov"]  =0.5*unafmalechrs.count(  r["chr"])  
    rfemale["ecov"]=0.5*unaffemalechrs.count(r["chr"])

    sexrangesUnaffectedMale.append(rmale)
    sexrangesUnaffectedFemale.append(rfemale)   #both contain { "chr", "start","end", "ecov": }


#sys.exit(0) 
#scale prior probabilities
if (unafmaleProb+unaffemaleProb)!=1:
    sys.stderr.write("\nWARNING: The probabilities for unaffected males and females do not add to 1, rescalling.\n") 
    s = (unafmaleProb+unaffemaleProb) 
    unafmaleProb   = unafmaleProb/s 
    unaffemaleProb = unaffemaleProb/s 



unafmaleProbScaled  = unafmaleProb 
unaffemaleProbScaled= unaffemaleProb 
conditionscpy=[] 

sys.stderr.write("\nScaling probabilities:") 
                     

    
    
#for each condition, adding a male and female condition for autosomal
#recomputing the unaffected male/female probabilities
for cd in range(0, len(conditions)):

    if conditions[cd]["type"] == "auto":#autosomes
       mpr = unafmaleProb   * conditions[cd]["priorprob"]  
       fpr = unaffemaleProb * conditions[cd]["priorprob"]  
       
       unafmaleProbScaled   -= mpr 
       unaffemaleProbScaled -= fpr

       malec    = copy.deepcopy(conditions[cd]) 
       femalec  = copy.deepcopy(conditions[cd]) 
       

       malec["namecond"]    = conditions[cd]["namecond"]+"_m" 
       femalec["namecond"]  = conditions[cd]["namecond"]+"_f" 

       malec["priorprob"]   = mpr 
       femalec["priorprob"] = fpr 

       conditionscpy.append(malec) 
       conditionscpy.append(femalec) 
       


    if conditions[cd]["type"] == "male":#male
        unafmaleProbScaled   -= conditions[cd]["priorprob"]          
        conditionscpy.append( conditions[cd] ) 
        

    if conditions[cd]["type"] == "female":#female
        unaffemaleProbScaled   -= conditions[cd]["priorprob"]  
        conditionscpy.append( conditions[cd] ) 
        



conditions = conditionscpy #It contains condition, type, prior probability, deletion or insertion, chromosome, start, end


conditions.insert(0, {"namecond":"unaffectedmale",
                    "type": "male",
                    "priorprob":unafmaleProbScaled,
                    "defChr":[]}) 

conditions.insert(0, {"namecond":"unaffectedfemale",
                    "type": "female",
                    "priorprob":unaffemaleProbScaled,
                    "defChr":[]}) 


sys.stderr.write("\nNew unaffected male probability: "+str(unafmaleProbScaled)) 
sys.stderr.write("\nNew unaffected female probability: "+str(unaffemaleProbScaled)) 

#It tells how many regions are affected because they have overlap, deletion or insertion
affectedrau = AffectedRegions(conditions, autoranges, sexrangesUnaffectedMale, sexrangesUnaffectedFemale)




#we parse the correspondance and find which bam files are we responsible for via bamindx

tfm=options.tmpfolder

datatoprocess=[] 
uniquebamfilesArray=[]

if options.file:
    if len(args) > 1:
        sys.stderr.write("ERROR: there must be only one file listing BAM/CRAM samples")
        sys.exit(0)
    else:
        infile = open(args[0], "r")
        for line in infile:
            uniquebamfilesArray.append(line[:-1])
        infile.close()
else:
    for i in range(len(args)):
        bamfile = str(args[i]) 
        bamfile = os.path.abspath( bamfile )  # resolve path

        uniquebamfilesArray.append(bamfile)





if True:
    b=uniquebamfilesArray[options.bamidx] 
    name_bam=b.split("/")[-1]
    name_bam = name_bam.split(".")[0]
    countfile = options.tmpfolder + "/counts_windows_" + name_bam
    
    if (not os.path.exists(countfile)):
        sys.stderr.write("\ERROR: The file {} with the read counts does not exist".format(countfile))
        sys.exit(1)

    outputfilek = open (results_path+'/'+ options.resultso+"_marthe_"+name_bam+".txt" , 'w' ) 
    outputfilek.write("#bamfile\t" + b + "\t") 
    outputfilek.write("\n") 



    ########################################################################
    #Read the counts file and create the corrected count file
    ########################################################################

    #merged_count_file = options.tmpfolder + "/merged_counts_"+name_bam
    filtered_count_file = options.tmpfolder + "/filtered_counts_"+name_bam
    
    #sys.stderr.write( "\nReading the file "+str(merged_count_file) +"\n" ) 
    #if(not os.path.exists(merged_count_file)):

        #If there is not merged count file, it creates it, it also creates countInfo and corrected count file
        #sys.stderr.write("Merged count file does not exist, generating it...")

        #countInfo = merge_files(options.tmpfolder, map_score_list, affectedrau, filtered_count_file)
        #THIS IS NOT FOR FINAL VERSION
        #First we read the tmp folder and store the read counts in the same file, then we do the subsampling reduction
        #Finally, we filter and produce the read count

    #fraction_list = [".5", ".05", ".005",".0005", ".00005", ".000005", ".0000005", ".00000009", ".00000008", ".00000007", ".00000006"]
    fraction_list = [".0005",  ".0004",  ".0003",  ".0002",  ".0001", ".00009", ".00008", ".00007", ".00006"]
    if options.subsample != 1:
        index_subsample = options.subsample-1
        final_coverage = fraction_list[index_subsample]
        subsample_value = get_subsample2(b, final_coverage, samtoolscmd)
    else:
        subsample_value = 1 #I don't subsample the read count
    
    #countInfo = subsample_file(subsample_value, options.tmpfolder, map_score_list, affectedrau, filtered_count_file)
    countInfo = subsample_file2(subsample_value, countfile, map_score_list, affectedrau, filtered_count_file)


    """else:
        #PROBABLY CORRECT FOR THE NEW VERSION
        #If there is merged count file, it creates the countInfo and corrected count file
        print("there is merged count file")
        #merged_count_file = ("/".join(pathofexecarray))+"/data/merged_counts_"+name_bam
        merged_count_file = options.tmpfolder + "/merged_counts_"+name_bam
        countInfo = []
        merged_count=open(merged_count_file, "r")
        
        for line in merged_count:
            elements= line.split("\t")
            rau=int(elements[0])
            mapped_reads=float(elements[1])
            map_score = float(map_score_list[rau])

            if map_score !=0:
                corrected_read = mapped_reads/map_score

            else:
                corrected_read = mapped_reads

            countInfo.append({"id": elements[0], "corrected_read":corrected_read, "count":mapped_reads ,"map_score":map_score})

        merged_count.close()
        #It creates the corrected count file from the merged counts file and the list of affected region
        filter_file(merged_count_file, filtered_count_file, affectedrau)"""



    #################################################
    #Filter the merged counts file and the gc file
    #Takes only the diploid chromosomes that will be used in the gc correction
    filtered_gcContent =  options.tmpfolder+"/hs37d5_filtered.gc"
    filter_file(gcContent, filtered_gcContent, affectedrau)


    path_to_script = ("/".join(pathofexecarray))+"/src/bam_loess.R"
    command = 'Rscript'
    args = [filtered_gcContent, filtered_count_file, results_path]  #Filtered files
    cmd = [command, path_to_script] + args
    print("Running the following command: \n" + command + path_to_script + " " + filtered_gcContent + " " + filtered_count_file + " " + results_path)
    outputGC=subprocess.check_output(cmd, universal_newlines = True)    #This has corrected mappability and GC biases


    
    outputGC=outputGC.split('\n')
    arrayGCcorr=[]
    j=0
    for i in range(0, len(outputGC)-1, 2):
        arrayGCcorr.append([outputGC[i], outputGC[i+1]])
        j+=1
    #Rscript /home/projects/marthe/karyokan_ang/src/bam_loess.R /home/projects/marthe/karyokan_ang/data/hs37d5.gc /home/projects/marthe/karyokan_ang/data/filtered_counts_nonAffectedMale_10.bam


    gcbinsDict = gcBins(arrayGCcorr)


    #####################################################
    #In this part I access the different GC bins
    #Then I check which gc values are inside each bin
    #I store the positions of those gc bins in a dictionary
    #I store the observed read counts of those gc bins in a dict
    #I store the observed read counts in a file to use with R



    gcbinsP_file = open(options.tmpfolder+"/observed_bins","w")
    
    tmp_pos = []
    tmp_read = []
    gcKeys = []
    for key in gcbinsDict:
        for row in range(len(gcranges)):
            #Select the correct gc bin
            if (float(gcranges[row]["gc"])>=float(key)) and (float(gcranges[row]["gc"])<float(key)+0.1):
                tmp_pos.append(row)
                tmp_read.append(str(int(countInfo[row]["corrected_read"])))
        
        if len(tmp_read) > 0:
            lambda_v = gcbinsDict[key]  #Save the value before creating the empty dictionary
            gcbinsDict[key] = {}
            gcbinsP_file.write("\t".join(tmp_read) + "\n")
            gcbinsDict[key]["lambda"] = lambda_v    #It is the predicted read count for that GC bin
            gcbinsDict[key]["pos"] = tmp_pos
            gcbinsDict[key]["read"] = tmp_read
            gcKeys.append(key)

    
            #Reset variables
            tmp_pos = []
            tmp_read = []
    gcbinsP_file.close()


    

    ########################################################
    #In this part I will have several arrays for the different gc bins
    #I will do the mle to estimate p for each bin
    #I will save the p for the different gc bins
    #Then, I need to have something where I save, which raus are included in each bin,
    #so I know which p I have to use for each rau


    path_to_script = "/".join(pathofexecarray)+"/src/p_estimation.R"
    gcbinsP_file = options.tmpfolder+"/observed_bins"
    command = 'Rscript'
    args = [gcbinsP_file] 
    cmd = [command, path_to_script] + args
    print("Running the following command: \n" + command + " " + path_to_script + " " + gcbinsP_file)
    p=subprocess.check_output(cmd, universal_newlines = True)
    p = p.split("\n")
    p = [float(item) if item != 'NA' else np.nan for item in p]

    #This p contains the p-estimation, mean and variance, so the list has to be divided in lists of 3 elements each
    info_list = []
    for i in range(0, len(p), 4):
        info_list.append(p[i:i+4])

    


    for i in range(len(gcKeys)):
        #The p estimations and the keys are in the same order
        gcbinsDict[gcKeys[i]]["p_est"] = info_list[i][0]
        gcbinsDict[gcKeys[i]]["r"] = info_list[i][1]
        gcbinsDict[gcKeys[i]]["mean"] = info_list[i][2]
        gcbinsDict[gcKeys[i]]["std"] = info_list[i][3]
    

    ##########################################################
    ##########################################################
    #It creates a dataframe with the gc bins, the p estimator, r, mean and variance.
    #This datafram will be stored in a csv file and read into R or spyder manually, or even in Excel
    list_lambda = []; list_p = []; list_r = []; list_std = []; list_gc = []
    for key, values in gcbinsDict.items():
        if gcbinsDict[key]!=-1:
            list_gc.append(key)
            list_lambda.append(str(gcbinsDict[key]["mean"]))
            list_p.append(str(gcbinsDict[key]["p_est"]))
            list_r.append(str(gcbinsDict[key]["r"]))
            list_std.append(str(gcbinsDict[key]["std"]))



    plot_file = open(options.tmpfolder + "/info_toplot.txt", "w")
    plot_file.write("\t".join(list_gc)+"\n")
    plot_file.write("\t".join(list_lambda)+"\n")
    plot_file.write("\t".join(list_std)+"\n")
    plot_file.write("\t".join(list_p)+"\n")
    plot_file.write("\t".join(list_r)+"\n")
    plot_file.close()



    #We call Rscript for plotting the staff
    path_to_script = "/".join(pathofexecarray)+"/src/plots.R"
    file_toplot = options.tmpfolder+"/info_toplot.txt"
    command = 'Rscript'
    args = [file_toplot, results_path] 
    cmd = [command, path_to_script] + args
    print("Running the following command:\n" + command + " " + path_to_script + " " + file_toplot)
    p=subprocess.check_output(cmd, universal_newlines = True)

    ##########################################################
    ##########################################################



    #################################################
    #we now have the GC correction    
    #################################################

    maxpostllik=-1.0*sys.float_info.max 
    maxpostllikCD=-1
    allpostllik=[] 
    allpostllikNoMax=[] 



    print ("Evaluating every condition")


    for cd in tqdm(range(0, len(conditions))): #for each condition
    #for cd in tqdm(range(0,1)):


        print ("\nEvaluating "+conditions[cd]["namecond"])
        #computing the lambda for non affected autosomal portions

        ########################################################################################################################
        # first compute the expected autosomal coverage for unaffected regions assuming the individual has condition "cd"
        ########################################################################################################################
        #compute GC bias and lambda on a per GC basis               

        #
        autoLambdaSum=0 
        autoLambdaN=0 
        lllikelihood=0 
        lllikelihood9 = 0



        #################################################
        #compute variance lambda parameters

        #precomputing the loglike for the poisson of most coverage at each GC                
        #currentgenlambda=autolambda 
        #GC * cov
        #this should be the probability of observing some coverage even if the segment was deleted
        #so a poisson using lambda mismappign rate


        #for each genomic range in the condition
        for rau in range(len(conditions[cd]["chrrange"])):    #conditions[cd]["chrrange"]={ "chr", "start", "end", "ecov":1 }    


            #the GC content is not defined
            if(gcranges[rau]["gc"] == -1):  #List containing dictionaries { "chr","start", "end", "gc": } 
                continue 

            #find the gc bin for that region
            gcround = f'{float(gcranges[rau]["gc"]):.2f}'
            #gcround_lower = float(gcround)
            #gcround_upper = gcround_lower + 0.01

            if gcbinsDict[gcround] == "NA":
                continue
            if( gcbinsDict[gcround] == -1 ): #GC bin was not defined by the Loess
                continue 
            
            
            if countInfo[rau]["map_score"]!=0:
                regionlambda = gcbinsDict[gcround]["lambda"]/countInfo[rau]["map_score"]
            else:
                regionlambda = gcbinsDict[gcround]["lambda"]
                #regionlambda = gcbinsDict[gcround]["lambda"]/0.01


            #print(regionlambda, int(countInfo[rau]["corrected_read"]))

            #the mismapping rate is the minimum rate at which we will see some coverage by mistake
            #computing the log likelihood for the given genomic window
            #This will be used to calculate the posterior probabilities

            #For unaffected male            
            if( conditions[cd]["chrrange"][rau]["ecov"] == 1): # if we have for this condition the normal expected diplod
                predicted = conditions[cd]["chrrange"][rau]["ecov"]*regionlambda

            else:   #if for this condition we don´t have the normal expected diploid
                if( conditions[cd]["chrrange"][rau]["ecov"] < mismappingrate):#if the expected coverage is so low that it cannot be distinguished from mismappings                   
                #this happens when ecov is 0
                    predicted = mismappingrate*regionlambda
                else:
                    predicted = conditions[cd]["chrrange"][rau]["ecov"]*regionlambda
            
            #Calculate variance 

            #I set the number of successes (n) asa high value for high coverage. For example n=50
            #For lower coverages, I will modify the number of successes
            #Then I calculate the probability of success using the predicted and the number of successes
            
            
            p = gcbinsDict[gcround]["p_est"]
            if p == 1:
                p = 0.99
            r = (predicted*p)/(1-p)

            llik = nbinom.logpmf(int(countInfo[rau]["corrected_read"]), r, p)
            #print(rau, gcranges[rau]["chr"], int(countInfo[rau]["corrected_read"]), int(predicted), conditions[cd]["chrrange"][rau]["ecov"],  llik)
            #print(rau, gcranges[rau]["chr"], int(countInfo[rau]["corrected_read"]), int(regionlambda), conditions[cd]["chrrange"][rau]["ecov"],  int(predicted), llik)
            #print(rau, gcranges[rau]["chr"], int(countInfo[rau]["corrected_read"]), int(predicted), llik)

        


            #if conditions[8]["chrrange"][rau]["ecov"] != conditions[cd]["chrrange"][rau]["ecov"]:
            #if llik != llik9:
                #print(rau, gcranges[rau]["chr"], countInfo[rau]["count"], round(predicted,2), round(predicted9,2), round(llik,2), round(llik9,2))


            #print(rau, gcranges[rau]["chr"], str(countInfo[rau]["corrected_read"]), str(predicted), p, llik)




            ########################################################################################
            #This is used for comparisons, remove when it is not needed
            ########################################################################################

            #42 klinefelter 
            #32 Smith Magenis_m
            #-1 XYY, Jacobs
            #36 22q_13_deletion_m


            
            if np.isnan(llik):	#treat NA values as 0, otherwise it will always be NA
                llik=0

            lllikelihood += llik


        #print(lllikelihood, lllikelihood9)

        #prior
        logprior=log(conditions[cd]["priorprob"]) 
        postllik=lllikelihood+logprior #prior*likelihood = log(prior)+log(likelihood)

        #postllik=postllik/len(conditions[cd]["chrrange"])
        
        #Posterior Distribution = Prior Distribution + Likelihood Function (“new evidence”)
        #The posterior probability is calculated like priorprobability + new evidence 
        #The new evidence is the likelihood function, which tells you what information is contained in your observed data               
        allpostllik.append(postllik) #all "posterior probs"
        

        

        #found a new most likely condition
        if( postllik > maxpostllik):
            maxpostllik   = postllik 
            maxpostllikCD = cd  #most likely condition


        print ("log-likelihood: "+str(lllikelihood)) 
        print ("log(prior): "+str(logprior)) 
        print ("lg(postllik): "+str(postllik)) 



    print("###################################################################################")
    print("\nThe condition with the max post is", conditions[maxpostllikCD]["namecond"])
    print("The original bam file was", name_bam)
    print("###################################################################################\n")
    

    #TODO REENABLE
    #end of loop for each condition
    #Here we need to compute the actual posterior probability
    #we need the probability of error

    #########################################################################
    #Calculate the error 
    #########################################################################


    #First method: error = best model/sum models
    #Second method: error = (sum models - best model)/sum models

    if(True):    
        #compute an array of values without the maximum value of the (likelihood*prior=posterior)
        for i in range(0,len(allpostllik)):
            if i != maxpostllikCD:
                allpostllikNoMax.append(allpostllik[i]) 

        if( len(allpostllik) < 2):
            sys.stderr.write("Error, not enough conditions\n") 
            sys.exit(1) 

        if( len(allpostllikNoMax) < 2):
            sys.stderr.write("Error, not enough conditions\n") 
            sys.exit(1) 

        for post in allpostllik:
            if post =='nan':
                sys.stderr.write("Some of the calculated posterior probabilities contain NAN values\nFurther analysis cannot be done")
                sys.exit()
        
        #If I arrive to this point is because there are no NAN values

        #Get the sum of the posteriors including the maximum
        sumpost =np.logaddexp(allpostllik[0], allpostllik[1])
        for i in range(2, len(allpostllik)):
            sumpost = np.logaddexp(sumpost, allpostllik[i]) 



        result_f = allpostllik[0]-sumpost    #female
        result_m = allpostllik[1]-sumpost    #male
        #result = allpostllik[2]-sumpost    #Down_m
        result = allpostllik[42]-sumpost    #Klinefelter
        #result = allpostllik[46]-sumpost    #Jacobs
        
        #result = allpostllik[43]-sumpost    #Turner

        posttoprint = [str(exp(allpostllik[i]-sumpost)) for i in range(len(allpostllik))]

        counter = 0
        for i in range(len(posttoprint)):
            if float(posttoprint[i])>0.9:
                print("The condition {} has a posterior probability of {} above 90%\nTherefore, we can say that it is the correct condition of the sample,\n".format(str(conditions[maxpostllikCD]["namecond"]), posttoprint[i]))
                max_post = posttoprint[i]
                counter += 1
            elif counter == len(posttoprint):
                print("None of the conditions has a posterior probability above 90%\n Therefore, we cannot know for sure which was the karyotype of the sample.\n")
            

        #Recover this part later
        """outputfilek.write("Conditions" + "\t" + "Prior prob" + "\t" + "Posterior prob" + "\n")

        for i in range(len(conditions)):
            outputfilek.write(conditions[i]["namecond"] + "\t" + str(conditions[i]["priorprob"]) + "\t" + posttoprint[i] + "\n")
        print("The result file with the prior and posterior probabilities for the conditions has been generated.\n")"""


        #print("the total is", sumpost)
        #print("the logarithm is", result)
        #print("the exponent is", math.exp(result))


        results_post = open( ("/".join(pathofexecarray))+"/data/results_post_mix", "a")
        results_post.write(str(max_post) + "\n" + str(conditions[maxpostllikCD]["namecond"]) + "\t" + name_bam + "\n")
        #results_post.write(str(math.exp(result_m)) + " " + str(math.exp(result_f)) + "\n" + str(conditions[maxpostllikCD]["namecond"]) + "\n" + name_bam +"\n")
        
        #results_post.write("The condition with the max post is " + conditions[maxpostllikCD]["namecond"] + "\n")
        #results_post.write("The original bam file was " + name_bam + "\n")
        #results_post.write("mean of reads: " + str(mean_count)+ "\n")

        results_post.close()



#end of loop for each bam file
print ("Printing results to: " + options.resultso)

outputfilek.close() 

print( "Recursively removing temp directory "+tfm)
#cmdrm = rmcmd+" -rfv "+tfm

#Remove files to save space in disk
os.remove(filtered_count_file)


sys.stderr.write("\nProgram finished succesfully\n")
