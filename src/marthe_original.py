#!/usr/bin/env python

import subprocess
from optparse import OptionParser
import sys, os, random
import time 
import numpy as np
import re
import copy
import tempfile 
from tqdm import tqdm
from math import exp, log
from scipy.stats import poisson
from function import *
import random
import optparse



mismappingrate=0.01


pathofexec  = os.path.abspath(sys.argv[0]) 
pathofexecarray=pathofexec.split('/')[:-2] 
pathofconfig=  ("/".join(pathofexecarray))+"/config/homosapiens_hg19.cfg" 
if(not os.path.exists(pathofconfig)):
    sys.stderr.write("\nERROR: The configuration file "+pathofconfig+" does not exist, re-clone the repository\n") 
    sys.exit(1) 
        

pathoffraglength=  ("/".join(pathofexecarray))+"/utils/baselengthFilter/src/fragLength" 
if(not os.path.exists(pathoffraglength)):
    sys.stderr.write("\nERROR: The executable file "+pathoffraglength+" does not exist, please type make in the utils directory\n") 
    sys.exit(1) 


pathofCountFragment=  ("/".join(pathofexecarray))+"/utils/countFrags/src/countFrags" 
if(not os.path.exists(pathoffraglength)):
    sys.stderr.write("\nERROR: The executable file "+pathoffraglength+" does not exist, please type make in the utils directory\n") 
    sys.exit(1) 


pathofparseCount=  ("/".join(pathofexecarray))+"/src/parseCount.py" 
if(not os.path.exists(pathofparseCount)):
    sys.stderr.write("\nERROR: The executable file "+pathofparseCount+" does not exist, please type make in the utils directory\n") 
    sys.exit(1) 

gcContent=  ("/".join(pathofexecarray))+"/data/hs37d5.gc" 
if(not os.path.exists(gcContent)):
    sys.stderr.write("\nERROR: The GC file "+gcContent+" does not exist\n") 
    sys.exit(1) 

    

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

"Bayesian inference of karyotype for low coverage data\n"+

"marthe [options] [bam file 1] [bam file 2]...") 


group_m = optparse.OptionGroup(parser, "Mandatory Options",
                    "These are the files that have to be used so the program can run")
group_m.add_option("-r", "--reference",      dest="ref",          help="Reference genome used for alignment",                 default=None,    type="string") 
group_m.add_option("--gc",                   dest="gcfile",       help="File containing the % of GC per window size (see above) default: "+gcContent+"",                  default=gcContent,    type="string") 
group_m.add_option("-c", "--conf",           dest="configfile",   help="Configuration for various conditions file to use for species/build default: "+str(pathofconfig)+"", default=pathofconfig, type="string") 
group_m.add_option("--winsize",              dest="winsize",      help="Size of genomic windows, default is 1Mbp",             default=1000000, type="int") 
group_m.add_option("--mismap",               dest="mismappingrate", help="Mismapping rate , default is "+str(mismappingrate)+"",            default=mismappingrate, type="float") 


parser.add_option_group(group_m)

group_i = optparse.OptionGroup(parser, "To improve options",
                                "These are the options in order to improve the performance")
group_i.add_option("--map",                  dest="mappability",     help="Use a mappability map in BED format (Recommended)",   default=None,    type="string") 
group_i.add_option("--minlen",               dest="minlen",       help="Minimum length of DNA fragment",             default=35, type="int") 
parser.add_option_group(group_i)


group_t = optparse.OptionGroup(parser, "Technical options",
                                "These are the options in order to improve the performance")
group_t.add_option("--tmpfolder",            dest="tmpfolder",    help="Temporary folder",                                    default="/tmp/",    type="string") 
group_t.add_option("-o","--outfile",         dest="resultso",     help="Output prefix to store the results",                  default="results",    type="string") 
group_t.add_option("--samtools",             dest="samtools",     help="Use this version of samtools",                        default="samtools",    type="string") 
group_t.add_option("--tabix",                dest="tabix",        help="Use this version of tabix",                           default="tabix",    type="string") 
group_t.add_option("--hpc",                  dest="hpc",          help="Use high-performance computing (for queueing systems ex: SGE)",          action="store_true") 
group_t.add_option("--resume",               dest="resume",       help="Resume by providing the temp directory used",                              type="string") 
group_t.add_option("--nice",                 dest="nice",         help="Nice the jobs",                                                          action="store_true") 
group_t.add_option("-t"  , "--threads",      dest="threads",      help="Number of threads to use if the local machine is used",  default=1,   type="int") 
group_t.add_option("--bamidx",               dest="bamidx",       help="Index of the BAM file to process",                              type="int")
group_t.add_option("-F",                     dest="file",         help="There is a list with names of BAM files",             action="store_true") 

#This has to be removed later
group_t.add_option("--subsample",           dest="subsample",   help="Subsample value for the simulations", default=1, type="int") 
parser.add_option_group(group_t)

#TODO ?
#parser.add_option(""  , "--branchl",     dest="branchlscale", help="Seq-gen branch scale, default is 0.00045",            default=0.00045, type="float") 
#parser.add_option(""  , "--chrlen",      dest="lengthchr",    help="Chromosome length, default is 10kb",              default=10000,   type="int") 
#parser.add_option("-c", "--numcont",     dest="numcont",      help="Number of present-day human contaminants, default is 2", default=2, type="int")


(options,args) = parser.parse_args()
if( options.resume == None):
    if( len(args)==0 ):
        sys.stderr.write("\nPlease use -h to see options.\n") 
        sys.exit(1) 

if( options.mappability == None):
    sys.stderr.write("\nWarning: Using a mappability is highly recommended.\n") 

else:
    if(not os.path.exists(options.mappability)):
        sys.stderr.write("\nERROR: The mappability does not exist.\n") 
        sys.exit(1) 

    if(not os.path.exists(options.mappability+".tbi")):
        sys.stderr.write("\nERROR: The tabix index for the mappability does not exist.\n") 
        sys.exit(1) 



mismappingrate=options.mismappingrate 

#detecting programs
samtoolscmd     = re.sub(b'\s+', b'', which(options.samtools))
samtoolscmd = samtoolscmd.decode()
tabixcmd        = re.sub(b'\s+', b'', which(options.tabix)) 
tabixcmd = tabixcmd.decode()
rmcmd     = "rm" 


if(options.ref ==  None):
    sys.stderr.write("Please specify the reference used for alignment, options -r or --reference") 
    sys.exit(1) 


if(not os.path.exists(options.ref+".fai")):
    sys.stderr.write("\nWarning index file:"+options.ref+".fai does not exist, we will create it, make sure you have permissions to write in that directory. If not, create a symlink.\n") 
    cmdindex = samtoolscmd + " faidx " + options.ref 
    subprocess.run(cmdindex, shell=True, universal_newlines=True)



chrnames   = {} 
chrranges  = [] 

#countfile = options.tmpfolder + "/counts_windows_" + str(options.winsize)

######################
#READING fai FILE #
######################

referencefai = open(options.ref+".fai", "r") 
#this code stores the start and end coordinates of the chromosomes found from samtools faidx 
#this info is stored in chrranges
listchr = list(range(1,23))+["X", "Y", "x", "y"]
listchr= [str(option) for option in listchr]

pathofregionsfile = options.tmpfolder + "/regions_windows"
regions_file = open(pathofregionsfile, "w")

#pathofregionsfile2 = options.tmpfolder + "/regions_windows_2"
#regions_file2 = open(pathofregionsfile2, "w")


for linefai in referencefai:

    faia=linefai.split("\t") 
    if faia[0] in listchr:
        #chrnames.append(   faia[0] ) 
        chrnames[ faia[0] ] = 0 

        #chrlengths.append( faia[1] ) 
        for c in range(1, (int(faia[1])-options.winsize), options.winsize):
            rangedict =	{ "chr": faia[0], "start": c,   "end": c+(options.winsize-1), "ecov":1 }     
            regions_file.write(str(faia[0])+":"+str(c)+"-"+str(c+(options.winsize-1))+"\n")
            #regions_file2.write("chr"+str(faia[0])+":"+str(c)+"-"+str(c+(options.winsize-1))+"\n")

            chrranges.append(rangedict) 
    
referencefai.close() 
regions_file.close()
#regions_file2.close()


#create the list of commands to run to evaluate coverage
#This has 2 blocks,     
#if( options.resume == None): is executed first to launch the coverage calculation
#else: #is executed after and launches parseCov for each bamfile


uniquebamfilesArray=[]
counterIndex=0

while True:

    if(len(args)==0):
        sys.stderr.write("ERROR: specify at least 1 BAM/CRAM file or file with list.\n") 
        sys.exit(0) 
    
    
    if options.file:
        print("hi")
        infile = open(args[0], "r")
        args_list = []
        for line in infile:
            line = checkLineJump(line)
            args_list.append(line)
        infile.close()
    
    else:
        args_list = args



    if( options.resume == None):
        if(options.tmpfolder == "/tmp/"):
            tfm=tempfile.mkdtemp(prefix = options.tmpfolder)
            sys.stderr.write("Creating temp directory "+tfm+"\n") 
        else:
            tfm = options.tmpfolder 
            if not os.path.exists( tfm ):
                os.mkdir( tfm , 0o755 ) 
                sys.stderr.write("Creating temp directory "+tfm+"\n") 
            else:
                sys.stderr.write( "\nTemp folder already exists: "+str(tfm) +"\n" )       
                     

        #iterate for each bam file
        myfilecounter=0 
        arraycmds=[] 
        fileHandleLC = open ( ""+tfm+"/listcommands1.txt",   'w' )  
        fileHandleIndex = open ( ""+tfm+"/listcommandsIndex.txt",   'w' )  
        

        for i in range(len(args_list)):

            bamfile = str(args_list[i]) 
            bamfile = os.path.abspath( bamfile )  # resolve path
            
            uniquebamfilesArray.append(bamfile)
            termination = bamfile.split(".")[-1]
            name_bam_complete = bamfile.split("/")[-1]
            name_bam = name_bam_complete.split(".")[0]
            
            #First check if the BAM or CRAM file exists
            if(not os.path.exists(bamfile)):
                sys.stderr.write("\nERROR File: " + bamfile + " does not exist, skipping\n") 
                continue 

                #sys.exit(0) 



            #Then it checks which BAM (or converted CRAM) do not have index
            #Write in a list and do in parallel
            if(not os.path.exists(bamfile+".bai")):
                sys.stderr.write("\nWarning Index file:" + bamfile + ".bai does not exist, we will create it. \nMake sure you have permissions to write in that directory, if not, create a symlink.\n") 
                cmdindex = samtoolscmd + " index " + bamfile 
                fileHandleIndex.write(cmdindex + "\n")
                counterIndex += 1

            #Check reference of bam file:
            #It goes to the second column of the second file and check if it is NS:1 or SN:chr1
            """cmdReference = samtoolscmd + " view -H " + bamfile + "  | sed -n ""2p"" | awk '{print $2}'"

            output=subprocess.check_output(cmdReference, shell=True, universal_newlines = True)
            output = checkLineJump(output)
            if output == "SN:1":
                pathofregionsfile = pathofregionsfile1
            elif output == "SN:chr1":
                pathofregionsfile = pathofregionsfile2
            else:
                sys.stderr.write("\nThe format of the BAM file must be SN:1 or SN:chr1 for the second column of the second line\n")"""

                
            command = ""
            if options.nice:
                command += "nice -19 "

            if options.mappability == None:
                command += pathofCountFragment + " -m " + str(options.minlen)                                + " " + bamfile + " " + pathofregionsfile + " |  gzip > " + options.tmpfolder + "/"+ str(i) + "_counts_windows_" + name_bam + ".gz"
            else:
                command += pathofCountFragment + " -m " + str(options.minlen) + " -L " + options.mappability + " " + bamfile + " " + pathofregionsfile + " |  gzip > " + options.tmpfolder + "/"+ str(i) + "_counts_windows_" + name_bam + ".gz"


            fileHandleLC.write(command +  "\n")



        fileHandleIndex.close()
        fileHandleLC.close()
        
        
        #Then, index the BAM (and converted CRAM) files
        if counterIndex > 0:
            cmdIndex="cat " + tfm + "/listcommandsIndex.txt | nice -19 parallel --slf  /home/projects/marthe/marthe/data/list_serv"
            #cmdIndex="cat " + tfm + "/listcommandsIndex.txt | parallel -j "+str(options.threads)

            print("Indexing BAM files")
            #subprocess.run(cmdIndex, shell=True, universal_newlines=True)

        #Count the reads for the BAM files and converted BAM files

        cmdtolaunch="cat " + tfm + "/listcommands1.txt | nice -19 parallel --slf  /home/projects/marthe/marthe/data/list_serv"
        #cmdIndex="cat " + tfm + "/listcommands1.txt | parallel -j "+str(options.threads)

        print("Counting the number of reads")
        if options.hpc:
            print ("Please run the commands manually either using:")
            print ("cat " + tfm + "/listcommands1.txt | parallel -j "+str(options.threads))
            print ("on the use a batch/queueing system to launch:")
            print ("cat " + tfm + "/listcommands1.txt | sbatch ...\n")
            print ("Once commands are done, rerun with:\n")
            if options.file:
                print ("python3 " + "/".join(pathofexecarray)+"/src/marthe.py   -r " + options.ref + "  -c " + str(options.configfile) + " -F --resume " + tfm + " --tmpfolder " + options.tmpfolder + " --hpc  -o " + options.resultso + " " + args[0]  + '\n')

            else:
                print ("python3 " + "/".join(pathofexecarray)+"/src/marthe.py   -r " + options.ref + "  -c " + str(options.configfile) + " --resume " + tfm + " --tmpfolder " + options.tmpfolder + " --hpc  -o " + options.resultso + " " + args[0]  + '\n')
            print ("exiting")
            break 
        else:
            print("Launching commands:\n" + cmdtolaunch)
            subprocess.run(cmdtolaunch, shell=True, universal_newlines=True)

            options.resume = tfm 
            print(options.resume)
            continue  #skip to next block and pretend it is a resume


    else: #resume, need to parse the count

        tfm=options.resume 
        print( "Generating listcommands2.txt")

        if options.mappability != None:
            fileHandleScore = open("" + tfm + "/listcommandsCount.txt", 'w' )  
            map_scores_file = tfm + "/map_score.txt"

            for r in chrranges:
                cmd_tabix = tabixcmd + " " + options.mappability + " " + str(r["chr"]) + ":" + str(r["start"]) + "-" + str(r["end"]) + " | "
                cmd_diff = "awk '!(NR%2){if ($3-$2 >=" + str(options.minlen) + ") {print int(($3-$2)/35)}}' |"
                cmd_sum =  "awk '{s+=$1} END {print (s*"+ str(options.minlen) + ")/" + str(options.winsize) + "}' >> " + map_scores_file
                command = cmd_tabix + cmd_diff + cmd_sum
                #subprocess.run(command, shell=True, universal_newlines=True)
            


        fileHandleLC = open ( "" + tfm + "/listcommands2.txt",   'w' ) 

        for i in range(len(args_list)):

            cmdtolaunch = "" 
            if(options.nice):
                cmdtolaunch = "nice -19  " 
            cmdtolaunch += "python3 " + pathofparseCount + " --bamidx " + str(i)

            #cmdtolaunch += " --resume " + options.tmpfolder + " " 
            for idx, arg in enumerate(sys.argv):

                if(idx > 0):
                    cmdtolaunch += " " + arg + " "
            fileHandleLC.write(cmdtolaunch + "\n") 

        fileHandleLC.close()
        
        cmdtolaunch="cat "+tfm+"/listcommands2.txt |  nice -19 parallel --slf /home/projects/marthe/marthe/data/list_serv"
        


        print("listcommands2.txt generated")



        if(options.hpc):

            print ("Please run the commands manually either using:")
            print ("cat " + tfm + "/listcommands2.txt | nice -19 parallel -j " + str(options.threads))
            print ("on the use a batch/queueing system to launch:")
            print ("cat " + tfm + "/listcommands2.txt | sbatch ...\n")
            print ("exiting")
            break 

        else:
            print("Launching commands:\n" + cmdtolaunch)
            subprocess.run(cmdtolaunch, shell=True, universal_newlines=True)

            options.resume = tfm 
            break  #skip to next block and pretend it is a resume
            #proceed as resume afterwards


        print ("Recursively removing temp directory " + tfm )
        cmdrm = rmcmd + " -rfv " + tfm


        break 


sys.exit(0)