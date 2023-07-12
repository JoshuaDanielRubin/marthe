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



mismappingrate=0.01




pathofexec  = os.path.abspath(sys.argv[0]) 
pathofexecarray=pathofexec.split('/')[:-2] 
pathofconfig=  ("/".join(pathofexecarray))+"/config/homosapiens_hg19.cfg" 
if(not os.path.exists(pathofconfig)):
    sys.stderr.write("\nERROR: The configuration file "+pathofconfig+" does not exist, re-clone the repository\n") 
    sys.exit(1) 
        

"""pathofsortuniq=  ("/".join(pathofexecarray))+"/utils/sortuniq/sortuniq" 
if(not os.path.exists(pathofsortuniq)):
    sys.stderr.write("\nERROR: The executable file "+pathofsortuniq+" does not exist, please type make in the utils directory\n") 
    sys.exit(1) """

pathoffraglength=  ("/".join(pathofexecarray))+"/utils/baselengthFilter/src/fragLength" 
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

    
#print(pathofconfig) 
parser = OptionParser(
"\n\n"+

" ##    ##    ###    ########  ##    ##  #######  ##    ##    ###    ##    ## \n"+
" ##   ##    ## ##   ##     ##  ##  ##  ##     ## ##   ##    ## ##   ###   ## \n"+
" ##  ##    ##   ##  ##     ##   ####   ##     ## ##  ##    ##   ##  ####  ## \n"+
" #####    ##     ## ########     ##    ##     ## #####    ##     ## ## ## ## \n"+
" ##  ##   ######### ##   ##      ##    ##     ## ##  ##   ######### ##  #### \n"+
" ##   ##  ##     ## ##    ##     ##    ##     ## ##   ##  ##     ## ##   ### \n"+
" ##    ## ##     ## ##     ##    ##     #######  ##    ## ##     ## ##    ## \n"+
"\n"+

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

#TODO ?
#parser.add_option(""  , "--branchl",     dest="branchlscale", help="Seq-gen branch scale, default is 0.00045",            default=0.00045, type="float") 
#parser.add_option(""  , "--chrlen",      dest="lengthchr",    help="Chromosome length, default is 10kb",              default=10000,   type="int") 
#parser.add_option("-c", "--numcont",     dest="numcont",      help="Number of present-day human contaminants, default is 2", default=2, type="int")
parser.add_option("--winsize",              dest="winsize",      help="Size of genomic windows, default is 1Mbp",             default=1000000, type="int") 
parser.add_option("--minlen",               dest="minlen",       help="Minimum length of DNA fragment",             default=35, type="int") 
parser.add_option("--gc",                   dest="gcfile",       help="File containing the % of GC per window size (see above) default: "+gcContent+"",                  default=gcContent,    type="string") 
parser.add_option("-c", "--conf",           dest="configfile",   help="Configuration for various conditions file to use for species/build default: "+str(pathofconfig)+"", default=pathofconfig, type="string") 

#This has to be removed later
parser.add_option("--subsample",           dest="subsample",   help="Subsample value for the simulations", default=1, type="int") 



(options,args) = parser.parse_args()
if( options.resume == None):
    if( len(args)==0 ):
        sys.stderr.write("\nPlease use -h to see options.\n") 
        sys.exit(1) 

if( options.map_mask == None):
    sys.stderr.write("\nWarning: Using a mappability is highly recommended.\n") 

else:
    if(not os.path.exists(options.map_mask)):
        sys.stderr.write("\nERROR: The mappability does not exist.\n") 
        sys.exit(1) 

    if(not os.path.exists(options.map_mask+".tbi")):
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
#TODO may be important
#chrlengths = [] 
chrranges  = [] 



######################
#READING fai FILE #
######################

referencefai = open(options.ref+".fai", "r") 
#this code stores the start and end coordinates of the chromosomes found from samtools faidx 
#this info is stored in chrranges
listchr = list(range(1,23))+["X", "Y", "x", "y"]
listchr= [str(option) for option in listchr]

regions_file = open(options.tmpfolder + "/regions_windows", "w")


for linefai in referencefai:

    faia=linefai.split("\t") 
    if faia[0] in listchr:
        #chrnames.append(   faia[0] ) 
        chrnames[ faia[0] ] = 0 

        #chrlengths.append( faia[1] ) 
        for c in range(1, (int(faia[1])-options.winsize), options.winsize):
            rangedict =	{ "chr": faia[0], "start": c,   "end": c+(options.winsize-1), "ecov":1 }     
            regions_file.write(str(faia[0])+":"+str(c)+"-"+str(c+(options.winsize-1)))

            chrranges.append(rangedict) 
    
referencefai.close() 
regions_file.close()


#create the list of commands to run to evaluate coverage
#This has 2 blocks,     
#if( options.resume == None): is executed first to launch the coverage calculation
#else: #is executed after and launches parseCov for each bamfile

while True:
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
        print("Generating listcommands1.txt and correspondance.dat")
        fileHandleLC = open ( ""+tfm+"/listcommands1.txt",   'w' )  
        fileHandleCR = open ( ""+tfm+"/correspondance.dat", 'w' )  

        if(len(args)==0):
            sys.stderr.write("ERROR: specify at least 1 BAM file.\n") 
            sys.exit(0) 

        for i in range(len(args)):

            bamfile = str(args[i]) 
            bamfile = os.path.abspath( bamfile )  # resolve path
            
            if(not os.path.exists(bamfile)):
                sys.stderr.write("\nERROR File: " + bamfile + " does not exist, skipping\n") 
                continue 
                #sys.exit(0) 

            #todo make it parallel
            if(not os.path.exists(bamfile+".bai")):
                sys.stderr.write("\nWarning Index file:" + bamfile + ".bai does not exist, we will create it. \nMake sure you have permissions to write in that directory, if not, create a symlink.\n") 
                cmdindex = samtoolscmd + " index " + bamfile 
                handle_job( cmdindex ) 
        
            #fraction_list = [".9", ".5", ".1", ".05", ".01", ".005", ".001", ".0005", ".0001", ".00005", ".00001", ".000005", ".000001", ".0000005", ".0000001", ".00000009", ".00000008", ".00000007", ".00000006"]
            fraction_list = [".5", ".05", ".005",".0005", ".00005", ".000005", ".0000005", ".00000009", ".00000008", ".00000007", ".00000006"]
            index_subsample = options.subsample-1
            final_coverage = fraction_list[index_subsample]
            #UNCOMMENT THIS LATER
            #subsample_value = get_subsample(bamfile, final_coverage)
            subsample_value = ".99"


            #print("The subsample fraction is", fraction_list[subsample])

            seed = random.randint(0, 100)
            #create read count jobs
            for r in chrranges:
                
                #create commands for the different options
                nice_section = "nice -19 "
                
                                
                map_unmap_section=samtoolscmd+" view  -u " + bamfile + "  "+str(r["chr"]) + ":"+str(r["start"])+"-"+str(r["end"])+" |  "+pathoffraglength+ " -q -l "+str(options.minlen)+" -u -  > "+tfm+"/"+str(myfilecounter)+".bam"
                #for subsampling
                mapped_section= samtoolscmd+" view  -u " + bamfile + " -s " + str(seed) + fraction_list[index_subsample] + " " + str(r["chr"])+":"+str(r["start"])+"-"+str(r["end"])+" |  "+pathoffraglength+ " -q -l "+str(options.minlen)+" -m -u -  > "+tfm+"/"+str(myfilecounter)+".bam"  
                #mapped_section= samtoolscmd+" view  -u " + bamfile + " -s " + str(seed) + subsample_value + " " + str(r["chr"])+":"+str(r["start"])+"-"+str(r["end"])+" |  "+pathoffraglength+ " -q -l "+str(options.minlen)+" -m -u -  > "+tfm+"/"+str(myfilecounter)+".bam"  

                #not subsampling
                #mapped_section= samtoolscmd+" view  -u " + bamfile  + " " + str(r["chr"])+":"+str(r["start"])+"-"+str(r["end"])+" |  "+pathoffraglength+ " -q -l "+str(options.minlen)+" -m -u -  > "+tfm+"/"+str(myfilecounter)+".bam"  
                index_section= " && "+samtoolscmd+" index " + tfm + "/" + str(myfilecounter)+".bam " 

                #after filter reads (mapped sections)
                reads_section=" && " + samtoolscmd +" view " +tfm+"/"+str(myfilecounter)+".bam  | awk '{print $3" + ' "\\t" ' +"$4" + ' "\\t" '+ "$4+length($10)-1}' > "+tfm+"/"+str(myfilecounter)+"_reads.bam"
                
                count_section1=" && "+samtoolscmd+" view  -c -F 0x80 "
                count_section2=" -T  "+options.ref+" "+tfm+"/"+str(myfilecounter)+".bam  "+str(r["chr"])+":"+str(r["start"])+"-"+str(r["end"])

                devnull_section = " 2> /dev/null "

                create_file_section= " > " +  tfm+"/"+str(myfilecounter)+".count "
                append_file_section=" | tee -a  " + tfm+"/"+str(myfilecounter)+".count "
                
                rm_section=" &&  rm -f "+tfm+"/"+str(myfilecounter)+".bam && rm -f "+tfm+"/"+str(myfilecounter)+".bam.bai"
                

                cmdtolaunch="" 
                if(options.nice):
                    cmdtolaunch=nice_section

                cmdtolaunch+= mapped_section
                cmdtolaunch+=index_section
                    
                #option_count=samtoolscmd+" view  -c -F 0x80 "  #-c count  and -F 0x80 to remove the second pair of reads
                cmdtolaunch+=count_section1
                #obtain mapped reads (unique+nonunique)
                cmdtolaunch+=count_section2
                
                if options.map_mask != None:
                    cmdtolaunch+=" -L   <( "+tabixcmd+" "+options.map_mask+" "+str(r["chr"])+":"+str(r["start"])+"-"+str(r["end"])+" ) " 


                cmdtolaunch+=devnull_section
                cmdtolaunch+=create_file_section

                cmdtolaunch+= rm_section #Remove BAM and bam.bai files


                arraycmds.append({"cmd":cmdtolaunch,"region":""+str(r["chr"])+":"+str(r["start"])+"-"+str(r["end"]),"bamfile":bamfile}) 
                fileHandleLC.write(cmdtolaunch+"\n") 

                fileHandleCR.write(str(myfilecounter)+"\t"+str(r["chr"])+":"+str(r["start"])+"-"+str(r["end"])+"\t"+bamfile+"\n") 
                myfilecounter+=1 

        #compute coverage threads
        fileHandleCR.close() 
        fileHandleLC.close() 
        print("listcommands1.txt and correspondance.dat files created")



        if(options.hpc):
            print ("Please run the commands manually either using:")
            print ("cat "+tfm+"/listcommands1.txt | parallel -j "+str(options.threads))
            print ("on the use a batch/queueing system to launch:")
            print ("cat "+tfm+"/listcommands1.txt | sbatch ...\n")
            print ("Once commands are done, rerun with:\n")
            print ("marthe.py   -r "+options.ref+"  -c "+str(options.configfile)+ " --resume "+tfm+" --tmpfolder "+options.tmpfolder+"  -o "+options.resultso + '\n')
            print ("exiting")
            break 
        else:

            #cmdtolaunch="cat "+tfm+"/listcommands1.txt | parallel  -j "+str(options.threads) 
            cmdtolaunch="cat " + tfm + "/listcommands1.txt | nice -19 parallel --slf " + "/".join(pathofexecarray)+"/data/list_serv"
            #if this option is taken and the file list_serv doesn't exist, create it first or send a warning
            print ("Launching commands:\n"+cmdtolaunch)
            subprocess.run(cmdtolaunch, shell=True, universal_newlines=True)

            options.resume=tfm 
            continue  #skip to next block and pretend it is a resume


    else: #resume, need to parse the count

        tfm=options.resume 
        print( "Generating listcommands2.txt")
        crrp= open(tfm+"/correspondance.dat", "r") 
        datatoprocess=[] 

        uniquebamfilesArray=[]

        crrp= open(tfm+"/correspondance.dat", "r");


        for linec in crrp: #linec = numer_file, chromosome:start-end, bamfile
            corr=linec.split() 
            datatoprocess.append({"myfilecounter":corr[0],"region":corr[1],"bamfile":corr[2]})   
            if corr[2] not in uniquebamfilesArray:  
                uniquebamfilesArray.append(corr[2])

        
        bamid=0 

        fileHandleLC = open ( "" + tfm + "/listcommands2.txt",   'w' )  

        #Content of listcommands2.txt
        #nice -19 python3 /home/projects/marthe/marthe/src/parseCount.py --bamidx 0  --nice  -r  
        #/home/databases/genomes/Homo_sapiens/hs37d5/hs37d5.fa  --tmpfolder  /home/projects/marthe/tmp  
        #--resume  /home/projects/marthe/tmp   
        #/home/projects/gabriel/simKaryo/emp/12samples/Labrana.hg19.flt.sort.rmdup.realign.md.bam



        for b in uniquebamfilesArray:
            cmdtolaunch = "" 
            if(options.nice):
                cmdtolaunch = "nice -19 python3 " 
            cmdtolaunch += pathofparseCount + " --bamidx " + str(bamid) + " " 
            #cmdtolaunch += " --resume " + options.tmpfolder + " " 
            for idx, arg in enumerate(sys.argv):
                if(idx > 0):
                    cmdtolaunch += " " + arg + " " 
            fileHandleLC.write(cmdtolaunch + "\n") 
            bamid+=1 

        fileHandleLC.close() 
        print("listcommands2.txt generated")



        if(options.hpc):

            print ("Please run the commands manually either using:")
            print ("cat " + tfm + "/listcommands2.txt | nice -19 parallel -j " + str(options.threads))
            print ("on the use a batch/queueing system to launch:")
            print ("cat " + tfm + "/listcommands2.txt | sbatch ...\n")
            print ("Once commands are done, rerun with:\n")
            print ("marthe.py   -r " + options.ref + "  -c " + str(options.configfile) + " --resume " + tfm + " --tmpfolder " + options.tmpfolder + "  -o " + options.resultso + "\n")
            print ("exiting")
            break 

        else:
            #cmdtolaunch="cat "+tfm+"/listcommands2.txt | parallel  -j "+str(options.threads) 
            cmdtolaunch="cat "+tfm+"/listcommands2.txt |  nice -19 parallel --slf " + "/".join(pathofexecarray)+"/data/list_serv"
            print ("Launching commands:\n" + cmdtolaunch)
            subprocess.run(cmdtolaunch, shell=True, universal_newlines=True)

            options.resume = tfm 
            break  #skip to next block and pretend it is a resume
            #proceed as resume afterwards

        print ("Recursively removing temp directory " + tfm )
        cmdrm = rmcmd + " -rfv " + tfm
        #handle_job( cmdrm ) 

        break 


sys.exit(0)
