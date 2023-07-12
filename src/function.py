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
from scipy.stats import nbinom
import gzip

def isFloat(stringToCheck):
    """
    It checks if the string provided can be converted to float
    """
    try:
        float(stringToCheck)
        return True
    except ValueError:
        return False


def isInt(stringToCheck):
    """
    It checks if the string provided can be converted to integer
    """
    try:
        int(stringToCheck)
        return True
    except ValueError:
        return False

def checkLineJump(string):
    if string.endswith("\n"):
        string = string[:-1]
    return string


def gcBins(data):  
    """
    It creates a dictionary of GC bins and stores the mean of the counts for each bin

    Input: list containing gc values and predicted read counts
    Output: dictionary with GC values as keys and read counts as values
    """

    bins = list(np.arange(0, 1.00, 0.01))
    formatted_bins = [ '%.2f' % elem for elem in bins ]
    dictofrepetitions = dict.fromkeys(formatted_bins, 0)
    dictofcounts = dict.fromkeys(formatted_bins, 0)


    for i in range(0, len(data)):
        value = f'{float(data[i][0]):.2f}'  #gc value is in the first column
        dictofcounts[str(value)]+= float(data[i][1])
        dictofrepetitions[str(value)]+=1

    for key, value in dictofcounts.items():

        if dictofrepetitions[key] != 0:
            dictofcounts[key] = round(dictofcounts[key]/dictofrepetitions[key],3)
            
        else:
            dictofcounts[key] = -1

    return dictofcounts


    

def which(program):
    """
    This function detect the program and check if it is installed
    Input: name of the program
    Output: output statment 
    """
    sys.stderr.write("Detecting program: "+program+"\n") 
    cjob = "type "+program
    sp = subprocess.Popen(["/bin/bash", "-i", "-c", cjob],
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)
    out, err = sp.communicate() 
    errcode  = sp.returncode 
    if(errcode != 0): #has finished but wrong code
        sys.stderr.write("Cannot find program: "+program+" please make sure it is installed\n") 
        sys.exit(1) 
    if(out.find(b"aliased") != -1): #was aliased
        return out[out.find(b" to ")+5:-2] 
    else:
        if(out.find(b" is ") != -1):
            return out[out.find(b" is ")+4:] 
        else:
            sys.stderr.write("Cannot seem to find program: "+program+" please make sure it is installed\n") 
            sys.exit(1) 

    return out




def readlineconfig(configfile, starting_string, options):
    """
    This function reads the config file and extracts the chromosomes

    Input:
        configfile: config file already opened
        starting_string: string indicating if we are extracting "autosomes" or "sexsomes"
        options: options provided by the user to the program

    Output:
        lineconfig: config file with some lines read
        somes: list of chromosomes (and probabities)
    """
    lineconfig = configfile.readline() #auto

    while(lineconfig.find("#")==0): #It is a commented line, skip
        configfile.readline() 

    if(lineconfig.find("#")!=-1):   #There is a comment
        lineconfig=lineconfig[lineconfig.find("#"):len(lineconfig)] 
    
    lineconfig = lineconfig.strip()  

    if(lineconfig.startswith(starting_string)):
        somesstr= lineconfig[len(starting_string):len(lineconfig)] 
        somes = somesstr.split() 
    
    else:
        if starting_string=='sex:':
            sys.stderr.write("\nParse error with config file:"+options.configfile+" The second line needs to define the sex chromosomes as such:\nsex: 1 2 ...") 

        elif starting_string=='autosomes:':
            sys.stderr.write("\nParse error with config file:"+options.configfile+" The first line needs to define the autosomes as such:\nautosomes: 1 2 ...")

        elif starting_string =='unaffectedfemale:':
            sys.stderr.write("\nParse error with config file:"+options.configfile+" The third line needs to define the karyotype of an unaffectedfemale as such:\nunaffectedfemale: [chr#1] [chr#2] [prob]\nex:unaffectedfemale: X X 0.5\nFound "+str(len(somes))+" fields instead of 3") 

        elif starting_string =='unaffectedmale:':
            sys.stderr.write("\nParse error with config file:"+options.configfile+" The third line needs to define the karyotype of an unaffectedmale as such:\nunaffectedmale: [chr#1] [chr#2] [prob]\nex:unaffectedmale: X X 0.5\nFound "+str(len(somes))+" fields instead of 3") 

        sys.exit(1) 


    return lineconfig, somes


def AffectedRegions(conditions, autoranges, sexrangesUnaffectedMale, sexrangesUnaffectedFemale):
    """
    This function determines which chromosomal regions are affected in each condition

    Input:
        conditions: list of dictionarys containing name of the condition, type, prior probability and defChr.
            defChr is another list of dictionaris containing type of modification (deletion/duplication/other), the chromosome, and the starting and end coordinades of the chromosome.
        autoranges: list of dictionaries containing for the autosomes: the chromosome name, start and end coordinades and the expected coverage of the chromosome.
        sexrangesUnaffectedMale: as autorranges but for the sexual chromosomes in unaffected male.
        sexrangesUnaffectedFemale: as autorranges but for the sexual chromosomes in unaffected female.

    Output:
        list of regions where a modification can occur considering all the conditions.

    """

    index_rau = set()
    for cd in range(0,len(conditions)): #for each condition
        print(conditions[cd]["namecond"]+" new prior probability:\t"+str(conditions[cd]["priorprob"] ))
        
        #find which ranges are affected non-affected
        chrrangescpy=copy.deepcopy(autoranges)  #first copy autosomes
        #autoranges contains { "chr", "start","end", "ecov": }

        if( conditions[cd]["type"] == "auto" ): #autosomal
            if(conditions[cd]["namecond"].endswith("_m")):
                chrrangescpy.extend(copy.deepcopy(sexrangesUnaffectedMale))  #copy male chromosomes

            elif(conditions[cd]["namecond"].endswith("_f")):
                chrrangescpy.extend(copy.deepcopy(sexrangesUnaffectedFemale))  #copy male chromosomes

            else:
                sys.stderr.write("\nERROR it seems that for condition "+str(conditions[cd]["namecond"])+" was not defined as male/female internally, contact programmer\n") 
                sys.exit(1)   
    
        #then copy sex chromosomes (unaffected)
        conditions[cd]["rangesautolength"] = len(autoranges) 

        if( conditions[cd]["type"] == "male" ):
            chrrangescpy.extend(copy.deepcopy(sexrangesUnaffectedMale))  #copy male chromosomes

        if( conditions[cd]["type"] == "female" ):
            chrrangescpy.extend(copy.deepcopy(sexrangesUnaffectedFemale))  #copy female chromosomes


        
        numberOfAffectedRegions=0 
        totalRegions=0 

        rau = 0
        for r in chrrangescpy: #for each genomic range
            
            #r contains { "chr", "start","end", "ecov": }
            overlap = 0 # 0 no overlap, 1 deleted  2 duplicated
            totalRegions += 1 
            
            for c in conditions[cd]["defChr"]: #for chromosome region defined for the condition
                doesoverlap = 0 
                if( r["chr"]  ==  c["chrname"]): #same chr
                    
                    if(c["coordst"] == -1): #does not have a coordinate, entire chromosome
                        doesoverlap = 1 
                    else: #it has coordinates
                        if(not(
                            r["end"] < c["coordst"] or  #has coordinates and overlaps
                            r["start"] > c["coorden"]
                        )):
                            doesoverlap = 1 


                    if doesoverlap:
                        index_rau.add(rau)
                        numberOfAffectedRegions += 1 
                        
                        if overlap != 0:
                            sys.stderr.write("\nERROR it seems that for condition "+str(conditions[cd])+" there are 2 genomic regions that overlap\n") 
                            sys.exit(1)     
                        if( c["type"] == "deletion"): 
                            r["ecov"] = r["ecov"]-0.5 
                            overlap=1 # 0 no overlap, 1 deleted  2 duplicated
                        if( c["type"] ==  "duplication"): 
                            r["ecov"] = r["ecov"]+0.5 
                            overlap=2 # 0 no overlap, 1 deleted  2 duplicated
 
            rau += 1 
        rau = 0  
        conditions[cd]["chrrange"] = chrrangescpy 
        print  ("{} affects {}/{} ({}%) chromosome windows".format(conditions[cd]["namecond"], str(numberOfAffectedRegions), str(totalRegions), str(round(100*(float(numberOfAffectedRegions)/float(totalRegions)),2)) ))

    return sorted(index_rau)


def filter_file (file1, file2, conditions_index):
    """
    It writes file2 reading file 1 and filtering by the regions that can be affected by a condition and the mappability score
    Input:
        file1: file to read
        file2: file to write
        conditions_index: index of the regions that may be affected by any condition in the config file.

    Ouput:
        None
    """

    infile = open(file1, "r")
    outfile = open (file2, "w")
    counter = 0
    for line in infile:
        if counter not in conditions_index:
            outfile.write(line)
        counter+=1
    infile.close()
    outfile.close()



def merge_files(directory_name, map_score_list, index, corrected_file):
    """
    It merges the files containing the number of counts. 
    It corrects the counts using the mappability score.
    It saves in a new file the read counts from the regions that are not affected by any condition and are highly unique.

    Input:
        directory_name: directory where the files containing the read counts are located.
        map_score_list: list with the mappability scores for the chromosomal regions.
        index: index of the regions that may be affected by any condition indicated in the config file.
        corrected_file: name of the file to write

    Output:
        List of dictionaries containing: the index of the chromosome, the corrected read count, the original read count and the mappability score for that region.
    
    """


    directory=os.path.abspath(directory_name) 
    countInfo = []
    outfile = open (corrected_file, "w")

    #Add a "/" if it is missing at the end of the path
    if directory_name[-1]!="/":
        directory_name+="/"


    entries = os.listdir(directory)
    entries = [int(entry.split(".")[0]) for entry in entries if entry.split(".")[-1]=="count"]  #This is only numbers
    test = list(range(0, 3088))

    #REMOVE ALL THIS PART LATER
    not_common_entries = [entry for entry in test if entry not in entries]
    if len(not_common_entries)>0:
        sys.stderr.write("There are missing files")
        sys.exit(1)


    entries.sort()

    for entry in entries:

        #Open the file and read the info (count of mapped reads)
        count_file = open(directory_name+str(entry)+".count", "r")
        line = count_file.readline()
        count_file.close()

        if len(line) == 0:
            sys.stderr.write("The file number {} does not seem to have content. Try generating the count files in the first step again".format(entry))
            sys.exit(1)


        rau = entry
        line = checkLineJump(line)
        isInt(line)
        isFloat(map_score_list[rau])
        mapped_reads = int(line)
        map_score = float(map_score_list[rau])
        if map_score !=0:
            corrected_read = mapped_reads/map_score

        else:
            corrected_read = mapped_reads

        #Store all the observed reads
        countInfo.append({"id": rau, "corrected_read": corrected_read, "count":mapped_reads, "map_score": map_score})


        #Save in a new file only the non affected chromosomes
        #rau is not in index (the list of affected chromosomes)
        #if rau not in index and map_score>=0.9:  
        if rau not in index:
            outfile.write(str(entry)+"\t"+str(mapped_reads)+"\n")



    outfile.close()
    

    print("The files have been merged")

    #for entry in entries:
        #os.remove(directory_name+str(entry)+".count")

    print("The count files have been removed")
    
    return countInfo


def get_subsample(bam, final_coverage):
    #Not for final version

    infile = open("/home/projects/marthe/tmp/coverage_empirical.txt", "r")
    dict_coverage = dict()
    for line in infile:
        bam_name = line.split()[0]
        coverage = line.split()[1]
        coverage = checkLineJump(coverage)
        dict_coverage[bam_name] = float(coverage)

    try:
        cov0 = dict_coverage[bam]
        print("The original coverage was", cov0)
    except:
        sys.stderr.write("The average coverage value of the sample is not included in the file\nGenerating it")
        cov0 = mean_coverage(bam, )
    cov1 = float("0"+final_coverage)
    print("The expected coverage is", cov1)
    if cov1>cov0:
        sys.stderr.write("The coverage provided is smaller than what you are looking for\n")
        sum_frac_return = ".99"
    else:
        sum_frac = cov1/cov0

        sum_frac_return =str(sum_frac)[1:]
        print("The needed subsample fraction is", sum_frac_return)

    return sum_frac_return

    
#THIS FUNCTION IS NOT IN THE FINAL VERSION
def get_subsample2(bam, final_coverage, samtoolscmd):
    #Not for final version
    #This is used when the we have produced the read counts for the original coverage
    #Then, we take the final coverage we are interested in and give back the fraction which needs to be used to get the specified coverage

    infile = open("/home/projects/marthe/tmp/coverage_empirical.txt", "r")
    dict_coverage = dict()
    for line in infile:
        bam_name = line.split()[0]
        coverage = line.split()[1]
        coverage = checkLineJump(coverage)
        dict_coverage[bam_name] = float(coverage)
    infile.close()

    try:
        cov0 = dict_coverage[bam]
        
    except:
        sys.stderr.write("The average coverage value of the sample is not included in the file\nGenerating it")
        cov0 = mean_coverage(bam, samtoolscmd)
        print(cov0)
        outfile = open("/home/projects/marthe/tmp/coverage_empirical.txt", "r")
        outfile.write(bam + str(cov0) + "\n")
        outfile.close()
        sys.stderr.write("\nThe average coverage of the sample has been added to the file")
    
    print("The original coverage was", cov0)
    cov1 = float("0"+final_coverage)
    print("The expected coverage is", cov1)

    if cov1>cov0:
        sys.stderr.write("The coverage provided is smaller than what you are looking for\n")
        sum_frac_return = 1
    elif cov1<cov0:
        sum_frac_return = cov1/cov0

    print("The needed subsample fraction is", sum_frac_return)

    return sum_frac_return


def mean_coverage(bam_file, samtoolscmd):
    command = samtoolscmd + " coverage -o tmp_coverage_file " + bam_file 
    print("Launching command: " +  command)
    subprocess.run(command, shell=True, universal_newlines = True)  


    infile = open("tmp_coverage_file", "r")
    listchr = list(range(1,23))+["X", "Y", "x", "y"]
    listchr= [str(option) for option in listchr]

    coverage_l = []
    for line in infile:
        coverage_l.append(float(line.split()[6]))
    infile.close()
    os.remove("tmp_coverage_file")

    mean_cov = np.mean(coverage_l)
    print(mean_cov)

    return mean_cov
    


def subsample_file(subsample_value, directory_name, map_score_list, index, corrected_file):
    directory=os.path.abspath(directory_name) 
    countInfo = []
    outfile = open (corrected_file, "w")

    #Add a "/" if it is missing at the end of the path
    if directory_name[-1]!="/":
        directory_name+="/"


    entries = os.listdir(directory)
    entries = [int(entry.split(".")[0]) for entry in entries if entry.split(".")[-1]=="count"]  #This is only numbers
    test = list(range(0, 3088))

    #REMOVE ALL THIS PART LATER
    not_common_entries = [entry for entry in test if entry not in entries]
    if len(not_common_entries)>0:
        sys.stderr.write("There are missing files")
        sys.exit(1)



    entries.sort()

    for entry in entries:

        #Open the file and read the info (count of mapped reads)
        count_file = open(directory_name+str(entry)+".count", "r")
        line = count_file.readline()
        count_file.close()

        if len(line) == 0:
            sys.stderr.write("The file number {} does not seem to have content. Try generating the count files in the first step again".format(entry))
            sys.exit(1)


        rau = entry
        line = checkLineJump(line)
        isInt(line)
        isFloat(map_score_list[rau])
        mapped_reads = int(line)
        map_score = float(map_score_list[rau])
        subsampled_mapped_reads = mapped_reads*subsample_value
        if map_score !=0:
            corrected_read = subsampled_mapped_reads/map_score

        else:
            corrected_read = subsampled_mapped_reads

        #Store all the observed reads
        countInfo.append({"id": rau, "corrected_read": corrected_read, "count":subsampled_mapped_reads, "map_score": map_score})


        #Save in a new file only the non affected chromosomes
        #rau is not in index (the list of affected chromosomes)
        #if rau not in index and map_score>=0.9:  
        if rau not in index:
            outfile.write(str(entry)+"\t"+str(subsampled_mapped_reads)+"\n")



    outfile.close()
    

    print("The files have been merged")

    #for entry in entries:
        #os.remove(directory_name+str(entry)+".count")

    print("The count files have been removed")
    
    return countInfo


def subsample_file2(subsample_value, infile,  index, corrected_file):

    outfile = open(corrected_file, "w")
    rau = 0
    countInfo = []

    listchr = list(range(1,23))+["X", "Y", "x", "y"]
    listchr= [str(option) for option in listchr]

    with gzip.open(infile, "r") as infile:
        for line in infile:
            line = line.decode("utf-8")
            elements = line.split()

            if len(elements)<2:
                sys.stderr.write("\nThe format of the count file should be: 1:1-1000000     136783\nOne column is missing")
                sys.exit(1)
                
            if elements[0][0] not in listchr:
                sys.stderr.write("\nThe format of the count file should be: 1:1-1000000     136783\nIt does not start with a chromosome")
                sys.exit(1)

            count = float(line.split()[1])
            rau+=1
            subsampled_reads = count*subsample_value

            #Store all the observed reads
            countInfo.append({"id": rau, "count":subsampled_reads})

            if rau not in index:
                outfile.write(str(rau)+"\t"+str(subsampled_reads)+"\n")

    #infile.close()
    outfile.close()

    return countInfo


def calc_map_score(map_file, window=1E6):

    #Ns_file = "/home/databases/genomes/Homo_sapiens/hs37d5/hs37d5.countNs.gz"
    #map_file = "/home/databases/genomes/Homo_sapiens/hs37d5/mappability/snpable99/mappable_99.bed.gz"
    f=gzip.open(map_file,'rb')
    file_content=f.read()
    lines= file_content.split(b"\n")
    map_bases = 0
    non_map_score_l = []
    counter = 1
    chr = b"1"
    listchr = list(range(1,23))+["X", "Y", "x", "y"]
    listchr= [str(option).encode("utf-8") for option in listchr]


    for line in lines:
        if len(line.split()) == 3:
            if line.split()[0] in listchr:
                start = line.split()[1]
                end = line.split()[2]
                start = start.decode("utf-8")
                start =int(start)
                end = end.decode("utf-8")
                end = int(end)
                diff = end-start

                if line.split()[0] == chr: #I am in the same chromosome
                    
                    if (start<=counter*window) and (end <= counter*window):
                        map_bases += diff
                    else:
                        map_score = map_bases/window
                        map_bases = 0
                        non_map_score_l.append(map_score)
                        counter += 1
                else:
                    chr = line.split()[0]
                    counter = 1

    return non_map_score_l