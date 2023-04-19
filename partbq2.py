#!/bin/env python

## SLURM variables

#SBATCH --account=mpcs56430
#SBATCH --job-name=porf
#SBATCH --output=%j_porf.out
#SBATCH --partition=broadwl

#SBATCH --cpus-per-task=1   # cores
#SBATCH --nodes=6            # number of nodes to run on       
#SBATCH --ntasks-per-node=1  # 
#SBATCH --ntasks=6         # total tasks to be launcedd

#SBATCH--exclusive
#SBATCH --time=00:05:00

import multiprocessing
from multiprocessing import Pool
from functools import partial
import time
import sys
import os

nucleic_codes = list("ACGTURYKMSWBDHVN-")
amino_codes = list("ABCDEFGHIKLMNPQRSTUVWYZX*-")
digits = list("0123456789")

def test_fasta(filename, isnucleic = False):
    f = open(filename, 'r')
    line = f.readline().strip() #https://www.guru99.com/python-file-readline.html
    while line:
        if line[0] == '>':
            if line[1] == ' ':
                return False
            #adjust based on type
        else:
            if isnucleic:
                for char in line:
                    if char not in digits and char.upper() not in nucleic_codes:
                        print("Error at: " + line)
                        return False
            else:
                for char in line:
                    if char not in digits and char.upper() not in amino_codes:
                        print("Error at: " + line)
                        return False
        line = f.readline().strip()
    f.close()
    return True #all lines are correct

def parse_fasta(filename):
    if (not test_fasta(filename, False)):
        print("Invalid file")
    f = open(filename, 'r')
    line = f.readline().strip() 
    desc = ""
    seq = ""
    while line:
        if line[0] == '>':
            if desc != "":
                yield (desc, seq)
                seq = ""
            desc = line[1:]
        else:
            seq += line
        line = f.readline().strip()
    yield (desc, seq)
    f.close()

def dnatorna(dnastring):
    rna = ""
    for letter in dnastring:
        if letter in valid:
            if (letter == 'T'):
                rna += "U"
            else:
                rna += letter
    return rna

valid = ['T', 'A', 'C', 'G']

codon_table = { "UUU": "F", "UUC": "F", "UUA":"L", "UUG":"L", "CUU": "L", "CUC": "L", "CUA":"L", "CUG":"L",
              "AUU": "I", "AUC": "I", "AUA":"I", "AUG":"M", "GUU": "V", "GUC": "V", "GUA":"V", "GUG":"V",
              "UCU":"S", "UCC":"S", "UCA":"S","UCG":"S","CCU":"P", "CCC":"P", "CCA":"P","CCG":"P",
              "ACU":"U", "ACC":"U", "ACA":"U","ACG":"U","GCU":"A", "GCC":"A", "GCA":"A","GCG":"A",
              "UAU":"Y","UAC":"Y","UAA":"Z","UAG":"Z","CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
              "AAU":"N","AAC":"N", "AAA":"K", "AAG":"K","GAU":"D","GAC":"D","GAA":"E","GAG":"E",
              "UGU":"C","UGC":"C","UGA":"Z","UGG":"W","CGU":"R","CGC":"R","CGA":"R","CGG":"R",
              "AGU":"S","AGC":"S","AGA":"R","AGG":"R","GGU":"G","GGC":"G","GGA":"G","GGG":"G"}

def singlereadingframe(offset, rnastring):
    starttime = time.perf_counter()
    if (offset < 0):
        rnastring = rnastring[::-1]
    loc = abs(offset) - 1
    code = ""
    start = -1
    end = -1
    count = 0
    while (loc + 2 < len(rnastring) and count < 2):
        codon = rnastring[loc:loc + 3]
        output = codon_table[codon]
        if (output == "M" and count == 0):
            count += 1
            start = loc
        if (output == "Z" and count == 1):
            count += 1
            end = loc
        if (count == 1):
            code += output
        loc += 3
    if (count == 2):
        #print("* " + offset + " | " + start + " | " + end + " | " + len(code) + " | " + code + " \n")
        print(f"* {offset} | {start} | {end} | {len(code)} | {code} \n")
    endtime = time.perf_counter()
    print("The subprogram took " + str(endtime -starttime) + " seconds.")

def preprocess(file_path):
    out = list(parse_fasta(file_path))
    dnatext = out[0][1]
    rnatext = dnatorna(dnatext)
    return rnatext


if __name__ == "__main__":
    # necessary to add cwd to path when script run 
    # by slurm (since it executes a copy)
    sys.path.append(os.getcwd()) 

    #job_id = os.environ["SLURM_JOB_ID"]
    nodes = int(os.environ["SLURM_NNODES"])
    #cpus_per_node = int(os.environ["SLURM_JOB_CPUS_PER_NODE"])
    #tasks_per_node =  int(os.environ["SLURM_NTASKS_PER_NODE"])
    #cpus_per_task =  int(os.environ["SLURM_CPUS_PER_TASK"])
    #nprocs =  int(os.environ["SLURM_NPROCS"])

    print("nodes = %d" % int(os.environ["SLURM_NNODES"]))
    #print("cpus_per_node = %s" % os.environ["SLURM_JOB_CPUS_PER_NODE"])
    #print("tasks_per_node = %d" %  int(os.environ["SLURM_NTASKS_PER_NODE"]))
    #print("cpus_per_task =  %d" % int(os.environ["SLURM_CPUS_PER_TASK"]))
    #print("nprocs =  %d" % int(os.environ["SLURM_NPROCS"]))

    print(os.environ)

    SLURM_NPROCS = 6#nodes * cpus_per_task

    print("Number of Processes: %d" % SLURM_NPROCS)
    start = time.time()
    print("Start time: %s" % start)
    file_path = "/home/ericwu4/dump/sequence.fasta"
    rnatext = preprocess(file_path)

    func = partial(singlereadingframe, rnastring = rnatext)

    # Set the number of processes to the number available to us; 
    # this isn't necessarily the best performing option.  You should 
    # test your application to determine this
    pool = multiprocessing.Pool(processes=SLURM_NPROCS)

    # The array of parameters that will be passed to your process worker.  
    # This creates array ["00","01"....]
    database_chunks = [-3,-2,-1,1,2,3]

    # Distribute the jobs to workers
    pool.map(func, database_chunks)

    print("Run time: %f" % (time.time() - start))