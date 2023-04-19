from multiprocessing import Pool
from functools import partial
import time

nucleic_codes = list("ACGTURYKMSWBDHVN-")
amino_codes = list("ABCDEFGHIKLMNPQRSTUVWYZX*-")
digits = list("0123456789")
codon_table = { "UUU": "F", "UUC": "F", "UUA":"L", "UUG":"L", "CUU": "L", "CUC": "L", "CUA":"L", "CUG":"L",
              "AUU": "I", "AUC": "I", "AUA":"I", "AUG":"M", "GUU": "V", "GUC": "V", "GUA":"V", "GUG":"V",
              "UCU":"S", "UCC":"S", "UCA":"S","UCG":"S","CCU":"P", "CCC":"P", "CCA":"P","CCG":"P",
              "ACU":"U", "ACC":"U", "ACA":"U","ACG":"U","GCU":"A", "GCC":"A", "GCA":"A","GCG":"A",
              "UAU":"Y","UAC":"Y","UAA":"Z","UAG":"Z","CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
              "AAU":"N","AAC":"N", "AAA":"K", "AAG":"K","GAU":"D","GAC":"D","GAA":"E","GAG":"E",
              "UGU":"C","UGC":"C","UGA":"Z","UGG":"W","CGU":"R","CGC":"R","CGA":"R","CGG":"R",
              "AGU":"S","AGC":"S","AGA":"R","AGG":"R","GGU":"G","GGC":"G","GGA":"G","GGG":"G"}


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
                        print(f"Error at: {line}")
                        return False
            else:
                for char in line:
                    if char not in digits and char.upper() not in amino_codes:
                        print(f"Error at: {line}")
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

valid = ['T', 'A', 'C', 'G']
def dnatorna(dnastring):
    rna = ""
    for letter in dnastring:
        if letter in valid:
            if (letter == 'T'):
                rna += "U"
            else:
                rna += letter
    return rna

def halfrnatoprotein(halfrnastring, isflipped = False):
    #print(halfrnastring + "!!!!\n")
    validframes = 0
    for i in range(3):
        loc = i
        code = ""
        start = -1
        end = -1
        count = 0
        #unfiltered = ""
        while (loc + 2 < len(halfrnastring) and count < 2):
            codon = halfrnastring[loc:loc + 3]
            output = codon_table[codon]
            #print(f"{codon} == {output}\n")
            if (output == "M" and count == 0):
                count += 1
                start = loc
            if (output == "Z" and count == 1):
                count += 1
                end = loc
            if (count == 1):
                code += output
            loc += 3
            #unfiltered += output
        if (count == 2):
            j = i + 1
            if (isflipped):
                j *= -1
            #print(f"* {j} | {start} | {end} | {len(code)} | {code} \n")
            validframes += 1
    return validframes
            
def fullrnatoprotein(rnastring):
    frames = halfrnatoprotein(rnastring)
    #https://www.w3schools.com/python/python_howto_reverse_string.asp
    fliprnatext = rnatext[::-1]
    frames += halfrnatoprotein(fliprnatext, True)
    return frames

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
    #if (count == 2):
        #print(f"* {offset} | {start} | {end} | {len(code)} | {code} \n")
    endtime = time.perf_counter()
    #print(f"The subprogram took {endtime - starttime} seconds.")

if __name__ == "__main__":
    #file_path = "C:/Users/ewu15/jup/Scientific Computing/week 1/data/data/dogma.fasta.txt"
    #file_path = "C:/Users/ewu15/jup/Scientific Computing/week 1/data/data/sars-cov2-genome-fasta-filtered.txt"
    file_path = "C:/Users/ewu15/jup/Scientific Computing/week 2/sequence.fasta"

    start = time.perf_counter()
    out = list(parse_fasta(file_path))
    dnatext = out[0][1]
    rnatext = dnatorna(dnatext)
    fullrnatoprotein(rnatext)
    end = time.perf_counter()
    print(f"The program took {end - start} seconds.")

    start = time.perf_counter()
    out = list(parse_fasta(file_path))
    dnatext = out[0][1]
    rnatext = dnatorna(dnatext)
    func = partial(singlereadingframe, rnastring = rnatext)
    offsets = [-3,-2,-1,1,2,3]
    pool = Pool(processes = 6)
    pool.map(func, offsets)
    end = time.perf_counter()
    print(f"The program took {end - start} seconds.")