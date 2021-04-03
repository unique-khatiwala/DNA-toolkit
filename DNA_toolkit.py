# bioinformatics DNA toolkit: youtuber rebelCoder
# collections library for counter
# f-string to display variable value for input variable in middle of string using f"string ... {variable}... string"
from collections import Counter
import random
class DNA_toolkit:

    nucleotidesDNA = ["A", "T", "C", "G"]
    nucleotidesRNA = ['A','U','C','G']
    compDictDNA = {"A":"T", "T":"A", "C":"G", "G":"C"}
    compDictRNA = {"A": "U", "U": "A", "C": "G", "G": "C"}
    codonDNA = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    codonRNA = { 'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
        'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',
    }


    def __init__(self,sequence='',type = 'DNA',length=20):
        self.sequence = sequence.upper() if sequence!='' else self.__generate_rnd(length)
        self.__seqVerifier()
        self.type = type.upper()
        self.codon = self.codonRNA if self.type =='RNA' else self.codonDNA
        self.compDict = self.compDictRNA if self.type == 'RNA' else self.codonDNA
        self.nucleotides = self.nucleotidesRNA if self.type =='RNA' else self.nucleotidesDNA

    def __generate_rnd(self,length = 20):
        sequence = [random.choice(self.nucleotides) for i in range(length)]


    def __seqVerifier(self):
        """ verifies the sequence for anomaly and standardizes its format """
        return set(self.nucleotides).issuperset(self.sequence)

    def ntCounter(self):
        """ counts the frequency of nucleotides"""
        counterDict = dict(Counter(self.sequence))
        return counterDict

    def revComplement(self,sequence=''):
        """ makes a reverse complement of DNA strand (double helix)
        another approach:
        mapped = sequence.translate(sequence.maketrans("ATGC", "TACG"))[::-1] """
        seq = list(self.sequence) if sequence == '' else list(sequence)
        complement="".join([self.compDict[bp] for bp in seq])
        return complement[::-1]

    def DNA2RNA(self):
        """ replaces T to U in base-pair sequence """
        if self.type == 'RNA':
            self.type  = 'DNA'
            return self.sequence.replace('U','T')

        elif self.type == 'DNA':
            self.type = 'RNA'
            return self.sequence.replace("T","U")


    def GCcontent(self):
        """ number of G/C base pairs out of total base pairs """
        content = dict(Counter(self.sequence))
        return round((content["G"]+content["C"])/len(self.sequence)*100,3)

    def protein_seq(self,sequence='',init = 0):
        """ returns a list of protein sequence """
        if sequence == '':
            sequence = self.sequence
        prot = [self.codon[sequence[codon:codon+3]] for codon in range(init,len(sequence)-2,3)]
        #prot = ("".join(prot)).replace('_','')
        return prot

    def aminoFreq(self, amino):
        """ Returns dictionary of codon sequences which match to single protein and its overall frequency"""
        tmp_list=[]
        for mark in range(0,len(self.sequence)-2,3):
            if self.codon[self.sequence[mark:mark+3]]==amino:
                tmp_list.append(self.sequence[mark:mark+3])

        freq = dict(Counter(tmp_list))
        freq = { key :round((value/(len(tmp_list)))*100,2) for (key,value) in freq.items()}
        return freq

    def frames(self,sequence=''):
        sequence = self.sequence if sequence =='' else sequence
        frame = []
        for i in range(0,3):
            frame.append(self.protein_seq(sequence,i))
            frame.append(self.protein_seq(self.revComplement(sequence),i))

        return frame

    def proteinLookup(self,frame):
        pt,flag = [],False
        for bp in frame:
            if (flag == True) and (bp == '_'):
                return "".join(pt)
            if (bp == 'M') or (flag == True):
                flag = True
                pt.append(bp)
        return "".join(pt) if pt!=[] else ''

    def all_frames_proteins(self, start = 0, stop = 0, ordered = False):
        frame = self.frames(self.sequence[start:stop] if stop>start else self.sequence)
        proteins = []

        for a_frame in frame:
            if self.proteinLookup(a_frame) != '':
                proteins.append(self.proteinLookup(a_frame))

        return sorted(proteins, key = len, reverse = ordered)


    def aminoCounter(self,sequence=''):
        ''' gives the amount of each amino acid is in each frame. And not the protein '''
        sequence = self.sequence if sequence == '' else sequence
        aa =[]
        for frame in self.frames(sequence):
            aa.append(dict(sorted(Counter(''.join(frame)).items())))
        return aa
    
    @staticmethod
    def hamming(seq1,seq2):
        z = zip(seq1,seq2) if len(seq1) == len(seq2) else []
        counter = 0
        differences = []
        for ele1,ele2 in z:
            counter+= 1 if ele1!=ele2 else 0
            differences.append((ele1,ele2)) if ele1!=ele2 else differences
        return counter

