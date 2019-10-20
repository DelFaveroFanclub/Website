
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 09:20:17 2018

@author: Michiel
"""

from collections.abc import Sequence
import reprlib
import ast
from json import dump
import numpy as np

try:
    import urllib3
except ModuleNotFoundError:
    print('urllib3 not found, BED.get_sequences_from_uscs() not available')

#Prepare codon library
RNA_BASES = 'UCAG'
DNA_BASES = 'TCAG'
AMINO_ACIDS = 'FLSYCWPHQRIMTNKVADEG'

CODONS = [a+b+c for a in RNA_BASES
                for b in RNA_BASES
                for c in RNA_BASES]
AMINO_ACID_CODE = list('FFLLSSSSTT') + ['STOP']*2 + list('CC') + ['STOP'] + list('WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
CODE = {}
for code, AA in zip(CODONS,AMINO_ACID_CODE):
    CODE[code] = AA

IUPAC_CODES = {'U': 'U', 'C': 'C', 'A': 'A', 'G': 'G', 'R': 'AG', 'Y': 'UC', 'K': 'UG', 'M': 'CA', 'S': 'CG', 'W': 'UA', 'B': 'UCG', 'D': 'UAG', 'H': 'UCA', 'V': 'CAG', 'N': 'UCAG'}
IUPAC_CODES_DNA = {key.replace('U','T'): value.replace('U','T') for key, value in IUPAC_CODES.items()}
IUPAC_REVERSE = {'T': 'T', 'C': 'C', 'A': 'A', 'G': 'G', 'AG': 'R', 'TC': 'Y', 'TG': 'K', 'CA': 'M', 'CG':'S', 'TA': 'W', 'TCG': 'B', 'TAG': 'D', 'TCA': 'H', 'CAG': 'V', 'TCAG': 'N'}

#Class voor DNA objects
class Gene(Sequence):
    BASES = DNA_BASES

    def __init__(self,seq='',file='',region=''):
        self.region = region
        if file:
            with open(file,'r') as f:
                first_line = f.readline()
                while first_line[0] in '>;':
                    first_line = f.readline()
                    continue
                self.seq = (first_line.strip() + f.read().replace('\n','')).upper()
        else:
            self.seq = seq.upper()

#Magic methods to use genes with generic python functions
    def __getitem__(self,index):
        return Gene(seq=self.seq[index])

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        n = type(self).__name__
        if self.region:
            s = self.region
        else:
            s = reprlib.repr(self.seq)
        return '{}({})'.format(n,s)

    def __repr__(self):
        return self.seq

    def __add__(self, other):
        if type(other) != type(self):
            raise TypeError('Operation not permitted between %s and %s' %(type(self), type(other)))
        return Gene(seq = self.seq + other.seq)

    def __contains__(self, item):
        if isinstance(item, (Gene, RNA, Protein, Motif)):
            return item.seq in self.seq
        else:
            return item in self.seq

    def toRNA(self):
        return RNA(self.seq.replace('T','U'))

    def complement(self):
        return Gene(self.seq.replace('ATGC','TACG')[::-1])

#Genius implementation to find startcodons
    def startcodon(self):
        return self.seq.find('A' + self.BASES[0] + 'G')

#Replace IUPAC codes by {bases}
    def expand(self):
        newseq = self.seq
        for letter in self.seq:
            if letter not in 'NRYWSKMBDHV':
                continue
            else:
                if isinstance(self,RNA):
                    ex = IUPAC_CODES[letter]
                else:
                    ex = IUPAC_CODES[letter].replace('U','T')
                fmt = '{' + '{}'.format(ex) + '}'
                newseq = newseq.replace(letter,fmt)

        return type(self)(seq=newseq)

#Classic translation
    def toProtein(self):
        mRNA = self.toRNA()
        coding_mRNA = mRNA[mRNA.startcodon():]
        number_of_codons = len(coding_mRNA) // 3
        aminos = []
        for i in range(number_of_codons):
            codon = coding_mRNA[3*i:3*i+3]
            if set(codon) & set('NRYWSKMBDHV'):
                a, b, c = (IUPAC_CODES[codon[i]] for i in range(3))
                pos_codons = [x + y + z for x in a for y in b for z in c]
                pos_AA = set(CODE[code] for code in pos_codons)
                new_amino = '{' + ''.join(pos_AA) + '}'
            else:
                new_amino = CODE[codon]
            if new_amino == 'STOP':
                break
            aminos.append(new_amino)
        return Protein(''.join(aminos))

    def toDNA(self):
        return self

#Create PWM from gene
    def PWM(self):
        pwm = np.zeros((len(self),4))
        for i,base in enumerate(self.seq):
            pos_bases = IUPAC_CODES_DNA[base]
            for b in pos_bases:
                pwm[i,DNA_BASES.index(b)] = 1/len(pos_bases)
        return pwm

#Class for RNA's, currently not used
class RNA(Gene):
    BASES = RNA_BASES

    def __init__(self,seq='',file=''):
        if file:
            with open(file,'r') as f:
                first_line = f.readline()
                while first_line[0] in '>;':
                    first_line = f.readline()
                    continue
                self.seq = (first_line.strip() + f.read().replace('\n','')).upper()
        else:
            self.seq = seq

    def toDNA(self):
        return Gene(self.seq.replace('U','T'))


#Class for Proteins, currently not used
class Protein(Sequence):
    AA = AMINO_ACIDS

    def __init__(self,seq='',file=''):
        if file:
            with open(file,'r') as f:
                first_line = f.readline()
                while first_line[0] in '>;':
                    first_line = f.readline()
                    continue
                self.seq = (first_line.strip() + f.read().replace('\n','')).upper()
        else:
            self.seq = seq

    def __getitem__(self,index):
        return self.seq[index]

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        s = reprlib.repr(self.seq)
        n = type(self).__name__
        return '{}({})'.format(n,s)

    def __repr__(self):
        return self.seq

    def toProtein(self):
        return self

    def toRNA(self):
        clist = []
        for aa in self.seq:
            clist.append(list(CODE.keys())[list(CODE.values()).index(aa)])
        return RNA(''.join(clist))

    def toDNA(self):
        return self.toRNA().toDNA()


#Class for motifs (short DNA stretches)
class Motif(Gene):

    def __init__(self,seq='', pwm = []):
        seq = seq.upper()
        pwm = np.array(pwm)

        if pwm.any():
            if pwm.ndim != 2 or np.shape(pwm)[-1] != 4:
                raise TypeError('PWM must be nx4 matrix')

            seq = ''
            for base in pwm:
                pos_bases = [DNA_BASES[i] for i in range(4) if base[i]]
                seq += IUPAC_REVERSE[''.join(pos_bases)]


        if seq and not pwm.any():
            pwm = np.zeros((len(seq),4))
            for i,base in enumerate(seq):
                pos_bases = IUPAC_CODES_DNA[base]
                for b in pos_bases:
                    pwm[i,DNA_BASES.index(b)] = 1/len(pos_bases)


        self.pwm = pwm
        self.seq = seq

#Generate all motifs corresponding to IUPAC motif or PWM
    def all_motifs(self):
        all_motifs = [self]
        motif_counter = 1
        unkown_bases = [x for x in self.seq if not x in 'TCAG']
        for base in unkown_bases:
            pos_bases = IUPAC_CODES_DNA[base]
            cur_motifs = all_motifs[-motif_counter:]
            for c in range(motif_counter):
                cur_motif = cur_motifs[c]
                for b in pos_bases:
                    ind = cur_motif.seq.index(base)
                    new_motif = cur_motif.seq[:ind] + b + cur_motif.seq[ind+1:]
                    all_motifs.append(Motif(seq=new_motif))
            motif_counter *= len(pos_bases)
        unique_motifs = all_motifs[-motif_counter:]

        return unique_motifs

#PWM-based overlap calculator
    def overlap(self,gene):
        pos_motif_list = [gene[i:i+len(self)] for i in range(len(gene)-len(self))]

        scorelist = []
        for gene in pos_motif_list:
            gene_pwm = gene.PWM()
            score = np.sum(gene_pwm * self.pwm)/len(self)
            scorelist.append(score)
        return max(scorelist)

#Is Motif present in gene?
    def is_present_in(self, gene, min_overlap = 1):
        result = False
        all_motifs = self.all_motifs()

        #Fast implementation for exact matches
        if min_overlap == 1:
            for motif in all_motifs:
                if motif in gene:
                    result = True
                    break

        #PWM-based overlap calculation
        else:
            if min_overlap <= self.overlap(gene):
                result = True

        return result

    def toDNA(self):
        return Gene(seq=self.seq)


#Class for reading and handling bed files
class BED():

    def __init__(self, bed_file=''):
        self._BEDFILE  = bed_file
        self.pos_list = []
        self.gene_list = []
        if self._BEDFILE:
            with open(self._BEDFILE, 'r') as f:
                header = True
                eof = False
                while header:
                    line = f.readline()
                    if line.startswith('chr'): header = False
                while not eof:
                    ch, start, stop = line.split()[:3]
                    self.pos_list.append('%s:%s-%s' % (ch,start,stop))
                    line = f.readline()
                    if not line: eof = True
            self.gene_list = self.pos_list

    def parse(self):
        return self.pos_list

    def __repr__(self):
        n = type(self).__name__
        l = repr(self.gene_list)
        return '<{}: {}>'.format(n,l)

    def __len__(self):
        return len(self.gene_list)

    def __getitem__(self,item):
        return self.gene_list[item]

    def __str__(self):
        return self._BEDFILE
    
    def get_sequences(self):
        pass

    def get_sequences_from_ucsc(self,ref_genome='hg19'):
        try:
            http = urllib3.connection_from_url('https://api.genome.ucsc.edu')
        except NameError:
            print('urllib3 not found, exit 1')
            return 1
        for i, pos in enumerate(self.pos_list):
            ch, dom = pos.split(':')
            start,stop = dom.split('-')
            url = '/getData/sequence?genome={};chrom={};start={};end={}'.format(ref_genome, ch,start,stop)
            try:
                response = http.request('GET',url)
            except urllib3.exceptions.MaxRetryError:
                print('Sequence could not be retrieved from ucsc (MaxRetryError), continuing')
                continue
            s = response.data.decode().strip()
            DNA = ast.literal_eval(s)['dna']
            self.gene_list[i] = Gene(seq=DNA,region=pos)
        return self.gene_list

    def save(self,filename='bed.json'):
        savedict = dict()
        for gene in self.gene_list:
            savedict[gene.region] = gene.seq
        with open(filename,'w') as f:
            dump(savedict, f)
            f.writelines('\n')
