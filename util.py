'''
Created on Dec 28, 2014
utilities for genome reconstruction
@author: Jim
'''
# used code from https://gist.github.com/philippbayer/1626642
import time, cProfile, pstats
import cPickle as pickle
from collections import OrderedDict

class Genome:
    def __init__(self, name, version, snps):
        self.name = name   # who the person is
        self.version = version # 36, 37
        self.chroms = []  # list of Chromo objects
        while len(snps) > 0:
            this_chrome = snps[0].get_chromosome()
            temp_snp_list = []
            while this_chrome == snps[0].get_chromosome():
                temp_snp_list.append(snps[0])
                snps.pop(0)
                if len(snps) == 0:
                    break
            self.chroms.append(Chromo(this_chrome, temp_snp_list))
        self.errors = 0  # add to errors any time something impossible happens in the snp calulation, maybe should be an error log, list of messages
    def __str__(self):
        return 'genome for ' + self.name + ' version ' + str(self.version) + ' number chromosomes ' + str(len(self.chroms))
    def get_chromosomes(self):
        return self.chroms  # return the list of all of the chromosomes
    def get_chromosome(self, name):
        for chrome in self.chroms:
            if chrome.get_name() == name:
                return chrome
        return None # if name isn't found in the chromosome list, return None
 
             
class Chromo:  
    '''
    
    '''  
    def __init__(self, name, snps=[]):
        '''
        name is a string, '1', '2', ..., 'Y', 'MT'
        snps is a list of Snp objects
        '''
        self.name = name
        self.snps = OrderedDict()  # store the snps in an OrderedDict to preserve the order of snps
        self.errors = 0  # add to errors any time something impossible happens in the snp calulation, maybe should be an error log, list of messages
        for snp in snps:
            self.snps[snp.get_name()] = snp 
    def __str__(self):
        return 'chromosome ' + self.name + ' contains ' + str(len(self.snps.keys())) + ' snps'   
    def get_name(self):
        return self.name   
    def add(self, snp):
        self.snps[snp.get_name()] = snp # add the snp, overwrites if it is already there   
    def get_snp_names(self):
        return self.snps.keys()
    def get_snp(self, name):
        try:
            return self.snps[name] 
        except KeyError:
            return False
    def percent_covered(self):
        '''
        for all of the snps in Genome, calculate percent coverage. 
        '''
        covered = 0
        for key in self.snps.keys():
            p_covered, x = self.snps[key].percent_covered()
            covered += p_covered
        covered = covered / len(self.snps.keys())
        return covered
    def update_geno(self, new_letter, name):
        '''
        updates snp with new letter. returns with an error if letter not ATGC or name not one of the snps
        '''
        if new_letter not in ['A', 'C', 'T', 'G']:
            return 'letter not a valid letter'
        if name not in self.snps.keys():
            return str(name) + ' not in genome'
        self.snps[name].update_geno(new_letter)  
    def update_snp(self, snp):
        try:
            current_snp = self.snps[snp.get_name()]
            new_letter = snp.get_uservariation().strip('-')
            current_snp.update_geno(new_letter)
        except KeyError:
             self.snps[snp.get_name()] = snp   
                  
class Snp:
    # Each SNP is one line in the original file
    def __init__(self, name, variation, chromosome, position, strand, uservariation, genotemp={}):
        self.name = name # e.g. RS1234
 
    # only deCODEme lists possible variations and strand
        if variation is None:
            self.variation = "Unknown"
        else:
            self.variation = variation # e.g. A/T
 
        self.chromosome = str(chromosome) # 1,2,3...
        self.position = int(position) # e.g. 213412
 
        if strand is None:
            self.strand = "Unknown"
        else:
            self.strand = strand # +,-
 
        self.uservariation = str(uservariation) # e.g. TT
        
        # genotemp holds a dictionary of what the snp might be, {'A':3, 'C':1}
        self.genotemp = genotemp
    def get_genotemp(self):
        return self.genotemp
    def get_name(self):
        return self.name
    def get_uservariation(self):
        return self.uservariation
    def get_chromosome(self):
        return self.chromosome
    def get_position(self):
        return self.position
    def get_strand(self):
        return self.strand
    def get_variation(self):
        return self.variation
    def set_uservariation(self, var):
        self.uservariation = var
    def update_geno(self, new_letter):
        if new_letter not in ['A', 'C', 'T', 'G']:
            return 'letter not a valid letter'
        if new_letter in self.genotemp.keys():
            self.genotemp[new_letter] += 1
        else:
            self.genotemp[new_letter] = 1
    def percent_covered(self):
        '''
        for this snp, how much of it is covered. 
        '''
        if len(self.genotemp.keys()) == 0:
            return 0.0, '--'
        if len(self.genotemp.keys()) == 1:
            if self.chromosome == 'Y':
                return 1.0, self.genotemp.keys()[0] + '-'
            if self.genotemp[self.genotemp.keys()[0]] < 2:
                return 0.5, self.genotemp.keys()[0] + '-'
            else:
                return 1.0, self.genotemp.keys()[0] + self.genotemp.keys()[0]
        if len(self.genotemp.keys()) == 2:
            return 1.0, self.genotemp.keys()[0] + self.genotemp.keys()[1] 
        pair = ''  
        total_tests = sum(self.genotemp.values()) # total number of trials
        for letter in self.genotemp.keys():
            if self.genotemp[letter]/total_tests >.4:
                pair += self.genotemp[letter]
                #pair.append(self.genotemp[letter])
        if len(pair) == 1 and self.genotemp[pair[0]]/total_tests > .8 and total_tests > 1:
            pair += pair  # if there is only one letter and it is the vast majority, make it homozygous for that letter
        covered = len(pair)/2.0
        print(pair)
        while len(pair) < 2:
            #pair.append('-')
            pair += '-'
        return covered, pair

class Phased_snp: 
    def __init__(self, name, mom, children, m_pattern='', d_pattern=''): 
        '''
        name, rsid string; mom, list of strings ['AG', ..] for each kid; f_dad the letter gotten from dad
        d_pattern, a pattern of inheritance to infer phasing, '0110'; f_mom the letter gotten from mom; m_pattern, string of mom's pattern
        '''
        self.name = name  # rs12...
        self.mom = mom    # 'AG'
        self.children = children # ['AG', 'AA', ...]
        self.f_dad = ''
        self.d_pattern = d_pattern
        self.f_mom = ''
        self.m_pattern = m_pattern
        self.m_dist = 0  # if either mom or dad's pattern changes, record the distance here
        self.d_dist = 0
        self.d_crossover = 0 # the child that the crossover occurred
        self.m_crossover = 0
        self._first_phase(m_pattern, d_pattern)

    def __str__(self):
        return self.name.ljust(14) + self.mom.ljust(5) + '.'.join(self.children).ljust(15) + self.f_dad.ljust(9) + \
                      self.d_pattern.ljust(8) + self.f_mom.ljust(9) + self.m_pattern.ljust(7) 
    def get_d_pattern(self):
        return self.d_pattern
    def get_m_pattern(self):
        return self.m_pattern
    def get_d_dist(self):
        return self.d_dist
    def get_m_dist(self):
        return self.m_dist
    def get_m_crossover(self):
        return self.m_crossover
    def get_d_crossover(self):
        return self.d_crossover
    def get_f_dad(self):
        return self.f_dad
    def get_f_mom(self):
        return self.f_mom
    def get_name(self):
        return self.name
    def phase(self, m_pattern, d_pattern):
        dad_snps = list(self.f_dad)
        mom_snps = list(self.f_mom)
        for times in range(2):
            if '-' not in dad_snps and '-' not in mom_snps: break  # may need to loop zero or one times
            d_zipped = zip(dad_snps, d_pattern)
            for index, x in enumerate(d_zipped): # make list like [('-', '0'), ('T', '1'), (), ()]
                if x[0] == '-':  # if the char is '-', look for a matching second char, '0' for example, that has a valid letter
                    for z in d_zipped:
                        if z[0] != '-' and z[1] == x[1]:
                            dad_snps[index] = z[0]
                            mom_snps[index] = self._other_letter(z[0], self.children[index])
            m_zipped = zip(mom_snps, m_pattern)
            for index, x in enumerate(m_zipped):
                if x[0] == '-':
                    for z in m_zipped:
                        if z[0] != '-' and z[1] == x[1]:
                            mom_snps[index] = z[0]
                            dad_snps[index] = self._other_letter(z[0], self.children[index])
        # now letters have been filled in. Maybe new patterns. Tricky to know what to do.
        self.f_mom = ''.join(mom_snps)
        self.f_dad = ''.join(dad_snps)
        self.m_pattern, self.m_dist, self.m_crossover = self._informative_pattern(self.f_mom, m_pattern)
        self.d_pattern, self.d_dist, self.d_crossover = self._informative_pattern(self.f_dad, d_pattern)
        
    def _first_phase(self, m_pattern='', d_pattern=''):
        self.f_mom, self.f_dad = self._pattern(self.mom, self.children)  # gets patterns like from mom '-TT-' and from dad 'CCCC'
        if d_pattern == '':
            self.d_pattern, self.d_dist, self.d_crossover = self._informative_pattern(self.f_dad)
        else:
            self.d_pattern, self.d_dist, self.d_crossover = self._informative_pattern(self.f_dad, d_pattern)
        if m_pattern == '':
            self.m_pattern, self.m_dist, self.m_crossover = self._informative_pattern(self.f_mom)
        else:
            self.m_pattern, self.m_dist, self.m_crossover = self._informative_pattern(self.f_mom, self.m_pattern)  
    def _other_letter(self, letter, snp):
        '''
        letter is A, T, G, or C; snp is AG, GC, GG, ...
        return the other letter, could actually be the same letter
        '''
        if letter == snp[0]:
            return snp[1]
        else:
            return snp[0]
    def _pattern(self, mom_snp, sib_list):
        '''
        mom_snp is a two letter string like 'AC'
        sib_list is a list of two letter strings like ['AC', ...]
        return two values, mom_pattern, dad_pattern
        the patterns are either None, or a string with the same number of letters as there are sibs
        '''
        mom_pattern = ['-'] * len(sib_list) # both mom and dad start with as many dashes as are siblings
        dad_pattern = mom_pattern[:]
        for index, sib in enumerate(sib_list):
            if sib[0] == sib[1]:
                dad_pattern[index] = sib[0] # homozygous, got the same letter from both parents
                if sib[0] not in mom_snp:
                    continue
                mom_pattern[index] = sib[0]
                continue
            if mom_snp[0] == mom_snp[1]:
                mom_pattern[index] = mom_snp[0] # if mom is homozygous, then the kid got that letter from mom
                if sib[0] != mom_snp[0] and sib[0] != '-':
                     # mom CC, sib CA, so C came from mom, A from Dad
                    dad_pattern[index] = sib[0]
                if sib[1] not in mom_snp and sib[1] != '-':
                    dad_pattern[index] = sib[1]
                continue
            if mom_snp == sib: # mom and sib are both AC for instance, can't tell anything
                continue
            if sib[0] in mom_snp and sib[1] not in mom_snp:
                mom_pattern[index] = sib[0]
                dad_pattern[index] = sib[1]
                continue
            if sib[1] in mom_snp and sib[0] not in mom_snp:
                mom_pattern[index] = sib[1]
                dad_pattern[index] = sib[0]
                continue
            #print('error, one of the above cases should have happened, mom_snp, sib', mom_snp, sib)
        return ''.join(mom_pattern), ''.join(dad_pattern)

    def _informative_pattern(self, pattern, c_pattern=''):
        '''
        take a pattern, and return empty string if not informative, or full string if informative
        c_pattern is the current pattern, like '0111' or '1100' or '' if there hasn't been a pattern yet
        informative are like 'AACC' or 'TCCC'
        there should only be two letters, but I don't test
        return '' if no informative pattern, or '0111' or some such informative pattern plus a Hamming distance and child
        if the Hamming distance is 0 or 1, all is well. If greater, then there's a 50/50 chance the resulting pattern is wrong
        child is the first position that has a change and is a number from 0 to n-1
        '''
        not_informative = ''
        if '-' in pattern:
            return not_informative, 0, 0
        if c_pattern == '':
            first_digit = '0'
            second_digit = '1'
        else:
            first_digit = c_pattern[0]
            if first_digit == '0':
                second_digit = '1'
            else:
                second_digit = '0'
        new_pattern = first_digit
        for letter in pattern[1:]:
            if letter == pattern[0]:
                new_pattern += first_digit
            else:
                new_pattern += second_digit
        if sum([int(i) for i in new_pattern]) == 0 or sum([int(i) for i in new_pattern]) == len(pattern):
            return not_informative, 0, 0
        if c_pattern == '':
            return new_pattern, 0, 0
        dist, child = hamming_distance(c_pattern, new_pattern)
        if dist <= len(pattern)/2:
            return new_pattern, dist, child
        #print('new_pattern, dist',  new_pattern, dist)
        final_pattern = ''
        for letter in new_pattern:
            if letter == '0':
                final_pattern += '1'
            else:
                final_pattern += '0'
        return final_pattern, hamming_distance(c_pattern, final_pattern)[0], hamming_distance(c_pattern, final_pattern)[1]
def test_phase():
    snp2 = Phased_snp('rs2', 'AC', ['AC', 'AA', 'AA', 'CC'])
    
def test_informative():
    snp2 = Phased_snp('rs2', 'AC', ['AC', 'AA', 'AA', 'CC'])
    assert snp2._informative_pattern('A-CC') == ('', 0, 0)
    assert snp2._informative_pattern('CCCC') == ('', 0, 0)
    assert snp2._informative_pattern('AACC') == ('0011', 0, 0)
    assert snp2._informative_pattern('AACC', '0011') == ('0011', 0, 0)
    assert snp2._informative_pattern('AACC', '0010') == ('0011', 1, 3)
    assert snp2._informative_pattern('CACC', '0011') == ('1011', 1, 0)
    assert snp2._informative_pattern('ACCC', '1011') == ('1000', 2, 2)
    assert snp2._informative_pattern('CCAC', '1110') == ('1101', 2, 2)

def test_Phased_snp():
    print('snp_name,   mom,     children, d_snps,   d_pat,  m_snps,    m_pat')
    snp2 = Phased_snp('rs2', 'AA', ['AA', 'AA', 'AG', '-G'])
    assert snp2.get_f_dad() == 'AAGG' and snp2.get_f_mom() == 'AAAA'
    #print(snp2)
    snp3 = Phased_snp('rs3', 'CT', ['CT', 'TT', 'TT', 'CT'], '', '0011')
    #
    assert snp3.get_f_dad() == '-TT-' and snp3.get_f_mom() == '-TT-'
    snp3.phase('', '0011')
    print(snp3)
    #assert snp3.get_m_pattern() == '0110'
    snp16 = Phased_snp('rs16', 'CC', ['CT', 'CC', 'CC', 'CC'])
    #print(snp16)
    #assert snp16.get_d_pattern() == '0111'
    snp19 = Phased_snp('rs19', 'GT', ['GT', 'GT', 'GT', 'GG'])
    #print(snp19)
    snp19.phase('0110', '0111') 
    #print(snp19) 
    #assert snp19.get_d_pattern() == '0111' and snp19.get_m_pattern() == '0110'
    
    snpx = Phased_snp('rsx', 'GG', ['GT', 'GT', 'GG', 'GG'], '', '1000')
    #print(snpx)
    #print(snpx.get_d_dist(), snpx.get_m_dist(), snpx.get_d_crossover(), snpx.get_m_crossover())
    #assert snpx.get_d_dist() == 1 and snpx.get_d_crossover() == 1
    snprs914994 = Phased_snp('rs914994', 'GT', ['GG', 'GT', 'GT', 'GT'], '0001', '1001')
    #print(snprs914994)
    snprs914994.phase('0001', '1001')
    #assert snprs914994.get_m_pattern() == '0001' and snprs914994.get_d_pattern() == '1001'
    
class Phase_chrom:
    '''
    phased chromosome for a family, a mom, children, holds the estimate for the father
    '''
    def __init__(self, name, phased_snps = [], m_start_pattern = '', d_start_pattern = ''):
        self.name = name  # '1', 'MT', 'X', ...
        self.phased_snps = phased_snps # a list of phased snps
        self.m_start_pattern =  m_start_pattern # the chromosome phasing pattern at the start of the chromosome, '0011' for mom
        self.d_start_pattern = d_start_pattern
        self.crossovers = [] # list of crossover objects
    def __str__(self):
        return 'name ' + self.name + ', number snps ' + str(len(self.phased_snps)) + ', number of crossovers ' + str(len(self.crossovers)) + \
            ', mom start pattern ' + self.m_start_pattern + ', dad start pattern ' + self.d_start_pattern
    def set_m_start_pattern(self, m_start_pattern):
        self.m_start_pattern = m_start_pattern
    def get_m_start_pattern(self):
        return self.m_start_pattern
    def set_d_start_pattern(self, d_start_pattern):
        self.d_start_pattern = d_start_pattern
    def get_d_start_pattern(self):
        return self.d_start_pattern
    def get_name(self):
        return self.name
    def add_snp(self, snp):
        self.phased_snps.append(snp)
    def get_snps(self):
        return self.phased_snps
    def add_crossovers(self, crossover):
        self.crossovers.append(crossover)
    def get_crossovers(self):
        return self.crossovers
              
class Crossover:
    '''
    when a crossover is found store relavant information here
    '''
    def __init__(self, snp_name, h_dist, child, parent):
        self.snp_name = snp_name # the name of the snp that the crossover was detected. the actual crossover will be earlier
        self.h_dist = h_dist  # hamming distance normally 1. if 2 or greater there will be a problem with phasing
        self.child = child    # 0 or 1 or ...
        self.parent = parent  # 0 or 1, mom or dad
    def __str__(self):
        return 'snp_name ' + self.snp_name + ', hamming distance ' + self.h_dist + ', child ' + self.child + ', parent ' + self.parent
def test_pattern():
    p = Phased_snp('name', 'CA', ['CC', 'AC', 'CT', 'TT'])
    assert p._pattern('CA', ['CC', 'AC', 'CT', 'TT']) == ('C-C-', 'C-TT')
    #assert Phased_snp._pattern('CA', ['CC', 'AC', 'CT', 'TT']) == ('C-C-', 'C-TT')
    #assert Phased_snp._pattern('CC', ['CC', 'AC', 'CT', 'TT']) == ('CCC-', 'CATT')   

def parseFTDna(filein):
    allSnps = []
    for line in filein:
        # skip first line
        if "RSID,CHROMOSOME,POSITION,RESULT" in line:
            continue
        lineArray = line.replace('"','').replace("\n","").split(",")
        allSnps.append(Snp(lineArray[0], None, lineArray[1], lineArray[2], None, lineArray[3]))
    return allSnps
def test_parseFTDna():
    #filein = open('c:\\documents\\Family Tree Maker\\DNA\\ftdna-37-JamesSmithRawData.csv')
    filein = open('c:\\Users\\Jim\\documents\\DNA\\ftdna-37-JamesSmithTestData3.csv', 'r')
    b_time = time.time()
    jim_genome = Genome('james_smith', 37, parseFTDna(filein))
    print('test_parseFTDna', (time.time() - b_time)/60, len(jim_genome.get_snp_names()))
def parse23andme(filein):
    allSnps = []
    for line in filein:
        if "#" in line: # skip comments
            continue
        lineArray = line.replace(" ", "").replace("\n","").replace("\r", "").split("\t")
        try:
            allSnps.append(Snp(lineArray[0], None, lineArray[1], lineArray[2], None, lineArray[3])) 
        except IndexError:
            print('lineArray exception in parse23andme', lineArray) 
    return allSnps
def test_parse23andme():
    #filein = open('c:\\Users\\Jim\\documents\\DNA\\23andme-37-CynthiaTest.txt')
    #filein = open('c:\\Users\\Jim\\documents\\DNA\\23andmeRay_Smith-37-Test2.txt')
    filein = open('/Users/jim/Documents/DNA/23andmeRay_Smith-37-Test2.txt')
    #filein = open('c:\\Users\\Jim\\documents\\DNA\\ftdna-37-JamesSmithTestData3.csv')
    b_time = time.time()
    #cindy_genome = Genome('cindy_vrazsity', 37, parse23andme(filein))
    ray_genome = Genome('ray_smith', 37, parse23andme(filein))
    #jim_genome = Genome('jim_smith', 37, parse23andme(filein))
    #print(cindy_genome)
    #print('test_parse23andme', (time.time() - b_time)/60, len(ray_genome.get_snp_names()))
    #print('len', len(jim_genome.get_snp_names()))
    #print(ray_genome.percent_covered())
    for chrom in ray_genome.get_chromosomes():
        for snp_name in chrom.get_snp_names():
            print('snp: ', chrom.get_snp(snp_name).get_uservariation())

def parse_files(f):
    '''
    test the file, f (a string), and call either parse23andme or parseftdna
    '''
    my_file = open(f, 'r')
    if '23andme' in f:
        return parse23andme(my_file)
    if 'ftdna' in f:
        return parseFTDna(my_file)
    return 'not a supported file type'
def test_parse_files():
    print(parse_files('/Users/jim/Documents/DNA/23andme-37-CarolynSmithTest.txt'))
def compute_father(mother, children = [], cousins = [], mitochondrial = [], ydna = [], father_genome = Genome('father', 37, [])):
    '''
    inputs; mother is a file name. There must be either sons or daughters file names in a list
    optionally stepchildren and cousins who could actually be niece, nephew, cousin, aunt, uncle
    output a genome file for a father with as much filled in as possible.
    Y chromosome would come from the sons if possible. I could add a separate Y chromosome entry if there were another person who shared the Y chromosome
    
    '''
    
    mother_genome = Genome('mother', 37, parse_files(mother))
    for child in children:  # I can't think of why to separate sons and daughters yet
        genome = Genome('child', 37, parse_files(child))
        for snp_name in genome.get_snp_names():
            if genome.get_snp(snp_name).get_chromosome() in ['X', 'MT']: # X & MT useless for creating father
                continue 
            if genome.get_snp(snp_name).get_chromosome() == 'Y':
                snp = genome.get_snp(snp_name).get_uservariation().strip('-')
                new_snp = Snp(genome.get_snp(snp_name).get_name(), genome.get_snp(snp_name).get_variation(), 'Y', genome.get_snp(snp_name).get_position(),
                               genome.get_snp(snp_name).get_strand(), snp, genotemp={snp:1})
            else:
                # make comparison of mother's snp to child's snp, to produce a letter in ['A', 'T', 'C', 'G', '-']
                new_snp = compare_parent_child(mother_genome.get_snp(snp_name), genome.get_snp(snp_name))
                father_genome.update_snp(new_snp)
#            if new_snp.get_name() in father_genome.get_snp_names(): ###
#                new_letter = new_snp.get_uservariation().strip('-')
#                if new_letter in ['A', 'C', 'T', 'G']:
#                    father_genome.update_geno(new_letter, new_snp.get_name()) # if the snp is already in the father's genome, merge it
#            else: # if the snp isn't in the father's genome, add it
#                father_genome.add(new_snp)
#            if snp_name == 'rs3094315':
#                print(snp_name, genome.get_snp(snp_name).get_uservariation(), mother_genome.get_snp(snp_name).get_uservariation())       
    return father_genome
def test_compute_father():
    #father = compute_father('c:\\Users\\Jim\\documents\\DNA\\ftdna-37-CarolynSmithTestData.csv', children = ['c:\\Users\\Jim\\documents\\DNA\\ftdna-37-JamesSmithTestData2.csv'])
    father = compute_father('/Users/Jim/documents/DNA/ftdna-37-CarolynSmithTestData.csv', children = ['/Users/Jim/documents/DNA/ftdna-37-JamesSmithTestData.csv'])
    print('father coverage from Jim, ftdna', father.percent_covered(), len(father.get_snp_names()))#start with Jim
#     father = compute_father('c:\\Users\\Jim\\documents\\DNA\\23andme-37-CarolynSmithTest2.txt', 
#                             children = ['c:\\Users\\Jim\\documents\\DNA\\23andmeRay_Smith-37-Test2.txt', 'c:\\Users\\Jim\\documents\\DNA\\23andme-37-Suzanne_BusjahnTest2.txt', 
#                                         'c:\\Users\\Jim\\documents\\DNA\\23andme-37-CynthiaTest2.txt'],
#                             father_genome = father)
#     print('father coverage from rest of siblings', father.percent_covered(), len(father.get_snp_names()))
#     for snp_name in father.get_snp_names():
#         print(snp_name, father.get_snp(snp_name).percent_covered()[1], father.get_snp(snp_name).get_chromosome(), father.get_snp(snp_name).get_genotemp())
#     save_object(father, 'c:\\Users\\Jim\\documents\\DNA\\father.pkl')
#     father = get_object('c:\\Users\\Jim\\documents\\DNA\\father.pkl')
#     print('recovered father', father)
#     for snp_name in father.get_snp_names():
#         print(snp_name, father.get_snp(snp_name).percent_covered()[1], father.get_snp(snp_name).get_chromosome(), father.get_snp(snp_name).get_genotemp())
    # add me
#    father = compute_father('c:\\Users\\Jim\\documents\\DNA\\ftdna-37-CarolynSmithTestData.csv', sons = ['c:\\Users\\Jim\\documents\\DNA\\ftdna-37-JamesSmithTestData.csv'])
#    print('father coverage from Ray', father.percent_covered(), len(father.get_snp_names()))
#    for snp_name in father.get_snp_names():
#        print(snp_name, father.get_snp(snp_name).percent_covered()[1], father.get_snp(snp_name).get_chromosome(), father.get_snp(snp_name).get_genotemp())    
    # add Sue
#    father = compute_father('c:\\Users\\Jim\\documents\\DNA\\ftdna-37-CarolynSmithTestData.csv', daughters = ['c:\\Users\\Jim\\documents\\DNA\\23andme-37-Suzanne_BusjahnTest2.txt', 'c:\\Users\\Jim\\documents\\DNA\\23andme-37-CynthiaTest2.txt'], 
#                            father_genome = father)
#    print('father coverage from Sue & Cindy', father.percent_covered())
#    for snp_name in father.get_snp_names():
#        print(snp_name, father.get_snp(snp_name).percent_covered()[1], father.get_snp(snp_name).get_chromosome(), 
#              father.get_snp(snp_name).get_genotemp())    
def compare_parent_child(parent, child):
    '''
    input two snps, parent and child. 
    I could use snp variation from Entrez to help with error checking http://danielecook.com/getting_snp_dat/
    '''
    
    snp = '--'
    return_snp = Snp(child.get_name(), child.get_variation(), child.get_chromosome(), child.get_position(), child.get_strand(), '-', genotemp={})
    if not parent:  # there is a child id but no parent, so check for homozygocity
        if child.get_uservariation()[0] == child.get_uservariation()[1]: # homozygous, so had to have gotten the letter from each parent
            snp = child.get_uservariation()[0]
            return_snp.set_uservariation(snp + '-' )
            return_snp.update_geno(snp)
            return return_snp  
        else:
            return return_snp
    else:
        if child.get_uservariation()[0] == child.get_uservariation()[1]: # homozygous, so had to have gotten the letter from each parent
            snp = child.get_uservariation()[0]
        elif child.get_uservariation()[0] not in parent.get_uservariation():
            snp = child.get_uservariation()[0]
        elif child.get_uservariation()[1] not in parent.get_uservariation():
            snp = child.get_uservariation()[1]
    if child.get_uservariation()[0] not in parent.get_uservariation() and child.get_uservariation()[1] not in parent.get_uservariation():
            #print('child snp not in parent. Mutation? Read Error? ', child.get_name())
            snp = '--'  # this shouldn't happen. Probably a read error, maybe a mutation
    if snp == '--': 
        return return_snp   
    else:
        return_snp.set_uservariation(snp + '-' )
        return_snp.update_geno(snp)
        #return_snp.genotemp[snp] = 1
        return return_snp  
    
def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, -1)
def get_object(filename):
    with open(filename, 'rb') as input:
        return pickle.load(input)
def hamming_distance(s1, s2):
    """
    Return the Hamming distance between equal-length sequences and first position with a mismatch
    if no mismatch, 0, 0; 
    """
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    dist = sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2)) 
    if dist == 0:
        return 0, 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            return dist, i
def runtests():
    test_Phased_snp()
    test_pattern()
    test_informative()
    #test_parse_files()
    #test_compute_father()
    #test_parseFTDna()
    #test_parse23andme()
    pass
    
runtests()
#cProfile.run('test_compute_father()', 'restats')
#p = pstats.Stats('restats')
#p.strip_dirs().sort_stats('time').print_stats(30)