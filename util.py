'''
Created on Dec 28, 2014
utilities for genome reconstruction
@author: Jim
'''
# used code from https://gist.github.com/philippbayer/1626642
import time, cProfile, pstats
import cPickle as pickle

class Genome:
    def __init__(self, name, version, snps):
        self.name = name   # who the person is
        self.version = version # 36, 37
        self.snps = {}  # store the snps in a dictionary
        self.errors = 0  # add to errors any time something impossible happens in the snp calulation, maybe should be an error log, list of messages
        for snp in snps:
            self.snps[snp.get_name()] = snp
    def __str__(self):
        x = self.snps.values()
        return 'genome for ' + self.name + ' version ' + str(self.version) + ' # snps ' + str(len(self.snps.values()))
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
        lineArray = line.replace(" ", "").replace("\n","").split("\t")
        allSnps.append(Snp(lineArray[0], None, lineArray[1], lineArray[2], None, lineArray[3]))  
    return allSnps
def test_parse23andme():
    #filein = open('c:\\Users\\Jim\\documents\\DNA\\23andme-37-CynthiaTest.txt')
    filein = open('c:\\Users\\Jim\\documents\\DNA\\23andmeRay_Smith-37-Test.txt')
    #filein = open('c:\\Users\\Jim\\documents\\DNA\\ftdna-37-JamesSmithTestData3.csv')
    b_time = time.time()
    #cindy_genome = Genome('cindy_vrazsity', 37, parse23andme(filein))
    ray_genome = Genome('ray_smith', 37, parse23andme(filein))
    #jim_genome = Genome('jim_smith', 37, parse23andme(filein))
    #print(cindy_genome)
    print('test_parse23andme', (time.time() - b_time)/60, len(ray_genome.get_snp_names()))
    #print('len', len(jim_genome.get_snp_names()))
    #print(ray_genome.percent_covered())

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
    father = compute_father('c:\\Users\\Jim\\documents\\DNA\\ftdna-37-CarolynSmithTestData.csv', children = ['c:\\Users\\Jim\\documents\\DNA\\ftdna-37-JamesSmithTestData2.csv'])
    print('father coverage from Jim, ftdna', father.percent_covered(), len(father.get_snp_names()))#start with Jim
    father = compute_father('c:\\Users\\Jim\\documents\\DNA\\23andme-37-CarolynSmithTest2.txt', 
                            children = ['c:\\Users\\Jim\\documents\\DNA\\23andmeRay_Smith-37-Test2.txt', 'c:\\Users\\Jim\\documents\\DNA\\23andme-37-Suzanne_BusjahnTest2.txt', 
                                        'c:\\Users\\Jim\\documents\\DNA\\23andme-37-CynthiaTest2.txt'],
                            father_genome = father)
    print('father coverage from rest of siblings', father.percent_covered(), len(father.get_snp_names()))
    for snp_name in father.get_snp_names():
        print(snp_name, father.get_snp(snp_name).percent_covered()[1], father.get_snp(snp_name).get_chromosome(), father.get_snp(snp_name).get_genotemp())
    save_object(father, 'c:\\Users\\Jim\\documents\\DNA\\father.pkl')
    father = get_object('c:\\Users\\Jim\\documents\\DNA\\father.pkl')
    print('recovered father', father)
    for snp_name in father.get_snp_names():
        print(snp_name, father.get_snp(snp_name).percent_covered()[1], father.get_snp(snp_name).get_chromosome(), father.get_snp(snp_name).get_genotemp())
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
     
def runtests():
    test_compute_father()
    #test_parseFTDna()
    #test_parse23andme()
    
runtests()
#cProfile.run('test_compute_father()', 'restats')
#p = pstats.Stats('restats')
#p.strip_dirs().sort_stats('time').print_stats(30)