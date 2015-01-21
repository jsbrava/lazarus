'''
Created on Jan 12, 2015

@author: Jim
'''
import util

def pattern(mom_snp, sib_list):
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
            if sib[0] in mom_snp:
                mom_pattern[index] = sib[0] # mom CC, sib CA, so C came from mom, A from Dad
                dad_pattern[index] = sib[1]
            else:
                mom_pattern[index] = sib[1] # or the other way around
                dad_pattern[index] = sib[0]
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
def test_pattern():
    
    assert pattern('CA', ['CC', 'AC', 'CT', 'TT']) == ('C-C-', 'C-TT')
    assert pattern('CC', ['CC', 'AC', 'CT', 'TT']) == ('CCC-', 'CATT')   

def informative_pattern(pattern):
    '''
    take a pattern, and return empty string if not informative, or full string if informative
    informative are like 'AACC' or 'TCCC'
    there should only be two letters, but I don't test
    '''
    not_informative = ''
    if '-' in pattern:
        return not_informative
    for letter in pattern[1:]:
        if pattern[0] != letter:
            return pattern
    return not_informative
def test_informative():
    assert informative_pattern('AACC') == 'AACC'
    assert informative_pattern('A-CC') == ''
    assert informative_pattern('CCCC') == ''

class Phase:
    def __init__(self, snp_id, mom_snp, sib_snp_list, mom_sibs='', mom_pat='', dad_sibs='', dad_pat=''):
        '''
        snp_id, string, 'rs234'
        mom_snp, string, 'AC'
        sib_snp_list, ['AC', ...]
        mom_sibs, string, the letter that mom gave to each of her children
        mom_pat, string, 'A--A' or '-C-C', useful pattern for phasing
        dad_sibs, string, the letter that dad gave to each of his children
        dad_pat, string, useful pattern for phasing
        based on Whit Athey's paper, http://www.jogg.info/62/files/Athey.pdf
        if you are missing mom, not dad, then just switch them
        '''
        self.snp_id = snp_id
        self.mom_snp = mom_snp
        self.sib_snp_list = sib_snp_list
        self.mom_sibs = mom_sibs
        self.mom_pat = mom_pat
        self.dad_sibs = dad_sibs
        self.dad_pat = dad_pat
        
def phase(mom, children):
    '''
    
    '''
    mom_genome = util.Genome('mom', 37, util.parse_files(mom))
    mom_chroms = mom_genome.get_chromosomes()
    children_genomes = [util.Genome('child', 37, util.parse_files(child)) for child in children]
    print(mom_genome)
    for chrom in mom_chroms:
        print(chrom)
        phase_chrom(chrom, [child.get_chromosome(chrom.get_name()) for child in children_genomes])
def test_phase():
    phase('/Users/jim/Documents/DNA/23andme-37-CarolynSmithTest2.txt', 
          ['/Users/jim/Documents/DNA/23andme-37-CynthiaTest2.txt',
            '/Users/jim/Documents/DNA/23andme-37-Suzanne_BusjahnTest2.txt',
            '/Users/jim/Documents/DNA/23andmeJames_Smith-37-Test2.txt',
            '/Users/jim/Documents/DNA/23andmeRay_Smith-37-Test2.txt']) 

def phase_chrom(chrome, children):
    '''
    chrome is the mother's (could be father's) chromosome as a util.Chromo object
    children is a list of Chrome objects of the children of the parent
    '''
    print('snp_name,         mom,                 children, d_snps,   d_pat,  m_snps,    m_pat')
    exclude_snps = []
    for snp_name in chrome.get_snp_names():
        
        children_snps = [] 
        for child in children:
            #print(snp_name, child.get_snp(snp_name), child.get_name())
            if child.get_snp(snp_name):
                children_snps.append(child.get_snp(snp_name).get_uservariation())
            else:
                children_snps.append('--')
            #children_snps += "."
        if no_snps(children_snps): exclude_snps.append(snp_name)
        mom_snp = chrome.get_snp(snp_name).get_uservariation()
        p = pattern(mom_snp, children_snps)
        if snp_name not in exclude_snps:
            print(snp_name.ljust(12), mom_snp, children_snps, p[1], informative_pattern(p[1]).ljust(6), p[0], informative_pattern(p[0]).ljust(6))
    #print('excluded: ', exclude_snps)
def no_snps(snp_list):
    for snp_str in snp_list:
        if snp_str != '--':
            return False
    return True

def run_tests():
    test_informative()
    test_phase()
    test_pattern()
    
run_tests()