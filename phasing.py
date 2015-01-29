'''
Created on Jan 12, 2015

@author: Jim
'''
import util
SNP_LIMIT = 500

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

#def informative_pattern(pattern, c_pattern=''):
#    '''
#    take a pattern, and return empty string if not informative, or full string if informative
#    c_pattern is the current pattern, like '0111' or '1100'
#    informative are like 'AACC' or 'TCCC'
#    there should only be two letters, but I don't test
#    return '' if no informative pattern, or '0111' or some such informative pattern plus a Hamming distance
#    if the Hamming distance is 0 or 1, all is well. If greater, then there's a 50/50 chance the resulting pattern is wrong
#    '''
#    not_informative = ''
#    if '-' in pattern:
#        return not_informative, 0
#    if c_pattern == '':
#        first_digit = '0'
#        second_digit = '1'
#    else:
#        first_digit = c_pattern[0]
#        if first_digit == '0':
#            second_digit = '1'
#        else:
#            second_digit = '0'
#    new_pattern = first_digit
#    for letter in pattern[1:]:
#        if letter == pattern[0]:
#            new_pattern += first_digit
#        else:
#            new_pattern += second_digit
#    if sum([int(i) for i in new_pattern]) == 0 or sum([int(i) for i in new_pattern]) == len(pattern):
#        return not_informative, 0
#    if c_pattern == '':
#        return new_pattern, 0
#    dist = util.hamming_distance(c_pattern, new_pattern)
#    if dist <= len(pattern)/2:
#        return new_pattern, dist
#    #print('new_pattern, dist',  new_pattern, dist)
#    final_pattern = ''
#    for letter in new_pattern:
#        if letter == '0':
#            final_pattern += '1'
#        else:
#            final_pattern += '0'
#    return final_pattern, util.hamming_distance(c_pattern, final_pattern)
#def test_informative():
#    
#    assert informative_pattern('A-CC') == ('', 0)
#    assert informative_pattern('CCCC') == ('', 0)
#    assert informative_pattern('AACC') == ('0011', 0)
#    print(informative_pattern('AACC', '0011'))
#    assert informative_pattern('AACC', '0011') == ('0011', 0)
#    assert informative_pattern('AACC', '0010') == ('0011', 1)
#    assert informative_pattern('CACC', '0011') == ('1011', 1)
#    assert informative_pattern('ACCC', '1011') == ('1000', 2)
#    assert informative_pattern('CCAC', '1110') == ('1101', 2)
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
    phase('c:\\Users\\jim\\Documents\\DNA\\23andme-37-CarolynSmithTest2.txt', 
          ['c:\\Users\\jim\\Documents\\DNA\\23andme-37-CynthiaTest2.txt',
            'c:\\Users\\jim\\Documents\\DNA\\23andme-37-Suzanne_BusjahnTest2.txt',
            'c:\\Users\\jim\\Documents\\DNA\\23andmeJames_Smith-37-Test2.txt',
            'c:\\Users\\jim\\Documents\\DNA\\23andmeRay_Smith-37-Test2.txt'])     
# for Mac
#    phase('/Users/jim/Documents/DNA/23andme-37-CarolynSmithTest2.txt', 
#          ['/Users/jim/Documents/DNA/23andme-37-CynthiaTest2.txt',
#            '/Users/jim/Documents/DNA/23andme-37-Suzanne_BusjahnTest2.txt',
#            '/Users/jim/Documents/DNA/23andmeJames_Smith-37-Test2.txt',
#            '/Users/jim/Documents/DNA/23andmeRay_Smith-37-Test2.txt']) 

def phase_chrom(chrome, children):
    '''
    chrome is the mother's (could be father's) chromosome as a util.Chromo object
    children is a list of Chrome objects of the children of the parent
    return a Phase_chrom object containing the phased snps
    '''

    p_chrom = util.Phase_chrom(chrome.get_name()) # create empty phased chromosome with just a name
    snp_count = 0
    include_snps = []
    list_of_snp_names = chrome.get_snp_names()
    m_start_pattern, d_start_pattern = look_ahead(chrome, children, 0)
    c_m_inform = m_start_pattern
    c_d_inform = d_start_pattern
    p_chrom.set_m_start_pattern(m_start_pattern)
    p_chrom.set_d_start_pattern(d_start_pattern)
    print('m_start_pattern', m_start_pattern, 'd_start_pattern', d_start_pattern)
    print('snp_name,         mom,                 children, d_snps,   d_pat,  m_snps,    m_pat')
    for snp_index, snp_name in enumerate(chrome.get_snp_names()): # do this for each snp in mom's chromosome
        children_snps = get_children_snps(children, snp_name)
        if not no_snps(children_snps): include_snps.append(snp_name)  # only include the snp names where the children aren't all '--'
        mom_snp = chrome.get_snp(snp_name).get_uservariation()
        if mom_snp == '--':
            continue
        #p = pattern(mom_snp, children_snps)
        p_snp = util.Phased_snp(snp_name, mom_snp, children_snps, c_d_inform, c_m_inform) # make a phased snp
        p_snp.phase(c_m_inform, c_d_inform)
        if p_snp.get_d_dist() > 0:
            p_chrom.add_crossovers(util.Crossover(snp_name, p_snp.get_d_dist(), p_snp.get_d_crossover(), 'dad'))    
            c_d_inform = p_snp.get_d_pattern()
        if p_snp.get_m_dist() > 0:
            p_chrom.add_crossovers(util.Crossover(snp_name, p_snp.get_m_dist(), p_snp.get_m_crossover(), 'mom'))    
            c_m_inform = p_snp.get_m_pattern()        
        if snp_count < SNP_LIMIT and p_snp.get_name() in include_snps:
            print(snp_name.ljust(12), mom_snp, children_snps, p_snp.get_f_dad(), p_snp.get_d_pattern().ljust(6), 
                  p_snp.get_f_mom(), p_snp.get_m_pattern().ljust(6))
        snp_count += 1
    '''
    At this point I've gone through the chromosome once. I've got a list of snps that can be calculated and,
    the starting pattern, and the snp names where I detected a break.
    '''
    print('mom init informative', p_chrom.get_m_start_pattern(), ' init dad informative', p_chrom.get_d_start_pattern())
    print('number of crossovers', len(p_chrom.get_crossovers()))
    #print('excluded: ', exclude_snps)
    return p_chrom

def look_ahead(chrome, children, snp_name_index):
    '''
    given a Phase_chrom object, an index into the snp names it contains and children, a list of Chrome objects for all of the children
    scan ahead in the list of snps to find a complete informative patterns for mom and dad
    return complete patterns if possible
    
    '''
    m_start_pattern = ''
    d_start_pattern = ''
    snp_names = chrome.get_snp_names()
    for time in ['first time', 'second time']:
        for snp_index in range(snp_name_index, len(snp_names)):
            snp_name = snp_names[snp_index]
            children_snps = get_children_snps(children, snp_name)
            #if not no_snps(children_snps): include_snps.append(snp_name)  # only include the snp names where the children aren't all '--'
            mom_snp = chrome.get_snp(snp_name).get_uservariation()
            #p = pattern(mom_snp, children_snps)
            p_snp = util.Phased_snp(snp_name, mom_snp, children_snps) # make a phased snp
            # first see check if dad's informative pattern has been set
            if d_start_pattern == '':
                if p_snp.get_d_pattern() != '':
                    d_start_pattern = p_snp.get_d_pattern()
                    break  # stop processing the chromosome and start over
                continue # don't process further than this until dad's informative pattern has been found
            # next look for mom's informative pattern
            if m_start_pattern == '':
                if p_snp.get_m_pattern() != '':
                    m_start_pattern = p_snp.get_m_pattern()
                    break # start over a second time, now with both mom and dad's starting informative pattern
                continue # don't process further until mom's patter is found
    return m_start_pattern, d_start_pattern

def get_children_snps(children, snp_name):
    '''
    given a snp_name string like 'rs123', and children, a list of Chrom objects
    put together the list of children snps for that name, ['AT', 'TT', ...]
    '''
    children_snps = [] 
    for index, child in enumerate(children):
        #print(snp_name, child.get_snp(snp_name), child.get_name())
        if child.get_snp(snp_name):
            children_snps.append(child.get_snp(snp_name).get_uservariation())
        else:
            children_snps.append('--')
    return children_snps

def no_snps(snp_list):
    for snp_str in snp_list:
        if snp_str != '--':
            return False
    return True
def pattern_change(pattern1, pattern2):
    if pattern1 == pattern2:
        return False
    else:
        return True
def run_tests():
    #test_informative()
    test_phase()
    test_pattern()
    
run_tests()