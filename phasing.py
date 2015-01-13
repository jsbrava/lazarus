'''
Created on Jan 12, 2015

@author: Jim
'''

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
    
    print(pattern('CA', ['CC', 'AC', 'CT', 'TT']))
    print(pattern('CC', ['CC', 'AC', 'CT', 'TT']))    

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
        
    
def run_tests():
    test_pattern()
    
run_tests()