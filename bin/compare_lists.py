#!/usr/bin/env python3
# Compare metagenome taxonomic abundance rankings.
# Fredrik Boulund 2016
# RankBiasedOverlap and AverageOverlap implementations by 
# Ritesh Agrawal <https://github.com/ragrawal/measures>
__author__ = "Fredrik Boulund"
__date__ = "2016"
__version__ = "0.9b"
__doc__ = "Compare metagenome taxonomic abundance rankings in kaijuReport format."

from sys import argv, exit
import argparse

def parse_args():
    desc = ("Compare metagenome taxonomic abundance rankings "
            "(expects lists in kaijuReport output format). "
            "{author} {year}. Version {version}.".format(author=__author__,
                year=__date__, version=__version__))
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("LIST1")
    parser.add_argument("LIST2")
    parser.add_argument("-p", default=0.98, type=float,
            help="p-value to use when comparing lists [%(default)s].")
    parser.add_argument("-t", default=25, type=int,
            help="Number of top list values to display alongside comparison [%(default)s].")
    
    if len(argv) < 2:
        parser.print_help()
        exit(1)
    
    return parser.parse_args()


def parse_list(fn):
    with open(fn) as f:
        f.readline() # Skip header lines
        f.readline() # Skip header lines
        for line in f:
            if line.startswith("-"):
                return
            yield " ".join(line.split()[2:]).strip()


def main(list1, list2, p=0.98, t=10):
    for rank, items in enumerate(zip(list1[:t], list2[:t]), start=1):
        if items[0] == items[1]:
            print("{:>4} | {:<60} | =".format(rank, items[0]))
        else:
            print("{:>4} | {:<60} | {}".format(rank, *items))
    print("Average Overlap score (only {} first entries):  {}".format(t, AverageOverlap(list1, list2, t)))
    print("Rank Biased Overlap score (entire list):        {}".format(RBO(list1, list2, p)))
    print("-"*20)


def RBO(l1, l2, p = 0.98):
    """Calculates Ranked Biased Overlap (RBO) score. 
        l1 -- Ranked List 1
        l2 -- Ranked List 2
    @author: Ritesh Agrawal
    @Date: 13 Feb 2013
    @Description: This is an implementation of rank biased overlap score 
    (Refererence: http://www.umiacs.umd.edu/~wew/papers/wmz10_tois.pdf). 
    This is a modified implementation of  https://github.com/maslinych/linis-scripts/blob/master/rbo_calc.py
    It is a linear implementation of the RBO and assumes there are no
    duplicates and doesn't handle for ties. 
    """
    if l1 == None: l1 = []
    if l2 == None: l2 = []
    
    sl,ll = sorted([(len(l1), l1),(len(l2),l2)])
    s, S = sl
    l, L = ll
    if s == 0: return 0

    # Calculate the overlaps at ranks 1 through l 
    # (the longer of the two lists)
    ss = set([]) # contains elements from the smaller list till depth i
    ls = set([]) # contains elements from the longer list till depth i
    x_d = {0: 0}
    sum1 = 0.0
    for i in range(l):
        x = L[i]
        y = S[i] if i < s else None
        d = i + 1
        
        # if two elements are same then 
        # we don't need to add to either of the set
        if x == y: 
            x_d[d] = x_d[d-1] + 1.0
        # else add items to respective list
        # and calculate overlap
        else: 
            ls.add(x) 
            if y != None: ss.add(y)
            x_d[d] = x_d[d-1] + (1.0 if x in ss else 0.0) + (1.0 if y in ls else 0.0)     
        #calculate average overlap
        sum1 += x_d[d]/d * pow(p, d)
        
    sum2 = 0.0
    for i in range(l-s):
        d = s+i+1
        sum2 += x_d[d]*(d-s)/(d*s)*pow(p,d)

    sum3 = ((x_d[l]-x_d[s])/l+x_d[s]/s)*pow(p,l)

    # Equation 32
    rbo_ext = (1-p)/p*(sum1+sum2)+sum3
    return rbo_ext


def AverageOverlap(l1, l2, depth = 10):
    """Calculates Average Overlap score. 
        l1 -- Ranked List 1
        l2 -- Ranked List 2
        depth -- depth

    @author: Ritesh Agrawal
    @Date: 13 Feb 2013
    @Description: This is an implementation of average overlap measure for 
    comparing two score 
    (Refererence: http://www.umiacs.umd.edu/~wew/papers/wmz10_tois.pdf). 
    This is a modified implementation of  https://github.com/maslinych/linis-scripts/blob/master/rbo_calc.py
    It is a linear implementation of the RBO and assumes there are no
    duplicates and doesn't handle for ties. 
    """
    if l1 == None: l1 = []
    if l2 == None: l2 = []

    sl, ll = sorted([(len(l1), l1),(len(l2),l2)])
    s, S = sl  # s = length of smaller list, S = Smaller List
    l, L = ll  # l = length of longer list, L = Longer list
    #sanity check
    if s == 0: return 0
    depth = depth if depth < l else l
    
    # Calculate fraction of overlap from rank  at ranks 1 through depth
    # (the longer of the two lists)
    ss = set([])
    ls = set([])
    overlap = {0: 0}  # overlap holds number of common elements at depth d 
    sum1 = 0.0  

    for i in range(depth):
        # get elements from the two list
        x = L[i]
        y = S[i] if i < s else None
        depth = i+1
        # if the two elements are same, then we don't need
        # to them to the list and just increment the 
        if x == y: 
            overlap[depth] = overlap[i] + 2
        #else add items to the two list
        else:
            ls.add(x)
            if y != None: ss.add(y)
            overlap[depth] = overlap[i] + (2 if x in ss else 0) + (2 if y in ls else 0) 
        sum1 = sum1 + float(overlap[depth])/(len(S[0:depth]) + depth)

    return sum1/depth
    
    

if __name__ == "__main__":

    options = parse_args()
    list1 = list(parse_list(options.LIST1))
    list2 = list(parse_list(options.LIST2))

    print("Comparing '{}' with '{}'".format(options.LIST1, options.LIST2))
    print("-"*20)
    main(list1, list2, options.p, options.t)
