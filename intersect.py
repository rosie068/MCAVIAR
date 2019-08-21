import sys
import argparse
import numpy as np 

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='MCAVIAR is a statistical framework that quantifies the probability of each variant '
        'to be causal while allowing with arbitrary number of causal variants.')
    parser.add_argument('-o', '--out', required=True, dest='output_file',
                        help='output file name')
    parser.add_argument('-e', '--set_eur', required=True, dest='set_eur',
                        help='set for eur')
    parser.add_argument('-a', '--set_asn', required=True, dest='set_asn',
                        help='set for asn')

    args = parser.parse_args()
    o_fn = args.output_file
    Eur = args.set_eur
    Asn = args.set_asn

    eur_arr = []
    asn_arr = []

    f = open(Eur, 'r')
    for line in f:
        line = line.strip()
        #array = line.split()
        eur_arr.append(line)
    
    s = open(Asn, 'r')
    for line in s:
    	line = line.strip()
    	#arr = line.split()
    	asn_arr.append(line)

    intersect = list(set(eur_arr)&set(asn_arr))

    z = open(o_fn, 'w')
    for i in range(len(intersect)):
    	z.write(str(intersect[i]))
    	z.write("\n")
    z.close()

