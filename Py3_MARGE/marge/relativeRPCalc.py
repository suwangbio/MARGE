import os,sys
import numpy as np

RP = sys.argv[1]
RPmedians = sys.argv[2]
outfilename = sys.argv[3]

   
def adjusted( RP, RPmedians, outfilename ):

    RPs = open(RP).readlines()
    rps = [ rp.strip().split( '\t' ) for rp in RPs ]

    rp = [ float( i[4] ) for i in rps ]

    histrps = open( RPmedians ).readlines()
    histrps = [ histrp.strip().split('\t') for histrp in histrps ]
    medianrp = [ float( i[0] ) for i in histrps ]

    adjrp = np.array( rp,dtype=np.float64 )/(np.array( medianrp,dtype=np.float64 ) + 1)
    adjrp = [ round(adjrp[i],3) for i in range(len(adjrp)) ]

    outf = open(outfilename,'w')
    for n in range(len(rps)):
        rpinfo = rps[n][0:5]
        rpinfo.append(str(adjrp[n]))
        info = rpinfo + rps[n][5:]
        outf.write('\t'.join(info) + '\n')
    outf.close()

def main():
    adjusted( RP, RPmedians, outfilename )

if __name__ == '__main__':
    main()
    
