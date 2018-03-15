import os,sys,time

UDHS = sys.argv[1]
name = sys.argv[2]
bwfile = sys.argv[3]
bwoverbed = sys.argv[4]


def K27ac_count(bwoverbed,bwfile,UDHS,name):
	signal_bigWigAverageOverBed=[]
	signal_bigWigAverageOverBed=os.popen("%s -sampleAroundCenter=1000 %s %s stdout | cut -f 4" % (bwoverbed, bwfile, UDHS)).read().split("\n")[:-1]
	n=1
	f=open(name+"_Strength.txt","w")
	for i in signal_bigWigAverageOverBed:
		f.write(str(n)+"\t"+str(i)+"\n")
		n=n+1
	f.close()


def main():
	K27ac_count(bwoverbed,bwfile,UDHS,name)

if __name__ == '__main__':
	main()

