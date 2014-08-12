# Finds Completely Free and Primitive Elements in an FiniteField Extension

from sage.all import *
from multiprocessing import Pool
load("./algorithmen.spyx")

filePath = "pcnsC_"

def runTest(data):
    p,n,totalCount = data
    r = 1
    q = p**r
    while q < n**4:
        #check if q has already processed
        isProcessed = False;
        if os.path.exists(filePath+str(n)):
            f = open(filePath+str(n),'r');
            for l in f.readlines():
                if re.search('^'+str(p)+"\t"+str(r),l) != None:
                    isProcessed = True;
                    print "(",p,",",r,") processed"
                    break;
            f.close();
        if not isProcessed:
            print "(",p,",",r,") not processed"
            x = findAnyPCN_internalC(GF(q,'a'),n)
            with open(filePath+str(n),'a') as f:
                if x == False:
                    f.write(str(p)+'\t'+str(r)+'\tFALSE\n')
                    print "\tFALSE"
                else:
                    f.write(str(p)+'\t'+str(r)+'\t'+str(x)+"\n"); 
                    print "\tfound"
            f.close();
        counter = sum(1 for line in open(filePath+str(n),'r'))
        print "\t", counter," / ", totalCount, " done -> "\
                , round(counter/totalCount*100,2),"%"

        r = r+1;
        q = p**r;



## Tests all q < n^4
def doTests(n):
    totalCount = 0
    for p in primes(n**4):
        r = 1
        q = p**r
        while q < n**4:
            totalCount += 1
            r += 1
            q = p**r

    pool = Pool();
    pool.imap_unordered(runTest, ([p,n,totalCount] for p in primes(n**4)))


## Tests all q < n^4
#def doTests(n):
    #for p in primes(n**4):
        #runTest([p,n])

