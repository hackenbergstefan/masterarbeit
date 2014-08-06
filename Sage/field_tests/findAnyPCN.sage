# Finds Completely Free and Primitive Elements in an FiniteField Extension

from sage.all import *
from multiprocessing import Pool
load("./algorithmen.spyx")

filePath = "pcns_"

def runTest(data):
    p = data[0]
    n = data[1]
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
            x = findPCNElement(GF(q,'a'),n)
            with open(filePath+str(n),'a') as f:
                if x == False:
                    f.write(str(p)+'\t'+str(r)+'\tFALSE\n')
                    print "\tFALSE"
                else:
                    f.write(str(p)+'\t'+str(r)+'\n'); #+'\t'+str(x.minpoly())+'\n')
                    print "\tfound"
            f.close();

        r = r+1;
        q = p**r;
    return "";



## Tests all q < n^4
#def doTests(n):
    #pool = Pool();
    #pool.map(runTest, itertools.izip(primes(n**4), itertools.repeat(n)));


# Tests all q < n^4
def doTests(n):
    for p in primes(n**4):
        runTest([p,n])

