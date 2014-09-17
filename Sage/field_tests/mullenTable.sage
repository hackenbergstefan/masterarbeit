from sage.all import *
from multiprocessing import Pool
import time
import datetime

load("./enumeratePCNs.spyx")


def enumeratePCNs_wrapper(qn):
    q = qn[0]
    n = qn[1]
    F = GF(Integer(q),'a');
    isProcessed = False;
    if os.path.exists(filePath):
        f = open(filePath,'r');
        for l in f.readlines():
            if re.search('^\(q,n\) = \('+str(q)+', '+str(n)+'\)',l) != None:
                isProcessed = True;
                print "(",q,",",n,") processed"
                return
        f.close();
    if not isProcessed:
        print "(",q,",",n,") not processed"
        ret = countCompleteSubmoduleGenerators(F,n,binaryPowers=True)
        print "(q,n) = ",(q,n)," -> (cn,pcn,submod_gens,time) = ", ret
        with open(filePath,'a') as f:
            f.write("(q,n) = "+str((q,n))+" -> (cn,pcn,submod_gens,time) = "+
                    str(ret)+"\n")
        f.close();

#------------------------------------------------------------------------------
# mullen setup
st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')
filePath = "mullenTableC_struct_"+st+".txt"
SETUP = \
[[2, xrange(2,19)], \
 [3, xrange(2,13)], \
 [4, xrange(2,10)], \
 [5, xrange(2,9)], \
 [7, xrange(2,7)], \
 [8, xrange(2,6)], \
 [9, xrange(2,6)], \
 ]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
filePath = "mullenTableC_struct_"+st+".txt"
SETUP = \
[[2, xrange(2,25)], \
 [3, xrange(2,19)], \
 [4, xrange(2,15)], \
 [5, xrange(2,13)], \
 [7, xrange(2,12)], \
 [8, xrange(2,10)], \
 [9, xrange(2,10)], \
 ]
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#SETUP = \
#[[3, xrange(18,21)]]
#------------------------------------------------------------------------------


GENLIST = []
for q, nlist in SETUP:
    for n in nlist:
        GENLIST += [[q,n]]

#------------------------------------------------------------------------------
#test n = 6, q = 2, 3, 4, ..., 32, ...
filePath = "mullenTableC_struct_n=3_"+st+".txt"
GENLIST = []
n = 3
for p in primes(100):
    for e in range(1,10):
        GENLIST += [[p**e,n]]

GENLIST = sorted(GENLIST)
##------------------------------------------------------------------------------



def main():
    pool = Pool(1);
    pool.imap_unordered(enumeratePCNs_wrapper, GENLIST)
