from sage.all import *
from multiprocessing import Pool
load("./algorithmen.spyx")

filePath = "mullenTableC.txt"

def process(qn):
    q = qn[0]
    n = qn[1]
    F = GF(Integer(q),'a');
    print "(q,n) = ",(q,n)," -> (cn,pcn,time) = ", countCNAndPCN(F,n)

def process_parallel(qn):
    q = qn[0]
    n = qn[1]
    F = GF(Integer(q),'a');
    print "(q,n) = ",(q,n)," -> (cn,pcn,time) = ", countCNAndPCN_parallel(F,n)

def process_submodules(qn):
    q = qn[0]
    n = qn[1]
    F = GF(Integer(q),'a');
    ret = countCompleteSubmoduleGenerators(F,n)
    print "(q,n) = ",(q,n)," -> (cn,pcn,submod_gens,time) = ", ret
    with open(filePath,'a') as f:
        f.write("(q,n) = "+str((q,n))+" -> (cn,pcn,submod_gens,time) = "+
                str(ret)+"\n")
    f.close();

def process_submodules_internalC(qn):
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
        ret = countCompleteSubmoduleGenerators_internalC(F,n)
        print "(q,n) = ",(q,n)," -> (cn,pcn,submod_gens,time) = ", ret
        with open(filePath,'a') as f:
            f.write("(q,n) = "+str((q,n))+" -> (cn,pcn,submod_gens,time) = "+
                    str(ret)+"\n")
        f.close();

#------------------------------------------------------------------------------
# mullen setup
#filePath = "mullenTableC_new.txt"
#SETUP = \
#[[2, xrange(2,19)], \
 #[3, xrange(2,13)], \
 #[4, xrange(2,10)], \
 #[5, xrange(2,9)], \
 #[7, xrange(2,7)], \
 #[8, xrange(2,6)], \
 #[9, xrange(2,6)], \
 #]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#SETUP = \
#[[3, xrange(18,21)]]
#------------------------------------------------------------------------------


#GENLIST = []
#for q, nlist in SETUP:
    #for n in nlist:
        #GENLIST += [[q,n]]

##------------------------------------------------------------------------------
# test n = 6, q = 2, 3, 4, ..., 32, ...
filePath = "mullenTableC_n=6.txt"
GENLIST = []
n = 6
for p in primes(100):
    for e in range(1,10):
        GENLIST += [[p**e,n]]

GENLIST = sorted(GENLIST)
##------------------------------------------------------------------------------

#def main():
    #pool = Pool();
    #pool.map(process, GENLIST)

#def main():
    #for qn in GENLIST:
        #process_parallel(qn)

def main():
    pool = Pool();
    pool.imap_unordered(process_submodules_internalC, GENLIST)
