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
        ret = countCompleteSubmoduleGenerators(F,n,binaryPowers=True,\
                testPrimitivity=TEST_PRIMITIVITY,\
                onlyNormal=TEST_ONLY_NORMAL)
        print "(q,n) = ",(q,n)," -> (cn,pcn,submod_gens,time) = ", ret
        with open(filePath,'a') as f:
            f.write("(q,n) = "+str((q,n))+" -> (cn,pcn,submod_gens,time) = "+
                    str(ret)+"\n")
        f.close();


###############################################################################
###############################################################################
SETUP_NUMBER = "pcn p5"
###############################################################################
###############################################################################
st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')

SETUP = []

#------------------------------------------------------------------------------
# mullen setup
if SETUP_NUMBER == "pcn p1":
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
    TEST_PRIMITIVITY = True
    TEST_ONLY_NORMAL = False
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
if SETUP_NUMBER == "pcn p2":
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
    TEST_PRIMITIVITY = True
    TEST_ONLY_NORMAL = False
#------------------------------------------------------------------------------

if SETUP_NUMBER == "pcn p3":
    filePath = "mullenTableC_struct_"+st+".txt"
    SETUP = \
    [[2, xrange(25,28)]]
    TEST_PRIMITIVITY = True
    TEST_ONLY_NORMAL = False

if SETUP_NUMBER == "pcn p4":
    filePath = "mullenTableC_struct_"+st+".txt"
    SETUP = \
    [[3, xrange(18,22)]]
    TEST_PRIMITIVITY = True
    TEST_ONLY_NORMAL = False


#------------------------------------------------------------------------------
# test higher primes
if SETUP_NUMBER == "pcn p5":
    filePath = "mullenTableC_struct_"+st+".txt"
    SETUP = \
    [[11, [2,5,7]], \
     [121, [2,5]], \
     [13, [2,5,7]], \
     [169, [2,5]], \
     [17, [2,5,7]], \
     [289, [2,5]], \
     [19, [2,5,7]], \
     [361, [2]], \
     [23, [2,5,7]], \
     [529, [2]], \
     [29, [2,5,7]], \
     [29**2, [2]], \
     [31, [2,5,7]], \
     [31**2, [2]], \
     [37, [2,5,7]], \
     [37**2, [2]], \
     [41, [2,5,7]], \
     [41**2, [2]], \
     [43, [2,5,7]], \
     [43**2, [2]], \
    ]
    TEST_PRIMITIVITY = True
    TEST_ONLY_NORMAL = False

###############################################################################
## Without Primitivity ########################################################
###############################################################################
#------------------------------------------------------------------------------
if SETUP_NUMBER == "cn p1":
    filePath = "mullenTableC_struct_noPrim_"+st+".txt"
    SETUP = \
    [[2, xrange(31,35)], \
     [3, xrange(21,30)], \
     [4, xrange(15,30)], \
     [5, xrange(13,30)], \
     [7, xrange(12,20)], \
     [8, xrange(10,20)], \
     [9, xrange(10,20)], \
     ]
    TEST_PRIMITIVITY = False
    TEST_ONLY_NORMAL = False
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
if SETUP_NUMBER == "cn p2":
    filePath = "mullenTableC_struct_noPrim_"+st+".txt"
    SETUP = \
    [[2, xrange(25,30)], \
     [3, xrange(19,30)], \
     [4, xrange(15,30)], \
     [5, xrange(13,30)], \
     [7, xrange(12,20)], \
     [8, xrange(10,20)], \
     [9, xrange(10,20)], \
     ]
    TEST_PRIMITIVITY = False
    TEST_ONLY_NORMAL = False
#------------------------------------------------------------------------------

###############################################################################
## Enumerate normals only P = ... #############################################
###############################################################################


#------------------------------------------------------------------------------
if SETUP_NUMBER == "pn p1":
    filePath = "mullenTableC_struct_norm_"+st+".txt"
    SETUP = \
    [[2, xrange(2,30)], \
     [3, xrange(2,30)], \
     [4, xrange(2,30)], \
     [5, xrange(2,30)], \
     [7, xrange(2,20)], \
     [8, xrange(2,20)], \
     [9, xrange(2,20)], \
     ]
    TEST_PRIMITIVITY = True
    TEST_ONLY_NORMAL = True
#------------------------------------------------------------------------------


GENLIST = []
for q, nlist in SETUP:
    p,e = list(factor(q))[0]
    for n in nlist:
        if not TEST_ONLY_NORMAL or not isCompletelyBasic(p,e,n):
            GENLIST += [[q,n]]

###############################################################################
## Enumerate N = ... ##########################################################
###############################################################################

#------------------------------------------------------------------------------
#test n = 3, q = 2, 3, 4, ..., 32, ...
if SETUP_NUMBER == "pcn n1":
    filePath = "mullenTableC_struct_n=3_"+st+".txt"
    GENLIST = []
    n = 3
    for p in primes(100):
        for e in range(1,10):
            GENLIST += [[p**e,n]]

    GENLIST = sorted(GENLIST)

    TEST_PRIMITIVITY = True

#------------------------------------------------------------------------------
#test n = 4, q = 2, 3, 4, ..., 32, ...
if SETUP_NUMBER == "pcn n2":
    filePath = "mullenTableC_struct_n=4_"+st+".txt"
    GENLIST = []
    n = 4
    for p in primes(100):
        for e in range(1,10):
            GENLIST += [[p**e,n]]

    GENLIST = sorted(GENLIST)
    TEST_PRIMITIVITY = True

#------------------------------------------------------------------------------
#test n = 6, q = 2, 3, 4, ..., 32, ...
if SETUP_NUMBER == "pcn n3":
    filePath = "mullenTableC_struct_n=6_"+st+".txt"
    GENLIST = []
    n = 6
    for p in primes(100):
        for e in range(1,10):
            GENLIST += [[p**e,n]]

    GENLIST = sorted(GENLIST)
    TEST_PRIMITIVITY = True


###############################################################################
## Without Primitivity ########################################################
###############################################################################

#------------------------------------------------------------------------------
#test n = 6, q = 43, ...
if SETUP_NUMBER == "cn n1":
    filePath = "mullenTableC_struct_n=6_noPrim_"+st+".txt"
    GENLIST = []
    n = 6
    for p in primes(100):
        for e in range(1,10):
            GENLIST += [[p**e,n]]

    GENLIST = sorted(GENLIST)[23:]
    TEST_PRIMITIVITY = False
##------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#test n = 4, q = 243, ...
if SETUP_NUMBER == "cn n2":
    filePath = "mullenTableC_struct_n=4_noPrim_"+st+".txt"
    GENLIST = []
    n = 4
    for p in primes(200):
        for e in range(1,20):
            GENLIST += [[p**e,n]]

    GENLIST = sorted(GENLIST)[243:]
    TEST_PRIMITIVITY = False
##------------------------------------------------------------------------------

def main():
    pool = Pool(2);
    pool.imap_unordered(enumeratePCNs_wrapper, GENLIST)
