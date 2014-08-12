from sage.all import *
from multiprocessing import Pool
import itertools
load("./algorithmen.spyx")

filePath = "extended_mullenTable.txt"

def process_submodules_internalC(pen):
    p,e,n = pen
    q = p**e
    F = GF(Integer(q),'a');
    ret = countCompleteSubmoduleGenerators_internalC(F,n, maxEta=60*60)
    if ret:
        print "(q,n) = ",(q,n)," -> (cn,pcn,submod_gens,time) = ", ret
        with open(filePath,'a') as f:
            f.write("(q,n) = "+str((q,n))+" -> (cn,pcn,submod_gens,time) = "+
                    str(ret)+"\n")
        f.close();
    else: 
        print "(q,n) = ",(q,n)," too long!"


def main():
    pool = Pool();
    pool.imap_unordered(process_submodules_internalC, \
            itertools.product(primes(1000),xrange(1,100),xrange(2,100)))

