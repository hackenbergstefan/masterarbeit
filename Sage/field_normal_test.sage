from sage.all import *

load("./algorithms_completely_free.spyx")

def test_all_complete_normal(F,n):
    E = F.extension(n,'a');
    return allCompletelyNormal(F,E);

def test_all_complete_normal_rnd():
    cmpNormals = [];
    notCmpNormals = [];
    for p in primes(100):
        if p == 2: continue;
        for q in map(lambda i: p**i, range(1,10)):
            for n in createFactorList(q+1,10):
                if q < 1000 and n < 500:
                    if test_all_complete_normal(GF(q,'a'), n):
                        cmpNormals = cmpNormals + [(q,n)]
                    else: notCmpNormals = notCmpNormals + [(q,n)]
                    print "cmpNormals", cmpNormals
                    print "notCmpNormals", notCmpNormals


def createFactorList(qp, l):
    fs = factor(qp);
    L = map( lambda (i,k): i, fs);
    L = filter(lambda i: i%2==1, L);
    for i in range(2,l+1):
        L = L + map(lambda j: j**i, L);
    return L
