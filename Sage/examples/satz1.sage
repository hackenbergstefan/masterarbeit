# Test the following theorem:
# m,s integers. gcd(m,s) = 1. f | Phi_m(x) monic, irreducible. It is:
# f(x^s) = product_(d|s) f_d(x) where f_d(x) | Phi_md(x) monic irreducible.
# Is u in F_(q^sm) with Ord_q(u) = f(x^s) then 
# Ord_q(Tr_(F_(q^sm), F_(q^m)) (u) ) = f_1(x)

from sage.all import *
import random as rnd

load("./algorithmen.spyx")

def mytest(maxtests):
    testcounter = 0;
    for p in primes:
        #if p == 2: continue;
        for r in exponents:
            for n in extensions:
                q = p**r
                # check conditions
                if p.divides(n): continue;
                s = ordn(squarefree(n),q)
                if gcd(n,s) != 1 or s == 1: continue;
                l = ordn(n,q);
                #if s == l: continue;
                sbar = s;
                while p.divides(sbar): sbar = sbar / p;
                if sbar == 1: continue;
                dotest(p,r,q,n,s,l,sbar);
    testcounter = testcounter + 1;
    if testcounter >= maxtests: return;

# finde Elemente mit ord_n(q) != ord_nu(n)(q)
def mytest2(maxtests):
    testcounter = 0;
    for p in primes:
        #if p == 2: continue;
        for r in exponents:
            for n in extensions:
                q = p**r
                # check conditions
                if p.divides(n): continue;
                s = ordn(squarefree(n),q)
                if gcd(n,s) != 1 or s == 1: continue;
                l = ordn(n,q);
                if s != l:
                    print "--- test p=",p, "; r=",r, "; q=",q, "; n=",n, \
                            "; s=",factor(s), "; l=", factor(l);
            testcounter = testcounter + 1;
            if testcounter >= maxtests: return;

def dotest(p,r,q,n,s,l,sbar):
    #create fields
    F = GF(q, 'a');
    Fx = PolynomialRing(F, 'x');
    x = Fx.gen();
    E = GF(q**n,'b');
    G = GF(q**(s*n), 'c');
    Phi = Fx.cyclotomic_polynomial(n).factor();
    Phixs = Fx.cyclotomic_polynomial(n)(x**s).factor();
    nFac = sum(map(lambda (i,j):j, Phi))
    nFacxs = sum(map(lambda (i,j):j, Phixs))

    print "----------------------------------------"
    print "--- test p=",p, "; r=",r, "; q=",q, "; n=",n, "; s=",s,\
            "; l=", l, "; sbar=",sbar
    #print "Phi_n = Phi_"+str(n)+" =", Phi,\
    print  "Phi_n      -> "+str(nFac)+" factors "+str(len(Phi))+" distinct"\
            " of deg "+str(euler_phi(n) / nFac)+" == ",l
    #print Phi_n(x^s) =", Phixs, \
    print  "Phi_n(x^s) -> "+str(nFacxs)+" factors "+str(len(Phixs))+" distinct"\
            " of deg "+str(euler_phi(n)*sbar/nFacxs)+" == ",l
    for (f,i) in Phi:
        print "\ttest f(x) =",f
        if i != 1: print "ERROR"
        fds = f(x**s).factor()
        nfds = sum(map(lambda (i,j):j, fds));
        print "\t=> f(x^s) = ",fds, \
                "\tcheck "+str(len(divisors(sbar)))\
                +" =? "+str(nfds)+" divisors:"
        for (fd,j) in fds:
            #if j != 1: print "ERROR";
            cnt=0;
            for d in divisors(sbar):
                cnt = cnt+1;
                Phind = (Fx.cyclotomic_polynomial(n*d))**j;
                if fd.divides(Phind):
                    #print "\t\tf_"+str(d)+" = ("+str(fd)+")^"+str(j)+" | "\
                            #+"Phi_"+str(d)+"n =",Phind
                    print "\t\tf_"+str(d)+" | (Phi_"+str(d)+"n)^"+str(j), \
                            " -> "+str(cnt)+"th factor"

        #find suitable element for testing tau-order
        u = findSuitableElement(F,s,n,f);
        print "\twe have Ord_q(u) = f(x^s) for u =",u
        print "\t    and Ord_q(Tr_(G|E)(u)) =", \
                tau_order(trace(G,E,u),F)



# find u in F_q^(sn) with Ord_q(u) = f(x^s)
def findSuitableElement(F,s,n,f):
    q = F.order();
    G = GF(q**(s*n), 'b');
    x = f.parent().gen();
    l = pim(n,q**s-1)
    Gx = PolynomialRing(G, 'x');
    for (u,i) in Gx.cyclotomic_polynomial(n*l).roots():
        if tau_order(u,F) == f(x**s): return u;
    #fac = []
    #for (f,i) in (x**n-1).factor():
        #fac = fac + [f]*i
    #counter = 0;
    #for u in G:
        #if u == 0 or u == 1: continue;
        #if tau_order(u,F,factors=fac) == f(x**s): return u;
        #counter = counter + 1;
        ##if counter%100 == 0: print counter," elements of ", G.order(), "tested"



primes = primes_first_n(100)
exponents = range(1,10)
extensions = range(2,100)

#mytest(10);
