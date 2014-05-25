from sage.all import *
import sys

load("./examples/satz1.sage")
load("./algorithmen.spyx")

def test1():
    #p= 2 ; r= 2 ; q= 4 ; n= 75 ; s= 2 ; l =  10;
    #p= 2 ; r= 4 ; q= 16 ; n= 99 ; s= 5 ; l =  15
    #p= 2 ; r= 2 ; q= 4 ; n= 49 ; s= 3 ; l =  21
    #p= 2 ; r= 1 ; q= 2 ; n= 7 ; s= 3 ; l =  3; sbar = 3;
    #p= 3 ; r= 1 ; q= 3 ; n= 22 ; s= 5 ; l =  5; sbar=5;
    #p= 3 ; r= 1 ; q= 3 ; n= 23 ; s= 11 ; l =  11; sbar = 11;
    #p= 3 ; r= 1 ; q= 3 ; n= 11 ; s= 5 ; l= 5 ; sbar= 5
    p= 3 ; r= 1 ; q= 3 ; n= 5 ; s= 4 ; l= 4 ; sbar= 4
    #p= 3 ; r= 1 ; q= 3 ; n= 37 ; s= 18 ; l= 18 ; sbar= 2;


    F = GF(q,'a');
    Fx = PolynomialRing(F,'x');
    e = ordn(n*sbar,q);
    E = F.extension(Integer(e),'c');
    Ex = PolynomialRing(E,'x');

    Phi = Fx.cyclotomic_polynomial(n);
    Phifac = Phi.factor();
    PhiE = Ex.cyclotomic_polynomial(n*s);
    z = - PhiE.factor()[0][0][0];
    PhiE = Ex.cyclotomic_polynomial(n*sbar);
    x = Ex.gen();
    f = (x - z**s)(x**s);
    
    dotest(p,r,q,n,s,l,sbar);

    print "----------------------------------------"
    print "--- test p=",p, "; r=",r, "; q=",q, "; n=",n, "; s=",s,\
            "; l=", l, "; sbar=",sbar


    h = Hom(F,E)[-1];
    print "take f(x) = (x-z^sbar) with z with ord = "+\
            str(z.multiplicative_order())\
            +" root of Phi_ns(x) in GF(q^e)"+\
            " where e = ord_(ns)(q)"
    print "f(x^s) = ", print_fac_in_zeta(factor_in_zeta(f,z));
    print " with orders       "+\
            join(map(lambda (i,_): str((z**i).multiplicative_order())\
            +"           ", factor_in_zeta(f,z)));
    print "check factors of Phi_n(x^s)";
    for (g,i) in Fx.cyclotomic_polynomial(n)(Fx.gen()**s).factor():
        for j in divisors(sbar):
            if g.divides(Fx.cyclotomic_polynomial(n*j)):
                print "g | Phi_n*"+str(j)+" for g=(",g,")^"+str(i)+" splits in";
        print "\t", print_fac_in_zeta(factor_in_zeta(Ex(g),z))
        print "\t          "+\
            join(map(lambda (i,_): str((z**i).multiplicative_order())\
            +"           ", factor_in_zeta(Ex(g),z)));

    for i in range(0,sbar):
        for d in divisors(sbar):
            if (x-z**(i*n+1)).divides(Ex.cyclotomic_polynomial(n*d)):
                print "i = ",i," belongs to Phi_n*"+str(d);

    #for (k,cs) in cosets(n,q):
        #product = prod(map(lambda (i,j): x - j**i, \
                #zip(cs,[zeta]*l)));
        #print "k = ",k
        #print "\t-> f= ",product
        #for (f,_) in Phifac:
            #if f.map_coefficients(h,E) == product:
                #print "\t",k, " ", f;


def cosets(n,q):
    L = []
    Ls = []
    l = ordn(n,q)
    for k in range(0,n):
        cs = coset_l(n,q,k);
        if not cs in Ls: #and len(cs) == l:
            Ls = Ls + [cs];
            L = L + [(k,cs)];
    return L

def coset_l(n,q,l):
    return sorted(list(set(map(lambda i: (l*q**i)%n, range(0,n)))));

def print_cosets(n,q):
    for (k,cs) in cosets(n,q):
        print k, cs


def factor_in_zeta(f,z):
    L = []
    for (g,j) in f.factor():
        for i in range(0,z.multiplicative_order()):
            if z**i == - g[0]:
                L = L + [(i,j)];
                break;
    return sorted(L)

def print_fac_in_zeta(L):
    return join(map(lambda (i,j): "(x - zeta^"+str(i)+")^"+str(j), L));




test1();
