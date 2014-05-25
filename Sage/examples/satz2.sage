from sage.all import *
import sys

load("./algorithmen.spyx")


def factor_in_zeta(f,h,z):
    L = []
    for (g,j) in f.factor():
        for i in range(0,z.multiplicative_order()+1):
            if z**i == - g.map_coefficients(h)[0]:
                L = L + [(i,j)];
                break;
    return L

def print_fac_in_zeta(L):
    return join(map(lambda (i,j): "(x-z^"+str(i)+")^"+str(j), L));

p= 3 ; r= 1 ; q= 3 ; n= 10 ; 

s = ordn(squarefree(n),q);
l = pim(n,q**s-1);
b = gcd(l,n);

F = GF(q,'a');
Fx = PolynomialRing(F,'x');
K = F.extension(Integer(s), 'b');
Kx = PolynomialRing(K,'x');
E = F.extension(Integer(n),'c');
Ex = PolynomialRing(E,'x');
G = F.extension(Integer(n*s), 'd');
Gx = PolynomialRing(G, 'x');

Phi = Fx.cyclotomic_polynomial(n);
Phifac = Phi.factor();
PhiE = Ex.cyclotomic_polynomial(n*l);

hK = Hom(K,G)[0];
hE = Hom(E,G)[0];
hF = Hom(F,G)[0];
hFK = Hom(F,K)[0];
idG = Hom(G,G)[0];



print "----------------------------------------"
print "--- test p=",p, "; r=",r, "; q=",q, "; n=",n, "; s=",s,\
        "; l=", l

print "--- Setting"
print "F = GF(q); K = GF(q^s); E = GF(q^n); G = GF(q^ns); "

print "take z primitive (ns)-th root of unity in G"
z = -Gx.cyclotomic_polynomial(n*s).factor()[0][0][0]

print "Phi_n(x) = ", Fx.cyclotomic_polynomial(n).factor(), " over F";
print "\t= ", print_fac_in_zeta(factor_in_zeta(Kx.cyclotomic_polynomial(n),hK,z)), " over K";
sys.stdout.write("\t=");
for (f,i) in Kx.cyclotomic_polynomial(n).factor():
    sys.stdout.write("["+\
            print_fac_in_zeta(factor_in_zeta(f,hK,z))+"] * ");
sys.stdout.write("over G\n");
print "\t= ", Ex.cyclotomic_polynomial(n).factor(), " over E";
sys.stdout.write("\t=");
for (f,i) in Ex.cyclotomic_polynomial(n).factor():
    sys.stdout.write("["+\
            print_fac_in_zeta(factor_in_zeta(f.map_coefficients(hE),idG,z))+"] * ");
sys.stdout.write("over G\n");

#print "Phi_ns(x) = ", Fx.cyclotomic_polynomial(n*s).factor(), " over F";
#sys.stdout.write("\t=");
#for (f,i) in Fx.cyclotomic_polynomial(n*s).factor():
    #sys.stdout.write("["+\
            #print_fac_in_zeta(factor_in_zeta(f.map_coefficients(hF),z))+"] * ");
#sys.stdout.write("over K\n");
#sys.stdout.write("\t=");
#for (f,i) in Fx.cyclotomic_polynomial(n*s).factor():
    #sys.stdout.write("{");
    #for (g,j) in f.map_coefficients(hFK).factor():
        #sys.stdout.write("["+\
                #print_fac_in_zeta(factor_in_zeta(f.map_coefficients(hK),z))+"] * ");
    #sys.stdout.write("} * ")
#sys.stdout.write("over G\n");

print "take u primitive (nl)-th root of unity"
u = -Gx.cyclotomic_polynomial(n*l).factor()[0][0][0]

print "calc tau-orders"
for i in range(0,s):
    print "-- i =",i
    print "Ord_(q^(ns))(u^(q^(n*"+str(i)+"))) = ", tau_order(u**(q**(n*i)),G), \
            "=",print_fac_in_zeta(factor_in_zeta(tau_order(u**(q**(n*i)),G),idG,z));
    print "Ord_(q^s)(u^(q^(n*"+str(i)+"))) = ",tau_order(u**(q**(n*i)),K),\
            "=",print_fac_in_zeta(factor_in_zeta(tau_order(u**(q**(n*i)),K),hK,z));
    print "Ord_q(u^(q^(n*"+str(i)+"))) = ", tau_order(u**(q**(n*i)), F).factor();
#print "Ord_q(Tr_(E'|E)(u)) = ", tau_order(trace(G,E,u),F);

print "calc tau-orders"
for i in range(0,s):
    print "-- i =",i
    print "Ord_(q^s)(u^(q^"+str(i)+")) = ",tau_order(u**(q**i),K);
    print "Ord_q(u^(q^"+str(i)+")) = ", tau_order(u**(q**i), F).factor();
    print "Ord_(q^n)(u^(q^"+str(i)+")) = ", tau_order(u**(q**i), E).factor();


#h = Hom(F,E)[-1];
#print "take f(x) = (x-z^sbar) with z with ord = "+\
        #str(z.multiplicative_order())\
        #+" root of Phi_ns(x) in GF(q^e)"+\
        #" where e = ord_(ns)(q)"
#print "f(x^s) = ", print_fac_in_zeta(factor_in_zeta(f,z));
#print " with orders       "+\
        #join(map(lambda (i,_): str((z**i).multiplicative_order())\
        #+"           ", factor_in_zeta(f,z)));
#print "check factors of Phi_n(x^s)";
#for (g,i) in Fx.cyclotomic_polynomial(n)(Fx.gen()**s).factor():
    #for j in divisors(sbar):
        #if g.divides(Fx.cyclotomic_polynomial(n*j)):
            #print "g | Phi_n*"+str(j)+" for g=(",g,")^"+str(i)+" splits in";
    #print "\t", print_fac_in_zeta(factor_in_zeta(Ex(g),z))
    #print "\t          "+\
        #join(map(lambda (i,_): str((z**i).multiplicative_order())\
        #+"           ", factor_in_zeta(Ex(g),z)));

#for i in range(0,sbar):
    #for d in divisors(sbar):
        #if (x-z**(i*n+1)).divides(Ex.cyclotomic_polynomial(n*d)):
            #print "i = ",i," belongs to Phi_n*"+str(d);
