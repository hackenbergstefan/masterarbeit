from sage.all import *
import time
load("./algorithmen.spyx")

TIME = time.clock();

q = 9;
n = 4;

#setup
n = Integer(n);
q = Integer(q);
F = GF(q,'a');
E = F.extension(n,'b');
facAll = dict();
prodsAll = dict();
fieldsAll = dict();
divs = divisors(n);
for d in divs:
    G = F.extension(d, 'c');
    Gx = PolynomialRing(G,'x');
    fieldsAll[d] = G;
    facAll[d] = list((Gx.gen()**(n/d)-1).factor());
    prodsAll[d] = dict();
    for i in xmrange(map(lambda (_,i): i+1, facAll[d])):
        pr = Gx.one();
        for jidx, j in enumerate(i):
            pr *= facAll[d][jidx][0]**j;
        prodsAll[d][str(i)] = pr;
#print facAll
#print prodsAll
#print fieldsAll
# do test
countCN = 0;
countPCN = 0;

perOld = 0;
for idx,x in enumerate(E):
#b = E.gen()
#x = b^4+b^3+b^2+b+1;
    isCn = True;
    #isCnBrute = True;
    for d in divs[:len(divs)-1]:
        if not isNormal(x,fieldsAll[d],facAll[d],prodsAll[d]):
        #if trace(E,fieldsAll[d],x) == 0:
            isCn = False;
        #if not isNormal_bruteforce(x,fieldsAll[d]):
            #isCnBrute = False;
    #if isCn != isCnBrute:
        #print "x = ",x
    if isCn: countCN += 1;
    if isCn and x.multiplicative_order() == E.order()-1:
        countPCN += 1;
    per = idx/E.order();
    if floor(100*per) > perOld: 
        print round(100*per),"% tested."
    perOld = floor(100*per);



print "CN = ",countCN;
print "PCN = ",countPCN

print "TIME = ",time.clock()-TIME,"s";




