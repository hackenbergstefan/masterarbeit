import ast
import re
import string
import subprocess
import itertools
import glob

load("./enumeratePCNs.spyx")


def pcns2CSV(filein):
    n = Integer(str.split(filein,"_")[1])
    fileout = "outputs/PCNS_"+str(n)+".csv"
    with open(filein) as f:
        with open(fileout,'a') as fout:
            fout.write("p,r,mipo\n")
            for line in f:
                split = str.split(line,"\t")
                p = Integer(split[0])
                r = Integer(split[1])
                el = map(Integer,list(ast.literal_eval(split[2])))
                E = GF(p**(r*n),'a')
                outstr = ""
                for i,fi in enumerate(E(el).minpoly()):
                    if fi != 0:
                        outstr += str(i)+":"+str(fi)+" "
                fout.write(str(p)+","+str(r)+","+outstr+"\n")
        fout.close()
    f.close()


def pcnsPolynom2CSV(filein):
    n = Integer(str.split(filein,"_")[2])
    fileout = "outputs/PCNS_weight_"+str(n)+".csv"
    with open(filein) as f:
        with open(fileout,'a') as fout:
            fout.write("p,r,mipo\n")
            for line in f:
                split = str.split(line,"\t")
                p = Integer(split[0])
                r = Integer(split[1])
                Px = PolynomialRing(GF(p),'x')
                outstr = ""
                for i,fi in enumerate(Px(split[2])):
                    if fi != 0:
                        outstr += str(i)+":"+str(fi)+" "
                fout.write(str(p)+","+str(r)+","+outstr+"\n")
        fout.close()
    f.close()

def enumsPCNP2csv(filein):
    fileout = "outputs/enumerationsPCN_P_"
    with open(filein) as f:
        for line in f:
            if all(c in string.whitespace for c in line): continue
            reMatch = re.compile("(\d+), (\d+)").search(line)
            reGroup = reMatch.groups()
            q = Integer(reGroup[0])
            n = Integer(reGroup[1])
            pr = tuple(factor(q)[0])
            p,r = pr

            reGroup = map(lambda s: s.replace(",",""),\
                    re.compile("(\d+), (\d+), {(.+)}")\
                    .search(line,reMatch.end()).groups())

            fileout_tmp = fileout+str(p)+".csv"
            alreadyDone = False
            if not os.path.exists(fileout_tmp):
                with open(fileout_tmp,'a') as fout:
                    fout.write("q,p,r,n,cn,pcn,gens\n")
                fout.close()
            else:
                with open(fileout_tmp) as fout:
                    if re.search(\
                            "\\n"+str(q)+","+str(p)+","+str(r)+","+str(n)\
                            ,fout.read()):
                        print "found ",p,q,r,n," in ",fileout_tmp
                        alreadyDone = True
                fout.close()
            if not alreadyDone:
                with open(fileout_tmp,'a') as fout:
                    fout.write(str(q)+","+str(p)+","+str(r)+","+\
                            str(n)+","+(",".join(reGroup))+"\n")
                fout.close()
                print "added ",str(q)+","+str(p)+","+str(r)+","+\
                            str(n)+","+(",".join(reGroup)),\
                            " to ",fileout_tmp


    f.close()

def enumsPNP2csv(filein):
    fileout = "outputs/enumerationsPN_P_"
    with open(filein) as f:
        for line in f:
            if all(c in string.whitespace for c in line): continue
            reMatch = re.compile("(\d+), (\d+)").search(line)
            reGroup = reMatch.groups()
            q = Integer(reGroup[0])
            n = Integer(reGroup[1])
            pr = tuple(factor(q)[0])
            p,r = pr

            reGroup = map(lambda s: s.replace(",",""),\
                    re.compile("(\d+), (\d+), {(.+)}")\
                    .search(line,reMatch.end()).groups())

            fileout_tmp = fileout+str(p)+".csv"
            alreadyDone = False
            if not os.path.exists(fileout_tmp):
                with open(fileout_tmp,'a') as fout:
                    fout.write("q,p,r,n,n,pn,gens\n")
                fout.close()
            else:
                with open(fileout_tmp) as fout:
                    if re.search(\
                            "\\n"+str(q)+","+str(p)+","+str(r)+","+str(n)\
                            ,fout.read()):
                        print "found ",p,q,r,n," in ",fileout_tmp
                        alreadyDone = True
                fout.close()
            if not alreadyDone:
                with open(fileout_tmp,'a') as fout:
                    fout.write(str(q)+","+str(p)+","+str(r)+","+\
                            str(n)+","+(",".join(reGroup))+"\n")
                fout.close()
                print "added ",str(q)+","+str(p)+","+str(r)+","+\
                            str(n)+","+(",".join(reGroup)),\
                            " to ",fileout_tmp


    f.close()


def enumsPCNN2csv(filein):
    fileoutN = "outputs/enumerationsPCN_N_"
    with open(filein) as f:
        for line in f:
            if all(c in string.whitespace for c in line): continue
            reMatch = re.compile("(\d+), (\d+)").search(line)
            reGroup = reMatch.groups()
            q = Integer(reGroup[0])
            n = Integer(reGroup[1])
            pr = tuple(factor(q)[0])
            p,r = pr

            reGroup = map(lambda s: s.replace(",",""),\
                    re.compile("(\d+), (\d+), {(.+)}")\
                    .search(line,reMatch.end()).groups())

            # check enumerationsN_
            fileout_tmp = fileoutN+str(n)+".csv"
            alreadyDone = False
            if not os.path.exists(fileout_tmp):
                with open(fileout_tmp,'a') as fout:
                    fout.write("q,p,r,n,cn,pcn,gens\n")
                fout.close()
            else:
                with open(fileout_tmp) as fout:
                    if re.search(\
                            "\\n"+str(q)+","+str(p)+","+str(r)+","+str(n)\
                            ,fout.read()):
                        print "found ",p,q,r,n," in ",fileout_tmp
                        alreadyDone = True
                fout.close()
            if not alreadyDone:
                with open(fileout_tmp,'a') as fout:
                    fout.write(str(q)+","+str(p)+","+str(r)+","+\
                            str(n)+","+(",".join(reGroup))+"\n")
                fout.close()
                print "added ",str(q)+","+str(p)+","+str(r)+","+\
                            str(n)+","+(",".join(reGroup)),\
                            " to ",fileout_tmp
    f.close()



def enumsPCN2Latex(filein):
    fileout = "../Latex/tables/"+(str.split(str.split(filein,"/")[-1],".")[0])+".tex"
    print fileout
    with open(fileout,'w') as fout:
        with open(filein) as f:
            for line in f:
                if all(c in string.whitespace for c in line): continue
                splits = str.split(line,",")
                try:
                    q = Integer(splits[0])
                    p = Integer(splits[1])
                    r = Integer(splits[2])
                    n = Integer(splits[3])
                    cn = Integer(splits[4])
                    pcn = Integer(splits[5])
                except:
                    continue
                regexGens = re.compile("\((\d+) (\d+) (\d+)\): (\d+)")\
                        .findall(splits[6])

                regexGens = sorted(regexGens,\
                        key=lambda t: Integer(t[0]))
                gensString = ""
                for idx, (k,t,pi,count) in enumerate(regexGens):
                    gensString += "$("+k+","+t+","+pi+")"
                    k = Integer(k)
                    t = Integer(t)
                    pi = Integer(pi)
                    if isRegular(p,r,k,t,pi):
                        gensString += "^\\dagger"
                    gensString += "$: "+count
                    if idx < len(regexGens)-1: gensString += ",\\ "

                isBasicString = ""
                if isCompletelyBasic(p,r,n):
                    isBasicString = "$^\\ast$"

                outString = str(q)\
                        +" & "+str(p)\
                        +" & "+str(r)\
                        +" & "+str(n)\
                        +" & "+str(cn)+isBasicString\
                        +" & "+str(pcn)\
                        +" & "+gensString\
                        +"\\\\"
                print outString
                fout.write(outString)

        f.close()
    fout.close()

def findPCN2Latex(fileins, border=lambda n: n**4):
    print "process ",fileins
    splitter = str.split(str.split(fileins[0],"/")[-1],"_")
    n = Integer(splitter[2])
    r = Integer(splitter[3])
    fileout = "../Latex/tables/pcns_"+str(n)+"_"+str(r)+"__"

    border = border(n)
    primeList = []
    p = 2
    while p**r < border:
        primeList += [p]
        p = next_prime(p)

    curFileNumber = 0
    fout = open(fileout+str(curFileNumber)+".tex", 'w')
    for filein in fileins:
        lineCount = Integer(\
                subprocess.check_output(["wc","-l",filein])\
                .split(" ")[0]) - 1
        with open(filein) as f:
            for curLineNum,line in enumerate(f):
                if all(c in string.whitespace for c in line): continue
                splits = str.split(line,"\t")
                try:
                    p = Integer(splits[0])
                    r = Integer(splits[1])
                    poly = splits[2].strip()
                except:
                    continue
                primeList.remove(p)
                #print "p=",p," r=",r," poly='",poly,"'"
                F = GF(p**r,'a')
                Fx = PolynomialRing(F,'x')
                polyF = Fx(poly)
                polyList = list(polyF)

                # check trinom
                outString = "\\textbf{"+str(p)+":} "
                if polyF.hamming_weight() == 3:
                    outString += "$"+str(polyList[0]).replace("*","")+"$"
                    # add (n-1)-th Coefficient if needed
                    if polyList[n-1] != F.one():
                        outString += ", $"\
                                +str(polyList[n-1]).replace("*","")\
                                +"$"
                else:
                    outString += ",\\ ".join([str(i)+":\\,$"\
                            +str(c).replace("*","")+"$" \
                            for i,c in enumerate(polyList[:-1]) if c != 0])

                if curLineNum < lineCount:
                    outString += ", "

                #print curLineNum, outString
                #COUNTER += 1
                #if COUNTER > 5000: break
                fout.write(outString+"\n")
                if curLineNum > 0 and curLineNum%1000 == 0:
                    fout.close()
                    curFileNumber += 1
                    fout = open(fileout+str(curFileNumber)+".tex", 'w')
            fout.close()
    f.close()
    print "max file: ", curFileNumber
    print "missing primes: ",primeList


def findPCN2Latex_wrapper(basePath,n):
    for r in itertools.count(1):
        globs = glob.glob(basePath+'pcns_trinom_'+str(n)+'_'+str(r)+'_*')
        if len(globs) == 0: return
        findPCN2Latex(globs)



###############################################################################
###############################################################################
###############################################################################
###############################################################################

def checkWrongComplBasics():
    for f in os.listdir("./outputs"):
        if not f.startswith('enumerations'): continue
        with open("./outputs/"+f,'r') as fin:
            for line in fin:
                if all(c in string.whitespace for c in line): continue
                splits = str.split(line,",")
                try:
                    q = Integer(splits[0])
                    p = Integer(splits[1])
                    r = Integer(splits[2])
                    n = Integer(splits[3])
                    cn = Integer(splits[4])
                    pcn = Integer(splits[5])
                except:
                    continue
                
                divs = get_not_completely_basic_divisors(p,r,n)
                divsWRONG = get_completely_basic_divisors_WRONG(p,r,n)
                if divs == divsWRONG:
                    continue
                print "unequal on p=",p," r=",r," n=",n," divs=",divs," divsWRONG=",divsWRONG


                regexGens = re.compile("\((\d+) (\d+) (\d+)\): (\d+)")\
                        .findall(splits[6])

                regexGens = sorted(regexGens,\
                        key=lambda t: Integer(t[0]))
                gensString = ""
                for idx, (k,t,pi,count) in enumerate(regexGens):
                    k = Integer(k)
                    t = Integer(t)
                    pi = Integer(pi)
                    if isRegular(p,r,k,t,pi):
                        continue

                    divsMod = divisors(get_module_character(k,t,pi))
                    if filter(lambda x: x in divs, divsMod) \
                            != filter(lambda x: x in divsWRONG, divsMod):
                        print "ERROR on: p=",p," r=",r," n=",n,\
                                " k=",k," t=",t," pi=",pi
                        print "\t",filter(lambda x: x in divs, divsMod)," vs wrong",\
                                filter(lambda x: x in divsWRONG, divsMod)

        fin.close()

