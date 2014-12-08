import ast
import re
import string
import subprocess
import itertools
import glob
import datetime

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


def enumsPCN2Latex_wrapper(basePath, onlyNList=[3,4,6]):
    globs = glob.glob(basePath+'mullenTableC_struct_*')
    allinfo = dict()
    
    for f in globs:
        with open(f) as fin:
            filename = f.split("/")[-1]
            isOnlyNormal = ( re.search("_norm_",filename) != None )
            
            dt = datetime.datetime.strptime(\
                    re.search("_(\d.*).txt",filename).groups()[0],\
                    "%Y-%m-%d_%H:%M:%S")
            print filename, " isOnlyNormal?", isOnlyNormal, " dt = ",dt
            for line in fin:
                if all(c in string.whitespace for c in line): continue
                try:
                    s = re.search("\((\d+), (\d+)\) -> .* = \((-?\d+), (-?\d+), \{(.*)\}", line)\
                            .groups()
                    q = Integer(s[0])
                    n = Integer(s[1])
                    cn = Integer(s[2])
                    pcn = Integer(s[3])
                    if pcn == 0: pcn = "--"
                except:
                    #print "FATAL ERROR on file ", f, " in line"
                    #print line
                    continue
                gens = re.findall("\((\d+), (\d+), (\d+)\): (\d+)", s[4])
                gens = sorted(gens,\
                        key=lambda t: Integer(t[0]))

                p,r = q.factor()[0]

                if cn < 0:
                    cn = prod(map(lambda gns: Integer(gns[3]), gens))

                key = (p,r,n)
                mustUpdate = True
                if allinfo.has_key(key):
                    mustUpdate = False
                    #print "already done: ", p,r,n
                    info = allinfo[key]
                    if isOnlyNormal and not info.has_key("norm"):
                        info["norm"] = cn
                        info["pnorm"] = pcn
                        #print "\tupdate norm!",info
                    elif not isOnlyNormal and not info.has_key("cn"):
                        info["cn"] = cn
                        info["pcn"] = pcn
                        #print "\tupdate cn!", info
                    elif dt > info["date"]: 
                        info["date"] = dt
                        if isOnlyNormal:
                            info["norm"] = cn
                            info["pnorm"] = pcn
                        else:
                            info["cn"] = cn
                            info["pcn"] = pcn
                            info["gens"] = gens
                if mustUpdate:
                    info = dict()
                    info["p"] = p
                    info["r"] = r
                    info["q"] = q
                    info["n"] = n
                    info["gens"] = gens
                    info["date"] = dt
                    if isOnlyNormal:
                        info["norm"] = cn
                        info["pnorm"] = pcn
                    else:
                        info["cn"] = cn
                        info["pcn"] = pcn
                    allinfo[key] = info
                
    #print sorted(allinfo.keys())
    # do outputs
    # first remove all old files
    purge("../Latex/tables/", "enumerations.*")
    purge("../Tables/Enumerations/", "enumerations.*")

    #do P outputs
    psDone = []
    pFilesCreated = []
    keysSorted = sorted(allinfo.keys())
    for idx,(p,r,n) in enumerate(keysSorted):
        #print "p=",p," r=",r," n=",n
        # continue if we have only one value for this p
        if not has_more_than(2, keysSorted, 0, p): continue
        info = allinfo[(p,r,n)]
        #print "info = ",info
        fileout = "../Latex/tables/enumerationsPCN_P_"+str(info["p"])+".tex"
        fileoutTables = "../Tables/Enumerations/enumerationsPCN_P_"+str(info["p"])+".csv"
        fileoutNorm = "../Latex/tables/enumerationsPN_P_"+str(info["p"])+".tex"
        fileoutTablesNorm = "../Tables/Enumerations/enumerationsPN_P_"+str(info["p"])+".csv"
        if info["p"] in psDone:
            fout = open(fileout,'a')
            foutTables = open(fileoutTables,'a')
            if info.has_key("norm"):
                foutNorm = open(fileoutNorm,'a')
                if not os.path.exists(fileoutTablesNorm):
                    foutTablesNorm = open(fileoutTablesNorm,'w')
                    foutTablesNorm.write("q, p, r, n, N, PN\n")
                else: foutTablesNorm = open(fileoutTablesNorm,'a')
        else:
            fout = open(fileout,'w')
            foutTables = open(fileoutTables,'w')
            foutTables.write("q, p, r, n, CN, PCN, gens\n")
            if info.has_key("norm"):
                foutNorm = open(fileoutNorm,'w')
                foutTablesNorm = open(fileoutTablesNorm,'w')
                foutTablesNorm.write("q, p, r, n, N, PN\n")
            psDone += [p]

        gensString = ""
        gensString_tables = ""
        for idx, (k,t,pi,count) in enumerate(info["gens"]):
            gensString += "$("+k+","+t+","+pi+")"
            gensString_tables += "("+k+" "+t+" "+pi+")"
            k = Integer(k)
            t = Integer(t)
            pi = Integer(pi)
            if isRegular(p,r,k,t,pi):
                gensString += "^\\dagger"
                gensString_tables += "*"
            gensString += "$: "+count
            gensString_tables += ": "+count
            if idx < len(info["gens"])-1: 
                gensString += ",\\ "
                gensString_tables += " "

        isBasicString = ""
        isBasicString_tables = ""
        if isCompletelyBasic(info["p"],info["r"],info["n"]):
            isBasicString = "$^\\ast$"
            isBasicString_tables = "*"

        outString = str(info["q"])\
                +" & "+str(info["p"])\
                +" & "+str(info["r"])\
                +" & "+str(info["n"])\
                +" & "+str(info["cn"])+isBasicString\
                +" & "+str(info["pcn"])\
                +" & "+gensString\
                +"\\\\"
        outString_tables = str(info["q"])\
                +", "+str(info["p"])\
                +", "+str(info["r"])\
                +", "+str(info["n"])\
                +", "+str(info["cn"])+isBasicString_tables\
                +", "+str(info["pcn"])\
                +", "+gensString_tables
        
        #print outString
        fout.write(outString+"\n")
        fout.close()
        foutTables.write(outString_tables+"\n")
        foutTables.close()

        if not p in pFilesCreated: pFilesCreated += [p]
        
        if info.has_key("norm"):
            outString = str(info["q"])\
                    +" & "+str(info["p"])\
                    +" & "+str(info["r"])\
                    +" & "+str(info["n"])\
                    +" & "+str(info["norm"])\
                    +" & "+str(info["pnorm"])\
                    +"\\\\"
            outString_tables = str(info["q"])\
                    +", "+str(info["p"])\
                    +", "+str(info["r"])\
                    +", "+str(info["n"])\
                    +", "+str(info["norm"])\
                    +", "+str(info["pnorm"])
            foutNorm.write(outString+"\n")
            foutNorm.close()
            foutTablesNorm.write(outString_tables+"\n")
            foutTablesNorm.close()
    print "pFilesCreated = ",pFilesCreated
    
    #do N outputs
    nsDone = []
    nFilesCreated = []
    keysSorted = sorted(allinfo.keys(), key=lambda prn: (prn[2],prn[0]**prn[1]))
    for idx,(p,r,n) in enumerate(keysSorted):
        #print "p=",p," r=",r," n=",n
        # continue if we have only one value for this n
        if ( not has_more_than(2, keysSorted, 2, n) ) \
                or ( onlyNList != None and len(onlyNList) > 0 and \
                       not n in onlyNList ): continue
        info = allinfo[(p,r,n)]
        #print "info = ",info
        fileout = "../Latex/tables/enumerationsPCN_N_"+str(info["n"])+".tex"
        fileoutNorm = "../Latex/tables/enumerationsPN_N_"+str(info["n"])+".tex"
        fileout_tables = "../Tables/Enumerations/enumerationsPCN_N_"+str(info["n"])+".csv"
        fileoutNorm_tables = "../Tables/Enumerations/enumerationsPN_N_"+str(info["n"])+".csv"
        if info["n"] in nsDone:
            fout = open(fileout,'a')
            fout_tables = open(fileout_tables,'a')
        else:
            fout = open(fileout,'w')
            fout_tables = open(fileout_tables,'w')
            fout_tables.write("q, p, r, n, CN, PCN, gens\n")
            nsDone += [n]
        gensString = ""
        gensString_tables = ""
        for idx, (k,t,pi,count) in enumerate(info["gens"]):
            gensString += "$("+k+","+t+","+pi+")"
            gensString_tables += "("+k+" "+t+" "+pi+")"
            k = Integer(k)
            t = Integer(t)
            pi = Integer(pi)
            if isRegular(p,r,k,t,pi):
                gensString += "^\\dagger"
                gensString_tables += "*"
            gensString += "$: "+count
            gensString_tables += ": "+count
            if idx < len(info["gens"])-1: 
                gensString += ",\\ "
                gensString_tables += " "

        isBasicString = ""
        isBasicString_tables = ""
        if isCompletelyBasic(info["p"],info["r"],info["n"]):
            isBasicString = "$^\\ast$"
            isBasicString_tables = "*"

        outString = str(info["q"])\
                +" & "+str(info["p"])\
                +" & "+str(info["r"])\
                +" & "+str(info["n"])\
                +" & "+str(info["cn"])+isBasicString\
                +" & "+str(info["pcn"])\
                +" & "+gensString\
                +"\\\\"
        outString_tables = str(info["q"])\
                +", "+str(info["p"])\
                +", "+str(info["r"])\
                +", "+str(info["n"])\
                +", "+str(info["cn"])+isBasicString_tables\
                +", "+str(info["pcn"])\
                +", "+gensString_tables
        
        #print outString
        fout.write(outString+"\n")
        fout.close()
        fout_tables.write(outString_tables+"\n")
        fout_tables.close()

        if not n in nFilesCreated: nFilesCreated += [n]
    
    print "nFilesCreated = ",nFilesCreated
    #return allinfo

def findPCN2CSV_allinone(basePath,n):
    fileins = []
    pairsToCheck = list(runGenerator(n**4))
    fileins += glob.glob(basePath+'pcns_trinom_'+str(n)+'_*')
    fileins += glob.glob(basePath+'pcns_additional_'+str(n)+'_*')


    processedPairs = dict()

    for filein in fileins:
        with open(filein) as f:
            dt = datetime.datetime.strptime(
                    re.search("_(\d+\-\d+\-\d+_\d+:\d+:\d)",filein).groups()[0],\
                    "%Y-%m-%d_%H:%M:%S")
            isAdditional = ( re.search("additional", filein) != None )
            for line in f:
                if all(c in string.whitespace for c in line): continue
                splits = str.split(line,"\t")
                try:
                    p = Integer(splits[0])
                    r = Integer(splits[1])
                    poly = splits[2].strip()
                except:
                    continue
                if not processedPairs.has_key((p,r)):
                    processedPairs[(p,r)] = (p,r,poly,dt,isAdditional)
                    pairsToCheck.remove((p,r))
                else:
                    # only update if newer and this file is not "additional"
                    if ( not isAdditional and processedPairs[(p,r)][4] )\
                            or ( isAdditional == processedPairs[(p,r)][4] and \
                            processedPairs[(p,r)][3] < dt ):
                        processedPairs[(p,r)] = (p,r,poly,dt,isAdditional)
        f.close()

    # do output
    # first clear old files
    purge("../Tables/", "pcns_"+str(n)+"_*")

    keysSorted = sorted(processedPairs.keys(), \
            key=lambda prpoly: (prpoly[1],prpoly[0]))
    for idx,key in enumerate(keysSorted):
        p,r,poly,dt,isAdditional = processedPairs[key]
        F = GF(p**r,'a')
        Fx = PolynomialRing(F,'x')
        polyF = Fx(poly)
        polyList = list(polyF)
        
        fout = open("../Tables/pcns_"+str(n)+"_"+str(r)+".csv",'a')
        # if new r or p write fileheader
        if idx == 0 or r != processedPairs[keysSorted[idx-1]][1]:
            if r == 1:
                fout.write("p,\tpoly\n")
            else:
                fout.write("p,\tpoly,\tmodulus\n")
        # write to ../Tables csv file
        isAdditionalString = ""
        if isAdditional: isAdditionalString += "!"
        if r == 1:
            fout.write(str(p)+isAdditionalString+",\t"+str(polyF)+"\n")
        else:
            fout.write(str(p)+isAdditionalString+",\t"+str(polyF)+",\t"\
                    +str(F.modulus().change_variable_name('a'))+"\n")

    fout.close()
    #print "max file: ", curFileNumber
    print "missing pairs: ", pairsToCheck
    with open(basePath+"missing_pairs_"+str(n),'w') as f:
        f.write(str(pairsToCheck))
    f.close()



def findPCN2Latex_wrapper(basePath,n):
    pairsToCheck = list(runGenerator(n**4))
    for r in itertools.count(1):
        globs = glob.glob(basePath+'pcns_trinom_'+str(n)+'_'+str(r)+'_*')
        if len(globs) == 0: break
        findPCN2Latex(globs,pairsToCheck=pairsToCheck)
    globs = glob.glob(basePath+'pcns_additional_'+str(n)+'_*')
    if len(globs) != 0: 
        findPCN2Latex(globs,pairsToCheck=pairsToCheck)

    print "missing pairs: ", pairsToCheck
    with open(basePath+"missing_pairs_"+str(n),'w') as f:
        f.write(str(pairsToCheck))
    f.close()


def runGenerator(border,pRange=None,rRange=None):
    if pRange == None:
        p = 1
    else:
        p = pRange[0]
    while p < border :
        p = next_prime(p)
        if pRange != None and p > pRange[1]: return
        # consider only rs in rRange
        if rRange != None:
            for r in xrange(rRange[0],rRange[1]+1):
                if p**r > border: break
                yield p,r
        # consider all rs
        else:
            r = 1
            q = p**r
            while q < border:
                yield p,r
                r += 1
                q = p**r

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# tests if l has more than n values of type l[..][idx] == value
def has_more_than(n, l, idx, value):
    return len(filter(lambda item: item[idx] == value, l)) > n


def purge(dir, pattern):
    for f in os.listdir(dir):
    	if re.search(pattern, f):
    		os.remove(os.path.join(dir, f))


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

