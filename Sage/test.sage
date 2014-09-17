import ast

def output2CSV(filein):
    n = Integer(str.split(filein,"_")[1])
    fileout = "outputs/PCNS_"+str(n)+".csv"
    with open(filein) as f:
        with open(fileout,'a') as fout:
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

