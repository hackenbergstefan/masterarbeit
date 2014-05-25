
# q = p = 7
# n = 9
# => n in Sq
# Zerfall von Phim für m teilt n
def example1():
    F = GF(7);
    P = PolynomialRing(F,'x');
    n = 9;
    for i in range(1,n+1):
        if 9%i == 0:
            print i, "\t",P.cyclotomic_polynomial(i), \
                "\t", P.cyclotomic_polynomial(i).factor();

example1();
