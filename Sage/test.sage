q = 2
n = 7
F = GF(q);
R.<x> = F['x'];
E.<w> = F.extension(n);
Q.<y> = E['y'];
f = 1
for i in range(1,n+1):
    if gcd(i,n) != 1:
        f = f * R.cyclotomic_polynomial(i)
print f.factor()
phi = Q.convert_map_from(R);
print phi
for i in f.factor():
    print phi(i[0]).roots()
