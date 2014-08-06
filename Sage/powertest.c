#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

void matmul(int *mat, int *vec, int *ret, int m, int charac){
    int i,j;
    for(i=0;i<m;i++){
        ret[i] = 0;
        for(j=0;j<m;j++){
            ret[i] += mat[i*m + j]*vec[j];
        }
        ret[i] %= charac;
    }
}

void initPoly(int *p, int m) {
    int i;
    for(i=0;i<m;i++)
        p[i]=0;
}

void printArr(int *arr, int l){
    int j;
    for(j=0;j<l;j++){
        printf(" %i", arr[j]);
    }
    printf("\n");
}

/* adds 2 polynomials, where p3 is the result
 * MUST have same length */
void addPoly(int *p1, int *p2, int *p3, int m, int charac) {
    int i;
    for(i=0;i<m;i++)
        p3[i]=(p1[i]+p2[i])%charac;
}

/* p1 - p2 , where p3 is the result
 * MUST have same length */
void subtrPoly(int *p1, int *p2, int *p3, int m, int charac) {
    int i;
    for(i=0;i<m;i++)
        p3[i]=(p1[i]-p2[i])%charac;
}

void multiplyPoly(int *p1, int m1, int *p2, int m2, int *p3, int m3, int charac) {
    int i,j;
    initPoly(p3,m3);
    for(i=0;i<m1;i++)
        for(j=0;j<m2;j++)
            p3[i+j] = (p3[i+j]+p1[i]*p2[j])%charac;
}


void moduloPoly(int *p1, int m1, int *mod, int m, int charac){
    int deg1, degmod;
    int i,j;
    int quo;
    //get degrees
    for(i=m1-1;i>=0;i--){
        if( p1[i] != 0){
            deg1 = i;
            break;
        }
    }
    for(i=m-1;i>=0;i--){
        if( mod[i] != 0){
            degmod = i;
            break;
        }
    }
    //printf("degMod=%i, degP1=%i\n", degmod,deg1);
    
    //printArr(p1,m1);

    //make polynomial division
    for(i=deg1-degmod; i>=0; i--){
        quo = (p1[i+degmod]*modInv(mod[degmod],charac))%charac;
        //printf("i=%i p1[i+degmod]=%i mod[degmod]=%i mod[degmod]^(-1)=%i ",i,p1[i+degmod],mod[degmod],modInv(mod[degmod],charac));
        //printf("quo=%i\n", quo);
        for(j=degmod;j>=0;j--){
            p1[i+j] = (p1[i+j] - mod[j]*quo)%charac;
        }
        //printArr(p1,m1);
    }
}




/*
 * calculates a^(-1) mod p
 */
int modInv(int a, int p){
    int s,t,r,old_s,old_t,old_r, quo, tmp;
    s = 0; old_s = 1;
    t = 1; old_t = 0;
    r = p; old_r = a;
    while(r != 0){
        quo = old_r / r;
        tmp = r; r = old_r - quo*tmp; old_r = tmp;
        tmp = s; s = old_s - quo*tmp; old_s = tmp;
        tmp = t; t = old_t - quo*tmp; old_t = tmp;
    }
    // bezout coeffs: old_s old_t
    // gcd: old_r
    // quotients for gcd: t s
    return abs(old_s);
}



/*
 * calculates g(sigma)(x) where g is a polynomial and sigma the frobenius
 * application of frobenius is given by mats, i.e.
 * mats is array of m x m matrices, x an array of m ints
 * g is an array of glen arrays of length m
 */
void applyFrob(int *x, int *x_mipo, int *g, int glen, int *mats, int *ret, int m, int charac){
    int mSize = m*m;
    int i,j;
    int *tmp = malloc(m*sizeof(int));
    int *tmp2 = malloc((m+glen)*sizeof(int));
    initPoly(ret,m);
    for(i=0;i<glen;i++){
        bool allZ = true;
        for(j=0;j<m;j++){
            if(g[i*glen+j] != 0){
                allZ = false;
                break;
            }
        }
        if(allZ == true)
            continue;
        //printf("i=%i",i);
        matmul(mats+i*mSize, x, tmp, m, charac);
        //printf(" mat*x = "); printArr(tmp, m);
        //printf("tmp*gi for \n\ttmp="); printArr(tmp,m);
        //printf("\tgi="); printArr(g+i*glen, m);
        multiplyPoly(tmp,m, g+i*glen, m, tmp2, 2*m, charac);
        //printf("\t=>tmp*gi = "); printArr(tmp2,2*m);
        //printf("tmp2 mod mipo for \n\ttmp2="); printArr(tmp2,2*m); 
        //printf("\tmipo="); printArr(x_mipo,m+1);
        moduloPoly(tmp2, 2*m, x_mipo, m+1, charac);
        //printf("\t=>tmp2 mod mipo = "); printArr(tmp2,m);
        for(j=0;j<m;j++){
            ret[j] += tmp2[j];
        }
    }
}
