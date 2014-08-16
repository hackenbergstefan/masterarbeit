#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>


unsigned long long ipow(int base, int exp)
{
    unsigned long long result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}


inline void matmul(int *mat, int *vec, int *ret, int m, int charac){
    int i,j;
    long tmp;
    for(i=0;i<m;i++){
        ret[i] = 0;
        for(j=0;j<m;j++){
            tmp = mat[i*m + j]*vec[j];
            ret[i] += (int)(tmp%charac);
        }
    }
}

inline void matmulShort(int *mat, int *vec, int *ret, int m, int charac){
    int i,j;
    for(i=0;i<m;i++){
        ret[i] = 0;
        for(j=0;j<m;j++){
            ret[i] += mat[i*m + j]*vec[j];
        }
    }
}

inline void initPoly(int *p, int m) {
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
inline void addPoly(int *p1, int *p2, int *p3, int m, int charac) {
    int i;
    for(i=0;i<m;i++)
        p3[i] = p1[i]+p2[i];
}

/* p1 - p2 , where p3 is the result
 * MUST have same length */
void subtrPoly(int *p1, int *p2, int *p3, int m, int charac) {
    int i;
    for(i=0;i<m;i++)
        p3[i] = p1[i]-p2[i];
}

inline void multiplyPoly(int *p1, int m1, int *p2, int m2, 
        int *p3, int m3, int charac) {
    int i,j;
    long tmp;
    initPoly(p3,m3);
    for(i=0;i<m1;i++){
        for(j=0;j<m2;j++){
            tmp = ((long)p1[i]*p2[j])%charac;
            p3[i+j] = (int)((p3[i+j]+tmp)%charac);
        }
    }
    /*for(i=0;i<m3;i++){*/
        /*p3[i] %= charac;*/
        /*if(p3[i] < 0) p3[i] += charac;*/
    /*}*/
}

inline void multiplyPolyShort(int *p1, int m1, int *p2, int m2, 
        int *p3, int m3, int charac) {
    int i,j;
    int deg1, deg2;
    for(i=m1-1;i>=0;i--){
        if(p1[i] != 0){
            deg1 = i;
            break;
        }
    }
    for(i=m2-1;i>=0;i--){
        if(p2[i] != 0){
            deg2 = i;
            break;
        }
    }
    initPoly(p3,m3);
    for(i=0;i<=deg1;i++){
        for(j=0;j<=deg2;j++){
            p3[i+j] += p1[i]*p2[j];
        }
    }
}

inline void multiplyPolyShortKnownDeg(int *p1, int m1, int *p2, int m2, int deg2,
        int *p3, int m3, int charac) {
    int i,j;
    int deg1, deg2Tmp;
    for(i=m1-1;i>=0;i--){
        if(p1[i] != 0){
            deg1 = i;
            break;
        }
    }
    for(i=m1-1;i>=0;i--){
        if(p2[i] != 0){
            deg2Tmp = i;
            break;
        }
    }
    initPoly(p3,m3);
    
    if(deg2 == -1) return;

    for(i=0;i<=deg1;i++){
        for(j=0;j<=deg2;j++){
            p3[i+j] += p1[i]*p2[j];
        }
    }
}



/*
 * calculates a^(-1) mod p
 */
inline int modInv(int a, int p){
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

inline void moduloPoly(int *p1, int m1, int *mod, int m, int charac){
    int deg1=0, degmod=0;
    int i=0,j=0;
    long quo=0;
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
    /*printf("poly = ");printArr(p1,m1);*/
    /*printf("mod = ");printArr(mod,m);*/
    /*printf("degMod=%i, degP1=%i\n", degmod,deg1);*/
    
    /*printArr(p1,m1);*/
    /*int *tmpArr = malloc(m1*sizeof(int));*/

    //make polynomial division
    int degmodInv = modInv(mod[degmod],charac);
    for(i=deg1-degmod; i>=0; i--){
        quo = (p1[i+degmod]*(long)degmodInv)%charac;
        /*printf("i=%i p1[i+degmod]=%i mod[degmod]=%i mod[degmod]^(-1)=%i ",i,p1[i+degmod],mod[degmod],modInv(mod[degmod],charac));*/
        /*printf("quo=%lu\n", quo);*/
        for(j=degmod;j>=0;j--){
            p1[i+j] = (int)((p1[i+j] - mod[j]*quo)%charac);
        }
        /*for(j=0;j<m1;j++){*/
            /*tmpArr[j] = p1[j];*/
            /*if(tmpArr[j] < 0) tmpArr[j] += charac;*/
        /*}*/

        /*printf("=> p1=");printArr(tmpArr,m1);*/
    }
    /*for(i=0;i<m1;i++){*/
        /*if(p1[i] < 0) p1[i] += charac;*/
    /*}*/
    /*free(tmpArr);*/
}

inline void moduloMonom(int *p1, int m1, int *mod, int m, int charac){
    int deg1, degmod=m-1;
    int quo;
    int i,j;
    //get degrees
    for(i=m1-1;i>=0;i--){
        if( p1[i] != 0){
            deg1 = i;
            break;
        }
    }
    //make polynomial division
    for(i=deg1-degmod; i>=0; i--){
        quo = p1[i+degmod]%charac;
        for(j=degmod;j>=0;j--){
            p1[i+j] = (p1[i+j] - mod[j]*quo)%charac;
        }
    }
}

inline void moduloMonomShort(int *p1, int m1, int *mod, int m, int charac){
    int deg1=0, degmod=m-1;
    int i=0,j=0;
    int quo=0;
    //get degrees
    for(i=m1-1;i>=0;i--){
        if( p1[i] != 0){
            deg1 = i;
            break;
        }
    }
    /*for(i=m-1;i>=0;i--){*/
        /*if( mod[i] != 0){*/
            /*degmod = i;*/
            /*break;*/
        /*}*/
    /*}*/

    //make polynomial division
    for(i=deg1-degmod; i>=0; i--){
        quo = p1[i+degmod];
        for(j=degmod;j>=0;j--){
            p1[i+j] -= mod[j]*quo;
        }
    }
}



/**
 * Square and multiply
 * x is modified!
 */
inline void powerPolyInt(int *x, int *x_mipo, int *ret, int m, 
        int power, int charac, int *tmp2){
    int i,j;
    initPoly(ret,m);
    ret[0] = 1;
    while(power > 0){
        if(power%2 == 1){
            multiplyPoly(ret,m, x,m, tmp2,2*m, charac);
            moduloPoly(tmp2,2*m,x_mipo,m+1,charac);
            for(i=0;i<m;i++)
                ret[i] = tmp2[i];
        }
        multiplyPoly(x,m, x,m, tmp2,2*m, charac);
        /*printf("powerPolyInt tmp2=");printArr(tmp2,2*m);*/
        moduloPoly(tmp2,2*m,x_mipo,m+1,charac);
        for(i=0;i<m;i++)
            x[i] = tmp2[i];
        power >>= 1;
    }
}
/**
 * Square and multiply
 * x is modified!
 */
inline void powerPoly(int *x, int *x_mipo, int *ret, int m, 
        char *power, int powerLen,
        int charac, int *tmp2){
    int i,j;
    initPoly(ret,m);
    ret[0] = 1;
    for(j=powerLen-1;j>=0;j--){
        if(power[j] == 1){
            multiplyPoly(ret,m, x,m, tmp2,2*m, charac);
            moduloMonom(tmp2,2*m,x_mipo,m+1,charac);
            for(i=0;i<m;i++)
                ret[i] = tmp2[i];
        }
        multiplyPoly(x,m, x,m, tmp2,2*m, charac);
        moduloMonom(tmp2,2*m,x_mipo,m+1,charac);
        for(i=0;i<m;i++)
            x[i] = tmp2[i];
    }
}


/**
 * Square and multiply
 * x is modified!
 */
inline void powerPolyShort(int *x, int *x_mipo, int *ret, int m, 
        char *power, int powerLen,
        int charac, int *tmp2){
    int i,j;
    initPoly(ret,m);
    ret[0] = 1;
    for(j=powerLen-1;j>=0;j--){
        if(power[j] == 1){
            multiplyPolyShort(ret,m, x,m, tmp2,2*m, charac);
            moduloMonomShort(tmp2,2*m,x_mipo,m+1,charac);
            for(i=0;i<m;i++)
                ret[i] = tmp2[i]%charac;
        }
        multiplyPolyShort(x,m, x,m, tmp2,2*m, charac);
        moduloMonomShort(tmp2,2*m,x_mipo,m+1,charac);
        for(i=0;i<m;i++)
            x[i] = tmp2[i]%charac;
    }
}

inline bool isOne(int *x, int m, int charac){
    int i=0;
    bool allZ = true;
    if(!(x[0] == 1 || x[0] == -charac+1))
        return false;
    for(i=1;i<m;i++){
        if(x[i] != 0){
            allZ = false;
            break;
        }
    }
    return allZ;
}

/**
 * multiplies 2 quadratic m x m matrices 
 */
void multMatrices(int *mat1, int *mat2, int *ret, int m, int charac){
    int i,j,k;
    long tmp;
    for(i=0;i<m;i++){
        for(j=0;j<m;j++){
            ret[i*m + j] = 0;
            for(k=0;k<m;k++){
                tmp = (long)mat1[i*m+k]*(long)mat2[k*m+j];
                ret[i*m + j] = (int)((ret[i*m + j] + tmp)%charac);
            }
        }
    }
}


void genMats(int *mipo, int m, int *mats, int maxPower, int charac, int q){
    int *tmp = malloc(m*sizeof(int));
    int *tmp2 = malloc(2*m*sizeof(int));
    int *x = malloc(m*sizeof(int));

    int i,j,k;
    int mSize = m*m;
    /*//put identity matrix at position 0*/
    for(i=0;i<m;i++)
        for(j=0;j<m;j++)
            if(j==i)
                /*matsInt[0][i][j] = 1;*/
                mats[0*mSize + i*m + j] = 1;
            else
                /*matsInt[0][i][j] = 0;*/
                mats[0*mSize + i*m + j] = 0;

    for(i=0;i<m;i++){
        initPoly(x,m);
        x[i] = 1;
        powerPolyInt(x, mipo, tmp, m, q, charac, tmp2);
        for(j=0;j<m;j++)
            /*matsInt[1][j][i] = tmp[j];*/
            mats[1*mSize + j*m + i] = tmp[j];
    }
    for(i=2;i<=maxPower;i++){
        /*multMatrices(matsInt[1],matsInt[i-1],matsInt[i], m, charac);*/
        multMatrices(mats+1*mSize, mats+(i-1)*mSize, mats+i*mSize, m, charac);
    }

    /*for(i=0;i<=maxPower;i++){*/
        /*printf("m^%i\n",i);*/
        /*for(j=0;j<m;j++)*/
            /*printArr(mats+mSize*i+j*m,m);*/
    /*}*/

    free(x);
    free(tmp2);
    free(tmp);
}



/**
 * Calc order of element
 * x is NOT modified!
 */
inline bool isPrimitive(int *x, int *x_mipo, int m, 
        char *barFactors, int* lenBarFactors, int countBarFactors, 
        int charac,
        int *tmp_x, int *tmp, int *firstPow, int *tmpRes, int *tmp2){
    int i,j, curPos;
    /*printf(" test for primitivity x=");printArr(x,m);*/
    //init firstPow
    for(j=0;j<m;j++) tmp_x[j] = x[j];
    powerPolyShort(tmp_x,x_mipo,tmp,m, 
            barFactors,
            lenBarFactors[0], charac, tmp2);
    if( isOne(tmp,m, charac) == true)
        return false;
    
    for(j=0;j<m;j++) firstPow[j] = tmp[j];
    curPos = lenBarFactors[0];

    for(i=1;i<countBarFactors;i++){
        // calc x^...
        for(j=0;j<m;j++)
            tmp_x[j] = x[j];
        powerPolyShort(tmp_x,x_mipo,tmp,m, 
                barFactors+curPos+lenBarFactors[2*i-1],
                lenBarFactors[2*i], charac, tmp2);
        // calc firstPow^...
        for(j=0;j<m;j++) tmp_x[j] = firstPow[j];
        powerPolyShort(tmp_x,x_mipo,tmpRes,m, 
                barFactors+curPos,
                lenBarFactors[2*i-1], charac, tmp2);

        multiplyPolyShort(tmpRes,m, tmp,m, tmp2, 2*m, charac);
        moduloMonomShort(tmp2,2*m, x_mipo,m+1,charac);
        for(j=0;j<m;j++) tmp2[j] %= charac;
        if( isOne(tmp2,m, charac) == true)
            return false;

        curPos += lenBarFactors[2*i-1]+lenBarFactors[2*i];
    }
    return true;
}


inline bool allZero(int *x, int m){
    int i=0;
    bool allZ = true;
    for(i=0;i<m;i++){
        if(x[i] != 0){
            allZ = false;
            break;
        }
    }
    return allZ;
}



/*
 * calculates g(sigma)(x) where g is a polynomial and sigma the frobenius
 * application of frobenius is given by mats, i.e.
 * mats is array of m x m matrices, x an array of m ints
 * g is an array of glen arrays of length m
 */
inline void applyFrob(int *x, int *x_mipo, int *g, int glen, int *mats, int frobpower,
        int *ret, int m, int charac, int *tmp, int *tmp2){
    //printf("x="); printArr(x,m);
    //printf("x_mipo="); printArr(x_mipo,m+1);
    //printf("g="); printArr(g,m*glen);
    int mSize = m*m;
    int i,j;
    initPoly(ret,m);
    for(i=0;i<glen;i++){
        if(allZero(g+i*m,m) == true){
            //printf("i=%i allZero = gi=",i); 
            //printArr(g+i*m,m);
            continue;
        }
        //printf("i=%i",i);
        matmul(mats+i*frobpower*mSize, x, tmp, m, charac);
        //printf(" mat*x = "); printArr(tmp, m);
        //printf("tmp*gi for \n\ttmp="); printArr(tmp,m);
        //printf("\tgi="); printArr(g+i*glen, m);
        multiplyPoly(tmp,m, g+i*m, m, tmp2, 2*m, charac);
        //printf("\t=>tmp*gi = "); printArr(tmp2,2*m);
        //printf("tmp2 mod mipo for \n\ttmp2="); printArr(tmp2,2*m); 
        //printf("\tmipo="); printArr(x_mipo,m+1);
        moduloPoly(tmp2, 2*m, x_mipo, m+1, charac);
        //printf("\t=>tmp2 mod mipo = "); printArr(tmp2,m);
        for(j=0;j<m;j++){
            ret[j] += tmp2[j];
        }
    }
    for(i=0;i<m;i++){
        ret[i] %= charac;
    }
}

/*
 * calculates g(sigma)(x) where g is a polynomial and sigma the frobenius
 * application of frobenius is given by mats, i.e.
 * mats is array of m x m matrices, x an array of m ints
 * g is an array of glen arrays of length m
 */
inline void applyFrobShort(int *x, int *x_mipo, 
        int *g, int glen, int *gCoeffDegs,
        int *mats, int frobpower,
        int *ret, int m, int charac, int *tmp, int *tmp2,
        int *matmulCache, bool *matmulCacheCalced){
    int mSize = m*m;
    int i,j,k,l;
    initPoly(ret,m);
    for(i=0;i<glen;i++){
        if(allZero(g+i*m,m) == true){
            continue;
        }
        if(matmulCacheCalced[i*frobpower] == 1){
            multiplyPolyShort(matmulCache+m*(i*frobpower),m,
                    g+i*m, m, tmp2, 2*m, charac);
        }else{
            matmulShort(mats+i*frobpower*mSize, x, tmp, m, charac);
            for(j=0;j<m;j++)
                matmulCache[m*(i*frobpower)+j] = tmp[j];
            matmulCacheCalced[i*frobpower] = 1;
            multiplyPolyShortKnownDeg(tmp,m, g+i*m, m, gCoeffDegs[i],
                    tmp2, 2*m, charac);
            /*multiplyPolyShort(tmp,m, g+i*m, m, */
                    /*tmp2, 2*m, charac);*/
        }
        moduloMonomShort(tmp2, 2*m, x_mipo, m+1, charac);
        for(j=0;j<m;j++){
            ret[j] += tmp2[j];
        }
    }
    for(i=0;i<m;i++){
        ret[i] %= charac;
    }
}

inline void testPolys(int *x, int *x_mipo, int decompCount,
        int *polys, int *polysLen, int *polysCount, bool *evalToZero,
        int *mats, int *frobPowers, 
        int *ret, int m, int charac, int *tmp, int *tmp2){
    int i,j;
    int curDecompPosition = 0;
    int curPolyPosition = 0;
    int lastZeroPoly = 0;
    int goodCounter = 0;
    /*printf("testPolys for x="), printArr(x,m);*/
    /*printf("polysCount="); printArr(polysCount,decompCount);*/

    for(i=0;i<decompCount;i++){
        /*printf("decomp: i=%i\n",i);*/
        goodCounter = 0;
        for(j=0;j<polysCount[i];j++){
            /*printf("\ttest poly j=%i",j); */
            applyFrob(x, x_mipo, 
                    polys+curPolyPosition, polysLen[curDecompPosition+j],
                    mats, frobPowers[curDecompPosition+j], 
                    ret, m, charac, tmp, tmp2);
            /*printf("\t=>ret=");printArr(ret,m);*/
            if( allZero(ret,m) == evalToZero[curDecompPosition+j] ){
                /*printf(" good\n");*/
                goodCounter += 1;
            }else{
                /*printf(" not good\n");*/
            }
            curPolyPosition += m*polysLen[curDecompPosition+j];
        }
        /*printf("\t all %i tested => goodCounter=%i\n", polysCount[i],goodCounter);*/
        if(goodCounter == polysCount[i]){
            ret[0] = i;
            return;
        }
        curDecompPosition += polysCount[i];
    }
    
    ret[0] = -1;
}

inline void testPolysShort(int *x, int *x_mipo, int decompCount,
        int *polys, int *polysLen, int *polysCoeffDegs, int *polysCount, 
        bool *evalToZero, int *mats, int *frobPowers, 
        int *ret, int m, int charac, int *tmp, int *tmp2,
        int *matmulCache, bool *matmulCacheCalced){
    int i,j;
    int curDecompPosition = 0;
    int curPolyPosition = 0;
    int curPolyCoeffPosition = 0;
    int lastZeroPoly = 0;
    int goodCounter = 0;

    for(i=0;i<decompCount;i++){
        goodCounter = 0;
        for(j=0;j<polysCount[i];j++){
            applyFrobShort(x, x_mipo, 
                    polys+curPolyPosition, polysLen[curDecompPosition+j],
                    polysCoeffDegs+curPolyCoeffPosition,
                    mats, frobPowers[curDecompPosition+j],
                    ret, m, charac, tmp, tmp2,
                    matmulCache, matmulCacheCalced);
            if( allZero(ret,m) == evalToZero[curDecompPosition+j] ){
                goodCounter += 1;
            }
            curPolyPosition += m*polysLen[curDecompPosition+j];
            curPolyCoeffPosition += polysLen[curDecompPosition+j];
        }
        if(goodCounter == polysCount[i]){
            ret[0] = i;
            return;
        }
        curDecompPosition += polysCount[i];
    }
    
    ret[0] = -1;
}

unsigned long long encodeArr(int *x, int m, int shiftSize){
    int i;
    unsigned long long ret = x[m-1];
    for(i=m-2;i>=0;i--){
        ret <<= shiftSize;
        ret += x[i];
    }
    return ret;
}


void decodeArr(unsigned long long enc, int *arr, int m, int shiftSize){
    int i;
    unsigned long long andOp = 1;
    for(i=1;i<shiftSize;i++){
        andOp <<= 1;
        andOp += 1;
    }
    /*printf("decodeArr: andOp=%i\n",andOp);*/
    for(i=0;i<m;i++){
        arr[i] = enc & andOp;
        enc >>= shiftSize;
    }
}

void decodeArrAdd(unsigned long long enc, int *arr, int m, int shiftSize){
    int i;
    unsigned long long andOp = 1;
    for(i=1;i<shiftSize;i++){
        andOp <<= 1;
        andOp += 1;
    }
    for(i=0;i<m;i++){
        arr[i] += (enc & andOp);
        enc >>= shiftSize;
    }
}

//make a linked list
struct Node {
    unsigned long long  x;
    struct Node * next;
};

/**
 * appends element to end of root, where element is copied to new array
 */
inline struct Node *appendToEnd(struct Node *root, int * element, int elLen, \
        int shiftSize){
    struct Node *nextNode = root;
    if( nextNode != 0){
        while(nextNode->next != 0){
            nextNode = nextNode->next;
        }
        if( nextNode->x != 0){
            nextNode->next = malloc( sizeof(struct Node) );
            nextNode = nextNode->next;
        }
        if( nextNode != 0){
            nextNode->next = 0;
            /*nextNode->x = malloc(elLen*sizeof(int));*/
            /*int i;*/
            /*for(i=0;i<elLen;i++){*/
                /*nextNode->x[i] = element[i];*/
            /*} */
            nextNode->x = encodeArr(element,elLen,shiftSize);
            return nextNode;
        }
    }
    return NULL;
}


inline void freeNode(struct Node* head){
    struct Node *next_n = NULL;
    struct Node *tmp_n = NULL;
    for(tmp_n=head; tmp_n !=NULL; ){
        next_n = tmp_n->next;
        /*printf("free node: %i ",n);*/
        /*free(tmp_n->x);*/
        /*printf(" data freed");*/
        free(tmp_n);
        /*printf(" tmp freed ");*/
        tmp_n = next_n;
        /*printf(" next=%i\n",n);*/
    }
    head = 0;
}

inline struct Node *appendNode(struct Node *curNode){
    curNode->next = malloc(sizeof(struct Node));
    curNode = curNode->next;
    curNode->x = 0;
    curNode->next = 0;
    return curNode;
}

unsigned long long processFFElements( int *x_mipo, int decompCount,
        int *polys, int *polysLen, int *polysCoeffDegs, int *polysCount, 
        bool *evalToZero, int *mats, int matLen, int *frobPowers, 
        int *genCounts, int m, int charac, int shiftSize,
        char *barFactors, int *lenBarFactors, int countBarFactors){
    time_t TIME = time(NULL);
    int i,j;
    
    int * x = malloc( m*sizeof(int) );
    int * ret = malloc( m*sizeof(int) );
    int * tmp = malloc( m*sizeof(int) );
    int * tmp_x = malloc( m*sizeof(int) );
    int * firstPow = malloc( m*sizeof(int) );
    int * tmpRes = malloc( m*sizeof(int) );
    int * tmp2 = malloc( 2*m*sizeof(int) );
    int * matmulCache = malloc(matLen*m*sizeof(int));
    bool * matmulCacheCalced = malloc(matLen*sizeof(bool));


    unsigned long long pcn = 0;

    initPoly(x,m);
    initPoly(genCounts,decompCount);

    struct Node **roots = malloc( decompCount*sizeof(struct Node) );
    struct Node **curRoots = malloc(decompCount*sizeof(struct Node));
    for(i=0;i<decompCount;i++){
        roots[i] = malloc( sizeof(struct Node) );
        roots[i]->x = 0;
        roots[i]->next = 0;
        curRoots[i] = roots[i];
    }

    
    while(1==1){
        for(i=0;i<matLen;i++) matmulCacheCalced[i] = 0;
        testPolysShort(x,x_mipo,decompCount,
                polys,polysLen,polysCoeffDegs, polysCount,evalToZero,
                mats, frobPowers,
                ret, m, charac, tmp, tmp2,
                matmulCache,matmulCacheCalced);
        /*for(i=0;i<matLen;i++)*/
            /*printf("%i ",matmulCacheCalced[i]);*/
        /*printf("\n");*/
        if( ret[0] != -1){
            genCounts[ret[0]] += 1;
            curRoots[ret[0]] = appendToEnd(curRoots[ret[0]], x, m, shiftSize);
            decodeArr(curRoots[ret[0]]->x,x,m,shiftSize);
        }
        //generate next element
        x[0] += 1;
        if( x[0] == charac ){
            for(i=0;i<m-1 && x[i]==charac;i++){
                x[i] = 0;
                x[i+1] += 1;
            }
            if( x[m-1]==charac ){
                break;
            }
        }
    }
    printf("CN time: %.2f\n", (double)(time(NULL)-TIME));

    for(i=0;i<decompCount;i++){
        curRoots[i] = roots[i];
    }
    while(1==1){
        decodeArr(curRoots[0]->x, x, m, shiftSize);
        for(i=1;i<decompCount;i++){
            decodeArrAdd(curRoots[i]->x, x, m, shiftSize);
        }
        for(i=0;i<m;i++)
            x[i] %= charac;
        moduloPoly(x,m,x_mipo,m+1,charac);

        //test primitivity
        if(isPrimitive(x,x_mipo,m,barFactors,lenBarFactors, countBarFactors,
                    charac,tmp_x,tmp,firstPow,tmpRes,tmp2) == true){
            pcn ++;
        }

        //nextElement
        curRoots[0] = curRoots[0]->next;
        if( curRoots[0] == 0 ){
            for(i=0;i<decompCount-1 && curRoots[i]==0;i++){
                curRoots[i] = roots[i];
                curRoots[i+1] = curRoots[i+1]->next;
            }
        }
        if( curRoots[decompCount-1] == 0){
            break;
        }
    }

    for(i=0;i<decompCount;i++){
        freeNode(roots[i]);
    }
    free(roots);
    free(curRoots);
    free(tmp2);
    free(tmp);
    free(tmpRes);
    free(firstPow);
    free(tmp_x);
    free(ret);
    free(matmulCache);
    free(matmulCacheCalced);
    
    printf("total time: %.2f\n", (double)(time(NULL)-TIME));
    return pcn;
}


double eta_processFFElements( int *x_mipo, int decompCount,
        int *polys, int *polysLen, int *polysCoeffDegs, int *polysCount, bool *evalToZero,
        int *mats, int matLen, int *frobPowers, 
        int *genCounts, int m, int charac, int shiftSize){
    struct timeval TIME1, TIME2;
    int i,j;
    
    
    int * x = malloc( m*sizeof(int) );
    int * ret = malloc( m*sizeof(int) );
    int * tmp = malloc( m*sizeof(int) );
    int * tmp2 = malloc( 2*m*sizeof(int) );
    int * matmulCache = malloc(matLen*m*sizeof(int));
    bool * matmulCacheCalced = malloc(matLen*sizeof(bool));

    initPoly(x,m);
    initPoly(genCounts,decompCount);

    struct Node **roots = malloc( decompCount*sizeof(struct Node) );
    struct Node **curRoots = malloc(decompCount*sizeof(struct Node));
    for(i=0;i<decompCount;i++){
        roots[i] = malloc( sizeof(struct Node) );
        roots[i]->x = 0;
        roots[i]->next = 0;
        curRoots[i] = roots[i];
    }

    int counter = 0;
    gettimeofday(&TIME1,NULL);
    
    while(1==1){
        for(i=0;i<matLen;i++) matmulCacheCalced[i] = 0;
        testPolysShort(x,x_mipo,decompCount,
                polys,polysLen,polysCoeffDegs, polysCount,evalToZero,
                mats, frobPowers,
                ret, m, charac, tmp, tmp2,
                matmulCache, matmulCacheCalced);
        if( ret[0] != -1){
            genCounts[ret[0]] += 1;
            curRoots[ret[0]] = appendToEnd(curRoots[ret[0]], x, m, shiftSize);
            decodeArr(curRoots[ret[0]]->x,x,m,shiftSize);
        }
        //generate next element
        x[0] += 1;
        if( x[0] == charac ){
            for(i=0;i<m-1 && x[i]==charac;i++){
                x[i] = 0;
                x[i+1] += 1;
            }
            if( x[m-1]==charac ){
                break;
            }
        }
        counter++;
        if(counter == 1000)
            break;
    }

    gettimeofday(&TIME2,NULL);
    
    for(i=0;i<decompCount;i++){
        /*printf("free root %i\n",i);*/
        freeNode(roots[i]);
    }
    free(roots);
    /*printf("roots freed!\n");*/
    free(curRoots);
    /*printf("curRoots freed!\n");*/
    free(tmp2);
    /*printf("tmp2 freed!\n");*/
    free(tmp);
    /*printf("tmp freed!\n");*/
    free(ret);
    /*printf("ret freed!\n");*/
    free(x);
    /*printf("x freed!\n");*/
    free(matmulCache);
    free(matmulCacheCalced);
    
    double timediff = (TIME2.tv_sec - TIME1.tv_sec +
         ((double)(TIME2.tv_usec - TIME1.tv_usec))/1000000.0);
    return 4*timediff *pow((double)charac,(double)m)/counter;
}


/**
 * returns the next CN element
 */
void findAnyPCN(int * x, int *x_mipo, 
        int *polys, int *polysLen, int *polysCount, bool *evalToZero,
        int *mats, int *frobPowers, 
        int m, int charac, 
        char *barFactors, int *lenBarFactors, int countBarFactors){
        /*int *ret, int *tmp, int * tmp2){*/
    time_t TIME = time(NULL);
    int i,j;
    
    
    int * ret = malloc( m*sizeof(int) );
    int * tmp = malloc( m*sizeof(int) );
    int * tmp_x = malloc( m*sizeof(int) );
    int * tmpRes = malloc( m*sizeof(int) );
    int * firstPow = malloc( m*sizeof(int) );
    int * tmp2 = malloc( 2*m*sizeof(int) );

    bool pcnFound = false;
    
    while(1==1){
        //generate next element
        x[0] += 1;
        if( x[0] == charac ){
            for(i=0;i<m-1 && x[i]==charac;i++){
                x[i] = 0;
                x[i+1] += 1;
            }
            if( x[m-1]==charac ){
                x[0] = -1;
                return;
            }
        }
        testPolys(x,x_mipo,1,
                polys,polysLen,polysCount,evalToZero,
                mats, frobPowers,
                ret, m, charac, tmp, tmp2);
        if( ret[0] != -1){
            if(isPrimitive(x,x_mipo,m,barFactors,lenBarFactors,countBarFactors,
                        charac, tmp_x,tmp,firstPow,tmpRes,tmp2) == true){
                pcnFound = true;
                break;
            }
        }
    }
    free(ret);
    free(tmp);
    free(tmpRes);
    free(firstPow);
    free(tmp_x);
    free(tmp2);

    if(pcnFound == false){
        x[0] = -1;
    }

    printf("C time: %.2f\n", (double)(time(NULL)-TIME));
}


bool findAnyPCN_useGen(int * x, int * generator, int *x_mipo, 
        int *polys, int *polysLen, int *polysCount, bool *evalToZero,
        int *mats, int *frobPowers, 
        int m, int charac, 
        int *powerTable, int lenPowerTable){
    time_t TIME = time(NULL);
    int i,j;
    
    
    int * ret = malloc( m*sizeof(int) );
    int * tmp = malloc( m*sizeof(int) );
    int * tmp_x = malloc( m*sizeof(int) );
    int * tmp2 = malloc( 2*m*sizeof(int) );

    bool pcnFound = false;

    // the generator
    for(j=0;j<m;j++) x[j] = generator[j];
    
    for(i=0;i<lenPowerTable;i++){
        if(powerTable[i] != 0){
            for(j=0;j<m;j++) tmp_x[j] = generator[j];
            powerPolyInt(tmp_x, x_mipo, ret, m, powerTable[i], charac, tmp2);

            multiplyPoly(x,m, ret,m, tmp2, 2*m, charac);
            moduloPoly(tmp2, 2*m, x_mipo, m+1, charac);
            for(j=0;j<m;j++) x[j] = tmp2[j];
        }

        testPolys(x,x_mipo,1,
                polys,polysLen,polysCount,evalToZero,
                mats, frobPowers,
                ret, m, charac, tmp, tmp2);
        if( ret[0] != -1){
            for(j=0;j<m;j++){
                if( x[j]<0 ) x[j] += charac;
            }
            pcnFound = true;
            break;
        }
    }
    free(ret);
    free(tmp);
    free(tmp_x);
    free(tmp2);

    printf("C time: %.2f\n", (double)(time(NULL)-TIME));
    return pcnFound;
}
int main(){
    int mats[] = 
        {1,0,0,0,0,0,0,0,0,0,
        0,1,0,0,0,0,0,0,0,0,
        0,0,1,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,0,0,0,
        0,0,0,0,1,0,0,0,0,0,
        0,0,0,0,0,1,0,0,0,0,
        0,0,0,0,0,0,1,0,0,0,
        0,0,0,0,0,0,0,1,0,0,
        0,0,0,0,0,0,0,0,1,0,
        0,0,0,0,0,0,0,0,0,1,
        1,0,0,0,0,1,0,1,0,0,
        0,0,0,0,0,1,0,1,1,1,
        0,1,0,0,0,1,1,1,0,1,
        0,0,0,0,0,1,1,1,0,0,
        0,0,1,0,0,0,1,1,0,1,
        0,0,0,0,0,1,1,0,1,0,
        0,0,0,1,0,1,0,0,0,1,
        0,0,0,0,0,0,1,1,1,0,
        0,0,0,0,1,0,1,0,0,0,
        0,0,0,0,0,0,0,1,1,1,
        1,0,0,0,0,0,0,0,0,0,
        0,0,0,0,1,1,1,0,1,1,
        0,0,0,1,0,1,0,1,0,1,
        0,0,0,1,0,0,0,1,0,1,
        0,1,0,1,0,0,0,1,0,1,
        0,0,0,1,1,0,0,0,1,1,
        0,0,0,0,0,0,0,0,0,1,
        0,0,0,1,1,1,0,1,1,1,
        0,0,1,1,0,1,1,1,0,0,
        0,0,0,0,1,0,0,0,0,1,
        1,0,0,0,0,1,0,1,0,0,
        0,0,1,1,1,0,1,0,0,1,
        0,0,0,0,0,0,1,1,1,1,
        0,0,0,0,0,1,0,1,0,1,
        0,0,0,0,0,0,0,0,1,0,
        0,0,1,0,1,1,1,1,1,0,
        0,0,0,0,0,0,0,1,1,1,
        0,0,1,0,1,0,1,0,1,0,
        0,1,0,1,0,0,0,1,0,0,
        0,0,1,0,0,0,1,0,1,0,
        1,0,0,0,0,0,0,0,0,0,
        0,1,1,1,0,1,1,0,1,0,
        0,0,0,1,1,1,0,0,0,0,
        0,0,0,0,0,1,0,0,1,1,
        0,0,0,0,1,0,1,0,0,0,
        0,1,1,1,1,1,1,1,0,1,
        0,0,0,0,1,0,0,0,0,1,
        0,1,1,1,1,0,1,0,0,1,
        0,0,0,0,0,0,0,1,0,1,
        0,1,0,1,1,0,0,1,0,0,
        1,0,0,0,0,1,0,1,0,0,
        0,1,0,1,1,1,0,1,0,1,
        0,0,1,0,0,0,1,0,1,1,
        0,0,0,0,1,1,0,1,0,1,
        0,0,1,1,0,1,1,1,0,0,
        0,1,1,1,0,1,1,0,0,1,
        0,0,1,0,0,0,1,0,1,0,
        0,1,1,1,0,0,1,1,0,1,
        0,0,0,0,0,0,1,0,0,1,
        0,0,1,0,0,0,1,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        0,0,1,0,0,1,0,1,0,1,
        0,1,0,1,1,0,0,0,1,1,
        0,0,1,0,0,1,1,1,1,0,
        0,1,0,1,0,0,0,1,0,0,
        0,1,0,1,0,1,1,0,1,0,
        0,1,0,1,1,0,0,1,0,0,
        0,1,0,1,0,0,1,1,1,0,
        0,0,0,1,0,1,0,1,1,0,
        0,1,0,1,0,0,1,1,0,0,
        1,0,0,0,0,1,0,1,0,0,
        0,1,0,0,0,0,1,1,1,0,
        0,0,1,0,1,0,1,0,0,1,
        0,1,0,1,1,1,0,0,0,0,
        0,0,0,0,0,0,0,1,0,1,
        0,0,0,1,1,0,1,0,0,0,
        0,0,1,0,0,0,1,0,0,0,
        0,0,0,1,1,1,1,1,0,0,
        0,0,0,0,1,0,0,0,0,0,
        0,0,0,1,0,1,0,1,0,0,
        1,0,0,0,0,0,0,0,0,0,
        0,0,0,1,1,0,0,0,0,0,
        0,1,1,1,0,0,0,1,1,0,
        0,0,1,0,0,1,1,1,0,0,
        0,0,0,0,0,0,1,0,0,1,
        0,0,1,1,0,0,0,0,0,0,
        0,1,0,1,0,0,1,1,0,0,
        0,0,1,1,0,1,0,1,0,0,
        0,0,1,0,0,0,1,1,0,1,
        0,0,0,0,0,0,1,0,0,0};
    int polys[] =  
        {1,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,1,0,1,0,1,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0,
        0,1,0,1,0,1,0,0,0,0,
        1,0,0,0,0,0,0,0,0,0};
    int polysLen[] = {3, 2, 2, 1, 9, 5, 5, 3, 3};
    bool evalToZero[] = {1, 0, 1, 0, 1, 0, 1, 0, 0};
    int frobPowers[] = {1, 1, 2, 2, 1, 1, 2, 2, 2};
    int polysCount[] = {4, 5};
    int xmipo[] = {1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1};
    int decompCount = 2;
    int m = 10;
    int charac = 2;
    int shiftSize = 1;
    char barFactors[] = {1, 0, 1, 0, 1, 0, 1, 0, 1, 
                         1, 0, 1, 1, 1, 0, 1,
                         1, 0, 0, 0, 0, 1}; //{341, 93, 33};
    int lenBarFactors[] = {9,7,6};
    int countBarFactors = 3;

    ////////////////////////////////////////////////////////////
    /*int *genCounts = malloc(decompCount*sizeof(int));*/

    /*unsigned long long pcn*/
        /*= processFFElements(xmipo, decompCount,*/
            /*polys, polysLen, polysCount, evalToZero,*/
            /*mats, frobPowers,*/
            /*genCounts, m, charac, shiftSize,*/
            /*barFactors,lenBarFactors,countBarFactors);*/
    /*free(genCounts);*/
    /*printf("pcn=%i",pcn);*/
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    /*int n = 10;*/
    /*int *mats2 = malloc((n+1)*m*m*sizeof(int));*/
    /*genMats(xmipo,m,mats2,n,charac,charac);*/
    /*free(mats2);*/
}
