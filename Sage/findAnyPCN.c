#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include "algs-base.h"

///////////////////////////////////////////////////////////////////////////////
// Setup FFElements ///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

struct FFElem{
    int *el;
    int *idcs;
    int len;
};

inline struct FFElem *mallocFFElem(int m){
    struct FFElem *ff = malloc(sizeof(struct FFElem));
    ff->el = malloc(m*sizeof(int));
    ff->idcs = malloc(m*sizeof(int));
    ff->len = 0;
    return ff;
}

inline void freeFFElem(struct FFElem *ff){
    if(ff==0) return;
    free(ff->el);
    free(ff->idcs);
    free(ff);
}

inline void freeFFElemMatrix(struct FFElem **mat, int len){
    if(mat==0) return;
    int i;
    for(i=0;i<len;i++) freeFFElem(mat[i]);
    free(mat);
}

/**
 * Copies the content of ff1 into ff2
 *
 * !! ff2 must be malloced!
 */
inline void copyFFElem(struct FFElem *ff1, struct FFElem *ff2){
    if(ff1 == ff2) return;
    int i;
    for(i=0;i < ff1->len;i++){
        ff2->idcs[i] = ff1->idcs[i];
        ff2->el[ ff1->idcs[i] ] = ff1->el[ ff1->idcs[i] ];
    }
    ff2->len = ff1->len;
}


void printFFElem(char *preName, struct FFElem *ff){
    printf("%s = ",preName);
    printf("\t{ el=");if(ff->len>0) printArr(ff->el,ff->idcs[0]+1);
    printf("\t  idcs=");printArr(ff->idcs,ff->len);
    printf("\t  len=%i}\n",ff->len);
}

void printFFElemShort(char *preName, struct FFElem *ff){
    printf("%s = ",preName);
    if(ff->len == 0){
        printf("0\n");
        return;
    }
    int i=0;
    int j=ff->len-1;
    while(j>=0){
        while(i<ff->idcs[j]){
            printf("0 ");
            i++;
        }
        printf("%i ",ff->el[ff->idcs[j]]);
        j--;
        i++;
    }
    printf("\n");
}

void printFFElemMatrix(struct FFElem **mat, int m){
    int row;
    for(row=0;row<m;row++){
        int i=0;
        int j=mat[row]->len-1;
        while(j>=0){
            while(i<mat[row]->idcs[j]){
                printf("0 ");
                i++;
            }
            printf("%i ",mat[row]->el[mat[row]->idcs[j]]);
            j--;
            i++;
        }
        while(i<m){
            printf("0 ");
            i++;
        }
        printf("\n");
    }
}

inline bool isOne(struct FFElem *ff){
    return (ff->len == 1 && ff->idcs[0] == 0 && ff->el[0] == 1);
}

inline bool isZero(struct FFElem *ff){
    return (ff->len == 0);
}

/**
 * Updates FFElem according to its ->el
 *
 * !! ->el must not have uninitialized or wrong values !!
 */
inline void updateFFElem(struct FFElem *ff,int m){
    int i,i2;
    i2 = 0;
    for(i=m-1;i>=0;i--){
        if(ff->el[i] != 0){
            ff->idcs[i2] = i;
            i2++;
        }
    }
    ff->len = i2;
}


///////////////////////////////////////////////////////////////////////////////
// Setup Polynomials //////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
struct FFPoly{
    struct FFElem **poly;
    int lenPoly;
};

inline struct FFPoly *mallocFFPoly(int m, int lenPoly){
    struct FFPoly *poly = malloc(lenPoly*sizeof(struct FFElem*));
    poly->lenPoly = lenPoly;
    int i;
    for(i=0;i<lenPoly;i++) poly->poly[i] = mallocFFElem(m);
    return poly;
}

inline void freeFFPoly(struct FFPoly *poly){
    int i;
    for(i=0;i<poly->lenPoly;i++) freeFFElem(poly->poly[i]);
    free(poly->poly);
    free(poly);
}

void printFFPoly(char *preName, struct FFPoly *poly){
    printf("%s: ",preName);
    int i;
    for(i=0;i<poly->lenPoly;i++)
        printFFElemShort("   ",poly->poly[i]);
    printf(" len=%i\n",poly->lenPoly);
}


///////////////////////////////////////////////////////////////////////////////
// Algorithms on FFElems //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/**
 * Adds two FFElems.
 *
 * !! ff1 may be same as ret !!
 * !! ff2 must not be same as ret !!
 */
inline void addFFElem(struct FFElem *ff1, struct FFElem *ff2,
        struct FFElem *ret,
        int *tmp,
        int *addTable){
    int i=0,j=0,k=0, i2;
    bool end = false;
    //handle trivial cases
    if(ff1->len == 0){
        copyFFElem(ff2,ret);
        return;
    }
    if(ff2->len == 0){
        copyFFElem(ff1,ret);
        return;
    }
    /*printf("addFFElem\n");*/
    /*printFFElem("ff1",ff1);*/
    /*printFFElem("ff2",ff2);*/
    copyArray(ff1->idcs,tmp,ff1->len);
    while( end == false ){
        while( tmp[i] != ff2->idcs[j] ){
            if( tmp[i] > ff2->idcs[j] ){
                /*printf("  take ff1 i=%i\n",i);*/
                ret->el[ tmp[i] ] = ff1->el[ tmp[i] ];
                ret->idcs[k] = tmp[i];
                i++; k++;
            }else if( tmp[i] < ff2->idcs[j] ){
                /*printf("  take ff2 j=%i\n",j);*/
                ret->el[ ff2->idcs[j] ] = ff2->el[ ff2->idcs[j] ];
                ret->idcs[k] = ff2->idcs[j];
                j++; k++;
            }
            if(i == ff1->len || j == ff2->len){
                end = true;
                break;
            }

        }
        if(end == true) break;
        /*printf("   same! i=%i j=%i\n",i,j);*/
        //tmp[i] == ff2->idcs[j]
        i2 = tmp[i];
        ret->el[i2] = addTable[ ff1->el[i2] + ff2->el[i2] ];
        if(ret->el[i2] != 0){
            ret->idcs[k] = i2;
            k++;
        }
        i++; j++;
        if(i == ff1->len || j == ff2->len) end = true;
    }
    /*printf("    take rest i=%i j=%i\n",i,j);*/
    //add rest of ff1 or ff2
    if(i != ff1->len ){
        while(i<ff1->len){
            /*printf("      take ff1 i=%i\n",i);*/
            ret->el[ tmp[i] ] = ff1->el[ tmp[i] ];
            ret->idcs[k] = tmp[i];
            i++; k++;
        }
    }else if(j != ff2->len){
        while(j<ff2->len){
            /*printf("      take ff2 j=%i\n",j);*/
            ret->el[ ff2->idcs[j] ] = ff2->el[ ff2->idcs[j] ];
            ret->idcs[k] = ff2->idcs[j];
            j++; k++;
        }
    }
    ret->len = k;
}

/**
 * Multiplies two FFElems and reduces the result by mipo.
 *
 * !! tmp must have at least length m!
 * !! ret must be malloced!
 */
inline void multiplyFFElem(struct FFElem *ff1, struct FFElem *ff2, 
        struct FFElem *ret, 
        struct FFElem *mipo, int *tmp, int m,
        int charac, int *addTable){
    /*printf("multiply ");printFFElemShort("ff1",ff1);*/
    /*printFFElemShort("         ff2",ff2);*/

    /* 
     * catch trivial cases
     */
    if(ff1->len == 0 || ff2->len == 0){
        ret->len = 0;
        return;
    }
    if(ff1->len == 1 && ff1->idcs[0] == 0 && ff1->el[0] == 1){
        copyFFElem(ff2, ret);
        return;
    }
    if(ff2->len == 1 && ff2->idcs[0] == 0 && ff2->el[0] == 1){
        copyFFElem(ff1, ret);
        return;
    }

    /*
     * Do multiplication
     */
    int maxlen = ff1->idcs[0] + ff2->idcs[0] + 1;
    int i,j,i2,j2,k;
    int max2 = maxlen;
    if( maxlen > m ){
        max2 = m;
        initPoly(tmp,maxlen-m);
    }
    initPoly(ret->el,max2);
    //multiply low part
    for(i=0;i<(ff1->len);i++){
        for(j=0;j<(ff2->len);j++){
            i2 = ff1->idcs[i];
            j2 = ff2->idcs[j];
            k = i2+j2;
            if(k<m){
                /*printf("mul i=%i j=%i ret[%i]=%i + ff1[%i]=%i * ff2[%i]=%i",*/
                        /*i,j,k,ret->el[k], i2, ff1->el[i2], j2, ff2->el[j2]);*/
                ret->el[k] = addTable[ ret->el[k] + 
                    (int)(((unsigned long)ff1->el[i2] * ff2->el[j2])%charac) ];
                /*printf("    =>ret[%i] = %i\n",k,ret->el[k]);*/
            }else{
                /*printf("k>m: i=%i j=%i k=%i tmp[%i]=%i",*/
                        /*i,j,k,k-m,tmp[k-m]);*/
                tmp[k-m] = addTable[ tmp[k-m] +
                    (int)(((unsigned long)ff1->el[i2] * ff2->el[j2])%charac) ];
                /*printf(" + %i * %i = %i\n",*/
                        /*ff1->el[i2],ff2->el[j2],tmp[k-m]);*/
            }
        }
    }

    /*printf("ret=");printArr(ret->el,m);*/
    /*printf("tmp=");printArr(tmp,m);*/

    /*
     * Reduce mod mipo
     */
    if(maxlen > m){
        int quo;
        for(i=maxlen-m-1;i>=0;i--){
            quo = tmp[i];
            if(quo == 0) continue;
            for(j=0;j<(mipo->len); j++){
                j2 = mipo->idcs[j];
                k = i+j2;
                if(k>=m){
                    tmp[k-m] = addTable[ tmp[k-m] - 
                        (int)(((unsigned long)mipo->el[j2]*quo)%charac) ];
                }else{
                    ret->el[k] = addTable[ ret->el[k] - 
                        (int)(((unsigned long)mipo->el[j2]*quo)%charac) ];
                }
            }
        }
    }

    /*
     * Recalc indices
     */
    i2 = 0;
    for(i=max2-1;i>=0;i--){
        if(ret->el[i] != 0){
            ret->idcs[i2] = i;
            i2++;
        }
    }
    ret->len = i2;
    /*printFFElemShort("            ",ret);*/
}

/**
 * Squares an FFElem
 *
 * !! ff is not modified !!
 * !! tmp must have at least length m !!
 */
inline void squareFFElem(struct FFElem *ff, struct FFElem *mipo,
        struct FFElem *ret, int *tmp, int m,
        int charac, int *addTable){
    /*printFFElemShort("square ff",ff);*/

    /* 
     * catch trivial cases
     */
    if(ff->len == 0){
        copyFFElem(ff,ret);
        return;
    }
    if(ff->len == 1 && ff->idcs[0] == 0 && ff->el[0] == 1){
        copyFFElem(ff,ret);
        return;
    }

    /*
     * Do multiplication
     */
    int maxlen = 2*ff->idcs[0] + 1;
    int i,j,i2,j2,k;
    int max2 = maxlen;
    if( maxlen > m ){
        max2 = m;
        initPoly(tmp,maxlen-m);
    }
    initPoly(ret->el,max2);
    for(i=0;i<(ff->len);i++){
        // same index must be squared
        i2 = ff->idcs[i];
        k = 2*i2;
        if(k<m){
            /*printf("square i=%i i2=%i k=%i ff[i2]=%i ret[k]=%i",*/
                    /*i,i2,k,ff->el[i2],ret->el[k]);*/
            ret->el[k] = addTable[ ret->el[k] + 
                (int)(((unsigned long)ff->el[i2]*ff->el[i2])%charac) ];
            /*printf("    =>ret[%i] = %i\n",k,ret->el[k]);*/
        }else{
            /*printf("square i=%i i2=%i k=%i ff[i2]=%i tmp[k-m]=%i",*/
                    /*i,i2,k,ff->el[i2],tmp[k-m]);*/
            tmp[k-m] = addTable[ tmp[k-m] +
                (int)(((unsigned long)ff->el[i2]*ff->el[i2])%charac) ];
            /*printf("    =>tmp[%i] = %i\n",k-m,tmp[k-m]);*/
        }
        // other indices only multipied and doubled
        for(j=i+1;j<(ff->len);j++){
            i2 = ff->idcs[i];
            j2 = ff->idcs[j];
            k = i2+j2;
            if(k<m){
                /*printf("mul i=%i i2=%i j=%i j2=%i k=%i ff[i2]=%i ff[j2]=%i ret[k]=%i",*/
                        /*i,i2,j,j2,k,ff->el[i2],ff->el[j2],ret->el[k]);*/
                ret->el[k] = addTable[ ret->el[k] + 
                    (int)((2 * (unsigned long)ff->el[i2] * ff->el[j2])%charac) ];
                /*printf("    =>ff[%i] = %i\n",k,ff->el[k]);*/
            }else{
                /*printf("mul i=%i i2=%i j=%i j2=%i k=%i ff[i2]=%i ff[j2]=%i tmp[k-m]=%i",*/
                        /*i,i2,j,j2,k,ff->el[i2],ff->el[j2],tmp[k-m]);*/
                tmp[k-m] = addTable[ tmp[k-m] +
                    (int)((2 * (unsigned long)ff->el[i2] * ff->el[j2])%charac) ];
                /*printf("    =>tmp[%i] = %i\n",k-m,tmp[k-m]);*/
            }
        }
    }

    /*printf("ret=");printArr(ret->el,m);*/
    /*printf("tmp=");printArr(tmp,m);*/

    /*
     * Reduce mod mipo
     */
    if(maxlen > m){
        int quo;
        for(i=maxlen-m-1;i>=0;i--){
            quo = tmp[i];
            if(quo == 0) continue;
            for(j=0;j<(mipo->len); j++){
                j2 = mipo->idcs[j];
                k = i+j2;
                if(k>=m){
                    tmp[k-m] = addTable[ tmp[k-m] - 
                        (int)(((unsigned long)mipo->el[j2]*quo)%charac) ];
                }else{
                    ret->el[k] = addTable[ ret->el[k] - 
                        (int)(((unsigned long)mipo->el[j2]*quo)%charac) ];
                }
            }
        }
    }

    /*
     * Recalc indices
     */
    i2 = 0;
    for(i=max2-1;i>=0;i--){
        if(ret->el[i] != 0){
            ret->idcs[i2] = i;
            i2++;
        }
    }
    ret->len = i2;
    /*printFFElemShort("            ",ret);*/
}

/**
 * Matrix multiplication
 *
 * !! tmp must have at least length m!
 */
inline void matmul(struct FFElem **mat, struct FFElem *ff,
        struct FFElem *ret, 
        int m, int charac, int *addTable){
    int i,j,i2, row;
    bool end;
    /*printf("matmul mat:\n");printFFElemMatrix(mat,m);*/
    /*printf(" * ");printFFElem("ff",ff);*/
    for(row=0;row<m;row++){
        /*printf("row=%i\n",row);*/
        ret->el[row] = 0;
        i=0; j=0;
        end = false;
        while(end == false){
            while(ff->idcs[i] != mat[row]->idcs[j]){
                if(ff->idcs[i] > mat[row]->idcs[j]) i++;
                else if(ff->idcs[i] < mat[row]->idcs[j]) j++;
                if(i == ff->len || j == mat[row]->len){
                    end = true;
                    break;
                }
            }
            /*printf("  i=%i j=%i end=%d\n",i,j,end);*/
            if(end == true) break;
            i2 = ff->idcs[i]; // == mat[row]->idcs[j]
            /*printf("  == i=%i j=%i i2=%i",i,j,i2);*/
            ret->el[row] = addTable[ ret->el[row] 
                + (int)(((unsigned long)mat[row]->el[i2]*ff->el[i2])%charac) ];
            i++;
            j++;
            if(i==ff->len || j==mat[row]->len) end = true;
        }
    }
    i2 = 0;
    for(i=m-1;i>=0;i--){
        if(ret->el[i] != 0){
            ret->idcs[i2] = i;
            i2++;
        }
    }
    ret->len = i2;
}

/**
 * Square and multiply
 *
 * ff is modified!
 */
inline void powerFFElemSqM(struct FFElem *ff, struct FFElem *mipo,
        struct FFElem *ret, 
        int m, int *power, int powerLen,
        int *tmp, struct FFElem *ffTmp,
        int charac, int *addTable){
    int i,j,k;
    int lenCurGap = 0;
    struct FFElem *ffSwitch = 0;
    struct FFElem *ffRetInt = ret;
    /*printFFElemShort("powerFFElemSqM ff",ff);*/
    // init ret to 1
    ffRetInt->el[0] = 1; ffRetInt->idcs[0] = 0; ffRetInt->len = 1;
    for(j=powerLen-1;j>=0;j--){
        /*printf("  power[%i]=%i\n",j,power[j]);*/
        if(power[j] == 1){
            multiplyFFElem(ffRetInt,ff,ffTmp, mipo,tmp,m,charac,addTable);
            /*printFFElemShort("     ffRetInt*ff",ffTmp);*/
            //switch ffTmp and ffRetInt
            ffSwitch = ffRetInt; ffRetInt = ffTmp; ffTmp = ffSwitch;
            /*copyFFElem(ffTmp,ffRetInt);*/
            /*printFFElemShort("     ffRetInt",ffRetInt);*/
        }
        if(j>0){
            squareFFElem(ff,mipo,ffTmp,tmp,m,charac,addTable);
            /*printFFElemShort("    ffTmp",ffTmp);*/
            /*printFFElemShort("    ff",ff);*/
            //switch ffTmp and ff
            ffSwitch = ff; ff = ffTmp; ffTmp = ffSwitch;
            /*printFFElem("     ff",ff);*/
        }
    }
    copyFFElem(ffRetInt,ret);
    /*printFFElemShort("   =>ret",ret);*/
}

/**
 * Square and multiply
 *
 * ff is modified!
 */
inline void powerFFElemSqMInt(struct FFElem *ff, struct FFElem *mipo,
        struct FFElem *ret, 
        int m, int power, 
        int *tmp, struct FFElem *ffTmp,
        int charac, int *addTable){
    int i,j,k;
    int lenCurGap = 0;
    struct FFElem *ffSwitch = 0;
    struct FFElem *ffRetInt = ret;
    /*printf("powerFFElemSqMInt power=%i",power); printFFElemShort(" ff",ff);*/
    // init ret to 1
    ffRetInt->el[0] = 1; ffRetInt->idcs[0] = 0; ffRetInt->len = 1;
    while( power > 0 ){
        /*printf("  power&1=%i\n",power&1);*/
        if(power & 1 == 1){
            multiplyFFElem(ffRetInt,ff,ffTmp, mipo,tmp,m,charac,addTable);
            /*printFFElemShort("     ffRetInt*ff",ffTmp);*/
            //switch ffTmp and ffRetInt
            ffSwitch = ffRetInt; ffRetInt = ffTmp; ffTmp = ffSwitch;
            /*printFFElemShort("     ffRetInt",ffRetInt);*/
        }

        if( power != 1){
            squareFFElem(ff,mipo,ffTmp,tmp,m,charac,addTable);
            /*printFFElemShort("    ffTmp",ffTmp);*/
            /*printFFElemShort("    ff",ff);*/
            //switch ffTmp and ff
            ffSwitch = ff; ff = ffTmp; ffTmp = ffSwitch;
            /*printFFElem("     ff",ff);*/
        }
        
        power >>= 1;
    }
    copyFFElem(ffRetInt,ret);
    /*printFFElemShort("   =>ret",ret);*/
}


/*
 * calculates g(sigma)(x) where g is a polynomial and sigma the frobenius
 * application of frobenius is given by mats
 */
inline void applyFrob(struct FFElem *ff, struct FFElem *mipo,
        struct FFPoly *poly,
        struct FFElem **mats,
        int frobPower, struct FFElem *ret, 
        int m, int *tmp, struct FFElem *ffTmp, struct FFElem *ffTmp2,
        struct FFElem **matmulCache, bool *matmulCacheCalced,
        int charac, int *addTable){
    int i,j;
    /*printFFElemShort("apply frob ff",ff);*/
    /*for(i=0;i<poly->lenPoly-1;i++){*/
        /*printf("mats+%i*%i*%i\n",i,frobPower,m);*/
        /*printFFElemMatrix(mats+i*frobPower,m);*/
    /*}*/
        
    ret->len = 0;
    for(i=0;i<poly->lenPoly;i++){
        /*printf("  apply i=%i",i);printFFElemShort(" poly[i]",poly->poly[i]);*/
        if(poly->poly[i]->len == 0) continue;
        j = i*frobPower-1;
        if(i>0 && matmulCacheCalced[j] == true){
            /*printf("  hasCache j=%i\n",j);*/
            multiplyFFElem(matmulCache[j],poly->poly[i],
                    ffTmp, mipo,
                    tmp,m,charac,addTable);
            addFFElem(ret,ffTmp,ret,tmp,addTable);
            /*printFFElemShort("  ->ffTmp",ffTmp);*/
        }else{
            /*printf("  noCache j=%i\n",j);*/
            if(i>0){
                matmul(mats+j*m, ff, ffTmp, m, charac,addTable);
                /*printFFElemShort("   matmul ffTmp",ffTmp);*/
                //update matmulCache
                copyFFElem(ffTmp, matmulCache[j]);
                matmulCacheCalced[j] = true;
            }else{
                copyFFElem(ff,ffTmp);
            }
            /*printFFElemShort("  ->ffTmp",ffTmp);*/
            //go on and multiply ffTmp with current coefficient
            multiplyFFElem(ffTmp, poly->poly[i],
                    ffTmp2, mipo,
                    tmp,m,charac,addTable);
            /*printFFElemShort("  ->ffTmp2",ffTmp2);*/
            addFFElem(ret,ffTmp2,ret,tmp,addTable);
            /*printFFElemShort("  =>ret",ret);*/
        }
    }
}


inline bool isCompletelyNormal(struct FFElem *ff, struct FFElem *mipo, 
        struct FFPoly **polys,
        int polysCount, 
        struct FFElem **mats, int *frobPowers,
        int m, int *tmp, 
        struct FFElem *ffTmp, struct FFElem *ffTmp2, struct FFElem *ffTmp3, 
        struct FFElem **matmulCache, bool *matmulCacheCalced,
        int charac, int *addTable){
    int i;
    int goodCounter = 0;
    /*printf("===test submod polysCount=%i ",polysCount);*/
    /*printFFElemShort("ff",ff);*/
    /*for(i=0;i<polysCount;i++)*/
        /*printFFPoly("   ",polys[i]);*/
    for(i=0;i<polysCount;i++){
        applyFrob(ff,mipo,
                polys[i],
                mats,frobPowers[i], ffTmp,
                m,tmp,ffTmp2,ffTmp3,
                matmulCache,matmulCacheCalced,
                charac,addTable);
        /*printf("   i=%i",i);printFFElemShort(" frobApply",ffTmp);*/
        if( isZero(ffTmp) ){
            /*printf(" false!\n===\n");*/
            return false;
        }
    }
    /*printf(" true!\n===\n");*/
    return true;
}

inline bool processFiniteField(struct FFElem *ff, struct FFElem *ffRet,
        int *nextPowers, int lenNextPowers, int maxPower,
        struct FFElem *mipo, 
        struct FFPoly **polys, int polysCount,
        struct FFElem **mats, int matLen, int *frobPowers,
        int m, int charac, int q,
        int *addTable){
    time_t TIME = time(NULL);
    int i,j;

    //setup temporary variables ----------------------------------------------
    int *tmp = malloc(m*sizeof(int));
    struct FFElem *fff = mallocFFElem(m);
    struct FFElem *ffLast = mallocFFElem(m);
    struct FFElem *ffTmp = mallocFFElem(m);
    struct FFElem *ffTmp2 = mallocFFElem(m);
    struct FFElem *ffTmp3 = mallocFFElem(m);
    struct FFElem *ffSwitch = 0;

    struct FFElem **powerCache = malloc((maxPower-1)*sizeof(struct FFElem*));
    for(i=0;i<maxPower-1;i++) powerCache[i] = 0;
    
    struct FFElem **matmulCache = malloc(matLen*sizeof(struct FFElem));
    for(i=0;i<matLen;i++) matmulCache[i] = mallocFFElem(m);
    bool *matmulCacheCalced = malloc(matLen*sizeof(bool));
    
    /*printf("nextPowers = ");printArr(nextPowers,lenNextPowers);*/
    /*printf("maxPower = %i\n",maxPower);*/
    /*printFFElemShort("ff",ff);*/

    // chase for pcn elements ------------------------------------------------
    bool done = false;
    // update ffLast to ffRet
    copyFFElem(ffRet,ffLast);
    for(i=0;i<lenNextPowers;i++){
        // calc current power of ff (== fff)
        if(nextPowers[i] == 1){  // no powering is required!
            /*printf("nextPowers[%i]==1\n",i);*/
            // multiply fff^current power with ffLast
            multiplyFFElem(ff,ffLast,ffTmp,
                    mipo,tmp,m,charac,addTable);
            //save result to ffLast
            ffSwitch = ffLast; ffLast = ffTmp; ffTmp = ffSwitch;
        }else if(powerCache[nextPowers[i]-2] != 0){ //powering is already done before
            /*printf("nextPowers[%i]==%i hasCache!\n",i,nextPowers[i]);*/
            /*printFFElemShort("    cache",powerCache[nextPowers[i]-2]);*/
            // multiply fff^current power with ffLast
            multiplyFFElem(powerCache[nextPowers[i]-2],ffLast,ffTmp,
                    mipo,tmp,m,charac,addTable);
            //save result to ffLast
            ffSwitch = ffLast; ffLast = ffTmp; ffTmp = ffSwitch;
        }else{ //powering was not done before; so do it now!
            /*printf("nextPowers[%i]==%i noCache!\n",i,nextPowers[i]);*/
            copyFFElem(ff,fff);
            //do powering
            powerFFElemSqMInt(fff,mipo,ffTmp,
                    m,nextPowers[i], tmp,ffTmp2,
                    charac,addTable);
            /*printFFElemShort("   fff^%i",ffTmp);*/
            //save result in cache
            j = nextPowers[i]-2;
            powerCache[j] = mallocFFElem(m);
            copyFFElem(ffTmp,powerCache[j]);
            // multiply fff^current power with ffLast
            multiplyFFElem(ffTmp,ffLast,ffTmp2,
                    mipo,tmp,m,charac,addTable);
            //save result to ffLast
            ffSwitch = ffLast; ffLast = ffTmp2; ffTmp2 = ffSwitch;
        }

        /*printFFElemShort("ffLast",ffLast);*/

        for(j=0;j<matLen;j++) matmulCacheCalced[j] = false;
        if(isCompletelyNormal(ffLast, mipo, polys, polysCount,
                    mats,frobPowers,
                    m,tmp,
                    ffTmp,ffTmp2,ffTmp3,
                    matmulCache,matmulCacheCalced,
                    charac,addTable)){
            /*printf("is CN!\n");*/
            done = true;
            break;
        }
        /*else{*/
            /*printf("not CN!\n");*/
        /*}*/
    }
    copyFFElem(ffLast,ffRet);
    
    //free temporary variables
    for(i=0;i<maxPower-1;i++) freeFFElem(powerCache[i]);
    free(powerCache);
    free(tmp);
    freeFFElem(fff);
    freeFFElem(ffLast);
    freeFFElem(ffTmp);
    freeFFElem(ffTmp2);
    freeFFElem(ffTmp3);
    for(i=0;i<matLen;i++) freeFFElem(matmulCache[i]);
    free(matmulCache);
    free(matmulCacheCalced);

    
    /*printf("total time: %.2f\n", (double)(time(NULL)-TIME));*/

    return done;
}



///////////////////////////////////////////////////////////////////////////////
// Helper for Creation of FFElems /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

inline void transposeMatrix(struct FFElem **mat, struct FFElem **matT, int m){
    int i,j,j2;
    for(i=0;i<m;i++) initPoly(matT[i]->el,m);
    for(i=0;i<m;i++){
        for(j=0;j<mat[i]->len;j++)
            matT[ mat[i]->idcs[j] ]->el[i] = mat[i]->el[ mat[i]->idcs[j] ];
    }
    for(i=0;i<m;i++){
        j2=0;
        for(j=m-1;j>=0;j--){
            if(matT[i]->el[j] != 0){
                matT[i]->idcs[j2] = j;
                j2++;
            }
        }
        matT[i]->len = j2;
    }
}

/**
 * Multiplies two matrices
 *
 * !! mat1, mat1, ret, matTmp must have same size and be quadratic !!
 */
inline void multiplyMatrices(struct FFElem **mat1, struct FFElem **mat2,
        struct FFElem **ret, 
        struct FFElem **matTmp,
        int m, int charac, int *addTable){
    int i,j,i2, row,col;
    bool end;
    /*printf("multiply matrices \n");*/
    /*printFFElemMatrix(mat1,m);printf(" * \n");*/
    /*printFFElemMatrix(mat2,m);*/
    transposeMatrix(mat2,matTmp,m);
    /*printf("matT = \n"); printFFElemMatrix(matTmp,m);*/
    for(row=0;row<m;row++){
        for(col=0;col<m;col++){
            ret[row]->el[col] = 0;
            i=0; j=0;
            end = false;
            while(end == false){
                while(mat1[row]->idcs[i] != matTmp[col]->idcs[j]){
                    if(mat1[row]->idcs[i] > matTmp[col]->idcs[j])  i++;
                    else if(mat1[row]->idcs[i] < matTmp[col]->idcs[j])  j++;
                    if(i == mat1[row]->len || j == matTmp[col]->len){
                        end = true;
                        break;
                    }
                }
                if(end == true) break;
                i2 = mat1[row]->idcs[i]; // == matTmp[col]->idcs[j]
                /*printf("row = %i col = %i i2 = %i ret[row]->el[col] = %i",*/
                        /*row,col,i2,ret[row]->el[col]);*/
                ret[row]->el[col] = addTable[ ret[row]->el[col] 
                    + (int)(((unsigned long)mat1[row]->el[i2]
                                *matTmp[col]->el[i2])%charac) ];
                /*printf(" + %i * %i = %i\n",mat1[row]->el[i2],matTmp[col]->el[i2],*/
                        /*ret[row]->el[col]);*/
                i++;
                j++;
                if(i==mat1[row]->len || j==matTmp[col]->len) end = true;
            }
        }
        i2 = 0;
        for(i=m-1;i>=0;i--){
            if(ret[row]->el[i] != 0){
                ret[row]->idcs[i2] = i;
                i2++;
            }
        }
        ret[row]->len = i2;
    }
}

/**
 * Generates Frobenius Matrices.
 *
 * !! mats are malloced here !!
 * !! all temporary variables are malloced inside !!
 */
struct FFElem **genFrobMats(struct FFElem *mipo,int m,int maxPower,int q,
        int charac, int *addTable){
    int i;
    int matSize = maxPower*m;
    struct FFElem **mats = malloc(matSize*sizeof(struct FFElem*));
    struct FFElem **matTmp = malloc(m*sizeof(struct FFElem*));
    struct FFElem *ffTmp = mallocFFElem(m);
    struct FFElem *ffTmp2 = mallocFFElem(m);
    int *tmp = malloc(m*sizeof(int));

    /*printf("genFrobMats m=%i q=%i maxPower=%i\n",m,q,maxPower);*/

    for(i=0;i<matSize;i++)
        mats[i] = mallocFFElem(m);
    for(i=0;i<m;i++)
        matTmp[i] = mallocFFElem(m);

    for(i=0;i<m;i++){
        ffTmp->el[i] = 1; ffTmp->idcs[0] = i; ffTmp->len = 1;
        /*printf("i=%i",i); printFFElemShort(" ffTmp",ffTmp);*/
        powerFFElemSqMInt(ffTmp, mipo, mats[i],
                m,q,tmp,ffTmp2,charac,addTable);
        /*printFFElemShort(" mats[i]",mats[i]);*/
    }
    //transpose first matrix
    transposeMatrix(mats,matTmp,m);
    for(i=0;i<m;i++){
        copyFFElem(matTmp[i],mats[i]);
    }
    /*printf("first matrix:\n");printFFElemMatrix(mats,m);*/
    for(i=1;i<maxPower;i++){
        multiplyMatrices(mats, mats+(i-1)*m, mats+i*m, matTmp,
                m,charac,addTable);
    }
    for(i=0;i<m;i++) freeFFElem(matTmp[i]);
    free(matTmp);
    free(tmp);
    free(ffTmp);
    free(ffTmp2);
    return mats;
}

inline int *genAddTable(int charac){
    int tmp = 2*(charac-1);
    int *addTable = malloc((2*tmp+1)*sizeof(int));
    addTable += tmp;
    int i;
    int cur=0;
    for(i=0;i<tmp+1;i++){
        addTable[i] = cur;
        cur++;
        if(cur == charac) cur = 0;
    }
    cur = charac-1;
    for(i=1;i<tmp+1;i++){
        addTable[-i] = cur;
        cur--;
        if(cur == -1) cur = charac-1;
    }
    
    /*printf("addTable p=%i ",charac);printArr(addTable-tmp,2*tmp+1);*/
    
    return addTable;
}


void main(){
    printf("int max: %i\n",INT_MAX);
    printf("uint max: %u\n",UINT_MAX);
    printf("ulong max: %lu\n",ULONG_MAX);
}
