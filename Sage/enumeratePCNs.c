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
// Setup Linked Lists /////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
struct Node {
    struct FFElem *ff;
    struct Node *next;
};

/**
 * appends element to end of root, where element is copied to new FFElem
 */
inline struct Node *appendToEnd(struct Node *root, struct FFElem *element,int m){
    struct Node *nextNode = root;
    if( nextNode != 0){
        while(nextNode->next != 0){
            nextNode = nextNode->next;
        }
        if( nextNode->ff != 0){
            nextNode->next = malloc( sizeof(struct Node) );
            nextNode = nextNode->next;
        }
        if( nextNode != 0){
            nextNode->next = 0;
            nextNode->ff = mallocFFElem(m);
            copyFFElem(element,nextNode->ff);
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
        freeFFElem(tmp_n->ff);
        /*printf(" data freed");*/
        free(tmp_n);
        /*printf(" tmp freed ");*/
        tmp_n = next_n;
        /*printf(" next=%i\n",n);*/
    }
    head = 0;
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
        int *multTable, int *addTable){
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
        int *multTable, int *addTable){
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
                    multTable[ ff1->el[i2] * ff2->el[j2] ] ];
                /*printf("    =>ret[%i] = %i\n",k,ret->el[k]);*/
            }else{
                /*printf("k>m: i=%i j=%i k=%i tmp[%i]=%i",*/
                        /*i,j,k,k-m,tmp[k-m]);*/
                tmp[k-m] = addTable[ tmp[k-m] +
                    multTable[ ff1->el[i2] * ff2->el[j2] ] ];
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
                        multTable[ mipo->el[j2]*quo ] ];
                }else{
                    ret->el[k] = addTable[ ret->el[k] - 
                        multTable[ mipo->el[j2]*quo ] ];
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
        int *multTable, int *addTable){
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
                multTable[ ff->el[i2]*ff->el[i2] ] ];
            /*printf("    =>ret[%i] = %i\n",k,ret->el[k]);*/
        }else{
            /*printf("square i=%i i2=%i k=%i ff[i2]=%i tmp[k-m]=%i",*/
                    /*i,i2,k,ff->el[i2],tmp[k-m]);*/
            tmp[k-m] = addTable[ tmp[k-m] +
                multTable[ ff->el[i2]*ff->el[i2] ] ];
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
                    multTable[ 2 * multTable[ ff->el[i2] * ff->el[j2] ] ] ];
                /*printf("    =>ff[%i] = %i\n",k,ff->el[k]);*/
            }else{
                /*printf("mul i=%i i2=%i j=%i j2=%i k=%i ff[i2]=%i ff[j2]=%i tmp[k-m]=%i",*/
                        /*i,i2,j,j2,k,ff->el[i2],ff->el[j2],tmp[k-m]);*/
                tmp[k-m] = addTable[ tmp[k-m] +
                    multTable[ 2 * multTable[ ff->el[i2] * ff->el[j2] ] ]];
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
                        multTable[ mipo->el[j2]*quo ] ];
                }else{
                    ret->el[k] = addTable[ ret->el[k] - 
                        multTable[ mipo->el[j2]*quo ] ];
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
        int m, int *multTable, int *addTable){
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
                + multTable[ mat[row]->el[i2]*ff->el[i2] ] ];
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
 * Square and multiply in charac
 * mat is powering by charac
 *
 * ff is modified!
 */
inline void powerFFElem(struct FFElem *ff, struct FFElem *mipo,
        struct FFElem *ret, 
        int m, int *power, int powerLen,
        struct FFElem **matCharac, int *tmp, struct FFElem *ffTmp,
        int *multTable, int *addTable){
    int i,j,k;
    int lenCurGap = 0;
    struct FFElem *ffSwitch = 0;
    struct FFElem *ffRetInt = ret;
    /*printFFElemShort("powerFFElem ff",ff);*/
    // init ret to 1
    ffRetInt->el[0] = 1; ffRetInt->idcs[0] = 0; ffRetInt->len = 1;
    for(j=powerLen-1;j>=0;j--){
        /*printf("  power[%i]=%i\n",j,power[j]);*/
        for(k=0;k<power[j];k++){
            multiplyFFElem(ffRetInt,ff,ffTmp, mipo,tmp,m,multTable,addTable);
            /*printFFElemShort("     ffRetInt*ff",ffTmp);*/
            ffSwitch = ffRetInt; ffRetInt = ffTmp; ffTmp = ffSwitch;
            /*printFFElemShort("     ffRetInt",ffRetInt);*/
        }
        if(j==0 || power[j-1] == 0){
            lenCurGap++;
            continue;
        }
        /*printf("     lenCurGap=%i",lenCurGap);*/
        matmul(matCharac+lenCurGap*m, ff, ffTmp, m, multTable,addTable);
        /*printFFElem("     mat*ff",ffTmp);*/
        ffSwitch = ff; ff = ffTmp; ffTmp = ffSwitch;
        /*printFFElem("     ff",ff);*/
        lenCurGap = 0;
    }
    copyFFElem(ffRetInt,ret);
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
        int *multTable, int *addTable){
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
            multiplyFFElem(ffRetInt,ff,ffTmp, mipo,tmp,m,multTable,addTable);
            /*printFFElemShort("     ffRetInt*ff",ffTmp);*/
            //switch ffTmp and ffRetInt
            ffSwitch = ffRetInt; ffRetInt = ffTmp; ffTmp = ffSwitch;
            /*copyFFElem(ffTmp,ffRetInt);*/
            /*printFFElemShort("     ffRetInt",ffRetInt);*/
        }
        if(j>0){
            squareFFElem(ff,mipo,ffTmp,tmp,m,multTable,addTable);
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

inline void powerFFElemInt(struct FFElem *ff, struct FFElem *mipo,
        struct FFElem *ret, 
        int m, int power, 
        int *tmp, struct FFElem *ffTmp,
        int *multTable, int *addTable){
    int i,j,k;
    // init ret to 1
    ret->el[0] = 1; ret->idcs[0] = 0; ret->len = 1;
    /*printf("powerFFElemInt power=%i",power); printFFElemShort(" ff",ff);*/
    while(power > 0){
        if(power & 1 == 1){
            multiplyFFElem(ret,ff,ffTmp, mipo,tmp,m,multTable,addTable);
            /*printFFElem("   ret*ff",ffTmp);*/
            copyFFElem(ffTmp,ret);
            /*printFFElem("   ret",ret);*/
        }
        multiplyFFElem(ff, ff, ffTmp, mipo, tmp,m,multTable,addTable);
        /*printFFElem("   mat*ff",ffTmp);*/
        copyFFElem(ffTmp,ff);
        /*printFFElem("   ff",ff);*/
        power >>= 1;
    }
}

/**
 * Calc order of element
 *
 * !! if matCharac is Zero, all powers are assumed as binary arrays !!
 *
 * !! fff,ffTmp,ffTmp2,ffTmp3,ffRet must be malloced !!
 * !! x is NOT modified !!
 */
inline bool isPrimitive(struct FFElem *ff, struct FFElem *mipo,
        int m,
        int *barFactors, int *lenBarFactors, int countBarFactors,
        int *commonBarFactor, int lenCommonBarFactor,
        int *commonBiggestBarFactor, int lenCommonBiggestBarFactor,
        struct FFElem **matCharac, 
        struct FFElem *fff, struct FFElem *ffff, struct FFElem *ffTmp,
        struct FFElem *ffTmp2, struct FFElem *ffRet, 
        int *tmp, int *multTable, int *addTable){
    int i;
    int curPos = 0;
    bool binarySqM = (matCharac == 0);
    struct FFElem *ffSwitch = 0;
    /*printFFElemShort("\n\ntestPrimitivity ff=",ff);*/

    copyFFElem(ff,fff);
    // all barFactors are power of commonBarFactor
    if(binarySqM)
        powerFFElemSqM(fff,mipo,ffTmp,
                m,commonBarFactor,lenCommonBarFactor,
                tmp,ffTmp2,
                multTable,addTable);
    else 
        powerFFElem(fff,mipo,ffTmp,
                m,commonBarFactor,lenCommonBarFactor,
                matCharac,tmp,ffTmp2,
                multTable,addTable);
    /*printf("commonBarFactor=");printArr(commonBarFactor,lenCommonBarFactor);*/
    /*printFFElem("   y",ffTmp);*/
    if(isOne(ffTmp)) return false;
    //switch ffTmp and fff
    ffSwitch = fff; fff = ffTmp; ffTmp = ffSwitch;
    copyFFElem(fff,ffff);
    //test first barFactor
    if(binarySqM)
        powerFFElemSqM(ffff,mipo,ffTmp,
                m,barFactors,lenBarFactors[0],
                tmp,ffTmp2,
                multTable,addTable);
    else
        powerFFElem(ffff,mipo,ffTmp,
                m,barFactors,lenBarFactors[0],
                matCharac,tmp,ffTmp2,
                multTable,addTable);
    /*printf("firstFactor=");printArr(barFactors,lenBarFactors[0]);*/
    /*printFFElem("   y^firstFactor",ffTmp);*/
    if(isOne(ffTmp)) return false;
    curPos += lenBarFactors[0];
    //test further factors which are powers of commonBiggestBarFactor
    //so first, calc y^commonBiggestBarFactor
    copyFFElem(fff, ffff);
    if(binarySqM)
        powerFFElemSqM(ffff,mipo,ffTmp,
                m,commonBiggestBarFactor,lenCommonBiggestBarFactor,
                tmp,ffTmp2,
                multTable,addTable);
    else
        powerFFElem(ffff,mipo,ffTmp,
                m,commonBiggestBarFactor,lenCommonBiggestBarFactor,
                matCharac,tmp,ffTmp2,
                multTable,addTable);
    if(isOne(ffTmp)) return false;
    ffSwitch = fff; fff = ffTmp; ffTmp = ffSwitch;
    /*printf("commonBiggestBarFactor=");printArr(commonBiggestBarFactor,lenCommonBiggestBarFactor);*/
    /*printFFElem("   z=y^commonBiggestBarFactor",fff);*/
    for(i=1;i<countBarFactors;i++){
        // copy z (fff) to ffff
        copyFFElem(fff,ffff);
        // *** ffff == fff == y^commonBiggestBarFactor
        if(binarySqM)
            powerFFElemSqM(ffff, mipo, ffTmp, 
                    m,barFactors+curPos, lenBarFactors[i],
                    tmp,ffTmp2,
                    multTable,addTable);
        else
            powerFFElem(ffff, mipo, ffTmp, 
                    m,barFactors+curPos, lenBarFactors[i],
                    matCharac,tmp,ffTmp2,
                    multTable,addTable);
        /*printf("barFactor[%i]=",i);printArr(barFactors+curPos,*/
                /*lenBarFactors[i]);*/
        /*printFFElem("   z^curFac",ffTmp);*/

        if(i>1){
            multiplyFFElem(ffRet,ffTmp,ffTmp2,mipo,
                    tmp,m,multTable,addTable);
        }else{
            ffSwitch = ffTmp2; ffTmp2 = ffTmp; ffTmp = ffSwitch;
        }
        /*printFFElem("ffTmp2",ffTmp2);*/
        if(isOne(ffTmp2)) return false;
        curPos += lenBarFactors[i];
        ffSwitch = ffRet; ffRet = ffTmp2; ffTmp2 = ffSwitch;
    }
    /*printf("isPrimitive!\n");*/
    return true;
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
        int *multTable, int *addTable){
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
                    tmp,m,multTable,addTable);
            addFFElem(ret,ffTmp,ret,tmp,multTable,addTable);
            /*printFFElemShort("  ->ffTmp",ffTmp);*/
        }else{
            /*printf("  noCache j=%i\n",j);*/
            if(i>0){
                matmul(mats+j*m, ff, ffTmp, m, multTable,addTable);
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
                    tmp,m,multTable,addTable);
            /*printFFElemShort("  ->ffTmp2",ffTmp2);*/
            addFFElem(ret,ffTmp2,ret,tmp,multTable,addTable);
            /*printFFElemShort("  =>ret",ret);*/
        }
    }
}

/*
 * calculates g(sigma)(x) where g is a polynomial and sigma the frobenius
 * application of frobenius is given by mats
 */
inline void applyFrob_noCache(struct FFElem *ff, struct FFElem *mipo,
        struct FFPoly *poly,
        struct FFElem **mats,
        int frobPower, struct FFElem *ret, 
        int m, int *tmp, struct FFElem *ffTmp, struct FFElem *ffTmp2,
        int *multTable, int *addTable){
    int i,j;
    /*printFFElemShort("apply frob noCache ff",ff);*/
    ret->len = 0;
    
    for(i=0;i<poly->lenPoly;i++){
        /*printf("  apply i=%i",i);printFFElemShort(" poly[i]=",poly->poly[i]);*/
        if(poly->poly[i]->len == 0) continue;
        if(i>0){
            j = i*frobPower-1;
            matmul(mats+j*m, ff, ffTmp, m, multTable,addTable);
        }else{
            copyFFElem(ff,ffTmp);
        }
        /*printFFElemShort("  ->ffTmp",ffTmp);*/
        multiplyFFElem(ffTmp, poly->poly[i],
                ffTmp2, mipo,
                tmp,m,multTable,addTable);
        /*printFFElemShort("  ->ffTmp2",ffTmp2);*/
        addFFElem(ret,ffTmp2,ret,tmp,multTable,addTable);
        /*printFFElemShort("  =>ret",ret);*/
    }
}


inline int testAllSubmods(struct FFElem *ff, struct FFElem *mipo, 
        int decompCount, struct FFPoly **polys,
        int *polysCountPerDecomp, bool *evalToZero, 
        struct FFElem **mats, int *frobPowers, bool *toTestIndicator,
        int m, int *tmp, 
        struct FFElem *ffTmp, struct FFElem *ffTmp2, struct FFElem *ffTmp3,
        struct FFElem **matmulCache, bool *matmulCacheCalced,
        int *multTable, int *addTable){
    if(ff->len == 0) return -1;
    int i,j,k;
    int goodCounter = 0;
    int curDecompPosition = 0;
    /*printf("=====\ntestAllSubmods ");printFFElemShort("ff",ff);*/
    /*printf("   polysCountPerDecomp");printArr(polysCountPerDecomp,decompCount);*/
    /*printf("   toTestIndicator    ");printBoolArr(toTestIndicator,decompCount);*/
    for(i=0;i<decompCount;i++){
        if(toTestIndicator != 0 && toTestIndicator[i] == false){
            curDecompPosition += polysCountPerDecomp[i];
            continue;
        }
        /*printf("   test gen=%i polysCount=%i\n",i,polysCountPerDecomp[i]);*/
        /*printf("        evalToZero=");printBoolArr(evalToZero+curDecompPosition,polysCountPerDecomp[i]);*/
        goodCounter = 0;
        for(j=0;j<polysCountPerDecomp[i];j++){
            /*printf("        j=%i",j);*/
            /*printf("        matmulCacheCalced=");printBoolArr(matmulCacheCalced,*/
                    /*polys[curDecompPosition+j]->lenPoly-1);*/
            /*printFFPoly("        poly",polys[curDecompPosition+j]);*/
            /*printf("        frobPower=%i evalToZero=%d\n",*/
                    /*frobPowers[curDecompPosition+j],evalToZero[curDecompPosition+j]);*/
            applyFrob(ff,mipo,
                    polys[curDecompPosition+j],
                    mats,frobPowers[curDecompPosition+j], ffTmp,
                    m,tmp,ffTmp2,ffTmp3,
                    matmulCache,matmulCacheCalced,
                    multTable,addTable);
            /*printFFElemShort("   applied",ffTmp);*/
            if( isZero(ffTmp) == evalToZero[curDecompPosition+j] ){
                goodCounter++;
            }else break;
        }
        if(goodCounter == polysCountPerDecomp[i]){
            /*printf("    return %i\n=====\n",i);*/
            return i;
        }
        curDecompPosition += polysCountPerDecomp[i];
    }
    /*printf("    return -1\n=====\n");*/
    return -1;
}

inline bool testSubmod(struct FFElem *ff, struct FFElem *mipo, 
        struct FFPoly **polys,
        int polysCount, bool *evalToZero, 
        struct FFElem **mats, int *frobPowers,
        int m, int *tmp, 
        struct FFElem *ffTmp, struct FFElem *ffTmp2, struct FFElem *ffTmp3, 
        struct FFElem **matmulCache, bool *matmulCacheCalced,
        int *multTable, int *addTable){
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
                multTable,addTable);
        /*printf("   i=%i",i);printFFElemShort(" frobApply",ffTmp);*/
        if( isZero(ffTmp) == evalToZero[i] ){
            /*printf("is good!\n");*/
            goodCounter++;
        }else{
            /*printf(" false!\n===\n");*/
            return false;
        }
    }
    if(goodCounter == polysCount){
        /*printf(" true!\n===\n");*/
        return true;
    }
    return false;
}


///////////////////////////////////////////////////////////////////////////////
// Central Algorithms for Enumeration /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
inline void calcSubmoduleElements(struct Node *root,
        struct FFElem *mipo,
        int maxLenPoly, 
        int *genCounts, int curGen,
        struct FFPoly **polys, int polysCount, bool *evalToZero,
        struct FFElem **mats, int matLen, int *frobPowers,
        struct FFElem **elementsF,
        int m, int q, int *tmp,
        struct FFElem *ffTmp, struct FFElem *ffTmp2, struct FFElem *ffTmp3, 
        struct FFElem *ffTmp4,
        struct FFElem **matmulCache, bool *matmulCacheCalced,
        int *multTable, int *addTable){
    int i,j;
    struct Node *curRoot = root;
    struct FFElem *ff = root->ff;
    int *curPoly = malloc( maxLenPoly*sizeof(int) );
    struct FFPoly *curFPoly = malloc( sizeof(struct FFPoly) );
    curFPoly->poly = malloc( maxLenPoly*sizeof(struct FFElem*) );
    curFPoly->lenPoly = 0;

    /*printf("========= calcSubmoduleElements curGen=%i maxLenPoly=%i\n",i,maxLenPoly);*/
    /*printFFElemShort("   ff",ff);*/
    
    initPoly(curPoly,maxLenPoly);
    curPoly[0] = 2;
    int curLenPoly = 1;
    if( q == 2 && maxLenPoly > 1){
        curLenPoly = 2;
        curPoly[0] = 0;
        curPoly[1] = 1;
    }
    /*printf("    curPoly");printArr(curPoly,maxLenPoly);*/
    if(q != 2 || maxLenPoly > 1){
        while(true){
            //setup curFPoly
            for(i=0;i<curLenPoly;i++)
                curFPoly->poly[i] = elementsF[curPoly[i]];
            curFPoly->lenPoly = curLenPoly;
            /*printFFPoly("   curFPoly",curFPoly);*/
            //apply Frobenius
            applyFrob_noCache(ff,mipo,
                    curFPoly,
                    mats,1, ffTmp, //return value
                    m,tmp,ffTmp2,ffTmp3,
                    multTable,addTable);
            /*printFFElemShort("   ->ffTmp",ffTmp);*/
            //test generated element
            for(i=0;i<matLen;i++) matmulCacheCalced[i] = false;
            if(testSubmod(ffTmp, mipo,
                    polys,polysCount,evalToZero,
                    mats,frobPowers,m,tmp,
                    ffTmp2,ffTmp3,ffTmp4,
                    matmulCache,matmulCacheCalced, multTable,addTable)){
                curRoot = appendToEnd(curRoot,ffTmp,m);
                genCounts[curGen]++;
            }
            //generate next element
            curPoly[0] += 1;
            if( curPoly[0] == q ){
                for(i=0;i<maxLenPoly-1 && curPoly[i]==q;i++){
                    curPoly[i] = 0;
                    curPoly[i+1] += 1;
                }
                if(i+1>curLenPoly)
                    curLenPoly = i+1;
                if( curPoly[maxLenPoly-1]==q){
                    break;
                }
            }
        }
    }
    free(curPoly);
    free(curFPoly->poly);
    free(curFPoly);
    /*printf("========\n");*/
}

/**
 * Processes last submodule as others before, but does not save generated
 * elements. Cycles through already generated elements and tests for 
 * primitivity:
 *
 * !! all temporary variables are generated inside !!
 */
unsigned long long processLastSubmoduleAndTestPrimitivity(struct Node **roots,
        struct FFElem *mipo, int decompCount,
        int maxLenPoly, 
        int *genCounts, 
        struct FFPoly **polys, int polysCount, bool *evalToZero,
        struct FFElem **mats, int matLen, int *frobPowers,
        struct FFElem **elementsF,
        int m, int q, 
        int *barFactors, int *lenBarFactors, int countBarFactors,
        int *commonBarFactor, int lenCommonBarFactor,
        int *commonBiggestBarFactor, int lenCommonBiggestBarFactor,
        struct FFElem **matCharac,
        struct FFElem **matmulCache, bool *matmulCacheCalced,
        int *multTable, int *addTable){
    /*printf("==============================================================\n");*/
    /*printf("processLast: start maxLenPoly=%i\n",maxLenPoly);*/
    //generate temporary variables
    struct FFElem *fff = mallocFFElem(m);
    struct FFElem *ffff = mallocFFElem(m);
    struct FFElem *ffTmp = mallocFFElem(m);
    struct FFElem *ffTmp2 = mallocFFElem(m);
    struct FFElem *ffTmp3 = mallocFFElem(m);
    struct FFElem *ffTmp4 = mallocFFElem(m);
    struct FFElem *ffTmp5 = mallocFFElem(m);
    int *tmp = malloc(m*sizeof(int));

    int i,j;
    int curGen = decompCount-1;
    struct Node **curRoots = malloc( decompCount*sizeof(struct Node*) );
    struct FFElem *ff = roots[curGen]->ff;
    int *curPoly = malloc( maxLenPoly*sizeof(int) );
    struct FFPoly *curFPoly = malloc( sizeof(struct FFPoly) );
    curFPoly->poly = malloc( maxLenPoly*sizeof(struct FFElem*) );
    curFPoly->lenPoly = 0;
    
    initPoly(curPoly,maxLenPoly);
    curPoly[0] = 1;
    int curLenPoly = 1;

    unsigned long long pcn = 0;
    while(true){
        //setup curFPoly
        for(i=0;i<curLenPoly;i++)
            curFPoly->poly[i] = elementsF[curPoly[i]];
        curFPoly->lenPoly = curLenPoly;
        /*printf("\ncurPoly=");printArr(curPoly,maxLenPoly);*/
        /*printFFPoly("curFPoly",curFPoly);*/
        //apply Frobenius
        applyFrob_noCache(ff,mipo,
                curFPoly,
                mats,1, fff, //return value
                m,tmp,ffTmp,ffTmp2,
                multTable,addTable);
        /*printFFElem("     =>f(x)",fff);*/
        //test generated element
        for(i=0;i<matLen;i++) matmulCacheCalced[i] = false;
        if(testSubmod(fff, mipo,
                polys,polysCount,evalToZero,
                mats,frobPowers,m,tmp,
                ffTmp,ffTmp2,ffTmp3,
                matmulCache,matmulCacheCalced, multTable,addTable)){
            /*printFFElemShort("found gen",fff);*/
            genCounts[curGen]++;
            // build element as sum of already calced Nodes and ffTmp
            for(i=0;i<decompCount;i++) curRoots[i] = roots[i];
            // cycle through Nodes, build element and test primitivity
            while(true){
                copyFFElem(fff,ffff);
                /*printf("   add others: \n");*/
                //build element
                for(i=0;i<decompCount-1;i++){
                    /*printf("     gen=%i",i);*/
                    /*printFFElemShort("",curRoots[i]->ff);*/
                    addFFElem(ffff, curRoots[i]->ff, ffff, tmp,
                            multTable,addTable);
                }
                /*printFFElemShort("   => got",ffff);*/
                /*printFFElemShort("test ff",ffff);*/
                //test primitivity
                if(isPrimitive(ffff, mipo,m,
                            barFactors,lenBarFactors,countBarFactors,
                            commonBarFactor,lenCommonBarFactor,
                            commonBiggestBarFactor,lenCommonBiggestBarFactor,
                            matCharac,
                            ffTmp,ffTmp2,ffTmp3,ffTmp4,ffTmp5,
                            tmp,multTable,addTable)){
                    /*printf("   isPrimitive!\n\n\n");*/
                    pcn++;
                }else{
                    /*printf("   is NOT primitive!\n\n\n");*/
                }

                //next element
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
        }
        //generate next element
        curPoly[0] += 1;
        if( curPoly[0] == q ){
            for(i=0;i<maxLenPoly-1 && curPoly[i]==q;i++){
                curPoly[i] = 0;
                curPoly[i+1] += 1;
            }
            if(i+1>curLenPoly)
                curLenPoly = i+1;
            if( curPoly[maxLenPoly-1]==q){
                break;
            }
        }
    }
    free(curPoly);
    free(curFPoly->poly);
    free(curFPoly);
    freeFFElem(fff);
    freeFFElem(ffff);
    freeFFElem(ffTmp);
    freeFFElem(ffTmp2);
    freeFFElem(ffTmp3);
    freeFFElem(ffTmp4);
    freeFFElem(ffTmp5);

    //we added first element twice
    genCounts[curGen]--;
    return pcn;
}


unsigned long long processFiniteField(struct FFElem *mipo, int decompCount,
        struct FFPoly **polys, int *polysCountPerDecomp,
        bool *evalToZero, 
        struct FFElem **mats, int matLen, int *frobPowers,
        int *genCounts, int m, int charac, int q,
        int *barFactors, int *lenBarFactors, int countBarFactors,
        int *commonBarFactor, int lenCommonBarFactor,
        int *commonBiggestBarFactor, int lenCommonBiggestBarFactor,
        struct FFElem **matCharac, struct FFElem **elementsF,
        int *multTable, int *addTable){
    time_t TIME = time(NULL);
    int i,j;

    //setup temporary variables ----------------------------------------------
    int *tmp = malloc(m*sizeof(int));
    struct FFElem *ff = mallocFFElem(m);
    initPoly(ff->el,m);
    struct FFElem *ffRet = mallocFFElem(m);
    struct FFElem *ffTmp = mallocFFElem(m);
    struct FFElem *ffTmp2 = mallocFFElem(m);
    struct FFElem *ffTmp3 = mallocFFElem(m);
    struct FFElem *ffTmp4 = mallocFFElem(m);
    
    struct FFElem **matmulCache = malloc(matLen*sizeof(struct FFElem));
    for(i=0;i<matLen;i++) matmulCache[i] = mallocFFElem(m);
    bool *matmulCacheCalced = malloc(matLen*sizeof(bool));
    
    bool *toTestIndicator = malloc(decompCount*sizeof(bool));
    struct Node **roots = malloc( decompCount*sizeof(struct Node) );
    struct Node **curRoots = malloc(decompCount*sizeof(struct Node*));
    for(i=0;i<decompCount;i++){
        roots[i] = malloc( sizeof(struct Node) );
        roots[i]->ff = 0;
        roots[i]->next = 0;
        curRoots[i] = roots[i];
        toTestIndicator[i] = true;
    }
    //------------------------------------------------------------------------
    
    int foundCounter = 0;
    initPoly(genCounts,decompCount);

    // chase for elements ----------------------------------------------------
    while(true){
        /*printFFElem("test element",ff);*/
        for(i=0;i<matLen;i++) matmulCacheCalced[i] = 0;
        int curGen = testAllSubmods(ff,mipo,decompCount,
                polys,polysCountPerDecomp,evalToZero,
                mats,frobPowers,toTestIndicator,
                m,tmp,ffTmp,ffTmp2,ffTmp3,
                matmulCache,matmulCacheCalced,
                multTable,addTable);
        if( curGen != -1 ){
            if(toTestIndicator[curGen] == true){
                genCounts[curGen]++;
                appendToEnd(roots[curGen], ff, m);
                foundCounter++;
                toTestIndicator[curGen] = false;
                /*printf("found curGen=%i",curGen);*/
                /*printFFElem("         ",ff);*/
            }
            if(foundCounter == decompCount) break;
        }
        //generate next element
        // (for sure there is a more efficient method)
        ff->el[0] += 1;
        if( ff->el[0] == charac ){
            for(i=0; i<m-1 && ff->el[i]==charac; i++){
                ff->el[i] = 0;
                ff->el[i+1] += 1;
            }
            if( ff->el[m-1] == charac )
                break;
        }
        updateFFElem(ff,m);
    }
    if( foundCounter != decompCount ){
        printf("BAAAD ERROR!!! foundCounter=%i < decompCount=%i\n",
                foundCounter,decompCount);
        exit(0);
    }
    printf("finding time: %.2f\n", (double)(time(NULL)-TIME));
    //------------------------------------------------------------------------
    


    // Process found elements ------------------------------------------------
    int curDecompPosition = 0;
    for(i=0;i<decompCount-1;i++){
        calcSubmoduleElements(roots[i], mipo,
                polys[curDecompPosition]->lenPoly-1, // *** == maxLenPoly
                genCounts,i, // *** i == curGen
                polys+curDecompPosition, polysCountPerDecomp[i],
                evalToZero+curDecompPosition,
                mats, matLen, frobPowers+curDecompPosition,
                elementsF,
                m,q,tmp,
                ffTmp,ffTmp2,ffTmp3,ffTmp4,
                matmulCache,matmulCacheCalced,
                multTable,addTable);
        curDecompPosition += polysCountPerDecomp[i];
    }
    printf("all not last time: %.2f\n", (double)(time(NULL)-TIME));
    //------------------------------------------------------------------------
    
    // Process last Decomposition --------------------------------------------
    int curGen = decompCount-1;
    unsigned long long pcn = 
        processLastSubmoduleAndTestPrimitivity(roots,mipo,decompCount,
            polys[curDecompPosition]->lenPoly-1, // *** == maxLenPoly
            genCounts,
            polys+curDecompPosition,polysCountPerDecomp[curGen],
            evalToZero+curDecompPosition,
            mats,matLen,frobPowers+curDecompPosition,
            elementsF,
            m,q,
            barFactors,lenBarFactors,countBarFactors,
            commonBarFactor,lenCommonBarFactor,
            commonBiggestBarFactor,lenCommonBiggestBarFactor,
            matCharac,
            matmulCache,matmulCacheCalced,
            multTable,addTable);
    //------------------------------------------------------------------------

    //free variables
    for(i=0;i<decompCount;i++)
        freeNode(roots[i]);
    free(roots); free(curRoots);

    //free temporary variables
    free(tmp);
    freeFFElem(ff);
    freeFFElem(ffRet);
    freeFFElem(ffTmp);
    freeFFElem(ffTmp2);
    freeFFElem(ffTmp3);
    freeFFElem(ffTmp4);
    for(i=0;i<matLen;i++) freeFFElem(matmulCache[i]);
    free(matmulCache);
    free(matmulCacheCalced);
    free(toTestIndicator);

    
    printf("total time: %.2f\n", (double)(time(NULL)-TIME));
    return pcn;
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
        int m, int *multTable, int *addTable){
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
                    + multTable[ mat1[row]->el[i2]*matTmp[col]->el[i2] ] ];
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
        int *multTable, int *addTable){
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
        powerFFElemInt(ffTmp, mipo, mats[i],
                m,q,tmp,ffTmp2,multTable,addTable);
        /*printf("i=%i",i); printFFElemShort(" ffTmp",ffTmp);*/
        /*printFFElemShort(" ffTmp2",ffTmp2);*/
    }
    //transpose first matrix
    transposeMatrix(mats,matTmp,m);
    for(i=0;i<m;i++){
        copyFFElem(matTmp[i],mats[i]);
    }
    /*printf("first matrix:\n");printFFElemMatrix(mats,m);*/
    for(i=1;i<maxPower;i++){
        multiplyMatrices(mats, mats+(i-1)*m, mats+i*m, matTmp,
                m,multTable,addTable);
    }
    for(i=0;i<m;i++) freeFFElem(matTmp[i]);
    free(matTmp);
    free(tmp);
    free(ffTmp);
    free(ffTmp2);
    return mats;
}


void main(){
    int multTableRaw[] = {2, 0, 1, 2, 0, 1, 2, 0, 1};
    int initialMultShift = 4;
    int *multTable = multTableRaw+initialMultShift;
    int addTableRaw[] = {2, 0, 1, 2, 0, 1, 2, 0, 1};
    int initialAddShift = 4;
    int *addTable = addTableRaw+initialAddShift;

    struct FFElem **matCharac = malloc(3*sizeof(struct FFElem*));
    matCharac[0] = malloc(sizeof(struct FFElem));
    matCharac[0]->el = (int[]){1,2,1};
    matCharac[0]->idcs = (int[]){2,1,0};
    matCharac[0]->len = 3;
    matCharac[1] = malloc(sizeof(struct FFElem));
    matCharac[1]->el = (int[]){0,1,1};
    matCharac[1]->idcs = (int[]){2,1,0};
    matCharac[1]->len = 2;
    matCharac[2] = malloc(sizeof(struct FFElem));
    matCharac[2]->el = (int[]){0,0,1};
    matCharac[2]->idcs = (int[]){2,0,0};
    matCharac[2]->len = 1;

    struct FFElem *mipo = malloc(sizeof(struct FFElem));
    mipo->el = (int[]){2, 2, 1, 0, 2, 0, 1};
    mipo->idcs = (int[]){6,4,2,1,0};
    mipo->len = 5;

    struct FFElem *x = malloc(sizeof(struct FFElem));
    x->el = (int[]){0,2,0,1,0,1};
    x->idcs = (int[]){5,3,1,0,0,0};
    x->len = 3;
    struct FFElem *x2 = malloc(sizeof(struct FFElem));
    x2->el = (int[]){0,2,0,1,0,1};
    x2->idcs = (int[]){5,3,1};
    x2->len = 3;

    int m = 6;
    int charac = 3;
    //////////////////////////////////////////////////////////////////////////

    int * tmp = malloc( m*sizeof(int) );

    struct FFElem *ret = mallocFFElem(m);
    
    // Mult test /////////////////////////////////////////////////////////////
    x->el = (int[]){1,0,0,0,0,0};
    x->idcs = (int[]){0,0,0,0,0,0};
    x->len = 1;
    printFFElem("x", x);
    printFFElem("mipo", mipo);
    multiplyFFElem(x,x,ret,mipo,tmp,m,multTable,addTable);
    printFFElem("ret", ret);
    squareFFElem(x,mipo,ret,tmp,m,multTable,addTable);
    printFFElem("x",x);
    
    // Matmul test ///////////////////////////////////////////////////////////
    /*int i;*/
    /*for(i=0;i<m;i++) printFFElem("matCharac",matCharac[i]);*/
    /*printFFElem("x",x);*/
    /*matmul(matCharac, x, ret, m, multTable, addTable);*/
    /*printFFElem("ret",ret);*/
    
    
    // Addition test /////////////////////////////////////////////////////////
    /*printFFElem("x",x);*/
    /*printFFElem("x2", x2);*/
    /*addFFElem(x,x2,ret,multTable,addTable);*/
    /*printFFElem("ret",ret);*/

    // Power test ////////////////////////////////////////////////////////////
    /*int power[] = {1,1,1}; //49*/
    /*int powerLen = 3;*/
    /*x->el = (int[]){0,0,1};*/
    /*x->idcs = (int[]){2,0,0};*/
    /*x->len = 1;*/
    /*int power[] = {1,0,0,1,1,1}; //39*/
    /*int powerLen = 6;*/
    /*printFFElem("x",x);*/
    /*struct FFElem *ffTmp = mallocFFElem(m);*/
    /*powerFFElem(x,mipo,ret,*/
            /*m,power,powerLen,*/
            /*matCharac,tmp,ffTmp,multTable,addTable);*/
    /*powerFFElemSqM(x,mipo,ret,*/
            /*m,power,powerLen,tmp,ffTmp,multTable,addTable);*/
    /*printFFElem("ret",ret);*/
    /*freeFFElem(ffTmp);*/
    
    // Primitivity test //////////////////////////////////////////////////////
    /*int barFactors[] = {2,1};*/
    /*int lenBarFactors[] = {1,1};*/
    /*int countBarFactors = 2;*/
    /*int biggestPrimeFactor[] = {1,1,1};*/
    /*int lenBiggestPrimeFactor = 3;*/
    /*struct FFElem *fff = mallocFFElem(m);*/
    /*struct FFElem *ffTmp = mallocFFElem(m);*/
    /*struct FFElem *ffTmp2 = mallocFFElem(m);*/
    /*struct FFElem *ffTmp3 = mallocFFElem(m);*/
    /*struct FFElem *ffTmp4 = mallocFFElem(m);*/
    /*x->el = (int[]){0,1,2};*/
    /*x->idcs = (int[]){2,1,0};*/
    /*x->len = 2;*/

    /*bool b = isPrimitive(x,mipo,m,*/
            /*barFactors,lenBarFactors,countBarFactors,*/
            /*biggestPrimeFactor,lenBiggestPrimeFactor,*/
            /*matCharac,*/
            /*fff,ffTmp,ffTmp2,ffTmp3,ffTmp4,*/
            /*tmp,multTable,addTable);*/
    /*printFFElem("x",x);*/
    /*printf(b?"is primitive!":"is not primitive!");*/
    /*freeFFElem(fff);*/
    /*freeFFElem(ffTmp);*/
    /*freeFFElem(ffTmp2);*/
    /*freeFFElem(ffTmp3);*/
    /*freeFFElem(ffTmp4);*/

    // Generate Matrices test ////////////////////////////////////////////////
    /*struct FFElem *ffTmp = mallocFFElem(m);*/
    /*struct FFElem *ffTmp2 = mallocFFElem(m);*/
    /*int maxPower = 4;*/
    /*struct FFElem **mats = genFrobMats(mipo,m, maxPower, 3, */
            /*multTable,addTable);*/
    /*int i;*/
    /*printf("mats = {\n");*/
    /*for(i=0;i<maxPower;i++){*/
        /*printFFElemMatrix(mats+i*m,m);*/
        /*printf("\n");*/
    /*}*/
    /*printf("}\n");*/
    /*freeFFElem(ffTmp);*/
    /*freeFFElem(ffTmp2);*/
    /*for(i=0;i<maxPower*m;i++)*/
        /*freeFFElem(mats[i]);*/
    /*free(mats);*/

    // Frobenius test ///////////////////////////////////////////////////////
    /*x->el = (int[]){2,0,2};*/
    /*x->idcs = (int[]){2,0,0};*/
    /*x->len = 2;*/
    /*struct FFPoly *poly = malloc(sizeof(struct FFPoly));*/
    /*poly->poly = malloc(3*sizeof(struct FFElem*));*/
    /*poly->poly[0] = malloc(sizeof(struct FFElem));*/
    /*poly->poly[0]->el = (int[]){1,0,0};*/
    /*poly->poly[0]->idcs = (int[]){0,0,0};*/
    /*poly->poly[0]->len = 1;*/
    /*poly->poly[1] = malloc(sizeof(struct FFElem));*/
    /*poly->poly[1]->el = (int[]){1,0,0};*/
    /*poly->poly[1]->idcs = (int[]){0,0,0};*/
    /*poly->poly[1]->len = 1;*/
    /*poly->poly[2] = malloc(sizeof(struct FFElem));*/
    /*poly->poly[2]->el = (int[]){1,1,0};*/
    /*poly->poly[2]->idcs = (int[]){1,0,0};*/
    /*poly->poly[2]->len = 2;*/
    /*poly->lenPoly = 3;*/
    /*int frobPower = 1;*/
    /*struct FFElem *fff = mallocFFElem(m);*/

    //cache
    /*struct FFElem **matmulCache = malloc(maxPower*sizeof(struct FFElem*));*/
    /*bool *matmulCacheCalced = malloc(maxPower*sizeof(bool));*/
    /*for(i=0;i<maxPower;i++){*/
        /*matmulCache[i] = mallocFFElem(m);*/
        /*matmulCacheCalced[i] = false;*/
    /*}*/

    /*applyFrob_noCache(x,mipo,poly,lenPoly,*/
            /*mats,frobPower,fff,m,tmp,ffTmp,ffTmp2,*/
            /*multTable,addTable);*/
    /*applyFrob(x,mipo,poly,*/
            /*mats,frobPower,fff,m,tmp,ffTmp,ffTmp2,*/
            /*matmulCache,matmulCacheCalced,*/
            /*multTable,addTable);*/
    /*printFFElem("fff",fff);*/

    /*for(i=0;i<maxPower;i++) freeFFElem(matmulCache[i]);*/
    /*free(matmulCache);*/
    /*free(matmulCacheCalced);*/

    /*freeFFElem(fff);*/
    /*freeFFElem(ffTmp);*/
    /*freeFFElem(ffTmp2);*/
    /*for(i=0;i<maxPower*m;i++)*/
        /*freeFFElem(mats[i]);*/
    /*free(mats);*/
    /*free(poly->poly[0]);*/
    /*free(poly->poly[1]);*/
    /*free(poly->poly[2]);*/
    /*free(poly->poly);*/
    /*free(poly);*/
    /////////////////////////////////////////////////////////////////////////
    free(tmp);
    free(x);
    free(x2);
    free(mipo);
    freeFFElem(ret);
    free(matCharac[0]); free(matCharac[1]); free(matCharac[2]); free(matCharac);
}
