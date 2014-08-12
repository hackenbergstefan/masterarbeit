#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

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
    for(i=0;i<m;i++){
        ret[i] = 0;
        for(j=0;j<m;j++){
            ret[i] += mat[i*m + j]*vec[j];
        }
        ret[i] %= charac;
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

inline void multiplyPoly(int *p1, int m1, int *p2, int m2, int *p3, int m3, int charac) {
    int i,j;
    initPoly(p3,m3);
    for(i=0;i<m1;i++)
        for(j=0;j<m2;j++)
            p3[i+j] = p3[i+j]+p1[i]*p2[j];
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
    int quo=0;
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
        quo = p1[i+degmod]*modInv(mod[degmod],charac);
        //printf("i=%i p1[i+degmod]=%i mod[degmod]=%i mod[degmod]^(-1)=%i ",i,p1[i+degmod],mod[degmod],modInv(mod[degmod],charac));
        //printf("quo=%i\n", quo);
        for(j=degmod;j>=0;j--){
            p1[i+j] = p1[i+j] - mod[j]*quo;
        }
        //printArr(p1,m1);
    }
    for(i=0;i<m1;i++){
        p1[i] %= charac;
    }
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

inline void testPolys(int *x, int *x_mipo, int decompCount,
        int *polys, int *polysLen, int *polysCount, bool *evalToZero,
        int *mats, int *frobPowers, 
        int *ret, int m, int charac, int *tmp, int *tmp2){
    int i,j;
    int curDecompPosition = 0;
    int curPolyPosition = 0;
    int lastZeroPoly = 0;
    int goodCounter = 0;
    /*printf("polysCount="); printArr(polysCount,decompCount);*/

    for(i=0;i<decompCount;i++){
        /*printf("decomp: i=%i\n",i);*/
        goodCounter = 0;
        for(j=0;j<polysCount[i];j++){
            /*printf("\ttest poly j=%i",j);*/
            applyFrob(x, x_mipo, 
                    polys+curPolyPosition, polysLen[curDecompPosition+j],
                    mats, frobPowers[curDecompPosition+j], 
                    ret, m, charac, tmp, tmp2);
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

char *processFFElements( int *x_mipo, int decompCount,
        int *polys, int *polysLen, int *polysCount, bool *evalToZero,
        int *mats, int *frobPowers, 
        int *genCounts, int m, int charac, int shiftSize){
    time_t TIME = time(NULL);
    int i,j;
    
    
    int * x = malloc( m*sizeof(int) );
    int * ret = malloc( m*sizeof(int) );
    int * tmp = malloc( m*sizeof(int) );
    int * tmp2 = malloc( 2*m*sizeof(int) );

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

    /*printf("begin testing %f\n", (double)(time(NULL)-TIME));*/
    
    while(1==1){
        testPolys(x,x_mipo,decompCount,
                polys,polysLen,polysCount,evalToZero,
                mats, frobPowers,
                ret, m, charac, tmp, tmp2);
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
    /*printf("allElements processed %f\n", (double)(time(NULL)-TIME));*/
    /*printf("genCounts = "); printArr(genCounts, decompCount);*/
    /*int prod = 1;*/
    /*for(i=0;i<decompCount;i++)*/
        /*prod *= genCounts[i];*/
    /*printf(" => CN = %i\n",prod);*/

    //combine elements
    /*struct Node *combinedElements = malloc(sizeof(struct Node));*/
    /*struct Node *curCombinedEnd = combinedElements;*/
    /*combinedElements->x = 0;*/
    /*combinedElements->next = 0;*/

    char *filepath = malloc(50*sizeof(char));
    sprintf(filepath,"tmp_cnSearch_%i",time(NULL));
    /*printf("filepath = %s",filepath);*/
    FILE * fp = fopen(filepath,"a");

    for(i=0;i<decompCount;i++){
        curRoots[i] = roots[i];
    }
    // debug ---
    /*int * y = malloc(m*sizeof(int));*/
    // --- debug
    while(1==1){
        /*printf("combine Element of:\n\t");*/
        decodeArr(curRoots[0]->x, x, m, shiftSize);
        /*decodeArr(curRoots[0]->x, y, m, shiftSize);*/
        /*printArr(y,m);*/
        for(i=1;i<decompCount;i++){
            decodeArrAdd(curRoots[i]->x, x, m, shiftSize);
            /*decodeArr(curRoots[i]->x, y, m, shiftSize);*/
            /*printf("\t"); printArr(y,m);*/
        }
        for(i=0;i<m;i++)
            x[i] %= charac;
        /*printf("    x = ");printArr(x,m);*/
        moduloPoly(x,m,x_mipo,m+1,charac);
        /*printf("append to end x="); printArr(x,m);*/
        /*curCombinedEnd = appendToEnd(curCombinedEnd,x,m,shiftSize);*/
        for(i=0;i<m;i++)
            fprintf(fp, "%i ", x[i]);
        fprintf(fp,"\n");


        //nextElement
        curRoots[0] = curRoots[0]->next;
        if( curRoots[0] == 0 ){
            for(i=0;i<decompCount-1 && curRoots[i]==0;i++){
                curRoots[i] = roots[i];
                curRoots[i+1] = curRoots[i+1]->next;
            }
        }
        if( curRoots[decompCount-1] == 0){
            /*printf("end of while reached\n");*/
            break;
        }
    }
    fclose(fp);
    // debug ---
    /*free(y);*/
    // --- debug

    /*printf("end reached!\n");*/

    
    /*for(i=0;i<decompCount;i++){*/
        /*printf("free curRoot %i\n",i);*/
        /*freeNode(curRoots[i]);*/
    /*}*/
    /*free(curRoots);*/
    /*free(curRoots);*/
    /*printf("curRoots freed!\n");*/

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
    
    printf("C time: %.2f\n", (double)(time(NULL)-TIME));
    return filepath;
}


double eta_processFFElements( int *x_mipo, int decompCount,
        int *polys, int *polysLen, int *polysCount, bool *evalToZero,
        int *mats, int *frobPowers, 
        int *genCounts, int m, int charac, int shiftSize){
    struct timeval TIME1, TIME2;
    int i,j;
    
    
    int * x = malloc( m*sizeof(int) );
    int * ret = malloc( m*sizeof(int) );
    int * tmp = malloc( m*sizeof(int) );
    int * tmp2 = malloc( 2*m*sizeof(int) );

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
        testPolys(x,x_mipo,decompCount,
                polys,polysLen,polysCount,evalToZero,
                mats, frobPowers,
                ret, m, charac, tmp, tmp2);
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
    
    double timediff = (TIME2.tv_sec - TIME1.tv_sec +
         ((double)(TIME2.tv_usec - TIME1.tv_usec))/1000000.0);
    return 2*timediff *pow((double)charac,(double)m)/counter;
}


/**
 * returns the next CN element
 */
void findNextCN(int * x, int *x_mipo, 
        int *polys, int *polysLen, int *polysCount, bool *evalToZero,
        int *mats, int *frobPowers, 
        int m, int charac,
        int *ret, int *tmp, int * tmp2){
    time_t TIME = time(NULL);
    int i,j;
    
    
    /*int * ret = malloc( m*sizeof(int) );*/
    /*int * tmp = malloc( m*sizeof(int) );*/
    /*int * tmp2 = malloc( 2*m*sizeof(int) );*/

    
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
            return;
        }
    }

    /*printf("C time: %.2f\n", (double)(time(NULL)-TIME));*/
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

    int *genCounts = malloc(decompCount*sizeof(int));

    /*char * filepath*/
        /*= processFFElements(xmipo, decompCount,*/
            /*polys, polysLen, polysCount, evalToZero,*/
            /*mats, frobPowers,*/
            /*genCounts, m, charac, shiftSize);*/
    free(genCounts);
    /*printf("filepath = %s",filepath);*/


}
