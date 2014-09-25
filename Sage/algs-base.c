#include <stdio.h>
#include <stdbool.h>

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
void printCharArr(char *arr, int l){
    int j;
    for(j=0;j<l;j++){
        printf(" %d", arr[j]);
    }
    printf("\n");
}
void printBoolArr(bool *arr, int l){
    int j;
    for(j=0;j<l;j++){
        printf(" %d", arr[j]);
    }
    printf("\n");
}
/**
 * copies arr1 into arr2
 */
inline void copyArray(int *arr1, int *arr2, int m){
    int i;
    for(i=0;i<m;i++) arr2[i] = arr1[i];
}

