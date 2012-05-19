#include "redis.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>
#include <float.h>


typedef struct sortedArrayKey {
  double x;
  double y;
  void* value;
} sortedArrayKey;

double* sorted_array_recursive_successor_search(double value, double* start, double* end){
  int numElements = ((int)(end - start) / sizeof(double));
  if (numElements < 16){
    while( !(*start >= value)){
      start++;
    }

    return start;
  }
  double* median = start + (numElements/2);

  if (value == *median){
    return median;
  }
 
  if (value < *median) {
    // recurse on left half of the array.
    double* ret = sorted_array_recursive_successor_search(value, start, median);
    if (*ret == DBL_MAX) {
      return median;
    }
    return ret;
  } else {
    double* ret = sorted_array_recursive_successor_search(value, median, end);
    return ret;
  }
}

//lexical comparison on x,y pairs, imitating the key concatenation thing.
int compareSortedArrayKeyByX (const void * a, const void * b)
{
  if  (((sortedArrayKey*)a)->x < ((sortedArrayKey*)b)->x ) return -1;
  else if (((sortedArrayKey*)a)->x > ((sortedArrayKey*)b)->x ) return 1;
  else return 0;
}

int compareSortedArrayKeyByY (const void * a, const void * b)
{
  if (((sortedArrayKey*)a)->y < ((sortedArrayKey*)b)->y ) return -1;
  else if (((sortedArrayKey*)a)->y > ((sortedArrayKey*)b)->y ) return 1;
  else return 0;
}

int resizeAmount = 10000;
int allElementsCount = 0;
int currentAllElementsArraySize = 0;

sortedArrayKey  *allElementsByX = NULL;
sortedArrayKey  *allElementsByY = NULL;
double *allX = NULL;
double *allY = NULL;


void sortedArray2DAddCommand(redisClient *c) {
  double x = strtod(c->argv[1]->ptr, NULL);
  double y = strtod(c->argv[2]->ptr, NULL);

  sortedArrayKey key;
  key.x = x;
  key.y = y;
  key.value = NULL; 

  if (allElementsCount == currentAllElementsArraySize) {
    currentAllElementsArraySize += resizeAmount;
    allElementsByX = (sortedArrayKey*) zrealloc((void*)allElementsByX, sizeof(sortedArrayKey)*currentAllElementsArraySize);
    allX = (double *) zrealloc((void*)allX, sizeof(double)*currentAllElementsArraySize);
    allElementsByY = (sortedArrayKey*) zrealloc((void*)allElementsByY, sizeof(sortedArrayKey)*currentAllElementsArraySize);
    allY = (double *) zrealloc((void*)allY, sizeof(double)*currentAllElementsArraySize);
  }

  allElementsByX[allElementsCount] = key;
  allElementsByY[allElementsCount] = key;
  allX[allElementsCount] = key.x;
  allY[allElementsCount] = key.y;

  allElementsCount++;
  addReplyLongLong(c, allElementsCount);
}


int compare(const void * a, const void * b){
  if (*(double *) a > *(double *) b) return 1;
  else if (*(double *) a < *(double *) b) return -1;
  else return 0;
}

void sortedArray2DBuildCommand(redisClient *c) {
  qsort(allElementsByX, allElementsCount, sizeof(sortedArrayKey), compareSortedArrayKeyByX);
  qsort(allX, allElementsCount, sizeof(double), compare);
  qsort(allElementsByY, allElementsCount, sizeof(sortedArrayKey), compareSortedArrayKeyByY);
  qsort(allY, allElementsCount, sizeof(double), compare);
  addReplyLongLong(c, 42);
}

void sortedArray2DStupidSearchCommand(redisClient *c) {
  double x1 = strtod(c->argv[1]->ptr, NULL);
  double x2 = strtod(c->argv[2]->ptr, NULL);
  double y1 = strtod(c->argv[3]->ptr, NULL);
  double y2 = strtod(c->argv[4]->ptr, NULL);

  int count = 0;
  int i = 0;
  for (i = 0; i < allElementsCount; i++){
    if (allElementsByX[i].x >= x1 && allElementsByX[i].x <= x2 && allElementsByX[i].y >= y1
	&& allElementsByX[i].y <= y2){
      count++;
    }
  }
  addReplyLongLong(c, count);
}

void sortedArray2DSearchCommand(redisClient *c) {
  double x1 = strtod(c->argv[1]->ptr, NULL);
  double x2 = strtod(c->argv[2]->ptr, NULL);
  double y1 = strtod(c->argv[3]->ptr, NULL);
  double y2 = strtod(c->argv[4]->ptr, NULL);

  double *start;
  int count = 0;

  //if X is narrower, then search by x first
  if (x2 - x1 < y2 - y1){
    printf("byX\n");
    start = sorted_array_recursive_successor_search(x1, allX, allX + allElementsCount);
    int i = start - allX;
    for ( ; i < allElementsCount && allElementsByX[i].x <= x2; i++){
      if (allElementsByX[i].y >= y1 && allElementsByX[i].y <= y2) count++;
    }
  } else {
    start = sorted_array_recursive_successor_search(y1, allY, allY + allElementsCount);
    int i = start - allY;
    for ( ; i < allElementsCount && allElementsByY[i].y <= y2; i++){
      if (allElementsByY[i].x >= x1 && allElementsByY[i].x <= x2) count++;
    }
  }

    addReplyLongLong(c, count);
}
