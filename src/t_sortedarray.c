#include "redis.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>
#include <float.h>

double tfkSortedArray_get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}
double sortedArray_totalQueryTime;
typedef struct sortedArrayKey {
  double x;
  double y;
  void* value;
} sortedArrayKey;


//lexical comparison on x,y pairs, imitating the key concatenation thing.
int compareSortedArrayKey (const void * a, const void * b)
{
  if  (((sortedArrayKey*)a)->x < ((sortedArrayKey*)b)->x ) return -1;
  else if (((sortedArrayKey*)a)->x > ((sortedArrayKey*)b)->x ) return 1;
  else if (((sortedArrayKey*)a)->y < ((sortedArrayKey*)b)->y ) return -1;
  else if (((sortedArrayKey*)a)->y > ((sortedArrayKey*)b)->y ) return 1;
  else return 0;
}


int resizeAmount = 10000;
int allElementsCount = 0;
int currentAllElementsArraySize = 0;
sortedArrayKey  *allElements = NULL;

void sortedArray2DAddCommand(redisClient *c) {
  double x = strtod(c->argv[1]->ptr, NULL);
  double y = strtod(c->argv[2]->ptr, NULL);

  sortedArrayKey key;
  key.x = x;
  key.y = y;
  key.value = NULL; 

  if (allElementsCount == currentAllElementsArraySize) {
    currentAllElementsArraySize += resizeAmount;
    allElements = (sortedArrayKey*) zrealloc((void*)allElements, sizeof(sortedArrayKey)*currentAllElementsArraySize);
  }

  allElements[allElementsCount] = key;
  allElementsCount++;
  addReplyLongLong(c, allElementsCount);
}

void sortedArray2DBuildCommand(redisClient *c) {
  qsort(allElements, allElementsCount, sizeof(sortedArrayKey), compareSortedArrayKey);
  addReplyLongLong(c, 42);
}

void sortedArray2DStupidSearchCommand(redisClient *c) {
  double x1 = strtod(c->argv[1]->ptr, NULL);
  double x2 = strtod(c->argv[2]->ptr, NULL);
  double y1 = strtod(c->argv[3]->ptr, NULL);
  double y2 = strtod(c->argv[4]->ptr, NULL);

  double start = tfkSortedArray_get_time(); 
  //this here does a stupid scan
  int count = 0;
  int i = 0;
  for (i = 0; i < allElementsCount; i++){
    if (allElements[i].x >= x1 && allElements[i].x <= x2 && allElements[i].y >= y1
	&& allElements[i].y <= y2){
      count++;
    }
  }
  double end = tfkSortedArray_get_time();
 
  sortedArray_totalQueryTime += end - start;
  printf("Total query time %f \n",sortedArray_totalQueryTime);
  addReplyLongLong(c, count);

}

void sortedArray2DSearchCommand(redisClient *c) {
  sortedArray2DStupidSearchCommand(c);
}
