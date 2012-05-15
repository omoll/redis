#include "redis.h"

#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>
// Each leaf in the tree contains this many elements.
// This must be greater than 2, otherwise splitting range
// using the median may result in infinite recursion.
static int QUAD_TREE_NODE_SIZE = 64;

double tfkLayeredRangeTree_get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

/*
  Type Defs
*/
typedef struct layeredRangeTreeKey{
  double x;
  double y;
  void* value;
} layeredRangeTreeKey;


// the inner level of the tree.
typedef struct layeredRangeTreeNodeLevel2 {
  // The rectangle corresponding to this node.
  double x1;
  double x2;
  double y1; 
  double y2;

  // a count of elments in the subtree for use in rebalancing.
  int elementCount;

  // TODO(tfkLayeredRangeTree): start is almost always 0.
  // start and end of range in elements array.
  int start;
  int end;
  layeredRangeTreeKey* elements; // undefined if not leaf.

  struct layeredRangeTreeNodeLevel2** children; // NULL if leaf
} layeredRangeTreeNodeLevel2;

// the top level of the tree.
typedef struct layeredRangeTreeNodeLevel1 {
  // The rectangle corresponding to this node.
  double x1;
  double x2;
  double y1; 
  double y2;

  // a count of elments in the subtree for use in rebalancing.
  int elementCount;

  // TODO(tfkLayeredRangeTree): start is almost always 0.
  // start and end of range in elements array.
  int start;
  int end;
  layeredRangeTreeKey* elements; // undefined if not leaf.

  struct layeredRangeTreeNodeLevel1** children; // NULL if leaf

  // a pointer to the second level of the tree containing the same points as in this subtree
  // sorted on the other coordinate.
  layeredRangeTreeNodeLevel2* secondLevelPointer;
} layeredRangeTreeNodeLevel1;


void layeredRangeTreeDeleteNodeLevel2(layeredRangeTreeNodeLevel2* n) {
  if (n->children == NULL) {
    zfree(n->elements);
    zfree(n);
  } else {
    layeredRangeTreeDeleteNodeLevel2(n->children[0]);
    layeredRangeTreeDeleteNodeLevel2(n->children[1]);
    zfree(n->children);
    zfree(n);
  }
}

// Delete node n and it's associated subree.
void layeredRangeTreeDeleteNodeLevel1(layeredRangeTreeNodeLevel1* n){
  if (n->children == NULL) {
    zfree(n->elements);
    zfree(n);
  } else {
    layeredRangeTreeDeleteNodeLevel1(n->children[0]);
    layeredRangeTreeDeleteNodeLevel1(n->children[1]);
    zfree(n->children);
    zfree(n);
  }

  if (n->secondLevelPointer != NULL){
    layeredRangeTreeDeleteNodeLevel2(n->secondLevelPointer);
  }
}

/*
  Comparison functions for median finding and sorting.
*/
int layeredRangeTreeElementCompareX(const void* a, const void* b){
  if (((layeredRangeTreeKey*) a)->x < ((layeredRangeTreeKey*) b)->x) {
    return -1;
  }
  if (((layeredRangeTreeKey*) a)->x > ((layeredRangeTreeKey*) b)->x) {
    return 1;
  }
  return 0;
}

int layeredRangeTreeElementCompareY(const void* a, const void* b){
  if (((layeredRangeTreeKey*) a)->y < ((layeredRangeTreeKey*) b)->y) {
    return -1;
  }
  if (((layeredRangeTreeKey*) a)->y > ((layeredRangeTreeKey*) b)->y) {
    return 1;
  }
  return 0;
}

void buildLayeredRangeTreeNodeLevel2 (layeredRangeTreeNodeLevel2* n, layeredRangeTreeKey* elements, int elementCount) {
  n->elementCount = elementCount;
  // Coarsen 
  if (n->elementCount <= QUAD_TREE_NODE_SIZE) {
    n->start = 0; 
    n->end = n->elementCount;
    n->elements = (layeredRangeTreeKey*) zmalloc(sizeof(layeredRangeTreeKey) * QUAD_TREE_NODE_SIZE);
    for (int i = 0; i < n->elementCount; i++) {
      n->elements[i] = elements[i];
    }
    return;
  }
 
  n->children = (layeredRangeTreeNodeLevel2**) zmalloc(sizeof(layeredRangeTreeNodeLevel2*)*2);

  qsort(elements, n->elementCount,
      sizeof(layeredRangeTreeKey), layeredRangeTreeElementCompareY);
  layeredRangeTreeKey medianY = elements[n->elementCount / 2];

  n->children[0] = (layeredRangeTreeNodeLevel2*) zmalloc(sizeof(layeredRangeTreeNodeLevel2));
  n->children[1] = (layeredRangeTreeNodeLevel2*) zmalloc(sizeof(layeredRangeTreeNodeLevel2));

  n->children[0]->elementCount = n->elementCount / 2 + n->elementCount % 2;
  n->children[1]->elementCount = n->elementCount - n->children[0]->elementCount; 

  n->children[0]->children = NULL;
  n->children[1]->children = NULL;

  n->children[0]->x1 = n->x1;
  n->children[0]->x2 = n->x2;
  n->children[0]->y1 = n->y1;
  n->children[0]->y2 = medianY.y;
  n->children[0]->start = 0;
  n->children[0]->end = 0;

  n->children[1]->x1 = n->x1;
  n->children[1]->x2 = n->x2;
  n->children[1]->y1 = medianY.y;
  n->children[1]->y2 = n->y2;
  n->children[1]->start = 0;
  n->children[1]->end = 0;

  buildLayeredRangeTreeNodeLevel2(n->children[0], elements, n->children[0]->elementCount);
  buildLayeredRangeTreeNodeLevel2(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount); 
}

// builds a second level of the quad tree.
void explodeTree(layeredRangeTreeNodeLevel1* n, layeredRangeTreeKey* elements, int elementCount){
  n->secondLevelPointer = NULL;
  if (n->elementCount <= QUAD_TREE_NODE_SIZE) {
    assert(elementCount == n->elementCount);
    // this is a leaf node, copy its elements into the array.
    memcpy(elements, n->elements, elementCount * sizeof(layeredRangeTreeKey));
    return;
  }
 
  explodeTree(n->children[0], elements, n->children[0]->elementCount); 
  explodeTree(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount);

  n->secondLevelPointer = (layeredRangeTreeNodeLevel2*) zmalloc(sizeof(layeredRangeTreeNodeLevel2));

  n->secondLevelPointer->x1 = n->x1;
  n->secondLevelPointer->x2 = n->x2;
  n->secondLevelPointer->y1 = n->y1;
  n->secondLevelPointer->y2 = n->y2;
  n->secondLevelPointer->children = NULL;
  n->secondLevelPointer->elements = NULL;
  n->secondLevelPointer->elementCount = n->elementCount;
  n->secondLevelPointer->start = 0;
  n->secondLevelPointer->end = 0;
  
  buildLayeredRangeTreeNodeLevel2(n->secondLevelPointer, elements, n->elementCount);
  assert(n->secondLevelPointer != NULL);
}

void buildLayeredRangeTreeNodeLevel1(layeredRangeTreeNodeLevel1* n, layeredRangeTreeKey* elements, int elementCount) {
  n->elementCount = elementCount;
  n->secondLevelPointer = NULL;

  // Coarsen 
  if (n->elementCount <= QUAD_TREE_NODE_SIZE) {
    n->start = 0; 
    n->end = n->elementCount;
    n->elements = (layeredRangeTreeKey*) zmalloc(sizeof(layeredRangeTreeKey) * QUAD_TREE_NODE_SIZE);
    for (int i = 0; i < n->elementCount; i++) {
      n->elements[i] = elements[i];
    }
    return;
  }
 
  n->children = (layeredRangeTreeNodeLevel1**) zmalloc(sizeof(layeredRangeTreeNodeLevel1*)*2);

  qsort(elements, n->elementCount,
      sizeof(layeredRangeTreeKey), layeredRangeTreeElementCompareX);

  layeredRangeTreeKey medianX = elements[n->elementCount / 2];

  n->children[0] = (layeredRangeTreeNodeLevel1*) zmalloc(sizeof(layeredRangeTreeNodeLevel1));
  n->children[1] = (layeredRangeTreeNodeLevel1*) zmalloc(sizeof(layeredRangeTreeNodeLevel1));

  n->children[0]->elementCount = n->elementCount / 2 + n->elementCount % 2;
  n->children[1]->elementCount = n->elementCount - n->children[0]->elementCount; 

  n->children[0]->children = NULL;
  n->children[1]->children = NULL;
  
  n->children[0]->x1 = n->x1;
  n->children[0]->x2 = medianX.x;
  n->children[0]->y1 = n->y1;
  n->children[0]->y2 = n->y2;
  n->children[0]->start = 0;
  n->children[0]->end = 0;

  n->children[1]->x1 = medianX.x;
  n->children[1]->x2 = n->x2;
  n->children[1]->y1 = n->y1;
  n->children[1]->y2 = n->y2;
  n->children[1]->start = 0;
  n->children[1]->end = 0;

  buildLayeredRangeTreeNodeLevel1(n->children[0], elements, n->children[0]->elementCount);
  buildLayeredRangeTreeNodeLevel1(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount); 
}

// returns true if a key falls within a given range.
bool layeredRangeTreeKeyInRange(layeredRangeTreeKey key, double x1, double x2, double y1, double y2) {
  if (key.x >= x1 && key.x < x2 && key.y > y1 && key.y <= y2) {
    return true;
  } else {
    return false;
  }
}

// returns true if a node intersects with the given range, false otherwise.
bool layeredRangeTreeNodeIntersectsRangeLevel1(layeredRangeTreeNodeLevel1* n, double x1, double x2, double y1, double y2) {
  // ranges lies to the left.
  if (x1 > n->x1 && x1 > n->x2 && x2 > n->x1 && x2 > n->x2) {
    return false;
  }
  if (x1 < n->x1 && x1 < n->x2 && x2 < n->x1 && x2 < n->x2) {
    return false;
  }
 
  if (y1 > n->y1 && y1 > n->y2 && y2 > n->y1 && y2 > n->y2) {
    return false;
  }
  if (y1 < n->y1 && y1 < n->y2 && y2 < n->y1 && y2 < n->y2) {
    return false;
  }
 
  // otherwise this range intersects with the range a.w. n.
  return true;
}

// returns true if a node intersects with the given range, false otherwise.
bool layeredRangeTreeNodeIntersectsRangeLevel2(layeredRangeTreeNodeLevel2* n, double x1, double x2, double y1, double y2) {
  // ranges lies to the left.
  if (x1 > n->x1 && x1 > n->x2 && x2 > n->x1 && x2 > n->x2) {
    return false;
  }
  if (x1 < n->x1 && x1 < n->x2 && x2 < n->x1 && x2 < n->x2) {
    return false;
  }
 
  if (y1 > n->y1 && y1 > n->y2 && y2 > n->y1 && y2 > n->y2) {
    return false;
  }
  if (y1 < n->y1 && y1 < n->y2 && y2 < n->y1 && y2 < n->y2) {
    return false;
  }
 
  // otherwise this range intersects with the range a.w. n.
  return true;
}


void layeredRangeTreeRangeSearchLevel2(layeredRangeTreeNodeLevel2 *n, double x1, double x2, double y1, double y2, int* count) {
  if (n->children == NULL) {
    assert(n->elementCount <= QUAD_TREE_NODE_SIZE);
    for (int i = n->start; i < n->end; i++) {
      if (layeredRangeTreeKeyInRange(n->elements[i], -200, 200, y1, y2)){
        (*count)++;
        //printf("Element (%f, %f) in range \n", n->elements[i].x, n->elements[i].y);
      }
    }
    return;
  }
 
  // call recursively on each child intersected by the range.
  if (layeredRangeTreeNodeIntersectsRangeLevel2(n->children[0], -200, 200, y1, y2)) {
    layeredRangeTreeRangeSearchLevel2(n->children[0], x1, x2, y1, y2, count);
  }

  if (layeredRangeTreeNodeIntersectsRangeLevel2(n->children[1], -200, 200, y1, y2)) {
    layeredRangeTreeRangeSearchLevel2(n->children[1], x1, x2, y1, y2, count);
  }
}

void layeredRangeTreeRangeSearchLevel1(layeredRangeTreeNodeLevel1 *n, double x1, double x2, double y1, double y2, int* count) {
  if (n->children == NULL) {
    assert(n->elementCount <= QUAD_TREE_NODE_SIZE);
    // if we get to a leaf just scan through the elements and report the ones in the range.
    for (int i = n->start; i < n->end; i++) {
      if (layeredRangeTreeKeyInRange(n->elements[i], x1, x2, y1, y2)){
        (*count)++;
        //printf("Element (%f, %f) in range \n", n->elements[i].x, n->elements[i].y);
      }
    }
    return;
  }

  assert(n->secondLevelPointer != NULL);

  // call the second level search if we're entirely within the x-range.
  if (n->x1 >= x1 && n->x2 <= x2) {
    layeredRangeTreeRangeSearchLevel2(n->secondLevelPointer, -200, 200, y1, y2, count);
    return;
  }

  // otherwise call recursively on each child intersected by the range.
  if (layeredRangeTreeNodeIntersectsRangeLevel1(n->children[0], x1, x2, -200, 200)) {
    layeredRangeTreeRangeSearchLevel1(n->children[0], x1, x2, y1, y2, count);
  }

  if (layeredRangeTreeNodeIntersectsRangeLevel1(n->children[1], x1, x2, -200, 200)) {
    layeredRangeTreeRangeSearchLevel1(n->children[1], x1, x2, y1, y2, count);
  }
}


// the root node
layeredRangeTreeNodeLevel1* rangeTreeRoot;

// collection of elements
layeredRangeTreeKey* allElements;
int allElementsCount;
int currentAllElementsArraySize;
double totalQueryTime;

void tfkLayeredRangeTree2DRangeSearchCommand(redisClient* c) {
  double x1 = strtod(c->argv[1]->ptr, NULL);
  double x2 = strtod(c->argv[2]->ptr, NULL);
  double y1 = strtod(c->argv[3]->ptr, NULL);
  double y2 = strtod(c->argv[4]->ptr, NULL);
  int count = 0;
  double start_time = tfkLayeredRangeTree_get_time();
  layeredRangeTreeRangeSearchLevel1(rangeTreeRoot, x1, x2, y1, y2, &count);
  double end_time = tfkLayeredRangeTree_get_time();
  totalQueryTime += end_time - start_time;
  printf("total query time %f rangeTreeCount: %d \n", totalQueryTime, count); 
  addReplyLongLong(c, count);
}

void tfkLayeredRangeTreeBuildTreeCommand(redisClient *c) {
  if (rangeTreeRoot != NULL){
    // delete the old quad tree.
    layeredRangeTreeDeleteNodeLevel1(rangeTreeRoot);
  }
  rangeTreeRoot = (layeredRangeTreeNodeLevel1*) zmalloc(sizeof(layeredRangeTreeNodeLevel1));

  rangeTreeRoot->x1 = -200;
  rangeTreeRoot->x2 = 200;
  rangeTreeRoot->y1 = -200;
  rangeTreeRoot->y2 = 200; 
  rangeTreeRoot->children = NULL;
  rangeTreeRoot->elementCount = allElementsCount;
  rangeTreeRoot->secondLevelPointer = NULL;
  // rebuild using the allElements array.

  // build the first level of the tree.
  buildLayeredRangeTreeNodeLevel1(rangeTreeRoot, allElements, allElementsCount);

  // build the second level of the tree.
  layeredRangeTreeKey* elementsCopy = (layeredRangeTreeKey*) zmalloc(sizeof(layeredRangeTreeKey) * rangeTreeRoot->elementCount);
  explodeTree(rangeTreeRoot, elementsCopy, allElementsCount);
  zfree(elementsCopy); 
  
  //walkLayeredRangeTree(rangeTreeRoot);

  addReplyLongLong(c, 42);
}

/* This generic command implements both ZADD and ZINCRBY. */
void tfkLayeredRangeTreeAddGenericCommand(redisClient *c, int incr) {
    double x = strtod(c->argv[1]->ptr, NULL);
    double y = strtod(c->argv[2]->ptr, NULL);
    
    int resizeAmount = 10000;
  
    // create a layeredRangeTreeKey
    layeredRangeTreeKey key;
    key.x = x;
    key.y = y;
    key.value = NULL; 
    if (allElementsCount == currentAllElementsArraySize) {
      // resize allElementsArray
      currentAllElementsArraySize += resizeAmount;
      allElements = (layeredRangeTreeKey*) zrealloc((void*)allElements, sizeof(layeredRangeTreeKey) * currentAllElementsArraySize);
    }
    allElements[allElementsCount] = key;
    allElementsCount++;
    addReplyLongLong(c, allElementsCount);
}

void layeredRangeTree2DAddCommand(redisClient *c) {
  tfkLayeredRangeTreeAddGenericCommand(c,0);
}

void layeredRangeTree2DBuildTreeCommand(redisClient *c) {
  tfkLayeredRangeTreeBuildTreeCommand(c);
}

void layeredRangeTree2DRangeSearchCommand(redisClient *c) {
  tfkLayeredRangeTree2DRangeSearchCommand(c);
}

