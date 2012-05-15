#include "redis.h"

#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <sys/time.h>
#include <sys/resource.h>
// Each leaf in the tree contains this many elements.
// This must be greater than 2, otherwise splitting range
// using the median may result in infinite recursion.
static int QUAD_TREE_NODE_SIZE = 3;
static double REBALANCE_FRACTION = 1.5;

double tfk_get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

/*
  Type Defs
*/
typedef struct quadTreeKey{
  double x;
  double y;
  void* value;
} quadTreeKey;

typedef struct quadTreeNode{
  // The rectangle corresponding to this node.
  double x1;
  double x2;
  double y1; 
  double y2;

  // a count of elments in the subtree for use in rebalancing.
  int elementCount;

  // TODO(tfk): start is almost always 0.
  // start and end of range in elements array.
  int start;
  int end;
  quadTreeKey* elements; // undefined if not leaf.

  /*  We either split horizontally or vertically.

      split == 1    split == 2
         Q1           Q1 Q2  
         Q2

      Bottom and Right sides are non-inclusive ranges.
  */
  int split;
  struct quadTreeNode** children; // NULL if leaf
} quadTreeNode;


// For debugging purposes.
quadTreeKey* getQuadTreeElements(quadTreeNode* a, quadTreeKey* elements){
  if (a->children == NULL) {
    for (int i = a->start; i < a->end; i++) {
      *elements = a->elements[i];
       elements++;
    }
    return elements;
  } else {
    elements = getQuadTreeElements(a->children[0], elements);
    return getQuadTreeElements(a->children[1], elements);
  }
}

// function declarations.
bool quadTreeInsert(quadTreeNode* n, quadTreeKey key);

// Delete node n and it's associated subree.
void quadTreeDeleteNode(quadTreeNode* n){
  if (n->children == NULL) {
    zfree(n->elements);
    zfree(n);
  } else {
    quadTreeDeleteNode(n->children[0]);
    quadTreeDeleteNode(n->children[1]);
    zfree(n->children);
    //zfree(n->elements);
    zfree(n);
  }
}

/*
  Comparison functions for median finding and sorting.
*/
int quadTreeElementCompareX(const void* a, const void* b){
  if (((quadTreeKey*) a)->x < ((quadTreeKey*) b)->x) {
    return -1;
  }
  if (((quadTreeKey*) a)->x > ((quadTreeKey*) b)->x) {
    return 1;
  }
  return 0;
}

int quadTreeElementCompareY(const void* a, const void* b){
  if (((quadTreeKey*) a)->y < ((quadTreeKey*) b)->y) {
    return -1;
  }
  if (((quadTreeKey*) a)->y > ((quadTreeKey*) b)->y) {
    return 1;
  }
  return 0;
}

/*
// Take a leaf and split it.
void splitQuadTreeNode(quadTreeNode* n){
  n->children = (quadTreeNode**) zmalloc(sizeof(quadTreeNode*)*2);

  // do a vertical split
  if (n->split == 1){
    // sort elements based on their y coordinate
    qsort(n->elements, n->end - n->start,
        sizeof(quadTreeKey), quadTreeElementCompareY);
    quadTreeKey medianY = n->elements[(n->end - n->start) / 2];
    n->children[0] = (quadTreeNode*) zmalloc(sizeof(quadTreeNode));
    n->children[1] = (quadTreeNode*) zmalloc(sizeof(quadTreeNode));
    
    n->children[0]->children = NULL;
    n->children[1]->children = NULL;
 
    n->children[0]->elements = (quadTreeKey*) zmalloc(sizeof(quadTreeKey)*QUAD_TREE_NODE_SIZE);
    n->children[1]->elements = (quadTreeKey*) zmalloc(sizeof(quadTreeKey)*QUAD_TREE_NODE_SIZE);

    n->children[0]->elementCount = 0;
    n->children[1]->elementCount = 0;

    n->children[0]->x1 = n->x1;
    n->children[0]->x2 = n->x2;
    n->children[0]->y1 = n->y1;
    n->children[0]->y2 = medianY.y;
    n->children[0]->split = 2;
    n->children[0]->start = 0;
    n->children[0]->end = 0;
    
    n->children[1]->x1 = n->x1;
    n->children[1]->x2 = n->x2;
    n->children[1]->y1 = medianY.y;
    n->children[1]->y2 = n->y2;
    n->children[1]->split = 2;
    n->children[1]->start = 0;
    n->children[1]->end = 0;
  }

  // do a horizontal split
  if (n->split == 2){
    // sort elements based on their y coordinate
    qsort(n->elements, n->end - n->start,
        sizeof(quadTreeKey), quadTreeElementCompareX);
    quadTreeKey medianX = n->elements[(n->end - n->start) / 2];
    n->children[0] = (quadTreeNode*) zmalloc(sizeof(quadTreeNode));
    n->children[1] = (quadTreeNode*) zmalloc(sizeof(quadTreeNode));
  
    n->children[0]->children = NULL;
    n->children[1]->children = NULL;
   

    n->children[0]->elements = (quadTreeKey*) zmalloc(sizeof(quadTreeKey)*QUAD_TREE_NODE_SIZE);
    n->children[1]->elements = (quadTreeKey*) zmalloc(sizeof(quadTreeKey)*QUAD_TREE_NODE_SIZE);

    n->children[0]->x1 = n->x1;
    n->children[0]->x2 = medianX.x;
    n->children[0]->y1 = n->y1;
    n->children[0]->y2 = n->y2;
    n->children[0]->split = 1;
    n->children[0]->start = 0;
    n->children[0]->end = 0;


    n->children[1]->x1 = medianX.x;
    n->children[1]->x2 = n->x2;
    n->children[1]->y1 = n->y2;
    n->children[1]->y2 = n->y2;
    n->children[1]->split = 1;
    n->children[1]->start = 0;
    n->children[1]->end = 0;
  }

  // insert all of this quad tree's elements into its children.
  for (int i = n->start; i < n->end; i++) {
    quadTreeInsert(n, n->elements[i]);
  }

  // now free these elements.
  zfree(n->elements);
  n->elements = NULL;
}
*/
void buildQuadTreeNode(quadTreeNode* n, quadTreeKey* elements, int elementCount, int split) {
  n->elementCount = elementCount;
  n->split = split;

  if (n->elementCount <= QUAD_TREE_NODE_SIZE) {
    n->start = 0; 
    n->end = n->elementCount;
    n->elements = (quadTreeKey*) zmalloc(sizeof(quadTreeKey) * QUAD_TREE_NODE_SIZE);
    for (int i = 0; i < n->elementCount; i++) {
      n->elements[i] = elements[i];
    }
    return;
  }
 
  n->children = (quadTreeNode**) zmalloc(sizeof(quadTreeNode*)*2);
  
  // find the median element of elements.
  if (n->split == 1){
    // we split vertically.
    qsort(elements, n->elementCount,
        sizeof(quadTreeKey), quadTreeElementCompareY);
    quadTreeKey medianY = elements[n->elementCount / 2];

    n->children[0] = (quadTreeNode*) zmalloc(sizeof(quadTreeNode));
    n->children[1] = (quadTreeNode*) zmalloc(sizeof(quadTreeNode));

    n->children[0]->elementCount = n->elementCount / 2 + n->elementCount % 2;
    n->children[1]->elementCount = n->elementCount - n->children[0]->elementCount; 

    n->children[0]->children = NULL;
    n->children[1]->children = NULL;

    //assert(n->y1 <= medianY.y);

    n->children[0]->x1 = n->x1;
    n->children[0]->x2 = n->x2;
    n->children[0]->y1 = n->y1;
    n->children[0]->y2 = medianY.y;
    n->children[0]->split = 2;
    n->children[0]->start = 0;
    n->children[0]->end = 0;

    //assert(medianY.y >= n->y2);
    n->children[1]->x1 = n->x1;
    n->children[1]->x2 = n->x2;
    n->children[1]->y1 = medianY.y;
    n->children[1]->y2 = n->y2;
    n->children[1]->split = 2;
    n->children[1]->start = 0;
    n->children[1]->end = 0;

    buildQuadTreeNode(n->children[0], elements, n->children[0]->elementCount, 2);
    buildQuadTreeNode(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount, 2); 
  } else {
    // we split horizontally.
    qsort(elements, n->elementCount,
        sizeof(quadTreeKey), quadTreeElementCompareX);
    quadTreeKey medianX = elements[n->elementCount / 2];

    n->children[0] = (quadTreeNode*) zmalloc(sizeof(quadTreeNode));
    n->children[1] = (quadTreeNode*) zmalloc(sizeof(quadTreeNode));

    n->children[0]->elementCount = n->elementCount / 2 + n->elementCount % 2;
    n->children[1]->elementCount = n->elementCount - n->children[0]->elementCount; 

    n->children[0]->children = NULL;
    n->children[1]->children = NULL;

    //assert(n->x1 <= medianX.x);
    n->children[0]->x1 = n->x1;
    n->children[0]->x2 = medianX.x;
    n->children[0]->y1 = n->y1;
    n->children[0]->y2 = n->y2;
    n->children[0]->split = 1;
    n->children[0]->start = 0;
    n->children[0]->end = 0;

    //assert(medianX.x <= n->x2);
    n->children[1]->x1 = medianX.x;
    n->children[1]->x2 = n->x2;
    n->children[1]->y1 = n->y1;
    n->children[1]->y2 = n->y2;
    n->children[1]->split = 1;
    n->children[1]->start = 0;
    n->children[1]->end = 0;

    buildQuadTreeNode(n->children[0], elements, n->children[0]->elementCount, 1);
    buildQuadTreeNode(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount, 1); 
  }
}

/*
void rebuildQuadTreeNode(quadTreeNode* n) {
  // allocate memory to get all the elements.
  quadTreeKey* elements = zmalloc(sizeof(quadTreeKey) * n->elementCount);
  getQuadTreeElements(n, elements);
  
  // now empty the quadtree
  quadTreeDeleteNode(n->children[0]); 
  quadTreeDeleteNode(n->children[1]);
  zfree(n->children);
  n->children = NULL;
  buildQuadTreeNode(n, elements);
  zfree(elements); 
}
*/
/*
bool quadTreeInsert(quadTreeNode* n, quadTreeKey key){
  // check if we are a leaf.
  if (n->children == NULL) {
    // check if key is a duplicate
    for (int i = n->start; i < n->end; i++) {
      if (key.x == n->elements[i].x && key.y == n->elements[i].y) {
        // this key already appears in this node.
        return false;
      }
    }
    if (n->end - n->start >= QUAD_TREE_NODE_SIZE) {
      splitQuadTreeNode(n);
    } else {
      n->elements[n->end] = key;
      n->end++;
      n->elementCount++;
      return true;
    }
  }
  bool added = false;
  if (n->split == 1) {
    // we split vertically.
    if (key.y > n->children[0]->y2) {
      // bottom half
      added = quadTreeInsert(n->children[0], key);
    } else {
      // top half
      added = quadTreeInsert(n->children[1], key);
    }
  }

  if (n->split == 2) {
    // we split horizontally
    if (key.x < n->children[0]->x2) {
      // left half
      added = quadTreeInsert(n->children[0], key);
    } else {
      // right half
      added = quadTreeInsert(n->children[1], key);
    }
  }
  
  if (added) {
    // add one to this node's elementCount
    n->elementCount++;
  }
 
  // check if we need to rebuild this subtree to restore balance.
  if (n->children[0]->elementCount > REBALANCE_FRACTION * n->children[1]->elementCount ||
      n->children[1]->elementCount > REBALANCE_FRACTION * n->children[0]->elementCount) {
    printf("Rebalancing because (%d, %d) \n", n->children[0]->elementCount, n->children[1]->elementCount);
    rebuildQuadTreeNode(n); 
  }

  return added;
}

// Deletes a key from the tree.
void quadTreeDelete(quadTreeNode* n, quadTreeKey key){
  // check if we are a leaf.
  if (n->children == NULL) {
    // check if the key is here.
    int scan = n->start;
    for (int i = n->start; i < n->end; i++) {
      n->elements[scan] = n->elements[i];
      if (key.x != n->elements[i].x || key.y != n->elements[i].y) {
        scan++;
      }
    }
    n->end = scan;
    return;
  }

  if (n->split == 1) {
    // we split vertically.
    if (key.y > n->children[0]->y2) {
      // bottom half
      quadTreeDelete(n->children[0], key);
    } else {
      // top half
      quadTreeDelete(n->children[1], key);
    }
  }

  if (n->split == 2) {
    // we split horizontally
    if (key.x < n->children[0]->x2) {
      // left half
      quadTreeDelete(n->children[0], key);
    } else {
      // right half
      quadTreeDelete(n->children[1], key);
    }
  }

  // check if everything can fit into this node.
  if (n->children[0]->end - n->children[0]->start +
      n->children[1]->end - n->children[1]->start <= QUAD_TREE_NODE_SIZE) {
    int j = n->start;
    n->elements = (quadTreeKey*) zmalloc(sizeof(quadTreeKey) * QUAD_TREE_NODE_SIZE);

    for (int i = n->children[0]->start; i < n->children[0]->end; i++) {
      n->elements[j] = n->children[0]->elements[i];
      j++;
    }

    for (int i = n->children[1]->start; i < n->children[1]->end; i++) {
      n->elements[j] = n->children[1]->elements[i];
      j++;
    }
    quadTreeDeleteNode(n->children[0]);
    quadTreeDeleteNode(n->children[1]);
    zfree(n->children);
    n->children = NULL;
  }
}
*/


// For debugging purposes.
void walkQuadTree(quadTreeNode* a){
  if (a->children == NULL) {
    /* for (int i = a->start; i < a->end; i++) { */
    /*   printf("element (%f, %f) \n", a->elements[i].x, a->elements[i].y);  */
    /* } */
    return;
  } else {
    walkQuadTree(a->children[0]);
    walkQuadTree(a->children[1]);
  }
}

// returns true if a key falls within a given range.
bool quadTreeKeyInRange(quadTreeKey key, double x1, double x2, double y1, double y2) {
  if (key.x >= x1 && key.x < x2 && key.y > y1 && key.y <= y2) {
    return true;
  } else {
    return false;
  }
}

// returns true if a node intersects with the given range, false otherwise.
bool quadTreeNodeIntersectsRange(quadTreeNode* n, double x1, double x2, double y1, double y2) {
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

void quadTreeRangeSearch(quadTreeNode *n, double x1, double x2, double y1, double y2, int* count) {
  if (n->children == NULL) {
    for (int i = n->start; i < n->end; i++) {
      if (quadTreeKeyInRange(n->elements[i], x1, x2, y1, y2)){
        (*count)++;
        /* printf("Element (%f, %f) in range \n", n->elements[i].x, n->elements[i].y); */
      }
    }
    return;
  }
  
  // call recursively on each child intersected by the range.
  if (quadTreeNodeIntersectsRange(n->children[0], x1, x2, y1, y2)) {
    quadTreeRangeSearch(n->children[0], x1, x2, y1, y2, count);
  }

  if (quadTreeNodeIntersectsRange(n->children[1], x1, x2, y1, y2)) {
    quadTreeRangeSearch(n->children[1], x1, x2, y1, y2, count);
  }
}


// the root node
quadTreeNode* rangeTreeRoot;

// collection of elements
quadTreeKey* allElements;
int allElementsCount;
int currentAllElementsArraySize;
double totalQueryTime;

void tfk2DRangeSearchCommand(redisClient* c) {
  double x1 = strtod(c->argv[1]->ptr, NULL);
  double x2 = strtod(c->argv[2]->ptr, NULL);
  double y1 = strtod(c->argv[3]->ptr, NULL);
  double y2 = strtod(c->argv[4]->ptr, NULL);
  int count = 0;
  double start_time = tfk_get_time();
  quadTreeRangeSearch(rangeTreeRoot, x1, x2, y1, y2, &count);
  double end_time = tfk_get_time();
/*
  for (int i = 0; i < allElementsCount; i++) {
    if (quadTreeKeyInRange(allElements[i], x1, x2, y1, y2)) {
      count++;
    }
  }
*/
  totalQueryTime += end_time - start_time;
  printf("total query time %f rangeTreeCount: %d \n", totalQueryTime, count); 
  addReplyLongLong(c, count);
}

void tfkBuildTreeCommand(redisClient *c) {
  if (rangeTreeRoot != NULL){
    // delete the old quad tree.
    quadTreeDeleteNode(rangeTreeRoot);
  }
  rangeTreeRoot = (quadTreeNode*) zmalloc(sizeof(quadTreeNode));

  rangeTreeRoot->x1 = -200;
  rangeTreeRoot->x2 = 200;
  rangeTreeRoot->y1 = -200;
  rangeTreeRoot->y2 = 200; 
  rangeTreeRoot->children = NULL;

  // rebuild using the allElements array.
  buildQuadTreeNode(rangeTreeRoot, allElements, allElementsCount, 1);
  walkQuadTree(rangeTreeRoot);

  addReplyLongLong(c, 42);
}

/* This generic command implements both ZADD and ZINCRBY. */
void tfkAddGenericCommand(redisClient *c, int incr) {
    double x = strtod(c->argv[1]->ptr, NULL);
    double y = strtod(c->argv[2]->ptr, NULL);
    
    int resizeAmount = 10000;
  
    // create a quadTreeKey
    quadTreeKey key;
    key.x = x;
    key.y = y;
    key.value = NULL; 
    if (allElementsCount == currentAllElementsArraySize) {
      // resize allElementsArray
      currentAllElementsArraySize += resizeAmount;
      allElements = (quadTreeKey*) zrealloc((void*)allElements, sizeof(quadTreeKey) * currentAllElementsArraySize);
    }
    allElements[allElementsCount] = key;
    allElementsCount++;
    addReplyLongLong(c, allElementsCount);
    /*
    if (a == NULL){ 
    a = (quadTreeNode*) zmalloc(sizeof(quadTreeNode));
    a->start = 0;
    a->end = 0;
    a->x1 = 0;
    a->x2 = 100;
    a->y1 = 0;
    a->y2 = 100; 
    a->children = NULL;
    a->split = 1;
    a->elements = (quadTreeKey*) zmalloc(sizeof(quadTreeKey)*QUAD_TREE_NODE_SIZE);
    a->elementCount = 0;
    } 
    quadTreeKey key1;
    key1.x = x;
    key1.y = y;
    quadTreeInsert(a, key1);
    walkQuadTree(a);*/ 
    /*quadTreeNode a;
    a.children = (quadTreeNode*) zmalloc(4*sizeof(quadTreeNode));
    a.start = 5;
    a.children[0].start = 4;*/ 
    //printf("Begin range search from x: 30 - 40 , y: 30 - 40 \n");
    //quadTreeRangeSearch(a, 30, 40, 30, 40); 
    //static char *nanerr = "resulting score is not a number (NaN)";
    //char* error = "ADD command not implemented nfnc;
    //addReplyError(c, error);
    //addReplyError(c, "nothing special");
    //addReplyLongLong(c, (long long)42);

/*    robj *key = c->argv[1];
    robj *ele;
    robj *zobj;
    robj *curobj;
*/
}

void zset2DAddCommand(redisClient *c) {
  tfkAddGenericCommand(c,0);
}

void zset2DBuildTreeCommand(redisClient *c) {
  tfkBuildTreeCommand(c);
}

void zset2DRangeSearchCommand(redisClient *c) {
  tfk2DRangeSearchCommand(c);
}

/*
void zcardCommand(redisClient *c) {
    robj *key = c->argv[1];
    robj *zobj;

    if ((zobj = lookupKeyReadOrReply(c,key,shared.czero)) == NULL ||
        checkType(c,zobj,REDIS_ZSET)) return;

    addReplyLongLong(c,zsetLength(zobj));
}
*/

