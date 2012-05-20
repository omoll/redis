#include "redis.h"

#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>
#include <float.h>
// Each leaf in the tree contains this many elements.
// This must be greater than 2, otherwise splitting range
// using the median may result in infinite recursion.
static int LAYERED_TREE_PLUS_OPT_NODE_SIZE = 8;

double tfkLayeredRangeTreePlusOpt_get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

int min_opt(int a, int b){
  return (a < b)? a : b;
}

typedef struct sortedArray2Key2 {
  double x;
  double y;
  void* value;
} sortedArray2Key2;

double* sorted_array2_recursive_successor_search(double value, double* start, double* end){
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
    double* ret = sorted_array2_recursive_successor_search(value, start, median);
    if (*ret == DBL_MAX) {
      return median;
    }
    return ret;
  } else {
    double* ret = sorted_array2_recursive_successor_search(value, median, end);
    return ret;
  }
}

//lexical comparison on x,y pairs, imitating the key concatenation thing.
int compare2SortedArray2KeyByX (const void * a, const void * b)
{
  if  (((sortedArray2Key2*)a)->x < ((sortedArray2Key2*)b)->x ) return -1;
  else if (((sortedArray2Key2*)a)->x > ((sortedArray2Key2*)b)->x ) return 1;
  else return 0;
}

int compare2SortedArray2KeyByY (const void * a, const void * b)
{
  if (((sortedArray2Key2*)a)->y < ((sortedArray2Key2*)b)->y ) return -1;
  else if (((sortedArray2Key2*)a)->y > ((sortedArray2Key2*)b)->y ) return 1;
  else return 0;
}

int resizeAmountOpt = 10000;
int allElementsCount_opt = 0;
int currentAllElementsArrayoptSize = 0;

sortedArray2Key2  *allElementsOptByX = NULL;
sortedArray2Key2  *allElementsOptByY = NULL;
double *allXopt = NULL;
double *allYopt = NULL;


void sortedArray22DAddCommand(redisClient *c) {
  double x = strtod(c->argv[1]->ptr, NULL);
  double y = strtod(c->argv[2]->ptr, NULL);

  sortedArray2Key2 key;
  key.x = x;
  key.y = y;
  key.value = NULL; 

  if (allElementsCount_opt == currentAllElementsArrayoptSize) {
    currentAllElementsArrayoptSize += resizeAmountOpt;
    allElementsOptByX = (sortedArray2Key2*) zrealloc((void*)allElementsOptByX, sizeof(sortedArray2Key2)*currentAllElementsArrayoptSize);
    allXopt = (double *) zrealloc((void*)allXopt, sizeof(double)*currentAllElementsArrayoptSize);
    allElementsOptByY = (sortedArray2Key2*) zrealloc((void*)allElementsOptByY, sizeof(sortedArray2Key2)*currentAllElementsArrayoptSize);
    allYopt = (double *) zrealloc((void*)allYopt, sizeof(double)*currentAllElementsArrayoptSize);
  }

  allElementsOptByX[allElementsCount_opt] = key;
  allElementsOptByY[allElementsCount_opt] = key;
  allXopt[allElementsCount_opt] = key.x;
  allYopt[allElementsCount_opt] = key.y;

  allElementsCount_opt++;
  //addReplyLongLong(c, allElementsCount_opt);
}


int compare2(const void * a, const void * b){
  if (*(double *) a > *(double *) b) return 1;
  else if (*(double *) a < *(double *) b) return -1;
  else return 0;
}

void sortedArray22DBuildCommand(redisClient *c) {
  qsort(allElementsOptByX, allElementsCount_opt, sizeof(sortedArray2Key2), compare2SortedArray2KeyByX);
  qsort(allXopt, allElementsCount_opt, sizeof(double), compare2);
  qsort(allElementsOptByY, allElementsCount_opt, sizeof(sortedArray2Key2), compare2SortedArray2KeyByY);
  qsort(allYopt, allElementsCount_opt, sizeof(double), compare2);
  //addReplyLongLong(c, 42);
}
typedef struct cascading_opt_node {
  int size; // array sizes
  double *ys;// eg ys
  double *xs; //eg x coords
  int *left_successor_opt; //positions for cascading_opt left (you keep track
  int *right_successor_opt; // positions for cascading_opt right

  int *left2_successor_opt;
  int *right2_successor_opt;

  struct cascading_opt_node *left;
  struct cascading_opt_node *left2;
  struct cascading_opt_node *right;
  struct cascading_opt_node *right2;

} cascading_opt_node_t;

void cascading_opt_node_print(cascading_opt_node_t *cn);

/*the value at *end is meant to be a sentinel MAX_INT
 *gives you the succesor position (smallest elt greater than or equal to value)
 */

double* linear_successor_opt_search(double value, double* start, double* end){
 while( !(*start >= value)){
   start++;
 }
 return start;
}

double* recursive_successor_opt_search(double value, double* start, double* end){
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
    double* ret = recursive_successor_opt_search(value, start, median);
    if (*ret == DBL_MAX) {
      return median;
    }
    return ret;
  } else {
    double* ret = recursive_successor_opt_search(value, median, end);
    return ret;
  }
}

/*CURRENTLY BUGGY, CHECK TEST CASES.
 */
double* binary_successor_opt_search(double value, double* start, double* end){

  double *initial_start = start;

  while(start < end){
    double * middle =  start + (end-start)/2;
    if ( value < *middle ) {
      end = middle;
    } else if ( value > *middle ){
      start = middle + 1;
    } else {

      // right now, do this stupidly
      while( *middle == value  && middle > start){
	middle--;
      }
      
      if (middle != initial_start){
	return middle + 1;
      } else {
	return middle;
      }

      return middle;
    }
  }

  return end;
}

/*
  binary searches for the index of the successor_opt (greater than or
  equal) to the value in the node array, 
  only meant to be used on the root.
*/
int arg_successor_opt_root(double yvalue, cascading_opt_node_t *root){
  //TODO:when debugged binary search and place here.
  return recursive_successor_opt_search(yvalue, root->ys, root->ys + root->size) - root->ys;
}

/*yindex is an index for on the parent
 */
int arg_successor_opt_left(int yindex, cascading_opt_node_t *parent){
  if (yindex < parent->size){
    return parent->left_successor_opt[yindex];
  } else {
    return parent->left->size;
  }
  /* int y = arg_successor_opt_root(parent->ys[yindex], parent->left); */
  /* assert(x == y); */
  /* return y; */
}

/*yindex is an index for on the parent
 */
int arg_successor_opt_left2(int yindex, cascading_opt_node_t *parent){
  if (yindex < parent->size){
    return parent->left2_successor_opt[yindex];
  } else {
    return parent->left2->size;
  }
  /* int y = arg_successor_opt_root(parent->ys[yindex], parent->left); */
  /* assert(x == y); */
  /* return y; */
}

int arg_successor_opt_right(int yindex, cascading_opt_node_t *parent){
  if (yindex < parent->size){
    return parent->right_successor_opt[yindex];
  } else {
    return parent->right->size;
  }
  /* int y = arg_successor_opt_root(parent->ys[yindex], parent->right); */
  /* assert(x== y); */
  /* return y; */
}

int arg_successor_opt_right2(int yindex, cascading_opt_node_t *parent){
  if (yindex < parent->size){
    return parent->right2_successor_opt[yindex];
  } else {
    return parent->right2->size;
  }
  /* int y = arg_successor_opt_root(parent->ys[yindex], parent->right); */
  /* assert(x== y); */
  /* return y; */
}

int cascading_opt_node_is_leaf(cascading_opt_node_t *cn){
  return cn->left == NULL;
}

void check_rep_cascading_opt_node(cascading_opt_node_t *c){
  assert(c->size > 0); 

  //leaf condition all or none should be NULL
  if (cascading_opt_node_is_leaf(c)){
    assert(c->left_successor_opt == NULL);
    assert(c->left == NULL); 
    assert(c->right_successor_opt == NULL);
    assert(c->right == NULL); 
  }
  // sentinel
  assert(c->ys[c->size] == DBL_MAX);
  int i;
  //no repeated ys, 
  for (i = 0; i < c->size; i++){
    assert(c->ys[i] <= c->ys[i+1]);
    assert(c->ys[i] < DBL_MAX);
  }

  if (!cascading_opt_node_is_leaf(c)){
    //sentinels inited to -1
    assert(*( c->left_successor_opt - 1 ) == -1);
    assert(*( c->right_successor_opt - 1 ) == -1);

      for (i = 0; i < c->size; i++){
	//check indices into left array all make sense
	assert( c->left_successor_opt[i] > -1 && c->left_successor_opt[i] >= c->left_successor_opt[i-1]
		&& c->left_successor_opt[i] <= c->left->size);

	//check indices into right child all make sense
	assert(c->right_successor_opt[i] > -1);
	assert(c->right_successor_opt[i] >= c->right_successor_opt[i-1]);
	assert(c->right_successor_opt[i] <= c->right->size);
      }


  }

}

void cascading_opt_node_leaf_init(cascading_opt_node_t *cn, double *x, double *y, int size){
  cn->size = size;
  cn->xs = x;
  cn->ys = (double*) zmalloc((size+1)*sizeof(double));
  memcpy(cn->ys, y, size*sizeof(double));
  cn->ys[size] = DBL_MAX;

  cn->left_successor_opt = NULL;
  cn->left2_successor_opt = NULL;
  cn->right_successor_opt = NULL;
  cn->right2_successor_opt = NULL;

  cn->left  = NULL;
  cn->left2 = NULL;
  cn->right  = NULL;
  cn->right2 = NULL;

  check_rep_cascading_opt_node(cn);
  assert(cascading_opt_node_is_leaf(cn));
}

//inits cascading_opt node cn,  given its two children and the values at a
//(normal)node 
void cascading_opt_node_parent_init(cascading_opt_node_t *cn, cascading_opt_node_t *left, cascading_opt_node_t* left2,
    cascading_opt_node_t *right, cascading_opt_node_t* right2){
  cn->size = left->size + left2->size + right->size + right2->size; //has all values of
						 //left all of right
						 //and the new values

  //allocates +1 size for sentinel value.
  cn->ys = (double *)zmalloc((cn->size+1)*sizeof(double));
  cn->ys[cn->size] = DBL_MAX;

  //does not need any sentinel value, simply the other coordinate
  cn->xs = (double *)zmalloc(cn->size*sizeof(double));

  //allocate an array of size +1 length, the -1 entry gets inited to 0
  cn->left_successor_opt = (int *) zmalloc((cn->size+1)*sizeof(int)) + 1;
  *(cn->left_successor_opt - 1) = -1;
  cn->left2_successor_opt = (int *) zmalloc((cn->size+1)*sizeof(int)) + 1;
  *(cn->left2_successor_opt - 1) = -1;

  cn->right_successor_opt = (int *) zmalloc((cn->size+1)*sizeof(int)) + 1;
  *(cn->right_successor_opt - 1) = -1;
  cn->right2_successor_opt = (int *) zmalloc((cn->size+1)*sizeof(int)) + 1;
  *(cn->right2_successor_opt - 1) = -1;

  cn->left = left;
  cn->left2 = left2;
  cn->right = right;
  cn->right2 = right2;

  int i, i2, j, j2;
  for (i = 0, i2 = 0, j = 0, j2 = 0;  i < left->size || i2 < left2->size || j < right->size || j2 < right2->size;){
    double* min_optchild = left->ys; 
    double min_opt_value = left->ys[i];
    
    if (min_opt_value > left2->ys[i2]) {
      min_opt_value = left2->ys[i2];
      min_optchild = left2->ys;
    }
    if (min_opt_value > right->ys[j]){
      min_opt_value = right->ys[j];
      min_optchild = right->ys;
    }
    if (min_opt_value > right2->ys[j2]){
      min_opt_value = right2->ys[j2];
      min_optchild = right2->ys;
    }

    if (min_optchild == left->ys){
      // left has the next smallest element
      cn->ys[i+i2+j+j2] = left->ys[i];
      cn->xs[i+i2+j+j2] = left->xs[i];
      i++;
    }
    if (min_optchild == left2->ys) {
      cn->ys[i+i2+j+j2] = left2->ys[i2];
      cn->xs[i+i2+j+j2] = left2->xs[i2];
      i2++;
    }
    if (min_optchild == right->ys){
      cn->ys[i+i2+j+j2] = right->ys[j];
      cn->xs[i+i2+j+j2] = right->xs[j];
      j++;
    }
    if (min_optchild == right2->ys){
      cn->ys[i+i2+j+j2] = right2->ys[j2];
      cn->xs[i+i2+j+j2] = right2->xs[j2];
      j2++;
    }
  }
  for (i = 0; i < cn->size; i++){
    cn->left_successor_opt[i] = arg_successor_opt_root(cn->ys[i], cn->left);
    cn->left2_successor_opt[i] = arg_successor_opt_root(cn->ys[i], cn->left2);
    cn->right_successor_opt[i] = arg_successor_opt_root(cn->ys[i], cn->right);
    cn->right2_successor_opt[i] = arg_successor_opt_root(cn->ys[i], cn->right2);
  }

  //  cascading_opt_node_print(cn);
  check_rep_cascading_opt_node(cn);
  assert(!cascading_opt_node_is_leaf(cn));
}

/*
  Type Defs
*/
typedef struct layeredRangeTreePlusOptKey{
  double x;
  double y;
  void* value;
} layeredRangeTreePlusOptKey;

// the inner level of the tree.
typedef struct layeredRangeTreePlusOptNodeLevel2 {
  // The rectangle corresponding to this node.
  double x1;
  double x2;
  double y1; 
  double y2;

  // a count of elments in the subtree for use in rebalancing.
  int elementCount;

  // TODO(tfkLayeredRangeTreePlusOpt): start is almost always 0.
  // start and end of range in elements array.
  int start;
  int end;
  layeredRangeTreePlusOptKey* elements; // undefined if not leaf.

  cascading_opt_node_t cnode;
  
  struct layeredRangeTreePlusOptNodeLevel2** children; // NULL if leaf
} layeredRangeTreePlusOptNodeLevel2;

// the top level of the tree.
typedef struct layeredRangeTreePlusOptNodeLevel1 {
  // The rectangle corresponding to this node.
  double x1;
  double x2;
  double y1; 
  double y2;

  // a count of elments in the subtree for use in rebalancing.
  int elementCount;

  // TODO(tfkLayeredRangeTreePlusOpt): start is almost always 0.
  // start and end of range in elements array.
  int start;
  int end;
  layeredRangeTreePlusOptKey* elements; // undefined if not leaf.

  struct layeredRangeTreePlusOptNodeLevel1** children; // NULL if leaf

  // a pointer to the second level of the tree containing the same points as in this subtree
  // sorted on the other coordinate.
  layeredRangeTreePlusOptNodeLevel2* secondLevelPointer;
} layeredRangeTreePlusOptNodeLevel1;





void layeredRangeTreePlusOptDeleteNodeLevel2(layeredRangeTreePlusOptNodeLevel2* n) {
  if (n->children == NULL) {
    zfree(n->elements);
    zfree(n);
  } else {
    layeredRangeTreePlusOptDeleteNodeLevel2(n->children[0]);
    layeredRangeTreePlusOptDeleteNodeLevel2(n->children[1]);
    layeredRangeTreePlusOptDeleteNodeLevel2(n->children[2]);
    layeredRangeTreePlusOptDeleteNodeLevel2(n->children[3]);
    zfree(n->children);
    zfree(n);
  }
}

// Delete node n and it's associated subree.
void layeredRangeTreePlusOptDeleteNodeLevel1(layeredRangeTreePlusOptNodeLevel1* n){
  if (n->children == NULL) {
    zfree(n->elements);
    zfree(n);
  } else {
    layeredRangeTreePlusOptDeleteNodeLevel1(n->children[0]);
    layeredRangeTreePlusOptDeleteNodeLevel1(n->children[1]);
    layeredRangeTreePlusOptDeleteNodeLevel1(n->children[2]);
    layeredRangeTreePlusOptDeleteNodeLevel1(n->children[3]);

    zfree(n->children);
    zfree(n);
  }

  if (n->secondLevelPointer != NULL){
    layeredRangeTreePlusOptDeleteNodeLevel2(n->secondLevelPointer);
  }
}

/*
  Comparison functions for median finding and sorting.
*/
int layeredRangeTreePlusOptElementCompareX(const void* a, const void* b){
  if (((layeredRangeTreePlusOptKey*) a)->x < ((layeredRangeTreePlusOptKey*) b)->x) {
    return -1;
  }
  if (((layeredRangeTreePlusOptKey*) a)->x > ((layeredRangeTreePlusOptKey*) b)->x) {
    return 1;
  }
  return 0;
}

int layeredRangeTreePlusOptElementCompareY(const void* a, const void* b){
  if (((layeredRangeTreePlusOptKey*) a)->y < ((layeredRangeTreePlusOptKey*) b)->y) {
    return -1;
  }
  if (((layeredRangeTreePlusOptKey*) a)->y > ((layeredRangeTreePlusOptKey*) b)->y) {
    return 1;
  }
  return 0;
}

void buildLayeredRangeTreePlusOptNodeLevel2 (layeredRangeTreePlusOptNodeLevel2* n, layeredRangeTreePlusOptKey* elements, int elementCount) {
  n->elementCount = elementCount;

  // Coarsen 
  if (n->elementCount <= LAYERED_TREE_PLUS_OPT_NODE_SIZE) {
    n->start = 0; 
    n->end = n->elementCount;
    n->elements = (layeredRangeTreePlusOptKey*) zmalloc(sizeof(layeredRangeTreePlusOptKey) * LAYERED_TREE_PLUS_OPT_NODE_SIZE);
    for (int i = 0; i < n->elementCount; i++) {
      n->elements[i] = elements[i];
    }
    return;
  }
 
  n->children = (layeredRangeTreePlusOptNodeLevel2**) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel2*)*4);

  qsort(elements, n->elementCount,
      sizeof(layeredRangeTreePlusOptKey), layeredRangeTreePlusOptElementCompareY);

  layeredRangeTreePlusOptKey medianY1 = elements[n->elementCount / 4];
  layeredRangeTreePlusOptKey medianY2 = elements[n->elementCount / 2];
  layeredRangeTreePlusOptKey medianY3 = elements[(3*n->elementCount) / 4];

  n->children[0] = (layeredRangeTreePlusOptNodeLevel2*) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel2));
  n->children[1] = (layeredRangeTreePlusOptNodeLevel2*) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel2));
  n->children[2] = (layeredRangeTreePlusOptNodeLevel2*) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel2));
  n->children[3] = (layeredRangeTreePlusOptNodeLevel2*) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel2));

  n->children[0]->elementCount = n->elementCount / 4 + (n->elementCount%4!=0);
  n->children[1]->elementCount = n->elementCount / 2 + (n->elementCount%2!=0) - n->children[0]->elementCount;
  n->children[2]->elementCount = (3*n->elementCount) / 4 + (3*n->elementCount)%4!=0
      - n->children[0]->elementCount - n->children[1]->elementCount;
  n->children[3]->elementCount = n->elementCount
      - n->children[0]->elementCount - n->children[1]->elementCount - n->children[2]->elementCount;


  n->children[0]->children = NULL;
  n->children[1]->children = NULL;
  n->children[2]->children = NULL;
  n->children[3]->children = NULL;

  n->children[0]->x1 = n->x1;
  n->children[0]->x2 = n->x2;
  n->children[0]->y1 = n->y1;
  n->children[0]->y2 = medianY1.y;
  n->children[0]->start = 0;
  n->children[0]->end = 0;

  n->children[1]->x1 = n->x1;
  n->children[1]->x2 = n->x2;
  n->children[1]->y1 = medianY1.y;
  n->children[1]->y2 = medianY2.y;
  n->children[1]->start = 0;
  n->children[1]->end = 0;

  n->children[2]->x1 = n->x1;
  n->children[2]->x2 = n->x2;
  n->children[2]->y1 = medianY2.y;
  n->children[2]->y2 = medianY3.y;
  n->children[2]->start = 0;
  n->children[2]->end = 0;

  n->children[3]->x1 = n->x1;
  n->children[3]->x2 = n->x2;
  n->children[3]->y1 = medianY3.y;
  n->children[3]->y2 = n->y2;
  n->children[3]->start = 0;
  n->children[3]->end = 0;

  buildLayeredRangeTreePlusOptNodeLevel2(n->children[0], elements, n->children[0]->elementCount);
  buildLayeredRangeTreePlusOptNodeLevel2(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount); 
  buildLayeredRangeTreePlusOptNodeLevel2(n->children[2], elements + n->children[0]->elementCount + n->children[1]->elementCount, n->children[1]->elementCount); 
  buildLayeredRangeTreePlusOptNodeLevel2(n->children[3],
      elements + n->children[0]->elementCount + n->children[1]->elementCount + n->children[2]->elementCount, n->children[2]->elementCount); 
}

// builds a second level of the quad tree.
void explodeTreePlusOpt(layeredRangeTreePlusOptNodeLevel1* n, layeredRangeTreePlusOptKey* elements, int elementCount){
  n->secondLevelPointer = NULL;
  if (n->elementCount <= LAYERED_TREE_PLUS_OPT_NODE_SIZE) {
    assert(elementCount == n->elementCount);
    n->secondLevelPointer = (layeredRangeTreePlusOptNodeLevel2*) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel2));

    // this is a leaf node, copy its elements into the array.
    memcpy(elements, n->elements, elementCount * sizeof(layeredRangeTreePlusOptKey));

    qsort(elements, n->elementCount,
        sizeof(layeredRangeTreePlusOptKey), layeredRangeTreePlusOptElementCompareY);

    double* xs = (double*) zmalloc(sizeof(double)*elementCount);  
    double* ys = (double*) zmalloc(sizeof(double)*elementCount);

    assert(xs != NULL && ys != NULL);

    for (int i = 0; i < elementCount; i++) {
      xs[i] = elements[i].x;
      ys[i] = elements[i].y;
    }

    cascading_opt_node_leaf_init(&n->secondLevelPointer->cnode, xs, ys, elementCount); 

    return;
  }
 
  explodeTreePlusOpt(n->children[0], elements, n->children[0]->elementCount); 
  explodeTreePlusOpt(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount);
  explodeTreePlusOpt(n->children[2],
      elements + n->children[0]->elementCount + n->children[1]->elementCount, n->children[2]->elementCount);
  explodeTreePlusOpt(n->children[3],
      elements + n->children[0]->elementCount + n->children[1]->elementCount + n->children[2]->elementCount, n->children[3]->elementCount);


  n->secondLevelPointer = (layeredRangeTreePlusOptNodeLevel2*) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel2));

  cascading_opt_node_parent_init(&n->secondLevelPointer->cnode, &n->children[0]->secondLevelPointer->cnode, &n->children[1]->secondLevelPointer->cnode,
      &n->children[2]->secondLevelPointer->cnode, &n->children[3]->secondLevelPointer->cnode);

  n->secondLevelPointer->x1 = n->x1;
  n->secondLevelPointer->x2 = n->x2;
  n->secondLevelPointer->y1 = n->y1;
  n->secondLevelPointer->y2 = n->y2;
  n->secondLevelPointer->children = NULL;
  n->secondLevelPointer->elements = NULL;
  n->secondLevelPointer->elementCount = n->elementCount;
  n->secondLevelPointer->start = 0;
  n->secondLevelPointer->end = 0;
  
  //buildLayeredRangeTreePlusOptNodeLevel2(n->secondLevelPointer, elements, n->elementCount);
  assert(n->secondLevelPointer != NULL);
}

void buildLayeredRangeTreePlusOptNodeLevel1(layeredRangeTreePlusOptNodeLevel1* n, layeredRangeTreePlusOptKey* elements, int elementCount) {
  n->elementCount = elementCount;
  n->secondLevelPointer = NULL;

  // Coarsen 
  if (n->elementCount <= LAYERED_TREE_PLUS_OPT_NODE_SIZE) {
    n->start = 0; 
    n->end = n->elementCount;
    n->elements = (layeredRangeTreePlusOptKey*) zmalloc(sizeof(layeredRangeTreePlusOptKey) * LAYERED_TREE_PLUS_OPT_NODE_SIZE);
    for (int i = 0; i < n->elementCount; i++) {
      n->elements[i] = elements[i];
    }
    return;
  }
 
  n->children = (layeredRangeTreePlusOptNodeLevel1**) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel1*)*4);

  qsort(elements, n->elementCount,
      sizeof(layeredRangeTreePlusOptKey), layeredRangeTreePlusOptElementCompareX);

  layeredRangeTreePlusOptKey medianX1 = elements[n->elementCount / 4];
  layeredRangeTreePlusOptKey medianX2 = elements[n->elementCount / 2];
  layeredRangeTreePlusOptKey medianX3 = elements[(3*n->elementCount) / 4];

  n->children[0] = (layeredRangeTreePlusOptNodeLevel1*) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel2));
  n->children[1] = (layeredRangeTreePlusOptNodeLevel1*) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel2));
  n->children[2] = (layeredRangeTreePlusOptNodeLevel1*) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel2));
  n->children[3] = (layeredRangeTreePlusOptNodeLevel1*) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel2));
/*
  n->children[0]->elementCount = n->elementCount / 4 + (n->elementCount%4!=0);
  n->children[1]->elementCount = n->elementCount / 2 + (n->elementCount%2!=0) - n->children[0]->elementCount;
  n->children[2]->elementCount = (3*n->elementCount) / 4 + (3*n->elementCount)%4!=0
      - n->children[0]->elementCount - n->children[1]->elementCount;
  n->children[3]->elementCount = n->elementCount
      - n->children[0]->elementCount - n->children[1]->elementCount - n->children[2]->elementCount;
*/
  int j = 1;
  for (int i = 0; i < n->elementCount; i++) {  
    if (i == n->elementCount/4){
      n->children[0]->elementCount = i;
    }
    if (i == n->elementCount/2) {
      n->children[1]->elementCount = i - n->children[0]->elementCount;
    }
    if (i == (3*n->elementCount)/4) {
      n->children[2]->elementCount = i - n->children[0]->elementCount - n->children[1]->elementCount;
    }
    j++;
  }
  n->children[3]->elementCount = n->elementCount - n->children[2]->elementCount - n->children[1]->elementCount - n->children[0]->elementCount;

  n->children[0]->children = NULL;
  n->children[1]->children = NULL;
  n->children[2]->children = NULL;
  n->children[3]->children = NULL;

  n->children[0]->x1 = n->x1;
  n->children[0]->x2 = medianX1.x;
  n->children[0]->y1 = n->y1;
  n->children[0]->y2 = n->y2;
  n->children[0]->start = 0;
  n->children[0]->end = 0;

  n->children[1]->x1 = medianX1.x;
  n->children[1]->x2 = medianX2.x;
  n->children[1]->y1 = n->y1;
  n->children[1]->y2 = n->y2;
  n->children[1]->start = 0;
  n->children[1]->end = 0;

  n->children[2]->x1 = medianX2.x;
  n->children[2]->x2 = medianX3.x;
  n->children[2]->y1 = n->y1;
  n->children[2]->y2 = n->y2;
  n->children[2]->start = 0;
  n->children[2]->end = 0;

  n->children[3]->x1 = medianX3.x;
  n->children[3]->x2 = n->x2;
  n->children[3]->y1 = n->y1;
  n->children[3]->y2 = n->y2;
  n->children[3]->start = 0;
  n->children[3]->end = 0;

  buildLayeredRangeTreePlusOptNodeLevel1(n->children[0], elements, n->children[0]->elementCount);
  buildLayeredRangeTreePlusOptNodeLevel1(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount);
  buildLayeredRangeTreePlusOptNodeLevel1(n->children[2],
      elements + n->children[0]->elementCount + n->children[1]->elementCount, n->children[2]->elementCount); 
  buildLayeredRangeTreePlusOptNodeLevel1(n->children[3],
      elements + n->children[0]->elementCount + n->children[1]->elementCount + n->children[2]->elementCount, n->children[3]->elementCount); 
}

// returns true if a key falls within a given range.
bool layeredRangeTreePlusOptKeyInRange(layeredRangeTreePlusOptKey key, double x1, double x2, double y1, double y2) {
  if (key.x >= x1 && key.x < x2 && key.y > y1 && key.y <= y2) {
    return true;
  } else {
    return false;
  }
}

// returns true if a node intersects with the given range, false otherwise.
bool layeredRangeTreePlusOptNodeIntersectsRangeLevel1(layeredRangeTreePlusOptNodeLevel1* n, double x1, double x2, double y1, double y2) {
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
bool layeredRangeTreePlusOptNodeIntersectsRangeLevel2(layeredRangeTreePlusOptNodeLevel2* n, double x1, double x2, double y1, double y2) {
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

void layeredRangeTreePlusOptRangeSearchLevel2(layeredRangeTreePlusOptNodeLevel2 *n, double x1, double x2, double y1, double y2, int* count) {
  if (n->children == NULL) {
    assert(n->elementCount <= LAYERED_TREE_PLUS_OPT_NODE_SIZE);
    for (int i = n->start; i < n->end; i++) {
      if (layeredRangeTreePlusOptKeyInRange(n->elements[i], -200, 200, y1, y2)){
        (*count)++;
        //printf("Element (%f, %f) in range \n", n->elements[i].x, n->elements[i].y);
      }
    }
    return;
  }
 
  // call recursively on each child intersected by the range.
  if (layeredRangeTreePlusOptNodeIntersectsRangeLevel2(n->children[0], -200, 200, y1, y2)) {
    layeredRangeTreePlusOptRangeSearchLevel2(n->children[0], x1, x2, y1, y2, count);
  }

  if (layeredRangeTreePlusOptNodeIntersectsRangeLevel2(n->children[1], -200, 200, y1, y2)) {
    layeredRangeTreePlusOptRangeSearchLevel2(n->children[1], x1, x2, y1, y2, count);
  }

  if (layeredRangeTreePlusOptNodeIntersectsRangeLevel2(n->children[2], -200, 200, y1, y2)) {
    layeredRangeTreePlusOptRangeSearchLevel2(n->children[2], x1, x2, y1, y2, count);
  }
  if (layeredRangeTreePlusOptNodeIntersectsRangeLevel2(n->children[3], -200, 200, y1, y2)) {
    layeredRangeTreePlusOptRangeSearchLevel2(n->children[3], x1, x2, y1, y2, count);
  }

}

void layeredRangeTreePlusOptRangeSearchLevel1(layeredRangeTreePlusOptNodeLevel1 *n, double x1, double x2, double y1, double y2, int* count, int offset) {
  if (offset == -1){
    // do the initial binary search.
    offset = arg_successor_opt_root(y1, &n->secondLevelPointer->cnode);
  }

  if (n->children == NULL) {
    assert(n->elementCount <= LAYERED_TREE_PLUS_OPT_NODE_SIZE);
    // if we get to a leaf just scan through the elements and report the ones in the range.
    for (int i = n->start; i < n->end; i++) {
      if (layeredRangeTreePlusOptKeyInRange(n->elements[i], x1, x2, y1, y2)){
        (*count)++;
        //printf("Element (%f, %f) in range \n", n->elements[i].x, n->elements[i].y);
      }
    }
    return;
  }

  assert(n->secondLevelPointer != NULL);

  // call the second level search if we're entirely within the x-range.
  if (n->x1 >= x1 && n->x2 <= x2) {
    for (;offset < n->secondLevelPointer->cnode.size && n->secondLevelPointer->cnode.ys[offset] <= y2; offset++ ){
            //printf("Element (%f, %f) in range \n", n->secondLevelPointer->cnode.xs[offset], n->secondLevelPointer->cnode.ys[offset]);
      (*count)++;
    }
    //layeredRangeTreePlusOptRangeSearchLevel2(n->secondLevelPointer, -200, 200, y1, y2, count);
    return;
  }

  // otherwise call recursively on each child intersected by the range.
  if (layeredRangeTreePlusOptNodeIntersectsRangeLevel1(n->children[0], x1, x2, -200, 200)) {
    layeredRangeTreePlusOptRangeSearchLevel1(n->children[0], x1, x2, y1, y2, count, 
        arg_successor_opt_left(offset, &n->secondLevelPointer->cnode));
  }

  if (layeredRangeTreePlusOptNodeIntersectsRangeLevel1(n->children[1], x1, x2, -200, 200)) {
    layeredRangeTreePlusOptRangeSearchLevel1(n->children[1], x1, x2, y1, y2, count,
        arg_successor_opt_left2(offset, &n->secondLevelPointer->cnode));
  }
  if (layeredRangeTreePlusOptNodeIntersectsRangeLevel1(n->children[2], x1, x2, -200, 200)) {
    layeredRangeTreePlusOptRangeSearchLevel1(n->children[2], x1, x2, y1, y2, count,
        arg_successor_opt_right(offset, &n->secondLevelPointer->cnode));
  }
  if (layeredRangeTreePlusOptNodeIntersectsRangeLevel1(n->children[3], x1, x2, -200, 200)) {
    layeredRangeTreePlusOptRangeSearchLevel1(n->children[3], x1, x2, y1, y2, count,
        arg_successor_opt_right2(offset, &n->secondLevelPointer->cnode));
  }

}


// the root node
layeredRangeTreePlusOptNodeLevel1* rangeTreePlusOptRoot;

// collection of elements
layeredRangeTreePlusOptKey* rtplus_opt_allElements;
int rtplus_opt_allElementsCount_opt;
int rtplus_opt_currentAllElementsArrayoptSize;
double rtplus_opt_totalQueryTime;

void tfkLayeredRangeTreePlusOpt2DRangeSearchCommand(redisClient* c) {
  double x1 = strtod(c->argv[1]->ptr, NULL);
  double x2 = strtod(c->argv[2]->ptr, NULL);
  double y1 = strtod(c->argv[3]->ptr, NULL);
  double y2 = strtod(c->argv[4]->ptr, NULL);
  int count = 0;

  // determine which data structure to use.
  double ratio = min_opt((x2-x1)/(y2-y1), (y2-y1)/(x2-x1));
  double density = (((double)rtplus_opt_allElementsCount_opt)/10000);
  double k = density*(y2-y1)*(x2-x1); // estimated number of points reported.
  double cutoff = rtplus_opt_allElementsCount_opt / (2*k);

  if (ratio < cutoff) {
  //if (true){
    // use the cascading tree.
    double start_time = tfkLayeredRangeTreePlusOpt_get_time();
    layeredRangeTreePlusOptRangeSearchLevel1(rangeTreePlusOptRoot, x1, x2, y1, y2, &count, -1);
    double end_time = tfkLayeredRangeTreePlusOpt_get_time();
    rtplus_opt_totalQueryTime += end_time - start_time;
  } else {
    // otherwise use filtering approach.
    double start_time = tfkLayeredRangeTreePlusOpt_get_time();
    double *scan_start;
    //if X is narrower, then search by x first
    if (x2 - x1 < y2 - y1){
      scan_start = sorted_array2_recursive_successor_search(x1, allXopt, allXopt + allElementsCount_opt);
      int i = scan_start - allXopt;
      for ( ; i < allElementsCount_opt && allElementsOptByX[i].x <= x2; i++){
        if (allElementsOptByX[i].y >= y1 && allElementsOptByX[i].y <= y2) count++;
      }
    } else {
      scan_start = sorted_array2_recursive_successor_search(y1, allYopt, allYopt + allElementsCount_opt);
      int i = scan_start - allYopt;
      for ( ; i < allElementsCount_opt && allElementsOptByY[i].y <= y2; i++){
        if (allElementsOptByY[i].x >= x1 && allElementsOptByY[i].x <= x2) count++;
      }
    }
    double end_time = tfkLayeredRangeTreePlusOpt_get_time();
    printf("All elemnts for sorted array %d \n", allElementsCount_opt);
    rtplus_opt_totalQueryTime += end_time - start_time;
  }
  printf("total query time %f rangeTreePlusOptCount: %d \n", rtplus_opt_totalQueryTime, count); 
  addReplyLongLong(c, count);
}

void tfkLayeredRangeTreePlusOptBuildTreePlusOptCommand(redisClient *c) {
  if (rangeTreePlusOptRoot != NULL){
    // delete the old quad tree.
    layeredRangeTreePlusOptDeleteNodeLevel1(rangeTreePlusOptRoot);
  }
  rangeTreePlusOptRoot = (layeredRangeTreePlusOptNodeLevel1*) zmalloc(sizeof(layeredRangeTreePlusOptNodeLevel1));

  rangeTreePlusOptRoot->x1 = -200;
  rangeTreePlusOptRoot->x2 = 200;
  rangeTreePlusOptRoot->y1 = -200;
  rangeTreePlusOptRoot->y2 = 200; 
  rangeTreePlusOptRoot->children = NULL;
  rangeTreePlusOptRoot->elementCount = rtplus_opt_allElementsCount_opt;
  rangeTreePlusOptRoot->secondLevelPointer = NULL;
  // rebuild using the rtplus_opt_allElements array.

  // build the first level of the tree.
  buildLayeredRangeTreePlusOptNodeLevel1(rangeTreePlusOptRoot, rtplus_opt_allElements, rtplus_opt_allElementsCount_opt);

  // build the second level of the tree.
  layeredRangeTreePlusOptKey* elementsCopy = (layeredRangeTreePlusOptKey*) zmalloc(sizeof(layeredRangeTreePlusOptKey) * rangeTreePlusOptRoot->elementCount);
  explodeTreePlusOpt(rangeTreePlusOptRoot, elementsCopy, rtplus_opt_allElementsCount_opt);
  zfree(elementsCopy); 
  
  //walkLayeredRangeTreePlusOpt(rangeTreePlusOptRoot);
  rtplus_opt_totalQueryTime = 0;
  addReplyLongLong(c, 42);
}

/* This generic command implements both ZADD and ZINCRBY. */
void tfkLayeredRangeTreePlusOptAddGenericCommand(redisClient *c, int incr) {
    double x = strtod(c->argv[1]->ptr, NULL);
    double y = strtod(c->argv[2]->ptr, NULL);
    
    int resizeAmountOpt = 10000;
  
    // create a layeredRangeTreePlusOptKey
    layeredRangeTreePlusOptKey key;
    key.x = x;
    key.y = y;
    key.value = NULL; 
    if (rtplus_opt_allElementsCount_opt == rtplus_opt_currentAllElementsArrayoptSize) {
      // resize rtplus_opt_allElementsArray
      rtplus_opt_currentAllElementsArrayoptSize += resizeAmountOpt;
      rtplus_opt_allElements = (layeredRangeTreePlusOptKey*) zrealloc((void*)rtplus_opt_allElements, sizeof(layeredRangeTreePlusOptKey) * rtplus_opt_currentAllElementsArrayoptSize);
    }
    rtplus_opt_allElements[rtplus_opt_allElementsCount_opt] = key;
    rtplus_opt_allElementsCount_opt++;
    addReplyLongLong(c, rtplus_opt_allElementsCount_opt);
}

void layeredRangeTreePlusOpt2DAddCommand(redisClient *c) {
  sortedArray22DAddCommand(c);
  tfkLayeredRangeTreePlusOptAddGenericCommand(c,0);
}

void layeredRangeTreePlusOpt2DBuildTreeCommand(redisClient *c) {
  sortedArray22DBuildCommand(c);
  tfkLayeredRangeTreePlusOptBuildTreePlusOptCommand(c);
}

void layeredRangeTreePlusOpt2DRangeSearchCommand(redisClient *c) {
  tfkLayeredRangeTreePlusOpt2DRangeSearchCommand(c);
}

