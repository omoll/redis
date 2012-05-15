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
static int QUAD_TREE_NODE_SIZE = 3;

double tfkLayeredRangeTreePlus_get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

/*
  Type Defs
*/
typedef struct layeredRangeTreePlusKey{
  double x;
  double y;
  void* value;
} layeredRangeTreePlusKey;


typedef struct cascading_node {
  int size; // array sizes
  double *ys;// eg ys
  double *xs; //eg x coords
  int *left_successor; //positions for cascading left (you keep track
  int *right_successor; // positions for cascading right

  struct cascading_node *left;
  struct cascading_node *right;
} cascading_node_t;

// the inner level of the tree.
typedef struct layeredRangeTreePlusNodeLevel2 {
  // The rectangle corresponding to this node.
  double x1;
  double x2;
  double y1; 
  double y2;

  // a count of elments in the subtree for use in rebalancing.
  int elementCount;

  // TODO(tfkLayeredRangeTreePlus): start is almost always 0.
  // start and end of range in elements array.
  int start;
  int end;
  layeredRangeTreePlusKey* elements; // undefined if not leaf.

  cascading_node_t cnode;
  
  struct layeredRangeTreePlusNodeLevel2** children; // NULL if leaf
} layeredRangeTreePlusNodeLevel2;

// the top level of the tree.
typedef struct layeredRangeTreePlusNodeLevel1 {
  // The rectangle corresponding to this node.
  double x1;
  double x2;
  double y1; 
  double y2;

  // a count of elments in the subtree for use in rebalancing.
  int elementCount;

  // TODO(tfkLayeredRangeTreePlus): start is almost always 0.
  // start and end of range in elements array.
  int start;
  int end;
  layeredRangeTreePlusKey* elements; // undefined if not leaf.

  struct layeredRangeTreePlusNodeLevel1** children; // NULL if leaf

  // a pointer to the second level of the tree containing the same points as in this subtree
  // sorted on the other coordinate.
  layeredRangeTreePlusNodeLevel2* secondLevelPointer;
} layeredRangeTreePlusNodeLevel1;

int min(int a, int b){
  return (a < b)? a : b;
}


/*the value at *end is meant to be a sentinel MAX_INT
 *gives you the succesor position (smallest elt greater than or equal to value)
 */
double* binary_successor_search(double value, double* start, double* end){

  double *initial_start = start;

  while(start < end){
    double * middle =  start + (end-start)/2;
    if ( value < *middle ) {
      end = middle;
    } else if ( value > *middle ){
      start = middle + 1;
    } else {
      while(middle >= initial_start && *middle == value  ){
	middle--;
      }
      
      if (middle != initial_start){
	return middle + 1;
      } else {
	return middle;
      }

    }
  }

  return end;
}

/*
  binary searches for the index of the successor (greater than or
  equal) to the value in the node array, 
  only meant to be used on the root.
*/
int arg_successor_root(double yvalue, cascading_node_t *root){
  return binary_successor_search(yvalue, root->ys, root->ys + root->size) - root->ys;
}

/*yindex is an index for on the parent
 */
int arg_successor_left(int yindex, cascading_node_t *parent){
  return parent->left_successor[yindex];
}

int arg_successor_right(int yindex, cascading_node_t *parent){
  return parent->right_successor[yindex];
}


int cascading_node_is_leaf(cascading_node_t *cn){
  return cn->left == NULL;
}

void check_rep_cascading_node(cascading_node_t *c){
  assert(c->size > 0); 

  //leaf condition all or none should be NULL
  if (cascading_node_is_leaf(c)){
    assert(c->left_successor == NULL);
    assert(c->left == NULL); 
    assert(c->right_successor == NULL);
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

  if (!cascading_node_is_leaf(c)){
    //sentinels inited to -1
    assert(*( c->left_successor - 1 ) == -1);
    assert(*( c->right_successor - 1 ) == -1);

      int max_size_left = 0; int max_size_right = 0;
      for (i = 0; i < c->size; i++){
	//check indices into left array all make sense
	assert( c->left_successor[i] > -1 && c->left_successor[i] >= c->left_successor[i-1]
		&& c->left_successor[i] <= c->left->size);
	if (c->left_successor[i] == c->left->size) max_size_left = 1;

	//check indices into right child all make sense
	assert(c->right_successor[i] > -1);
	assert(c->right_successor[i] >= c->right_successor[i-1]);
	assert(c->right_successor[i] <= c->right->size);

	if (c->right_successor[i] == c->right->size) max_size_right = 1;
      }

      assert(max_size_right ^ max_size_left);
  }




}

void cascading_node_leaf_init(cascading_node_t *cn, double *x, double *y, int size){
  cn->size = size;
  cn->xs = x;
  cn->ys = (double*) zmalloc((size+1)*sizeof(double));
  memcpy(cn->ys, y, size*sizeof(double));
  cn->ys[size] = DBL_MAX;

  cn->left_successor = NULL;
  cn->right_successor = NULL;
  cn->left  = NULL;
  cn->right  = NULL;

  check_rep_cascading_node(cn);
  assert(cascading_node_is_leaf(cn));
}

//inits cascading node cn,  given its two children and the values at a
//(normal)node 
void cascading_node_parent_init(cascading_node_t *cn, cascading_node_t *left, cascading_node_t *right){
  cn->size = left->size + right->size; //has all values of
						 //left all of right
						 //and the new values

  //allocates +1 size for sentinel value.
  cn->ys = (double *)zmalloc((cn->size+1)*sizeof(double));
  cn->ys[cn->size] = DBL_MAX;

  //does not need any sentinel value, simply the other coordinate
  cn->xs = (double *)zmalloc(cn->size*sizeof(double));

  //allocate an array of size +1 length, the -1 entry gets inited to 0
  cn->left_successor = (int *)zmalloc((cn->size+1)*sizeof(int)) + 1;
  *(cn->left_successor - 1) = -1;
  cn->right_successor = (int *)zmalloc((cn->size+1)*sizeof(int)) + 1;
  *(cn->right_successor - 1) = -1;

  cn->left = left;
  cn->right = right;

  
  
  

  int i, j;
  for (i = 0, j = 0;  i < left->size || j < right->size; ){
    if (left->ys[i] < right->ys[j]){ //left has next smallest elt.
      cn->ys[i+j] = left->ys[i];
      cn->xs[i+j] = left->xs[i];

      /* cn->left_successor[i+j] = i; */

      /* cn->right_successor[i+j] = cn->right_successor[i+j-1] + 1; */
      /* //this is wrong. too. */

      i++; 
    } else { //right has next smallest
      cn->ys[i+j] = right->ys[j];
      cn->xs[i+j] = right->xs[j];

      /* cn->left_successor[i+j] = cn->left_successor[i+j-1] + 1; */
      /* //this may be wrong,  because it keeps increasing the other one,  */
      /* // but  when it maxes out it shoudl no longer be increased. */

      /* cn->right_successor[i+j] = j; */
      j++;
    }
  }

  for (i = 0; i < cn->size; i++){
    cn->left_successor[i] = arg_successor_root(cn->ys[i], cn->left);
    cn->right_successor[i] = arg_successor_root(cn->ys[i], cn->right);
  }

  check_rep_cascading_node(cn);
  assert(!cascading_node_is_leaf(cn));
}



/*
 *@a and b are  elements in the two siblings.
 */

void cascading_node_print(cascading_node_t *cn){
  printf("size: %d\n", cn->size);

  printf("ys:");
  int i;
  for (i = 0; i < cn->size; i ++){
    printf("%1.6f,", cn->ys[i]);
  }
  printf("\n");

  printf("xs:");
  for (i = 0; i < cn->size; i ++){
    printf("%1.6f,", cn->xs[i]);
  }
  printf("\n");

  if (cascading_node_is_leaf(cn)){
    printf("(leaf ends)\n");
  } else {
    printf("left indices:");
    for (i = 0; i < cn->size; i ++){
      printf("%d,", cn->left_successor[i]);
    }

    printf("right indices:");
    for (i = 0; i < cn->size; i ++){
      printf("%d,", cn->right_successor[i]);
    }
  }
}




void layeredRangeTreePlusDeleteNodeLevel2(layeredRangeTreePlusNodeLevel2* n) {
  if (n->children == NULL) {
    zfree(n->elements);
    zfree(n);
  } else {
    layeredRangeTreePlusDeleteNodeLevel2(n->children[0]);
    layeredRangeTreePlusDeleteNodeLevel2(n->children[1]);
    zfree(n->children);
    zfree(n);
  }
}

// Delete node n and it's associated subree.
void layeredRangeTreePlusDeleteNodeLevel1(layeredRangeTreePlusNodeLevel1* n){
  if (n->children == NULL) {
    zfree(n->elements);
    zfree(n);
  } else {
    layeredRangeTreePlusDeleteNodeLevel1(n->children[0]);
    layeredRangeTreePlusDeleteNodeLevel1(n->children[1]);
    zfree(n->children);
    zfree(n);
  }

  if (n->secondLevelPointer != NULL){
    layeredRangeTreePlusDeleteNodeLevel2(n->secondLevelPointer);
  }
}

/*
  Comparison functions for median finding and sorting.
*/
int layeredRangeTreePlusElementCompareX(const void* a, const void* b){
  if (((layeredRangeTreePlusKey*) a)->x < ((layeredRangeTreePlusKey*) b)->x) {
    return -1;
  }
  if (((layeredRangeTreePlusKey*) a)->x > ((layeredRangeTreePlusKey*) b)->x) {
    return 1;
  }
  return 0;
}

int layeredRangeTreePlusElementCompareY(const void* a, const void* b){
  if (((layeredRangeTreePlusKey*) a)->y < ((layeredRangeTreePlusKey*) b)->y) {
    return -1;
  }
  if (((layeredRangeTreePlusKey*) a)->y > ((layeredRangeTreePlusKey*) b)->y) {
    return 1;
  }
  return 0;
}

void buildLayeredRangeTreePlusNodeLevel2 (layeredRangeTreePlusNodeLevel2* n, layeredRangeTreePlusKey* elements, int elementCount) {
  n->elementCount = elementCount;

  // Coarsen 
  if (n->elementCount <= QUAD_TREE_NODE_SIZE) {
    n->start = 0; 
    n->end = n->elementCount;
    n->elements = (layeredRangeTreePlusKey*) zmalloc(sizeof(layeredRangeTreePlusKey) * QUAD_TREE_NODE_SIZE);
    for (int i = 0; i < n->elementCount; i++) {
      n->elements[i] = elements[i];
    }
    return;
  }
 
  n->children = (layeredRangeTreePlusNodeLevel2**) zmalloc(sizeof(layeredRangeTreePlusNodeLevel2*)*2);

  qsort(elements, n->elementCount,
      sizeof(layeredRangeTreePlusKey), layeredRangeTreePlusElementCompareY);
  layeredRangeTreePlusKey medianY = elements[n->elementCount / 2];

  n->children[0] = (layeredRangeTreePlusNodeLevel2*) zmalloc(sizeof(layeredRangeTreePlusNodeLevel2));
  n->children[1] = (layeredRangeTreePlusNodeLevel2*) zmalloc(sizeof(layeredRangeTreePlusNodeLevel2));

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

  buildLayeredRangeTreePlusNodeLevel2(n->children[0], elements, n->children[0]->elementCount);
  buildLayeredRangeTreePlusNodeLevel2(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount); 
}

// builds a second level of the quad tree.
void explodeTreePlus(layeredRangeTreePlusNodeLevel1* n, layeredRangeTreePlusKey* elements, int elementCount){
  n->secondLevelPointer = NULL;
  if (n->elementCount <= QUAD_TREE_NODE_SIZE) {
    assert(elementCount == n->elementCount);
    n->secondLevelPointer = (layeredRangeTreePlusNodeLevel2*) zmalloc(sizeof(layeredRangeTreePlusNodeLevel2));

    // this is a leaf node, copy its elements into the array.
    memcpy(elements, n->elements, elementCount * sizeof(layeredRangeTreePlusKey));

    qsort(elements, n->elementCount,
        sizeof(layeredRangeTreePlusKey), layeredRangeTreePlusElementCompareY);

    double* xs = (double*) zmalloc(sizeof(double)*elementCount);  
    double* ys = (double*) zmalloc(sizeof(double)*elementCount);

    assert(xs != NULL && ys != NULL);

    for (int i = 0; i < elementCount; i++) {
      xs[i] = elements[i].x;
      ys[i] = elements[i].y;
    }

    cascading_node_leaf_init(&n->secondLevelPointer->cnode, xs, ys, elementCount); 

    return;
  }
 
  explodeTreePlus(n->children[0], elements, n->children[0]->elementCount); 
  explodeTreePlus(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount);

  n->secondLevelPointer = (layeredRangeTreePlusNodeLevel2*) zmalloc(sizeof(layeredRangeTreePlusNodeLevel2));

  cascading_node_parent_init(&n->secondLevelPointer->cnode, &n->children[0]->secondLevelPointer->cnode, &n->children[1]->secondLevelPointer->cnode);

  n->secondLevelPointer->x1 = n->x1;
  n->secondLevelPointer->x2 = n->x2;
  n->secondLevelPointer->y1 = n->y1;
  n->secondLevelPointer->y2 = n->y2;
  n->secondLevelPointer->children = NULL;
  n->secondLevelPointer->elements = NULL;
  n->secondLevelPointer->elementCount = n->elementCount;
  n->secondLevelPointer->start = 0;
  n->secondLevelPointer->end = 0;
  
  //buildLayeredRangeTreePlusNodeLevel2(n->secondLevelPointer, elements, n->elementCount);
  assert(n->secondLevelPointer != NULL);
}

void buildLayeredRangeTreePlusNodeLevel1(layeredRangeTreePlusNodeLevel1* n, layeredRangeTreePlusKey* elements, int elementCount) {
  n->elementCount = elementCount;
  n->secondLevelPointer = NULL;

  // Coarsen 
  if (n->elementCount <= QUAD_TREE_NODE_SIZE) {
    n->start = 0; 
    n->end = n->elementCount;
    n->elements = (layeredRangeTreePlusKey*) zmalloc(sizeof(layeredRangeTreePlusKey) * QUAD_TREE_NODE_SIZE);
    for (int i = 0; i < n->elementCount; i++) {
      n->elements[i] = elements[i];
    }
    return;
  }
 
  n->children = (layeredRangeTreePlusNodeLevel1**) zmalloc(sizeof(layeredRangeTreePlusNodeLevel1*)*2);

  qsort(elements, n->elementCount,
      sizeof(layeredRangeTreePlusKey), layeredRangeTreePlusElementCompareX);

  layeredRangeTreePlusKey medianX = elements[n->elementCount / 2];

  n->children[0] = (layeredRangeTreePlusNodeLevel1*) zmalloc(sizeof(layeredRangeTreePlusNodeLevel1));
  n->children[1] = (layeredRangeTreePlusNodeLevel1*) zmalloc(sizeof(layeredRangeTreePlusNodeLevel1));

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

  buildLayeredRangeTreePlusNodeLevel1(n->children[0], elements, n->children[0]->elementCount);
  buildLayeredRangeTreePlusNodeLevel1(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount); 
}

// returns true if a key falls within a given range.
bool layeredRangeTreePlusKeyInRange(layeredRangeTreePlusKey key, double x1, double x2, double y1, double y2) {
  if (key.x >= x1 && key.x < x2 && key.y > y1 && key.y <= y2) {
    return true;
  } else {
    return false;
  }
}

// returns true if a node intersects with the given range, false otherwise.
bool layeredRangeTreePlusNodeIntersectsRangeLevel1(layeredRangeTreePlusNodeLevel1* n, double x1, double x2, double y1, double y2) {
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
bool layeredRangeTreePlusNodeIntersectsRangeLevel2(layeredRangeTreePlusNodeLevel2* n, double x1, double x2, double y1, double y2) {
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


void layeredRangeTreePlusRangeSearchLevel2(layeredRangeTreePlusNodeLevel2 *n, double x1, double x2, double y1, double y2, int* count) {
  if (n->children == NULL) {
    assert(n->elementCount <= QUAD_TREE_NODE_SIZE);
    for (int i = n->start; i < n->end; i++) {
      if (layeredRangeTreePlusKeyInRange(n->elements[i], -200, 200, y1, y2)){
        (*count)++;
        //printf("Element (%f, %f) in range \n", n->elements[i].x, n->elements[i].y);
      }
    }
    return;
  }
 
  // call recursively on each child intersected by the range.
  if (layeredRangeTreePlusNodeIntersectsRangeLevel2(n->children[0], -200, 200, y1, y2)) {
    layeredRangeTreePlusRangeSearchLevel2(n->children[0], x1, x2, y1, y2, count);
  }

  if (layeredRangeTreePlusNodeIntersectsRangeLevel2(n->children[1], -200, 200, y1, y2)) {
    layeredRangeTreePlusRangeSearchLevel2(n->children[1], x1, x2, y1, y2, count);
  }
}

void layeredRangeTreePlusRangeSearchLevel1(layeredRangeTreePlusNodeLevel1 *n, double x1, double x2, double y1, double y2, int* count, int offset) {
  if (offset == -1){
    // do the initial binary search.
    offset = arg_successor_root(y1, &n->secondLevelPointer->cnode);
  }

  if (n->children == NULL) {
    assert(n->elementCount <= QUAD_TREE_NODE_SIZE);
    // if we get to a leaf just scan through the elements and report the ones in the range.
    for (int i = n->start; i < n->end; i++) {
      if (layeredRangeTreePlusKeyInRange(n->elements[i], x1, x2, y1, y2)){
        (*count)++;
        //printf("Element (%f, %f) in range \n", n->elements[i].x, n->elements[i].y);
      }
    }
    return;
  }

  assert(n->secondLevelPointer != NULL);

  // call the second level search if we're entirely within the x-range.
  if (n->x1 >= x1 && n->x2 <= x2) {
    for (;offset < n->secondLevelPointer->cnode.size; offset++ ){
      printf("Element (%f, %f) in range \n", n->secondLevelPointer->cnode.xs[offset], n->secondLevelPointer->cnode.ys[offset]);
      (*count)++;
    }
    //layeredRangeTreePlusRangeSearchLevel2(n->secondLevelPointer, -200, 200, y1, y2, count);
    return;
  }

  // otherwise call recursively on each child intersected by the range.
  if (layeredRangeTreePlusNodeIntersectsRangeLevel1(n->children[0], x1, x2, -200, 200)) {
    layeredRangeTreePlusRangeSearchLevel1(n->children[0], x1, x2, y1, y2, count, 
        arg_successor_left(offset, &n->secondLevelPointer->cnode));
  }

  if (layeredRangeTreePlusNodeIntersectsRangeLevel1(n->children[1], x1, x2, -200, 200)) {
    layeredRangeTreePlusRangeSearchLevel1(n->children[1], x1, x2, y1, y2, count,
        arg_successor_right(offset, &n->secondLevelPointer->cnode));
  }
}


// the root node
layeredRangeTreePlusNodeLevel1* rangeTreePlusRoot;

// collection of elements
layeredRangeTreePlusKey* allElements;
int allElementsCount;
int currentAllElementsArraySize;
double totalQueryTime;

void tfkLayeredRangeTreePlus2DRangeSearchCommand(redisClient* c) {
  double x1 = strtod(c->argv[1]->ptr, NULL);
  double x2 = strtod(c->argv[2]->ptr, NULL);
  double y1 = strtod(c->argv[3]->ptr, NULL);
  double y2 = strtod(c->argv[4]->ptr, NULL);
  int count = 0;
  double start_time = tfkLayeredRangeTreePlus_get_time();
  layeredRangeTreePlusRangeSearchLevel1(rangeTreePlusRoot, x1, x2, y1, y2, &count, -1);
  double end_time = tfkLayeredRangeTreePlus_get_time();
  totalQueryTime += end_time - start_time;
  printf("total query time %f rangeTreePlusCount: %d \n", totalQueryTime, count); 
  addReplyLongLong(c, count);
}

void tfkLayeredRangeTreePlusBuildTreePlusCommand(redisClient *c) {
  if (rangeTreePlusRoot != NULL){
    // delete the old quad tree.
    layeredRangeTreePlusDeleteNodeLevel1(rangeTreePlusRoot);
  }
  rangeTreePlusRoot = (layeredRangeTreePlusNodeLevel1*) zmalloc(sizeof(layeredRangeTreePlusNodeLevel1));

  rangeTreePlusRoot->x1 = -200;
  rangeTreePlusRoot->x2 = 200;
  rangeTreePlusRoot->y1 = -200;
  rangeTreePlusRoot->y2 = 200; 
  rangeTreePlusRoot->children = NULL;
  rangeTreePlusRoot->elementCount = allElementsCount;
  rangeTreePlusRoot->secondLevelPointer = NULL;
  // rebuild using the allElements array.

  // build the first level of the tree.
  buildLayeredRangeTreePlusNodeLevel1(rangeTreePlusRoot, allElements, allElementsCount);

  // build the second level of the tree.
  layeredRangeTreePlusKey* elementsCopy = (layeredRangeTreePlusKey*) zmalloc(sizeof(layeredRangeTreePlusKey) * rangeTreePlusRoot->elementCount);
  explodeTreePlus(rangeTreePlusRoot, elementsCopy, allElementsCount);
  zfree(elementsCopy); 
  
  //walkLayeredRangeTreePlus(rangeTreePlusRoot);

  addReplyLongLong(c, 42);
}

/* This generic command implements both ZADD and ZINCRBY. */
void tfkLayeredRangeTreePlusAddGenericCommand(redisClient *c, int incr) {
    double x = strtod(c->argv[1]->ptr, NULL);
    double y = strtod(c->argv[2]->ptr, NULL);
    
    int resizeAmount = 10000;
  
    // create a layeredRangeTreePlusKey
    layeredRangeTreePlusKey key;
    key.x = x;
    key.y = y;
    key.value = NULL; 
    if (allElementsCount == currentAllElementsArraySize) {
      // resize allElementsArray
      currentAllElementsArraySize += resizeAmount;
      allElements = (layeredRangeTreePlusKey*) zrealloc((void*)allElements, sizeof(layeredRangeTreePlusKey) * currentAllElementsArraySize);
    }
    allElements[allElementsCount] = key;
    allElementsCount++;
    addReplyLongLong(c, allElementsCount);
}

void layeredRangeTreePlus2DAddCommand(redisClient *c) {
  tfkLayeredRangeTreePlusAddGenericCommand(c,0);
}

void layeredRangeTreePlus2DBuildTreeCommand(redisClient *c) {
  tfkLayeredRangeTreePlusBuildTreePlusCommand(c);
}

void layeredRangeTreePlus2DRangeSearchCommand(redisClient *c) {
  tfkLayeredRangeTreePlus2DRangeSearchCommand(c);
}

