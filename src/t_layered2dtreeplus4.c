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
static int LAYERED_TREE_PLUS4_NODE_SIZE = 64;

double tfkLayeredRangeTreePlus4_get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

int min4(int a, int b){
  return (a < b)? a : b;
}

typedef struct cascading4_node {
  int size; // array sizes
  double *ys;// eg ys
  double *xs; //eg x coords
  int *left_successor4; //positions for cascading4 left (you keep track
  int *right_successor4; // positions for cascading4 right

  int *left2_successor4;
  int *right2_successor4;

  struct cascading4_node *left;
  struct cascading4_node *left2;
  struct cascading4_node *right;
  struct cascading4_node *right2;

} cascading4_node_t;

void cascading4_node_print(cascading4_node_t *cn);

/*the value at *end is meant to be a sentinel MAX_INT
 *gives you the succesor position (smallest elt greater than or equal to value)
 */

double* linear_successor4_search4(double value, double* start, double* end){
 while( !(*start >= value)){
   start++;
 }
 return start;
}

double* recursive_successor4_search4(double value, double* start, double* end){
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
    double* ret = recursive_successor4_search4(value, start, median);
    if (*ret == DBL_MAX) {
      return median;
    }
    return ret;
  } else {
    double* ret = recursive_successor4_search4(value, median, end);
    return ret;
  }
}

/*CURRENTLY BUGGY, CHECK TEST CASES.
 */
double* binary_successor4_search4(double value, double* start, double* end){

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
  binary search4es for the index of the successor4 (greater than or
  equal) to the value in the node array, 
  only meant to be used on the root.
*/
int arg_successor4_root(double yvalue, cascading4_node_t *root){
  //TODO:when debugged binary search4 and place here.
  return recursive_successor4_search4(yvalue, root->ys, root->ys + root->size) - root->ys;
}

/*yindex is an index for on the parent
 */
int arg_successor4_left(int yindex, cascading4_node_t *parent){
  if (yindex < parent->size){
    return parent->left_successor4[yindex];
  } else {
    return parent->left->size;
  }
  /* int y = arg_successor4_root(parent->ys[yindex], parent->left); */
  /* assert(x == y); */
  /* return y; */
}

/*yindex is an index for on the parent
 */
int arg_successor4_left2(int yindex, cascading4_node_t *parent){
  if (yindex < parent->size){
    return parent->left2_successor4[yindex];
  } else {
    return parent->left2->size;
  }
  /* int y = arg_successor4_root(parent->ys[yindex], parent->left); */
  /* assert(x == y); */
  /* return y; */
}

int arg_successor4_right(int yindex, cascading4_node_t *parent){
  if (yindex < parent->size){
    return parent->right_successor4[yindex];
  } else {
    return parent->right->size;
  }
  /* int y = arg_successor4_root(parent->ys[yindex], parent->right); */
  /* assert(x== y); */
  /* return y; */
}

int arg_successor4_right2(int yindex, cascading4_node_t *parent){
  if (yindex < parent->size){
    return parent->right2_successor4[yindex];
  } else {
    return parent->right2->size;
  }
  /* int y = arg_successor4_root(parent->ys[yindex], parent->right); */
  /* assert(x== y); */
  /* return y; */
}

int cascading4_node_is_leaf(cascading4_node_t *cn){
  return cn->left == NULL;
}

void check_rep_cascading4_node(cascading4_node_t *c){
  assert(c->size > 0); 

  //leaf condition all or none should be NULL
  if (cascading4_node_is_leaf(c)){
    assert(c->left_successor4 == NULL);
    assert(c->left == NULL); 
    assert(c->right_successor4 == NULL);
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

  if (!cascading4_node_is_leaf(c)){
    //sentinels inited to -1
    assert(*( c->left_successor4 - 1 ) == -1);
    assert(*( c->right_successor4 - 1 ) == -1);

      for (i = 0; i < c->size; i++){
	//check indices into left array all make sense
	assert( c->left_successor4[i] > -1 && c->left_successor4[i] >= c->left_successor4[i-1]
		&& c->left_successor4[i] <= c->left->size);

	//check indices into right child all make sense
	assert(c->right_successor4[i] > -1);
	assert(c->right_successor4[i] >= c->right_successor4[i-1]);
	assert(c->right_successor4[i] <= c->right->size);
      }


  }

}

void cascading4_node_leaf_init(cascading4_node_t *cn, double *x, double *y, int size){
  cn->size = size;
  cn->xs = x;
  cn->ys = (double*) zmalloc((size+1)*sizeof(double));
  memcpy(cn->ys, y, size*sizeof(double));
  cn->ys[size] = DBL_MAX;

  cn->left_successor4 = NULL;
  cn->left2_successor4 = NULL;
  cn->right_successor4 = NULL;
  cn->right2_successor4 = NULL;

  cn->left  = NULL;
  cn->left2 = NULL;
  cn->right  = NULL;
  cn->right2 = NULL;

  check_rep_cascading4_node(cn);
  assert(cascading4_node_is_leaf(cn));
}

//inits cascading4 node cn,  given its two children and the values at a
//(normal)node 
void cascading4_node_parent_init(cascading4_node_t *cn, cascading4_node_t *left, cascading4_node_t* left2,
    cascading4_node_t *right, cascading4_node_t* right2){
  cn->size = left->size + left2->size + right->size + right2->size; //has all values of
						 //left all of right
						 //and the new values

  //allocates +1 size for sentinel value.
  cn->ys = (double *)zmalloc((cn->size+1)*sizeof(double));
  cn->ys[cn->size] = DBL_MAX;

  //does not need any sentinel value, simply the other coordinate
  cn->xs = (double *)zmalloc(cn->size*sizeof(double));

  //allocate an array of size +1 length, the -1 entry gets inited to 0
  cn->left_successor4 = (int *) zmalloc((cn->size+1)*sizeof(int)) + 1;
  *(cn->left_successor4 - 1) = -1;
  cn->left2_successor4 = (int *) zmalloc((cn->size+1)*sizeof(int)) + 1;
  *(cn->left2_successor4 - 1) = -1;

  cn->right_successor4 = (int *) zmalloc((cn->size+1)*sizeof(int)) + 1;
  *(cn->right_successor4 - 1) = -1;
  cn->right2_successor4 = (int *) zmalloc((cn->size+1)*sizeof(int)) + 1;
  *(cn->right2_successor4 - 1) = -1;

  cn->left = left;
  cn->left2 = left2;
  cn->right = right;
  cn->right2 = right2;

  int i, i2, j, j2;
  for (i = 0, i2 = 0, j = 0, j2 = 0;  i < left->size || i2 < left2->size || j < right->size || j2 < right2->size;){
    double* min4child = left->ys; 
    double min4_value = left->ys[i];
    
    if (min4_value > left2->ys[i2]) {
      min4_value = left2->ys[i2];
      min4child = left2->ys;
    }
    if (min4_value > right->ys[j]){
      min4_value = right->ys[j];
      min4child = right->ys;
    }
    if (min4_value > right2->ys[j2]){
      min4_value = right2->ys[j2];
      min4child = right2->ys;
    }

    if (min4child == left->ys){
      // left has the next smallest element
      cn->ys[i+i2+j+j2] = left->ys[i];
      cn->xs[i+i2+j+j2] = left->xs[i];
      i++;
    }
    if (min4child == left2->ys) {
      cn->ys[i+i2+j+j2] = left2->ys[i2];
      cn->xs[i+i2+j+j2] = left2->xs[i2];
      i2++;
    }
    if (min4child == right->ys){
      cn->ys[i+i2+j+j2] = right->ys[j];
      cn->xs[i+i2+j+j2] = right->xs[j];
      j++;
    }
    if (min4child == right2->ys){
      cn->ys[i+i2+j+j2] = right2->ys[j2];
      cn->xs[i+i2+j+j2] = right2->xs[j2];
      j2++;
    }
  }
  for (i = 0; i < cn->size; i++){
    cn->left_successor4[i] = arg_successor4_root(cn->ys[i], cn->left);
    cn->left2_successor4[i] = arg_successor4_root(cn->ys[i], cn->left2);
    cn->right_successor4[i] = arg_successor4_root(cn->ys[i], cn->right);
    cn->right2_successor4[i] = arg_successor4_root(cn->ys[i], cn->right2);
  }

  //  cascading4_node_print(cn);
  check_rep_cascading4_node(cn);
  assert(!cascading4_node_is_leaf(cn));
}

/*
  Type Defs
*/
typedef struct layeredRangeTreePlus4Key{
  double x;
  double y;
  void* value;
} layeredRangeTreePlus4Key;

// the inner level of the tree.
typedef struct layeredRangeTreePlus4NodeLevel2 {
  // The rectangle corresponding to this node.
  double x1;
  double x2;
  double y1; 
  double y2;

  // a count of elments in the subtree for use in rebalancing.
  int elementCount;

  // TODO(tfkLayeredRangeTreePlus4): start is almost always 0.
  // start and end of range in elements array.
  int start;
  int end;
  layeredRangeTreePlus4Key* elements; // undefined if not leaf.

  cascading4_node_t cnode;
  
  struct layeredRangeTreePlus4NodeLevel2** children; // NULL if leaf
} layeredRangeTreePlus4NodeLevel2;

// the top level of the tree.
typedef struct layeredRangeTreePlus4NodeLevel1 {
  // The rectangle corresponding to this node.
  double x1;
  double x2;
  double y1; 
  double y2;

  // a count of elments in the subtree for use in rebalancing.
  int elementCount;

  // TODO(tfkLayeredRangeTreePlus4): start is almost always 0.
  // start and end of range in elements array.
  int start;
  int end;
  layeredRangeTreePlus4Key* elements; // undefined if not leaf.

  struct layeredRangeTreePlus4NodeLevel1** children; // NULL if leaf

  // a pointer to the second level of the tree containing the same points as in this subtree
  // sorted on the other coordinate.
  layeredRangeTreePlus4NodeLevel2* secondLevelPointer;
} layeredRangeTreePlus4NodeLevel1;





void layeredRangeTreePlus4DeleteNodeLevel2(layeredRangeTreePlus4NodeLevel2* n) {
  if (n->children == NULL) {
    zfree(n->elements);
    zfree(n);
  } else {
    layeredRangeTreePlus4DeleteNodeLevel2(n->children[0]);
    layeredRangeTreePlus4DeleteNodeLevel2(n->children[1]);
    layeredRangeTreePlus4DeleteNodeLevel2(n->children[2]);
    layeredRangeTreePlus4DeleteNodeLevel2(n->children[3]);
    zfree(n->children);
    zfree(n);
  }
}

// Delete node n and it's associated subree.
void layeredRangeTreePlus4DeleteNodeLevel1(layeredRangeTreePlus4NodeLevel1* n){
  if (n->children == NULL) {
    zfree(n->elements);
    zfree(n);
  } else {
    layeredRangeTreePlus4DeleteNodeLevel1(n->children[0]);
    layeredRangeTreePlus4DeleteNodeLevel1(n->children[1]);
    layeredRangeTreePlus4DeleteNodeLevel1(n->children[2]);
    layeredRangeTreePlus4DeleteNodeLevel1(n->children[3]);

    zfree(n->children);
    zfree(n);
  }

  if (n->secondLevelPointer != NULL){
    layeredRangeTreePlus4DeleteNodeLevel2(n->secondLevelPointer);
  }
}

/*
  Comparison functions for median finding and sorting.
*/
int layeredRangeTreePlus4ElementCompareX(const void* a, const void* b){
  if (((layeredRangeTreePlus4Key*) a)->x < ((layeredRangeTreePlus4Key*) b)->x) {
    return -1;
  }
  if (((layeredRangeTreePlus4Key*) a)->x > ((layeredRangeTreePlus4Key*) b)->x) {
    return 1;
  }
  return 0;
}

int layeredRangeTreePlus4ElementCompareY(const void* a, const void* b){
  if (((layeredRangeTreePlus4Key*) a)->y < ((layeredRangeTreePlus4Key*) b)->y) {
    return -1;
  }
  if (((layeredRangeTreePlus4Key*) a)->y > ((layeredRangeTreePlus4Key*) b)->y) {
    return 1;
  }
  return 0;
}

void buildLayeredRangeTreePlus4NodeLevel2 (layeredRangeTreePlus4NodeLevel2* n, layeredRangeTreePlus4Key* elements, int elementCount) {
  n->elementCount = elementCount;

  // Coarsen 
  if (n->elementCount <= LAYERED_TREE_PLUS4_NODE_SIZE) {
    n->start = 0; 
    n->end = n->elementCount;
    n->elements = (layeredRangeTreePlus4Key*) zmalloc(sizeof(layeredRangeTreePlus4Key) * LAYERED_TREE_PLUS4_NODE_SIZE);
    for (int i = 0; i < n->elementCount; i++) {
      n->elements[i] = elements[i];
    }
    return;
  }
 
  n->children = (layeredRangeTreePlus4NodeLevel2**) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel2*)*4);

  qsort(elements, n->elementCount,
      sizeof(layeredRangeTreePlus4Key), layeredRangeTreePlus4ElementCompareY);

  layeredRangeTreePlus4Key medianY1 = elements[n->elementCount / 4];
  layeredRangeTreePlus4Key medianY2 = elements[n->elementCount / 2];
  layeredRangeTreePlus4Key medianY3 = elements[(3*n->elementCount) / 4];

  n->children[0] = (layeredRangeTreePlus4NodeLevel2*) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel2));
  n->children[1] = (layeredRangeTreePlus4NodeLevel2*) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel2));
  n->children[2] = (layeredRangeTreePlus4NodeLevel2*) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel2));
  n->children[3] = (layeredRangeTreePlus4NodeLevel2*) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel2));

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

  buildLayeredRangeTreePlus4NodeLevel2(n->children[0], elements, n->children[0]->elementCount);
  buildLayeredRangeTreePlus4NodeLevel2(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount); 
  buildLayeredRangeTreePlus4NodeLevel2(n->children[2], elements + n->children[0]->elementCount + n->children[1]->elementCount, n->children[1]->elementCount); 
  buildLayeredRangeTreePlus4NodeLevel2(n->children[3],
      elements + n->children[0]->elementCount + n->children[1]->elementCount + n->children[2]->elementCount, n->children[2]->elementCount); 
}

// builds a second level of the quad tree.
void explodeTreePlus4(layeredRangeTreePlus4NodeLevel1* n, layeredRangeTreePlus4Key* elements, int elementCount){
  n->secondLevelPointer = NULL;
  if (n->elementCount <= LAYERED_TREE_PLUS4_NODE_SIZE) {
    assert(elementCount == n->elementCount);
    n->secondLevelPointer = (layeredRangeTreePlus4NodeLevel2*) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel2));

    // this is a leaf node, copy its elements into the array.
    memcpy(elements, n->elements, elementCount * sizeof(layeredRangeTreePlus4Key));

    qsort(elements, n->elementCount,
        sizeof(layeredRangeTreePlus4Key), layeredRangeTreePlus4ElementCompareY);

    double* xs = (double*) zmalloc(sizeof(double)*elementCount);  
    double* ys = (double*) zmalloc(sizeof(double)*elementCount);

    assert(xs != NULL && ys != NULL);

    for (int i = 0; i < elementCount; i++) {
      xs[i] = elements[i].x;
      ys[i] = elements[i].y;
    }

    cascading4_node_leaf_init(&n->secondLevelPointer->cnode, xs, ys, elementCount); 

    return;
  }
 
  explodeTreePlus4(n->children[0], elements, n->children[0]->elementCount); 
  explodeTreePlus4(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount);
  explodeTreePlus4(n->children[2],
      elements + n->children[0]->elementCount + n->children[1]->elementCount, n->children[2]->elementCount);
  explodeTreePlus4(n->children[3],
      elements + n->children[0]->elementCount + n->children[1]->elementCount + n->children[2]->elementCount, n->children[3]->elementCount);


  n->secondLevelPointer = (layeredRangeTreePlus4NodeLevel2*) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel2));

  cascading4_node_parent_init(&n->secondLevelPointer->cnode, &n->children[0]->secondLevelPointer->cnode, &n->children[1]->secondLevelPointer->cnode,
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
  
  //buildLayeredRangeTreePlus4NodeLevel2(n->secondLevelPointer, elements, n->elementCount);
  assert(n->secondLevelPointer != NULL);
}

void buildLayeredRangeTreePlus4NodeLevel1(layeredRangeTreePlus4NodeLevel1* n, layeredRangeTreePlus4Key* elements, int elementCount) {
  n->elementCount = elementCount;
  n->secondLevelPointer = NULL;

  // Coarsen 
  if (n->elementCount <= LAYERED_TREE_PLUS4_NODE_SIZE) {
    n->start = 0; 
    n->end = n->elementCount;
    n->elements = (layeredRangeTreePlus4Key*) zmalloc(sizeof(layeredRangeTreePlus4Key) * LAYERED_TREE_PLUS4_NODE_SIZE);
    for (int i = 0; i < n->elementCount; i++) {
      n->elements[i] = elements[i];
    }
    return;
  }
 
  n->children = (layeredRangeTreePlus4NodeLevel1**) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel1*)*4);

  qsort(elements, n->elementCount,
      sizeof(layeredRangeTreePlus4Key), layeredRangeTreePlus4ElementCompareX);

  layeredRangeTreePlus4Key medianX1 = elements[n->elementCount / 4];
  layeredRangeTreePlus4Key medianX2 = elements[n->elementCount / 2];
  layeredRangeTreePlus4Key medianX3 = elements[(3*n->elementCount) / 4];

  n->children[0] = (layeredRangeTreePlus4NodeLevel1*) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel2));
  n->children[1] = (layeredRangeTreePlus4NodeLevel1*) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel2));
  n->children[2] = (layeredRangeTreePlus4NodeLevel1*) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel2));
  n->children[3] = (layeredRangeTreePlus4NodeLevel1*) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel2));
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

  buildLayeredRangeTreePlus4NodeLevel1(n->children[0], elements, n->children[0]->elementCount);
  buildLayeredRangeTreePlus4NodeLevel1(n->children[1], elements + n->children[0]->elementCount, n->children[1]->elementCount);
  buildLayeredRangeTreePlus4NodeLevel1(n->children[2],
      elements + n->children[0]->elementCount + n->children[1]->elementCount, n->children[2]->elementCount); 
  buildLayeredRangeTreePlus4NodeLevel1(n->children[3],
      elements + n->children[0]->elementCount + n->children[1]->elementCount + n->children[2]->elementCount, n->children[3]->elementCount); 
}

// returns true if a key falls within a given range.
bool layeredRangeTreePlus4KeyInRange(layeredRangeTreePlus4Key key, double x1, double x2, double y1, double y2) {
  if (key.x >= x1 && key.x < x2 && key.y > y1 && key.y <= y2) {
    return true;
  } else {
    return false;
  }
}

// returns true if a node intersects with the given range, false otherwise.
bool layeredRangeTreePlus4NodeIntersectsRangeLevel1(layeredRangeTreePlus4NodeLevel1* n, double x1, double x2, double y1, double y2) {
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
bool layeredRangeTreePlus4NodeIntersectsRangeLevel2(layeredRangeTreePlus4NodeLevel2* n, double x1, double x2, double y1, double y2) {
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

void layeredRangeTreePlus4RangeSearchLevel2(layeredRangeTreePlus4NodeLevel2 *n, double x1, double x2, double y1, double y2, int* count) {
  if (n->children == NULL) {
    assert(n->elementCount <= LAYERED_TREE_PLUS4_NODE_SIZE);
    for (int i = n->start; i < n->end; i++) {
      if (layeredRangeTreePlus4KeyInRange(n->elements[i], -200, 200, y1, y2)){
        (*count)++;
        //printf("Element (%f, %f) in range \n", n->elements[i].x, n->elements[i].y);
      }
    }
    return;
  }
 
  // call recursively on each child intersected by the range.
  if (layeredRangeTreePlus4NodeIntersectsRangeLevel2(n->children[0], -200, 200, y1, y2)) {
    layeredRangeTreePlus4RangeSearchLevel2(n->children[0], x1, x2, y1, y2, count);
  }

  if (layeredRangeTreePlus4NodeIntersectsRangeLevel2(n->children[1], -200, 200, y1, y2)) {
    layeredRangeTreePlus4RangeSearchLevel2(n->children[1], x1, x2, y1, y2, count);
  }

  if (layeredRangeTreePlus4NodeIntersectsRangeLevel2(n->children[2], -200, 200, y1, y2)) {
    layeredRangeTreePlus4RangeSearchLevel2(n->children[2], x1, x2, y1, y2, count);
  }
  if (layeredRangeTreePlus4NodeIntersectsRangeLevel2(n->children[3], -200, 200, y1, y2)) {
    layeredRangeTreePlus4RangeSearchLevel2(n->children[3], x1, x2, y1, y2, count);
  }

}

void layeredRangeTreePlus4RangeSearchLevel1(layeredRangeTreePlus4NodeLevel1 *n, double x1, double x2, double y1, double y2, int* count, int offset) {
  if (offset == -1){
    // do the initial binary search4.
    offset = arg_successor4_root(y1, &n->secondLevelPointer->cnode);
  }

  if (n->children == NULL) {
    assert(n->elementCount <= LAYERED_TREE_PLUS4_NODE_SIZE);
    // if we get to a leaf just scan through the elements and report the ones in the range.
    for (int i = n->start; i < n->end; i++) {
      if (layeredRangeTreePlus4KeyInRange(n->elements[i], x1, x2, y1, y2)){
        (*count)++;
        //printf("Element (%f, %f) in range \n", n->elements[i].x, n->elements[i].y);
      }
    }
    return;
  }

  assert(n->secondLevelPointer != NULL);

  // call the second level search4 if we're entirely within the x-range.
  if (n->x1 >= x1 && n->x2 <= x2) {
    for (;offset < n->secondLevelPointer->cnode.size && n->secondLevelPointer->cnode.ys[offset] <= y2; offset++ ){
            //printf("Element (%f, %f) in range \n", n->secondLevelPointer->cnode.xs[offset], n->secondLevelPointer->cnode.ys[offset]);
      (*count)++;
    }
    //layeredRangeTreePlus4RangeSearchLevel2(n->secondLevelPointer, -200, 200, y1, y2, count);
    return;
  }

  // otherwise call recursively on each child intersected by the range.
  if (layeredRangeTreePlus4NodeIntersectsRangeLevel1(n->children[0], x1, x2, -200, 200)) {
    layeredRangeTreePlus4RangeSearchLevel1(n->children[0], x1, x2, y1, y2, count, 
        arg_successor4_left(offset, &n->secondLevelPointer->cnode));
  }

  if (layeredRangeTreePlus4NodeIntersectsRangeLevel1(n->children[1], x1, x2, -200, 200)) {
    layeredRangeTreePlus4RangeSearchLevel1(n->children[1], x1, x2, y1, y2, count,
        arg_successor4_left2(offset, &n->secondLevelPointer->cnode));
  }
  if (layeredRangeTreePlus4NodeIntersectsRangeLevel1(n->children[2], x1, x2, -200, 200)) {
    layeredRangeTreePlus4RangeSearchLevel1(n->children[2], x1, x2, y1, y2, count,
        arg_successor4_right(offset, &n->secondLevelPointer->cnode));
  }
  if (layeredRangeTreePlus4NodeIntersectsRangeLevel1(n->children[3], x1, x2, -200, 200)) {
    layeredRangeTreePlus4RangeSearchLevel1(n->children[3], x1, x2, y1, y2, count,
        arg_successor4_right2(offset, &n->secondLevelPointer->cnode));
  }

}


// the root node
layeredRangeTreePlus4NodeLevel1* rangeTreePlus4Root;

// collection of elements
layeredRangeTreePlus4Key* rtplus4_allElements;
int rtplus4_allElementsCount;
int rtplus4_currentAllElementsArraySize;
double rtplus4_totalQueryTime;

void tfkLayeredRangeTreePlus42DRangeSearchCommand(redisClient* c) {
  double x1 = strtod(c->argv[1]->ptr, NULL);
  double x2 = strtod(c->argv[2]->ptr, NULL);
  double y1 = strtod(c->argv[3]->ptr, NULL);
  double y2 = strtod(c->argv[4]->ptr, NULL);
  int count = 0;
  double start_time = tfkLayeredRangeTreePlus4_get_time();
  layeredRangeTreePlus4RangeSearchLevel1(rangeTreePlus4Root, x1, x2, y1, y2, &count, -1);
  double end_time = tfkLayeredRangeTreePlus4_get_time();
  rtplus4_totalQueryTime += end_time - start_time;
  printf("total query time %f rangeTreePlus4Count: %d \n", rtplus4_totalQueryTime, count); 
  addReplyLongLong(c, count);
}

void tfkLayeredRangeTreePlus4BuildTreePlus4Command(redisClient *c) {
  if (rangeTreePlus4Root != NULL){
    // delete the old quad tree.
    layeredRangeTreePlus4DeleteNodeLevel1(rangeTreePlus4Root);
  }
  rangeTreePlus4Root = (layeredRangeTreePlus4NodeLevel1*) zmalloc(sizeof(layeredRangeTreePlus4NodeLevel1));

  rangeTreePlus4Root->x1 = -200;
  rangeTreePlus4Root->x2 = 200;
  rangeTreePlus4Root->y1 = -200;
  rangeTreePlus4Root->y2 = 200; 
  rangeTreePlus4Root->children = NULL;
  rangeTreePlus4Root->elementCount = rtplus4_allElementsCount;
  rangeTreePlus4Root->secondLevelPointer = NULL;
  // rebuild using the rtplus4_allElements array.

  // build the first level of the tree.
  buildLayeredRangeTreePlus4NodeLevel1(rangeTreePlus4Root, rtplus4_allElements, rtplus4_allElementsCount);

  // build the second level of the tree.
  layeredRangeTreePlus4Key* elementsCopy = (layeredRangeTreePlus4Key*) zmalloc(sizeof(layeredRangeTreePlus4Key) * rangeTreePlus4Root->elementCount);
  explodeTreePlus4(rangeTreePlus4Root, elementsCopy, rtplus4_allElementsCount);
  zfree(elementsCopy); 
  
  //walkLayeredRangeTreePlus4(rangeTreePlus4Root);
  rtplus4_totalQueryTime = 0;
  addReplyLongLong(c, 42);
}

/* This generic command implements both ZADD and ZINCRBY. */
void tfkLayeredRangeTreePlus4AddGenericCommand(redisClient *c, int incr) {
    double x = strtod(c->argv[1]->ptr, NULL);
    double y = strtod(c->argv[2]->ptr, NULL);
    
    int resizeAmount = 10000;
  
    // create a layeredRangeTreePlus4Key
    layeredRangeTreePlus4Key key;
    key.x = x;
    key.y = y;
    key.value = NULL; 
    if (rtplus4_allElementsCount == rtplus4_currentAllElementsArraySize) {
      // resize rtplus4_allElementsArray
      rtplus4_currentAllElementsArraySize += resizeAmount;
      rtplus4_allElements = (layeredRangeTreePlus4Key*) zrealloc((void*)rtplus4_allElements, sizeof(layeredRangeTreePlus4Key) * rtplus4_currentAllElementsArraySize);
    }
    rtplus4_allElements[rtplus4_allElementsCount] = key;
    rtplus4_allElementsCount++;
    addReplyLongLong(c, rtplus4_allElementsCount);
}

void layeredRangeTreePlus42DAddCommand(redisClient *c) {
  tfkLayeredRangeTreePlus4AddGenericCommand(c,0);
}

void layeredRangeTreePlus42DBuildTreeCommand(redisClient *c) {
  tfkLayeredRangeTreePlus4BuildTreePlus4Command(c);
}

void layeredRangeTreePlus42DRangeSearchCommand(redisClient *c) {
  tfkLayeredRangeTreePlus42DRangeSearchCommand(c);
}

