#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

int min(int a, int b){
  return (a < b)? a : b;
}

typedef struct cascading_node {
  int size; // array sizes
  double *ys;// eg ys
  double *xs; //eg x coords
  int *left_successor; //positions for cascading left (you keep track
  int *right_successor; // positions for cascading right

  struct cascading_node *left;
  struct cascading_node *right;
} cascading_node_t;

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
    assert(c->ys[i] < c->ys[i+1]);
    // assert(c->ys[i] < DBL_MAX);
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
  cn->ys = (double*) malloc((size+1)*sizeof(double));
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
  cn->ys = (double *)malloc((cn->size+1)*sizeof(double));
  cn->ys[cn->size] = DBL_MAX;

  //does not need any sentinel value, simply the other coordinate
  cn->xs = (double *)malloc(cn->size*sizeof(double));

  //allocate an array of size +1 length, the -1 entry gets inited to 0
  cn->left_successor = (int *) malloc((cn->size+1)*sizeof(int)) + 1;
  *(cn->left_successor - 1) = -1;
  cn->right_successor = (int *) malloc((cn->size+1)*sizeof(int)) + 1;
  *(cn->right_successor - 1) = -1;

  cn->left = left;
  cn->right = right;


  int i, j;
  for (i = 0, j = 0;  i < left->size || j < right->size; ){
    if (left->ys[i] < right->ys[j]){ //left has next smallest elt.
      cn->ys[i+j] = left->ys[i];
      cn->xs[i+j] = left->xs[i];

      cn->left_successor[i+j] = i;
      cn->right_successor[i+j] = cn->right_successor[i+j-1] + 1;
      i++; 
    } else { //right has next smallest
      cn->ys[i+j] = right->ys[j];
      cn->xs[i+j] = right->xs[j];

      cn->left_successor[i+j] = cn->left_successor[i+j-1] + 1;//copy previous
      cn->right_successor[i+j] = j;
      j++;
    }
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

/*the value at *end is meant to be a sentinel MAX_INT
 *gives you the succesor position (smallest elt greater than or equal to value)
 */
double* binary_successor_search(double value, double* start, double* end){

  while(start < end){
    double * middle =  start + (end-start)/2;
    if ( value < *middle ) {
      end = middle;
    } else if ( value > *middle ){
      start = middle + 1;
    } else {
      return middle;
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



void test1(){
 double avalx[] = {1, 3, 5};
 double avaly[] = {100.0, 300.0, 500.0};
 cascading_node_t anode;
 cascading_node_leaf_init(&anode, avalx, avaly, 3);

 double bvalx[] = {2, 4};
 double bvaly[] = {200, 400};
 cascading_node_t bnode;
 cascading_node_leaf_init(&bnode, bvalx, bvaly, 2);

 cascading_node_t cnode;
 cascading_node_parent_init(&cnode, &anode, &bnode);

 /* printf("anode:\n"); */
 /* cascading_node_print(&anode); */
 /* printf("bnode:\n"); */
 /* cascading_node_print(&bnode); */
 /* printf("cnode:\n"); */
 /* cascading_node_print(&cnode); */
 /* printf("done\n"); */

 double expected_ys[] = {100, 200, 300, 400, 500};
 double expected_xs[] = {1, 2, 3, 4, 5};
 int expected_left[] = {0, 1, 1, 2,  2};
 int expected_right[] = {0, 0, 1, 1, 2};
 int expected_size = 5;
 assert(cnode.size == expected_size);

 int i = 0;
 for (i = 0; i < cnode.size; i++){
   assert(cnode.ys[i] == expected_ys[i]);
   assert(cnode.xs[i] == expected_xs[i]);
   assert(cnode.left_successor[i] == expected_left[i]);
   assert(cnode.right_successor[i] == expected_right[i]);
 }

 /*write test case later*/
 /* int ays[] = {0, 1, 2, 3, 4}; */
 /* int axs[] = {10, 11, 12, 13, 14}; */
 /* int anode2 = cascading_node */
 
 /* int bys[] = {5, 6, 7}; */
 /* int bxs[] = {15, 16, 17}; */

 int index;
 assert(arg_successor_root(100, &cnode) == 0);
 assert(arg_successor_root(150, &cnode) == 1);
 index = arg_successor_root(150, &cnode);
 assert(arg_successor_left(index, &cnode) == arg_successor_root(150, &anode));
 assert(arg_successor_right(index, &cnode) == arg_successor_root(150, &bnode));

}


void testBinarySuccessorSearch(){
  double a1[] = {0.0, 2.0, 4.0, 6.0, DBL_MAX};
  assert(*binary_successor_search(-10, a1, a1+4) == 0);
  assert(*binary_successor_search(0, a1, a1+4) == 0);
  assert(*binary_successor_search(1, a1, a1+4) == 2);
  assert(*binary_successor_search(2, a1, a1+4) == 2);
  assert(*binary_successor_search(3, a1, a1+4) == 4);
  assert(*binary_successor_search(4, a1, a1+4) == 4);
  assert(*binary_successor_search(10, a1, a1+4) == DBL_MAX);
  assert(*binary_successor_search(100, a1, a1+4) == DBL_MAX);
}

int main(){
  testBinarySuccessorSearch();
  test1();
}
