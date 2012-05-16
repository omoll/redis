1) need to use binary search on the root so that its efficient (right now linear search) 
   right now my binary search is flawed. 
binary search, and also, this binary search is meant to find the successor, so even if a number is
not int he array, we cannot just return -1 or something, it needs to be the next thing.

plus, it needs to be 'duplicate resistant' binary search, 
2) can optimize by stopping the  first tree traversal as soon as the y tree is empty already (right now
we wait until we have gone down all the way in the X-tree, and then check the y, but it may be the y index already tells use there is no point. we can use the offset to figure out if there is no point)

--------
bugs found:
     -binary search itself (still buggy), replaced with linear search
     -dealing with duplicates during binary search
     -dealing with duplicates at merge time
     -not stopping when y was too big, instead stopping until the end of array
     -not dealing with the case where the offset for y is already the max it can be
(ie no soluton) and feeding that in to the child offset calculation again.     
