import csv
from random import uniform
from sys import argv
import sys

(DIST2D, DISTDIAG) = disttypes = ["dist2d", "distdiag"]
(QUERIES, POINTS) = tasktypes = ["queries", "points"]

LIMIT = 50
RECSIZE = 20

FMT = "%1.6f"

if len(argv) == 1:
    print "usage:", "queries|points" "num", "distType:", disttypes
    sys.exit(1)

mode = argv[1]
argv = argv[1:]

if mode ==  POINTS:
    numpoints = int(argv[1])
    assert argv[2] in disttypes, "invalid distribution"

    writer = csv.writer(sys.stdout, delimiter=",")
    for i in range(numpoints):
        sample = (i, FMT % uniform(-LIMIT, LIMIT), FMT % uniform(-LIMIT, LIMIT))
        writer.writerow(sample)
    
    sys.stdout.flush()
    exit(0)

elif mode == QUERIES:
    numqueries = int(argv[1])
    writer = csv.writer(sys.stdout, delimiter="|")

    for i in range(numqueries):
        (minx,miny) = (uniform(-LIMIT, LIMIT-RECSIZE), uniform(-LIMIT, LIMIT-RECSIZE))
        maxy = miny + RECSIZE
        maxx = minx + RECSIZE

        query = (",".join([str(i),str(i+1)]), ",".join([FMT % minx, FMT % maxx]), \
                     ",".join([FMT % miny, FMT % maxy]))
        
        writer.writerow(query)

    sys.stdout.flush()
    exit(0)
