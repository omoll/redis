import csv
from random import uniform
from sys import argv
from math import sqrt
import sys

(DIST2D, DISTDIAG) = disttypes = ["dist2d", "distdiag"]
(QUERIES, POINTS) = tasktypes = ["queries", "points"]

LIMIT = 50
RECSIZE = 10

FMT = "%1.6f"

if len(argv) == 1:
    print "usage:", "queries|points" "num", "distType", "ratio"
    print "distType can be: ", disttypes
    sys.exit(1)

mode = argv[1]
assert mode in tasktypes, "invalid mode"
argv = argv[1:]

dist = argv[2]

if mode == QUERIES and dist == DIST2D:
    ratio = float(argv[3])
    horizontal = sqrt(float(RECSIZE**2)/ratio)
    vertical = ratio*horizontal
    assert (horizontal*vertical - RECSIZE**2)/float(RECSIZE) < 0.01
    assert (vertical/float(horizontal))/ratio  - 1 < 0.01


if mode ==  POINTS:
    numpoints = int(argv[1])

    assert dist in disttypes, "invalid distribution"
    writer = csv.writer(sys.stdout, delimiter=",")

    if dist == DIST2D: #points in a 2D plane
        for i in range(numpoints):
            sample = (i, FMT % uniform(-LIMIT, LIMIT), FMT % uniform(-LIMIT, LIMIT))
            writer.writerow(sample)

    elif dist == DISTDIAG: #only points on diagonal
        assert False, "check code before running this"
        for i in range(numpoints):
            coord = FMT % uniform(-LIMIT, LIMIT)
            sample = (i, coord, coord)
            writer.writerow(sample)

    sys.stdout.flush()
    exit(0)

elif mode == QUERIES:
    numqueries = int(argv[1])
    writer = csv.writer(sys.stdout, delimiter="|")

    if dist == DIST2D: #query anywhere, fully within the range
        for i in range(numqueries):
            (minx,miny) = (uniform(-LIMIT, LIMIT-horizontal), uniform(-LIMIT, LIMIT-vertical))
            maxy = miny + vertical
            maxx = minx + horizontal

            query = (",".join([str(i),str(i+1)]), ",".join([FMT % minx, FMT % maxx]), \
                     ",".join([FMT % miny, FMT % maxy]))
        
            writer.writerow(query)
    elif dist == DISTDIAG: #query squares only on diag
        assert False, "check code before running this"
        for i in range(numqueries):
            minx = miny = uniform(-LIMIT, LIMIT-RECSIZE)
            maxy = maxx = miny + RECSIZE

            query = (",".join([str(i),str(i+1)]), ",".join([FMT % minx, FMT % maxx]), \
                     ",".join([FMT % miny, FMT % maxy]))
        
            writer.writerow(query)

    sys.stdout.flush()
    exit(0)
