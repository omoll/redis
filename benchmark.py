import redis
import sys
from sys import argv
import csv


if len(argv) == 1:
    print "usage: ", argv[0], "datafile", "queryfile", "key", "port", "mode"
    sys.exit(1)

datafile = argv[1]
benchmarkfile = argv[2]

key = argv[3]
portnum = int(argv[4])
implementation = argv[5]

implementations =  ["filtered", "kdtree"]
(FILTERED, KDTREE) = implementations
assert implementation in implementations

#kd tree commands:
command1 = ["add2d"] #x y
command2 = ["range2d"] #xmin xmax ymin ymax
command3 = ["build2d"]

#skiplist and filter commands
writecommand = ["zadd2d", key] #key x y
readcommand = ["zrangebyscore", key] #key xmin xmax ymin ymax

client = redis.StrictRedis(host='localhost', port=portnum, db=0)

"""
The files staring with "i" are inserts of format: timestamp, lat, lon
The files staring with "q" are queries of format: timestart, timeend | minlat, maxlat | minlon, maxlon
"""
#r.execute_command(*["ZADD2D", "test1", 2.0, 3.0, 'twothree'])
print "starting inserts"
insertreader = csv.reader(open(datafile, 'r'), delimiter = ',')
counter = 0
for line in insertreader:
    (timestamp, lat, lon) = line
    try:
        # redis syntax: zadd key score member. it seems it won't accept repeated member values.

        if implementation == FILTERED:
            insert2dcommand = writecommand + [lat, lon, str(counter)]
        elif implementation == KDTREE:
            insert2dcommand = command1 + [lat, lon]

        counter += 1
       # need to add lon for 2D range queries
        if counter == 1:
            print "first command:", " ".join(insert2dcommand)
            
        client.execute_command(*insert2dcommand)

    except Exception as e:
        print "Exception. insert command", " ".join(insert2dCommand), "failed:", e
        sys.exit(1)

print "done with inserts"

if implementation == KDTREE:
    print "building..."
    client.execute_command(*["build2d"])
    print "done building..."

print "starting read benchmark"

benchreader = csv.reader(open(benchmarkfile, 'r'),  delimiter='|')

numquery = 0
for line in benchreader:
    (minlat, maxlat) = line[1].split(",") 
    (minlon, maxlon) = line[2].split(",")
    try:

        if implementation == FILTERED:
            query = readcommand + [minlat, maxlat, minlon, maxlon]
        elif implementation == KDTREE:
            query = command2 + [minlat, maxlat, minlon, maxlon]

        answer = client.execute_command(*query)
        print answer
        numquery += 1
        # if numquery % 100 == 0:
        #     print "queries done: ", numquery
    except Exception as e:
        print "Exception. call: ", query, "failed:", e
        sys.exit(1)

print "benchmark done."
