import redis
import sys
from sys import argv
import csv


implementations =  ["kdtree", "layered", "cascading", "cascading4", "array"]
(KDTREE, LAYERED, CASCADING, CASCADING4, ARRAY) = implementations

if len(argv) < 5:
    print "usage: ", argv[0], "datafile", "queryfile", "port", "implementation"
    print "implementation can be any of: ", implementations
    sys.exit(1)

datafile = argv[1]
benchmarkfile = argv[2]
#key = argv[3] no longer needed
portnum = int(argv[3])
implementation = argv[4]

assert implementation in implementations

# FILTERED: 
# {
# 'add':'zadd2d',
# 'build':'zcard',
# 'range': 'zrangebyscore'
# }, 


#zset commands:
commands = {

#kd tree commands:
KDTREE: 
{'add':'add2dx',
 'build':'build2dx',
 'range':'range2dx'
},

LAYERED: 
{
'add':'add2dlayered',
'build':'build2dlayered',
'range':'range2dlayered'
},

CASCADING: 
{
'add':'add2dplus',
'build':'build2dplus',
'range':'range2dplus'
},
CASCADING4: 
{
'add':'add2dplus4',
'build':'build2dplus4',
'range':'range2dplus4'
},
ARRAY: 
{
'add':'add2darray',
'build':'build2darray',
'range':'range2dstupid'
}
}
# #after implementing, may want 'range2darray'

# writeCmd = ["add2d"] #x y
# rangeCmd = ["range2d"] #xmin xmax ymin ymax
# buildCmd = ["build2"]#build2d

# #skiplist and filter commands
# writecommand = ["zadd2d", key] #key x y
# readcommand = ["zrangebyscore", key] #key xmin xmax ymin ymax

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
        insert2dcommand = [commands[implementation]['add']] + [lat, lon, str(counter)]
        counter += 1
        if counter == 1:
            print "first command:", " ".join(insert2dcommand)
        client.execute_command(*insert2dcommand)

    except Exception as e:
        print "Exception. insert command", " ".join(insert2dcommand), "failed:", e
        sys.exit(1)

print "done with inserts"
print "building..."
client.execute_command(*[commands[implementation]['build']])
print "done building..."

print "starting read benchmark"
benchreader = csv.reader(open(benchmarkfile, 'r'),  delimiter='|')

numquery = 0
for line in benchreader:
    (minlat, maxlat) = line[1].split(",") 
    (minlon, maxlon) = line[2].split(",")
    try:

        query = [commands[implementation]['range']] + [minlat, maxlat, minlon, maxlon]
        answer = client.execute_command(*query)
        print answer
        numquery += 1
        # if numquery % 100 == 0:
        #     print "queries done: ", numquery
    except Exception as e:
        print "Exception. call: ", query, "failed:", e
        sys.exit(1)

print "benchmark done."
