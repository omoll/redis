import redis
import sys
from sys import argv
import csv


if len(argv) == 1:
    print "usage: ", argv[0], "datafile", "queryfile", "key", "port"
    sys.exit(1)

datafile = argv[1]
benchmarkfile = argv[2]

key = argv[3]
portnum = int(argv[4])

writecommand = "zadd2d"
readcommand = "zrangebyscore" # need to change to new zrange command

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
        insert2dcommand = [writecommand, key, lat, lon, str(counter)]
        counter += 1
        
       # need to add lon for 2D range queries
        if counter == 1:
            print " ".join(insert2dcommand)
            
        client.execute_command(*insert2dcommand)

    except Exception as e:
        print "Exception. insert command", "".join(insert), "failed:", e
        sys.exit(1)

print "done with inserts"
sys.exit(0)

print "starting read benchmark"

benchreader = csv.reader(open(benchmarkfile, 'r'),  delimiter='|')
for line in benchreader:
    (minlat, maxlat) = line[1].split(",") 
    (minlon, maxlon) = line[2].split(",")
    try:
        # redis syntax: zrangebyscore key minscore maxscore. 
        query = prefix + readcommand + [minlat, maxlat]
        # will need to add minlon maxlon for 2d range queries
        print query
        check_call(query)
    except Exception as e:
        print "Exception. call: ", query, "failed:", e
        sys.exit(1)

print "benchmark done."
