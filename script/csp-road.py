#!/bin/python

from subprocess import Popen

case_l = ["ny", "fla", "nw", "ne", "lks"]
log1 = open("result1.txt", "a")

for case in case_l:
    topo = "../../CSP/inputs/" + str(case) + ".graph"
    for i in range(1, 3):
        demand = "../../CSP/inputs/rquery_" + str(case) + "_s" + str(i)
        config = open("../inputs/config", 'w')
        print(topo)
        print(demand)
        input_file = topo
        demand_file = demand
        config.write("2\n")
        config.write(input_file + "\n")
        config.write(demand_file + "\n")
        config.write("1\n")
        config.close()

        process = Popen(["./vine"], stdout=log1)
        process.wait()
        log1.flush()
