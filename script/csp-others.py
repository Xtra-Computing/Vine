#!/bin/python

from subprocess import Popen

case_l = ["tech-internet-as.mtx.txt", "tech-p2p-Gnutella09.mtx.txt"]
case_l += ["soc-gplus.mtx.txt", "soc-slashdot.mtx.txt"]

log3 = open("cpu.txt", "a")

for case in case_l:
    topo = "../../CSP/inputs/" + str(case)
    demand = "../../CSP/inputs/demand.csv"
    config = open("./inputs/config", 'w')
    input_file = topo
    demand_file = demand
    config.write("1\n")
    config.write(input_file + "\n")
    config.write(demand_file + "\n")
    config.write("1\n")
    config.close()

    process = Popen(["./vine"], stdout=log3)
    process.wait()
    log3.flush()
