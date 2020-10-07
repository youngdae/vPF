#!/usr/bin/env python

from write_gdx import *
from write_matpower import *
import sys
import os

def usage(progname):
    print("Usage: " + progname + " godir")
    print("  godir: directory containing raw, rop, inl and con files")

def get_next_section_start(lines, q, inc):
    while not lines[q].startswith("0"):
        q += inc
    return q

def read_raw_data(filename):
    f = open(filename, "r")
    lines = f.readlines()
    data = {}
    p = 0

    # Read baseMVA.
    data["baseMVA"] = lines[0].split(',')[1].strip()

    p = 3
    keys = ["raw_bus", "raw_load", "raw_fshunt", "raw_gen", "raw_branch"]
    for k in keys:
        q = get_next_section_start(lines, p, 1)
        data[k] = [[x.strip() for x in lines[j].split(",")] for j in range(p,q)]
        p = q + 1

    # Each transfomer data consists of 4 lines.
    q = get_next_section_start(lines, p, 4)
    data["raw_trans"] = [[x.strip() for x in ','.join(lines[j:j+4]).split(",")] for j in range(p,q,4)]
    p = q + 1

    # Skip sections.
    for i in range(0,10):
        q = get_next_section_start(lines, p, 1)
        p = q + 1

    q = get_next_section_start(lines, p, 1)
    data["raw_sshunt"] = [[x.strip() for x in lines[j].split(",")] for j in range(p,q)]
    p = q + 1
    f.close()
    return data

def arrange_bus(raw):
    baseMVA = float(raw["baseMVA"])
    busdict = {v[0]:i for i,v in enumerate(raw["raw_bus"])}
    # By default, bus type is set to PQ-bus.
    raw["bus"] = [[v[0],1,0,0,0,0,v[4],v[7],v[8],v[2],v[5],v[9],v[10],0,0,0] for i,v in enumerate(raw["raw_bus"])]

    # Add load data, Pd and Qd.
    for e in raw["raw_load"]:
        if int(e[2]) == 1:    # if status == 1
            raw["bus"][busdict[e[0]]][2] += float(e[5])
            raw["bus"][busdict[e[0]]][3] += float(e[6])

    # Add fixed shunt data, Gs and Bs.
    for e in raw["raw_fshunt"]:
        if int(e[2]) == 1:    # if status == 1
            raw["bus"][busdict[e[0]]][4] += (float(e[3]) / baseMVA)
            raw["bus"][busdict[e[0]]][5] += (float(e[4]) / baseMVA)

    # Add switched shunt data.
    for e in raw["raw_sshunt"]:
        if int(e[3]) == 1:    # if status == 1
            binit = float(e[9])/baseMVA
            bl = [(float(e[i])*float(e[i+1]))/baseMVA for i in range(10,26,2)]
            bcs_up = 0
            bcs_lo = 0
            nbl = -1
            for i in range(0,8):
                if bl[i] == 0:
                    break
                nbl += 1
            if nbl != -1:
                bcs_up = sum(max(0, bl[l]) for l in range(0,nbl+1))
                bcs_lo = sum(min(0, bl[l]) for l in range(0,nbl+1))
            raw["bus"][busdict[e[0]]][-3] = binit
            raw["bus"][busdict[e[0]]][-2] = bcs_up
            raw["bus"][busdict[e[0]]][-1] = bcs_lo

    # Set type of bus.
    pmax = -1
    pmax_idx = -1
    for e in raw["raw_gen"]:
        if int(e[14]) != 0:
            raw["bus"][busdict[e[0]]][1] = 2   # Set bus type to PV-bus.
            if float(e[16]) > pmax:
                pmax = float(e[16])
                pmax_idx = e[0]
    raw["bus"][busdict[pmax_idx]][1] = 3   # Set the slack bus.

def arrange_gen(raw):
    raw["gen"] = [[v[0],v[2],v[3],v[4],v[5],0,0,v[14],v[16],v[17],0,0,0,0,0,0,0,0,0,0,0] for i,v in enumerate(raw["raw_gen"])]

def arrange_branch(raw):
    baseMVA = float(raw["baseMVA"])
    branch = [[v[0],v[1],v[3],v[4],v[5],v[6],v[7],v[8],0,0,v[13],-360,360,0,0] for i,v in enumerate(raw["raw_branch"])]
    trans = [[v[0],v[1],v[21],v[22],0,v[27],v[28],v[29],float(v[24])/float(v[41]),v[26],v[11],-360,360,v[7],v[8]] for i,v in enumerate(raw["raw_trans"])]
    raw["branch"] = branch + trans

def main(argv):
    raw = read_raw_data(argv[1] + "/case.raw")
    realpath = os.path.realpath(argv[1]).split('/')
    case = realpath[-2] + "_" + realpath[-1]

    arrange_bus(raw)
    arrange_gen(raw)
    arrange_branch(raw)
    merge_duplicate_generators(raw, False)

    print(" ** Statistics of case.raw")
    print("    baseMVA         : {0:8.2f}".format(float(raw["baseMVA"])))
    print("  # buses           : {0:8d}".format(len(raw["raw_bus"])))
    print("  # fixed shunts    : {0:8d}".format(len(raw["raw_fshunt"])))
    print("  # generators      : {0:8d}".format(len(raw["raw_gen"])))
    print("  # branches        : {0:8d}".format(len(raw["raw_branch"])))
    print("  # transformers    : {0:8d}".format(len(raw["raw_trans"])))
    print("  # switched shunts : {0:8d}".format(len(raw["raw_sshunt"])))

    # write_gdx() should be called first since write_matpower() may change raw.
    write_gdx(raw, case)
    write_matpower(raw, case, True, True)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        usage(sys.argv[0])
        sys.exit(-1)
    main(sys.argv)
