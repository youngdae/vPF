#/usr/bin/env python

from write_gdx import *
from write_matpower import *
import sys
import os

def usage(progname):
    print("Usage: " + progname + " matpower_filename")

def next_key(lines, pos):
    key = ""
    while pos < len(lines):
        terms = lines[pos].split()
        if terms and terms[0] in ["mpc.baseMVA", "mpc.bus", "mpc.gen", "mpc.branch", "mpc.gencost"]:
            key = terms[0].split(".")[1]
            break
        else:
            pos += 1

    return key, pos

def read_raw_data(filename):
    # Read data from the MATPOWER case file.
    f = open(filename, "r")
    lines = f.readlines()
    p = 0
    data = {}
    while p < len(lines):
        key, p = next_key(lines, p)
        p += 1

        if key == "":
            break
        elif key == "baseMVA":
            data["baseMVA"] = lines[p-1].split("=")[1].strip().replace(";","")
        else:
            q = p
            while not lines[q].startswith("];"):
                q += 1
            data[key] = [[x.strip().replace(";","") for x in lines[j].split()] for j in range(p,q)]
    f.close()
    merge_duplicate_generators(data)

    # Add gij_m and bij_m values. Matpower does not have them so we just attach 0.
    for i,v in enumerate(data["branch"]):
        data["branch"][i] = v[0:13] + [0,0]

    for i,v in enumerate(data["bus"]):
        data["bus"][i] = v[0:13] + [0,0,0]

    return data

def main(argv):
    data = read_raw_data(argv[1])
    case = os.path.splitext(os.path.basename(argv[1]))[0]
    print(" ** Statistics of", case)
    print("  # buses     : ", len(data["bus"]))
    print("  # generators: ", len(data["gen"]))
    print("  # branches  : ", len(data["branch"]))
    print("  # gencost   : ", len(data["gencost"]))

    write_gdx(data, case)
    write_matpower(data, case)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        usage(sys.argv[0])
        exit(-1)
    main(sys.argv)
