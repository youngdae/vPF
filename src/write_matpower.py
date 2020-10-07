#!/usr/bin/env python

from gams import *
import os

def write_matpower(data, case, gencost_dummy=False, cs=False):
    baseMVA = float(data["baseMVA"])
    f = open(case + "_merged.m", "w")
    if case in ["Network_06O-124_scenario_17", "Network_09O-195_scenario_27", \
                "Network_12O-050_scenario_2", "Network_13O-103_scenario_15",  \
                "Network_20O-100_scenario_13", "Network_25R-060_scenario_10", \
                "Network_25O-060_scenario_10", "Network_30O-048_scenario_10"]:
        busdict = {v[0]:i for i,v in enumerate(data["bus"])}
        bus_ref = [v[0] for v in data["bus"] if v[1] == 3]
        ws = GamsWorkspace()
        db = ws.add_database_from_gdx(os.getcwd() + "/" + case + "_start.gdx")
        for rec in db["PG"]:
            data["gen_aggr"][int(rec.keys[0])-1][1] = rec.level*baseMVA
        for rec in db["CS"]:
            data["bus"][busdict[rec.keys[0]]][13] = rec.level
        for rec in db["VM"]:
            if rec.keys[0] in bus_ref:
                data["bus"][busdict[rec.keys[0]]][7] = rec.level
    # Identify buses having controllable generators.
    pv_dict = {}
    for g in data["gen_aggr"]:
        if float(g[3]) != float(g[4]) and not (g[0] in pv_dict):
            pv_dict[g[0]] = 1
    f.write("function mpc = " + case + "_merged\n\n")
    f.write("mpc.version = '2';\n\n")
    f.write("mpc.baseMVA = " + data["baseMVA"] + ";\n\n")
    f.write("mpc.bus = [\n")
    for b in data["bus"]:
        bcopy = b[:]
        if int(bcopy[1]) == 2 and not (bcopy[0] in pv_dict):
            bcopy[1] = 1
        if cs:
            bcopy[5] += (bcopy[13]*baseMVA)
        f.write('\t' + '\t'.join(str(e) for e in bcopy[:13]) + ";\n")
    f.write("];\n\n")
    f.write("mpc.gen = [\n")
    for g in data["gen_aggr"]:
        f.write('\t' + '\t'.join(str(e) for e in g[:21]) + ";\n")
    f.write("];\n\n")
    f.write("mpc.branch = [\n")
    for l in data["branch"]:
        f.write('\t' + '\t'.join(str(e) for e in l[:13]) + ";\n")
    f.write("];\n\n")

    f.write("mpc.gencost =[\n")
    if gencost_dummy:
        for g in data["gen_aggr"]:
            f.write('\t' + "2\t0\t0\t3\t0\t1\t0;\n")
    else:
        for c in data["gencost_aggr"]:
            f.write('\t' + '\t'.join(str(e) for e in c[:7]) + ";\n")
    f.write("];\n")
    f.close()


