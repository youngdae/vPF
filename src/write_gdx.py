#!/usr/bin/env python

from gams import *
import sys
import os
import numpy as np

GENSTR = "gen_aggr"
GENCOST_STR = "gencost_aggr"

def merge_duplicate_generators(data, gencost=True):
    # Merge generators that share the same bus by aggregating their min/max.
    gen_bus = {}
    gen_bus_count = {}
    ngen = len(data["gen"])
    data[GENSTR] = []
    if gencost:
        data[GENCOST_STR] = []
    j = 0
    for i in range(ngen):
        if int(data["gen"][i][7]) != 0:
            bus_i = int(data["gen"][i][0])
            if bus_i in gen_bus:
                # Check if the signs on min/max are consistent.
                # Otherwise, they cannot be aggregated.
                idx = gen_bus[bus_i]
#                sign_diff = sum([1 for x in [3, 4, 8, 9] if float(data[GENSTR][idx][x])*float(data["gen"][i][x]) < 0])
#                if sign_diff > 0:
#                    print("Error: cannot aggregate generator " + str(i+1) + " due to sign conflicts.")
#                    sys.exit(-1)
                for x in [1, 2, 3, 4, 8, 9]:
                    data[GENSTR][idx][x] = float(data[GENSTR][idx][x]) + float(data["gen"][i][x])
                gen_bus_count[bus_i] += 1
            else:
                gen_bus[bus_i] = j
                gen_bus_count[bus_i] = 1
                data[GENSTR].append(data["gen"][i])
                if gencost:
                    data[GENCOST_STR].append(data["gencost"][i])
                j += 1
    ndup = sum([1 for k,v in gen_bus_count.items() if v > 1])
    if ndup > 0:
        print("Warning: there are " + str(ndup) + " buses having more than 1 generator.")


def add_bus(db, data):
    b = db.add_set("b", 1, "bus index")
    nbus = len(data["bus"])
    busref = []
    for i in range(nbus):
        b.add_record(data["bus"][i][0])
        if int(data["bus"][i][1]) == 3:
            busref.append(data["bus"][i][0])
    if not busref:
        print("Error: reference bus was not found.")
        sys.exit(-1)
    b_ref = db.add_set_dc("bus_ref", [db["b"]], "reference bus index")
    data["busref"] = []
    for r in busref:
        b_ref.add_record(r)
        data["busref"].append(int(r))

    bf_name = [ "bus_i", "type", "Pd", "Qd", "Gs", "Bs", "area",
                "Vm", "Va", "baseKV", "zone", "Vmax", "Vmin", "Cs", "Csmax", "Csmin" ]
    bf = db.add_set("bf", 1, "bus field name")
    for i in bf_name:
        bf.add_record(i)
    bus = db.add_parameter_dc("bus", [b, bf], "bus data")
    b_cs = db.add_set_dc("bus_cs", [db["b"]], "bus switched shunt index")
    for i in range(nbus):
        for k, v in enumerate(bf_name):
            if v == "Va":
                bus.add_record((data["bus"][i][0], v)).value = float(data["bus"][i][k]) * (np.pi/180)
            else:
                bus.add_record((data["bus"][i][0], v)).value = float(data["bus"][i][k])
        if float(data["bus"][i][-1]) != 0 or float(data["bus"][i][-2]) != 0:
            b_cs.add_record(data["bus"][i][0])

def add_gen(db, data):
    g = db.add_set("g", 1, "generator index")
    ngen = len(data[GENSTR])
    j = 0
    for i in range(ngen):
        if int(data[GENSTR][i][7]) != 0:
            j += 1
            g.add_record(str(j))

    gf_name = [ "bus", "Pg", "Qg", "Qmax", "Qmin", "Vg", "mBase",
                "status", "Pmax", "Pmin", "Pc1", "Pc2", "Qc1min",
                "Qc1max", "Qc2min", "Qc2max", "ramp_agc", "ramp_10",
                "ramp_30", "ramp_q", "apf" ]
    gf = db.add_set("gf", 1, "generator field name")
    for i in gf_name:
        gf.add_record(i)
    gen = db.add_parameter_dc("gen", [g, gf], "generator data")
    j = 0
    for i in range(ngen):
        if int(data[GENSTR][i][7]) != 0:
            j += 1
            for k, v in enumerate(gf_name):
                if v in ["Pg", "Qg", "Qmax", "Qmin", "Pmax", "Pmin"]:
                    gen.add_record((str(j), v)).value = float(data[GENSTR][i][k]) / float(data["baseMVA"])
                else:
                    gen.add_record((str(j), v)).value = float(data[GENSTR][i][k])

def add_branch(db, data):
    l = db.add_set("l", 1, "branch index")
    nbr = len(data["branch"])
    j = 0
    for i in range(nbr):
        if int(data["branch"][i][10]) != 0:
            j += 1
            l.add_record(str(j))
    lf_name = [ "fbus", "tbus", "r", "x", "b", "rateA", "rateB",
                "rateC", "ratio", "angle", "status", "angmin", "angmax", "mag1", "mag2" ]
    lf = db.add_set("lf", 1, "branch field name")
    for i in lf_name:
        lf.add_record(i)
    linelim = db.add_set_dc("linelim", [db["l"]], "line indices having line limit")
    br = db.add_parameter_dc("br", [l, lf], "branch data")
    j = 0
    for i in range(nbr):
        if int(data["branch"][i][10]) != 0:
            j += 1
            for k, v in enumerate(lf_name):
                br.add_record((str(j), v)).value = float(data["branch"][i][k])
            rateA = float(data["branch"][i][5])
            if rateA > 0 and rateA < 1e10:
                linelim.add_record(str(j))

def add_join_bus_br(db, data):
    tij_reg = db.add_set_dc("tij_reg", [db["l"]], "is tap a regulator")
    tij_reg_bus = db.add_set_dc("tij_reg_bus", [db["l"], db["b"]],
                              "a pair of line and bus where tij controls the bus")
    fr_bus = db.add_set_dc("fr_bus", [db["l"], db["b"]], "line indices from bus")
    to_bus = db.add_set_dc("to_bus", [db["l"], db["b"]], "line indices from bus")
    reg_inserted = {}
    nbr = len(data["branch"])
    j = 0
    for i in range(nbr):
        if int(data["branch"][i][10]) != 0:
            j += 1
            fr = data["branch"][i][0]
            to = data["branch"][i][1]
            fr_bus.add_record((str(j), fr))
            to_bus.add_record((str(j), to))
            # At most one bus per line is added, and at most one line
            # can control each bus.
            if float(data["branch"][i][8]) != 0.0:
                if fr not in reg_inserted:
                    tij_reg.add_record(str(j))
                    tij_reg_bus.add_record((str(j), fr))
                    reg_inserted[fr] = 1
                elif to not in reg_inserted:
                    tij_reg.add_record(str(j))
                    tij_reg_bus.add_record((str(j), to))
                    reg_inserted[to] = 1

def add_join_bus_gen(db, data):
    bus_reg = db.add_set_dc("bus_reg", [db["b"]], "is a regulated bus")
    gen_bus = db.add_set_dc("gen_bus", [db["g"], db["b"]], "generator indices attached to bus")
    gen_ref = db.add_set_dc("gen_ref", [db["g"]], "generator indices attached to a slack bus")
    ngen = len(data[GENSTR])
    pv_dict = {}
    j = 0

    for i in range(ngen):
        if int(data[GENSTR][i][7]) != 0:
            j += 1
            gen_bus.add_record((str(j), data[GENSTR][i][0]))
            if int(data[GENSTR][i][0]) in data["busref"]:
                gen_ref.add_record(str(j))
            else:
                # Bus is regulated only if we could use reactive power,
                # i.e., Qmin != Qmax.
                if float(data[GENSTR][i][3]) != float(data[GENSTR][i][4]):
                    if not (data[GENSTR][i][0] in pv_dict):
                        pv_dict[data[GENSTR][i][0]] = 1
                        bus_reg.add_record(data[GENSTR][i][0])
                    else:
                        pv_dict[data[GENSTR][i][0]] += 1

    nbus_gen = sum([1 for k,v in pv_dict.items() if v > 1])
    if nbus_gen > 0:
        print("Error: there are " + str(nbus_gen) + " regulated buses having more than 1 generator.")
        print("       execute merge routine() to merge them into one generator.")
        sys.exit(-1)

def add_admittance_shunt(db, data):
    gij = db.add_parameter_dc("gij", [db["l"]], "series conductance of branch l")
    bij = db.add_parameter_dc("bij", [db["l"]], "series susceptance of branch l")
    bsh_ij = db.add_parameter_dc("bsh_ij", [db["l"]], "shunt susceptance of branch l")
    tij = db.add_parameter_dc("tij", [db["l"]], "ratio of transformer l")
    phi_ij = db.add_parameter_dc("phi_ij", [db["l"]], "angle of transformer l")
    gij_m = db.add_parameter_dc("gij_m", [db["l"]], "magnetizing conductance of branch l")
    bij_m = db.add_parameter_dc("bij_m", [db["l"]], "magnetizing susceptance of branch l")
    br = db["br"]
    for rec in db["l"]:
        r = br.find_record((rec.keys[0], "r")).value
        x = br.find_record((rec.keys[0], "x")).value
        b = br.find_record((rec.keys[0], "b")).value
        ratio = br.find_record((rec.keys[0], "ratio")).value
        angle = br.find_record((rec.keys[0], "angle")).value
        mag1 = br.find_record((rec.keys[0], "mag1")).value
        mag2 = br.find_record((rec.keys[0], "mag2")).value

        Ys = 1 / (r + x*1j)
        tap = 1.0
        phi = 0.0
        if ratio != 0.0:
            tap = ratio
        if angle != 0.0:
            phi = angle * (np.pi/180)

        gij.add_record(rec.keys[0]).value = Ys.real
        bij.add_record(rec.keys[0]).value = Ys.imag
        bsh_ij.add_record(rec.keys[0]).value = b
        tij.add_record(rec.keys[0]).value = tap
        phi_ij.add_record(rec.keys[0]).value = phi
        gij_m.add_record(rec.keys[0]).value = mag1
        bij_m.add_record(rec.keys[0]).value = mag2

    gs = db.add_parameter_dc("gs", [db["b"]], "shunt conductance at bus i")
    bs = db.add_parameter_dc("bs", [db["b"]], "shunt susceptance at bus i")
    bus = db["bus"]
    for rec in db["b"]:
        Gs = bus.find_record((rec.keys[0], "Gs")).value
        Bs = bus.find_record((rec.keys[0], "Bs")).value
        gs.add_record(rec.keys[0]).value = Gs / float(data["baseMVA"])
        bs.add_record(rec.keys[0]).value = Bs / float(data["baseMVA"])

def add_ybus(db, data):
    YffR = db.add_parameter_dc("YffR", [db["l"]], "real value of Yff")
    YffI = db.add_parameter_dc("YffI", [db["l"]], "imaginary value of Yff")
    YttR = db.add_parameter_dc("YttR", [db["l"]], "real value of Ytt")
    YttI = db.add_parameter_dc("YttI", [db["l"]], "imaginary value of Ytt")
    YftR = db.add_parameter_dc("YftR", [db["l"]], "real value of Yft")
    YftI = db.add_parameter_dc("YftI", [db["l"]], "imaginary value of Yft")
    YtfR = db.add_parameter_dc("YtfR", [db["l"]], "real value of Ytf")
    YtfI = db.add_parameter_dc("YtfI", [db["l"]], "imaginary value of Ytf")
    br = db["br"]
    for rec in db["l"]:
        r = br.find_record((rec.keys[0], "r")).value
        x = br.find_record((rec.keys[0], "x")).value
        b = br.find_record((rec.keys[0], "b")).value
        ratio = br.find_record((rec.keys[0], "ratio")).value
        angle = br.find_record((rec.keys[0], "angle")).value

        Ys = 1 / (r + x*1j)
        if ratio == 0.0:
            tap = 1.0
        else:
            tap = ratio
        tap *= np.exp(angle * (np.pi/180) * 1j)
        Ytt = Ys + (b/2)*1j
        Yff = Ytt / (tap * np.conj(tap))
        Yft = -Ys / np.conj(tap)
        Ytf = -Ys / tap

        YffR.add_record(rec.keys[0]).value = Yff.real
        YffI.add_record(rec.keys[0]).value = Yff.imag
        YttR.add_record(rec.keys[0]).value = Ytt.real
        YttI.add_record(rec.keys[0]).value = Ytt.imag
        YftR.add_record(rec.keys[0]).value = Yft.real
        YftI.add_record(rec.keys[0]).value = Yft.imag
        YtfR.add_record(rec.keys[0]).value = Ytf.real
        YtfI.add_record(rec.keys[0]).value = Ytf.imag

    YshR = db.add_parameter_dc("YshR", [db["b"]], "real value of Ysh")
    YshI = db.add_parameter_dc("YshI", [db["b"]], "imaginary value of Ysh")
    bus = db["bus"]
    for rec in db["b"]:
        Gs = bus.find_record((rec.keys[0], "Gs")).value
        Bs = bus.find_record((rec.keys[0], "Bs")).value
        YshR.add_record(rec.keys[0]).value = Gs / float(data["baseMVA"])
        YshI.add_record(rec.keys[0]).value = Bs / float(data["baseMVA"])

def write_gdx(data, case):
    # Write into a GDX file.
    ws = GamsWorkspace()
    db = ws.add_database()

    baseMVA = db.add_parameter("baseMVA", 0, "")
    baseMVA.add_record().value = float(data["baseMVA"])

    # Add raw data into db.
    add_bus(db, data)
    add_gen(db, data)
    add_branch(db, data)

    # Join bus-generator and bus-branch.
    add_join_bus_gen(db, data)
    add_join_bus_br(db, data)

    # Add a nodal admittance and shunt values.
    add_admittance_shunt(db, data)

    db.export(os.getcwd() + "/" + case + ".gdx")
