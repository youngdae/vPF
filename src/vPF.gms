$ontext

We formulate a voltage-controlled power flow equation as a mixed
complementarity problem (MCP) and solve it using the PATH solver.

Contributor: Youngdae Kim (05.31.2020)

$offtext

$if not set case $set case case9
$if not set pf $set pf 0
$if not set qlim $set qlim 0
$if not set mpec_ineq $set mpec_ineq 0
$if not set mpec_ineq_noerr $set mpec_ineq_noerr 0
$if not set mpec_fb $set mpec_fb 0
$if not set fb_func $set fb_func 0
$if not set cs $set cs 0
$if not set cs_dev $set cs_dev 0.01
$if not set feas $set feas 0
$if not set relax_pg $set relax_pg 0
$if not set load_start $set load_start 0
$if not set tau $set tau 0
$if not set tau_dev $set tau_dev 0.01

sets
    b(*)             bus indices,
    bf(*)            bus fields,
    bus_ref(b)       slack bus index,
    bus_reg(b)       is a regulated bus,
    bus_cs(b)        bus index containing switched shunts,
    g(*)             generator indices,
    gf(*)            generator fields,
    gen_ref(g)       is generator attached to a slack bus
    gen_bus(g,b)     gen-bus mapping,
    l(*)             line indices,
    fr_bus(l,b)      branch-frombus mapping,
    to_bus(l,b)      branch-tobus mapping,
    tij_reg(l)       is tap at branch l a regulator,
    tij_reg_bus(l,b) tap l regulates bus b
    ;

parameters
    baseMVA      baseMVA,
    bus(b,bf)    bus table,
    gen(g,gf)    generator table,
    gij(l)       series conductance of branch l,
    bij(l)       series susceptance of branch l,
    bsh_ij(l)    shunt susceptance of branch l,
    gij_m(l)     magnetizing conductance,
    bij_m(l)     magnetizing susceptance,
    tij(l)       ratio of transformer l,
    phi_ij(l)    angle of transformer l,
    gs(b)        shunt conductance at bus b,
    bs(b)        shunt susceptance at bug s
    ;

$gdxin %case%
$load b bf bus_ref bus_reg bus_cs g gf gen_ref gen_bus l fr_bus to_bus tij_reg tij_reg_bus
$load baseMVA bus gen gij bij bsh_ij gij_m bij_m tij phi_ij gs bs
$gdxin

alias(b,bb);

variables
    VM(b)    voltage magnitude,
    VA(b)    voltage angle,
    PG(g)    real power,
    QG(g)    reative power,
    CS(b)    switched shunt,
    TAU(l)   transformer ratio,
    OBJ      objective variable
    ;

equations
    real_balance(b)         real power balance,
    reactive_balance(b)     reactive power balance,
    reactive_balance_pq(b)  reactive power balance at PQ buses,
    vc_qlimit(g)            voltage control using Q,
    defobj                  auxiliary objective
    ;

equations
    vc_sshunt_upv_p(b)       voltage control using switched shunts,
    vc_sshunt_upv_m(b)       voltage control using switched shunts,
    vc_sshunt_vpv_p(b)       voltage control using switched shunts,
    vc_sshunt_vpv_m(b)       voltage control using switched shunts,
    vc_sshunt_upv(b)         voltage control using switched shunts,
    real_balance_tau(b)      real power balance where tau is a variable,
    reactive_balance_tau(b)  reactive power balance where tau is a variable,
    vc_tau_upv_p(l)          voltage control using transformer,
    vc_tau_upv_m(l)          voltage control using transformer,
    vc_tau_vpv_p(l)          voltage control using transformer,
    vc_tau_vpv_m(l)          voltage control using transformer,
    vc_tau_upv(l)            voltage control using transformer
    ;

positive variables
    VM_CS_PLUS(b)   PV bus value goes over set-point,
    VM_CS_MINUS(b)  PV bus value goes below set-point,
    CS_PLUS(b)      switched shunt goes over set-point,
    CS_MINUS(b)     switched shunt goes below set-point,
    VM_TAU_PLUS(l)  PV bus value goes over set-point,
    VM_TAU_MINUS(l) PV bus value goes below set-point,
    TAU_PLUS(l)     ratio goes over set-point,
    TAU_MINUS(l)    ratio goes below set-point
    ;

* ----- Variables and equations for mpec_ineq -----
variables
    ERR1(b)          real power mismatch error,
    ERR2(b)          reactive power mismatch error,
    OBJ_MPEC_INEQ    objective variable for MPEC with inequalities
    ;

positive variables
    VM_Q_PLUS(b)     voltage magnitude plus,
    VM_Q_MINUS(b)    voltage magnitude minus
    ;

equations
    real_balance_err1(b)      real power allowing mismatch,
    reactive_balance_err2(b)  reactive power allowing mismatch,
    vm_q_plus_mpec_ineq(b)    complementarity voltage magnitude plus,
    vm_q_minus_mpec_ineq(b)   complementarity voltage magnitude minus,
    vm_setpoint_mpec(b)       setpoint deviation of regulated buses,
    defobj_mpec_ineq          objective equation for MPEC with inequalities
    ;
* ----- END of mpec_ineq -----

* ----- Fischer-Burmeister function -----
parameters
    alpha_q    weight on phi_q                  / 1e-4 /,
    alpha_r    weight on phi_r                  / 1e-1 /,
    alpha_u    weight on phi_u                  / 1e-1 /,
    alpha_va   weight on phi_va                 / 1e-3 /,
    alpha_d    weight on phi_d                  / 1e-3 /,
    alpha_y    weight on phi_y                  / 1e-3 /,
    alpha_z    weight on phi_z                  / 1e-3 /,
    delta_v    factor to normalize deviation    / 2e-1 /,
    xi         smoothing parameter              / 1e-8 /
    ;

variables
    OBJ_MPEC_FB    objective variable of MPEC with FB
    ;

equations
    vm_q_plus_mpec_fb(b)     complementarity voltage magnitude plus,
    vm_q_minus_mpec_fb(b)    complementarity voltage magnitude minus,
    defobj_mpec_fb           objective equation for MPEC with FB
    ;
* ----- END of Fischer-Burmeister function -----

parameters
    vm_setpoint(b)   set-point of voltage magnitude,
    vm_lo(b)         lower bound of voltage magnitude,
    vm_up(b)         upper bound of voltage magnitude,
    tau_lo(l)        lower bound of ratio,
    tau_up(l)        upper bound of ratio,
    cs_setpoint(b)   set-point of switched shunts,
    cs_lo(b)         lower bound of switched shunts,
    cs_up(b)         upper bound of switched shunts
    ;

real_balance(b)$(not bus_ref(b))..
    (sum(g$gen_bus(g,b), baseMVA*PG(g)) - bus(b,'Pd')) / baseMVA
    =E=
    (sum(l$fr_bus(l,b), gij(l)/sqr(tij(l)) + gij_m(l)) + sum(l$to_bus(l,b), gij(l)) + gs(b))*sqr(VM(b))
    - sum((l,bb)$(fr_bus(l,b) and to_bus(l,bb)), VM(b)*VM(bb)*((gij(l)/tij(l))*cos(VA(b) - VA(bb) - phi_ij(l)) + (bij(l)/tij(l))*sin(VA(b) - VA(bb) - phi_ij(l))))
    - sum((l,bb)$(fr_bus(l,bb) and to_bus(l,b)), VM(b)*VM(bb)*((gij(l)/tij(l))*cos(VA(b) - VA(bb) + phi_ij(l)) + (bij(l)/tij(l))*sin(VA(b) - VA(bb) + phi_ij(l))));

reactive_balance_pq(b)$(not bus_ref(b) and not bus_reg(b))..
    (-bus(b,'Qd')) / baseMVA
    =E=
    (sum(l$fr_bus(l,b), -((bij(l) + bsh_ij(l)/2)/sqr(tij(l))) + bij_m(l)) + sum(l$to_bus(l,b), -(bij(l)+bsh_ij(l)/2)) - bs(b) - CS(b))*sqr(VM(b))
    + sum((l,bb)$(fr_bus(l,b) and to_bus(l,bb)), VM(b)*VM(bb)*((bij(l)/tij(l))*cos(VA(b) - VA(bb) - phi_ij(l)) - (gij(l)/tij(l))*sin(VA(b) - VA(bb) - phi_ij(l))))
    + sum((l,bb)$(fr_bus(l,bb) and to_bus(l,b)), VM(b)*VM(bb)*((bij(l)/tij(l))*cos(VA(b) - VA(bb) + phi_ij(l)) - (gij(l)/tij(l))*sin(VA(b) - VA(bb) + phi_ij(l))));

model m_pf / real_balance.VA, reactive_balance_pq.VM /;

* ----- Voltage set-point regulation using reactive power -----
reactive_balance(b)$(not bus_ref(b))..
    (sum(g$gen_bus(g,b), baseMVA*QG(g)) - bus(b,'Qd')) / baseMVA
    =E=
    (sum(l$fr_bus(l,b), -((bij(l) + bsh_ij(l)/2)/sqr(tij(l))) + bij_m(l)) + sum(l$to_bus(l,b), -(bij(l)+bsh_ij(l)/2)) - bs(b) - CS(b))*sqr(VM(b))
    + sum((l,bb)$(fr_bus(l,b) and to_bus(l,bb)), VM(b)*VM(bb)*((bij(l)/tij(l))*cos(VA(b) - VA(bb) - phi_ij(l)) - (gij(l)/tij(l))*sin(VA(b) - VA(bb) - phi_ij(l))))
    + sum((l,bb)$(fr_bus(l,bb) and to_bus(l,b)), VM(b)*VM(bb)*((bij(l)/tij(l))*cos(VA(b) - VA(bb) + phi_ij(l)) - (gij(l)/tij(l))*sin(VA(b) - VA(bb) + phi_ij(l))));

vc_qlimit(g)..
    sum(b$gen_bus(g,b), VM(b) - vm_setpoint(b)) =E= 0;

model m_qlim / real_balance.VA, reactive_balance.VM, vc_qlimit.QG /;

* ----- Feasibility problem -----
defobj..
    OBJ =E= 0;

model m_feas / real_balance, reactive_balance, defobj /;

* ----- Voltage band regulation using switched shunts -----
vc_sshunt_upv_p(b)$bus_cs(b)..
    VM(b) - vm_lo(b) + VM_CS_MINUS(b) =G= 0;

vc_sshunt_upv_m(b)$bus_cs(b)..
    vm_up(b) - VM(b) + VM_CS_PLUS(b) =G= 0;

vc_sshunt_vpv_p(b)$bus_cs(b)..
    CS(b) - cs_lo(b) =G= 0;

vc_sshunt_vpv_m(b)$bus_cs(b)..
    cs_up(b) - CS(b) =G= 0;

vc_sshunt_upv(b)$bus_cs(b)..
    CS(b) - cs_setpoint(b) - CS_PLUS(b) + CS_MINUS(b) =E= 0;

model m_sshunt / real_balance.VA, reactive_balance.VM, vc_qlimit.QG,
                 vc_sshunt_upv_p.CS_PLUS, vc_sshunt_upv_m.CS_MINUS,
                 vc_sshunt_vpv_p.VM_CS_PLUS, vc_sshunt_vpv_m.VM_CS_MINUS,
                 vc_sshunt_upv.CS /;

* ----- Voltage band regulation using transformers -----
real_balance_tau(b)$(not bus_ref(b))..
    (sum(g$gen_bus(g,b), baseMVA*PG(g)) - bus(b,'Pd')) / baseMVA
    =E=
    (sum(l$fr_bus(l,b), gij(l)/sqr(TAU(l)) + gij_m(l)) + sum(l$to_bus(l,b), gij(l)) + gs(b))*sqr(VM(b))
    - sum((l,bb)$(fr_bus(l,b) and to_bus(l,bb)), VM(b)*VM(bb)*((gij(l)/TAU(l))*cos(VA(b) - VA(bb) - phi_ij(l)) + (bij(l)/TAU(l))*sin(VA(b) - VA(bb) - phi_ij(l))))
    - sum((l,bb)$(fr_bus(l,bb) and to_bus(l,b)), VM(b)*VM(bb)*((gij(l)/TAU(l))*cos(VA(b) - VA(bb) + phi_ij(l)) + (bij(l)/TAU(l))*sin(VA(b) - VA(bb) + phi_ij(l))));

reactive_balance_tau(b)$(not bus_ref(b))..
    (sum(g$gen_bus(g,b), baseMVA*QG(g)) - bus(b,'Qd')) / baseMVA
    =E=
    (sum(l$fr_bus(l,b), -((bij(l) + bsh_ij(l)/2)/sqr(TAU(l))) + bij_m(l)) + sum(l$to_bus(l,b), -(bij(l)+bsh_ij(l)/2)) - bs(b) - CS(b))*sqr(VM(b))
    + sum((l,bb)$(fr_bus(l,b) and to_bus(l,bb)), VM(b)*VM(bb)*((bij(l)/TAU(l))*cos(VA(b) - VA(bb) - phi_ij(l)) - (gij(l)/TAU(l))*sin(VA(b) - VA(bb) - phi_ij(l))))
    + sum((l,bb)$(fr_bus(l,bb) and to_bus(l,b)), VM(b)*VM(bb)*((bij(l)/TAU(l))*cos(VA(b) - VA(bb) + phi_ij(l)) - (gij(l)/TAU(l))*sin(VA(b) - VA(bb) + phi_ij(l))));

vc_tau_upv_p(l)$tij_reg(l)..
    sum(b$tij_reg_bus(l,b), VM(b) - vm_lo(b)) + VM_TAU_MINUS(l) =G= 0;

vc_tau_upv_m(l)$tij_reg(l)..
    sum(b$tij_reg_bus(l,b), vm_up(b) - VM(b)) + VM_TAU_PLUS(l) =G= 0;

vc_tau_vpv_p(l)$tij_reg(l)..
    TAU(l) - tau_lo(l) =G= 0;

vc_tau_vpv_m(l)$tij_reg(l)..
    tau_up(l) - TAU(l) =G= 0;

vc_tau_upv(l)$tij_reg(l)..
    TAU(l) - tij(l) - TAU_PLUS(l) + TAU_MINUS(l) =E= 0;

model m_tau / real_balance_tau.VA, reactive_balance_tau.VM, vc_qlimit.QG,
              vc_tau_upv_p.TAU_PLUS, vc_tau_upv_m.TAU_MINUS,
              vc_tau_vpv_p.VM_TAU_PLUS, vc_tau_vpv_m.VM_TAU_MINUS,
              vc_tau_upv.TAU /;

* ----- Voltage band regulation using transformers and switched shunts -----
model m_sshunt_tau / real_balance_tau.VA, reactive_balance_tau.VM, vc_qlimit.QG,
                     vc_sshunt_upv_p.CS_PLUS, vc_sshunt_upv_m.CS_MINUS,
                     vc_sshunt_vpv_p.VM_CS_PLUS, vc_sshunt_vpv_m.VM_CS_MINUS,
                     vc_sshunt_upv.CS,
                     vc_tau_upv_p.TAU_PLUS, vc_tau_upv_m.TAU_MINUS,
                     vc_tau_vpv_p.VM_TAU_PLUS, vc_tau_vpv_m.VM_TAU_MINUS,
                     vc_tau_upv.TAU /;

* ----- MPEC reformulation using inequalities -----
real_balance_err1(b)$(not bus_ref(b))..
    (sum(g$gen_bus(g,b), baseMVA*PG(g)) - bus(b,'Pd')) / baseMVA - ERR1(b)
    =E=
    (sum(l$fr_bus(l,b), gij(l)/sqr(tij(l)) + gij_m(l)) + sum(l$to_bus(l,b), gij(l)) + gs(b))*sqr(VM(b))
    - sum((l,bb)$(fr_bus(l,b) and to_bus(l,bb)), VM(b)*VM(bb)*((gij(l)/tij(l))*cos(VA(b) - VA(bb) - phi_ij(l)) + (bij(l)/tij(l))*sin(VA(b) - VA(bb) - phi_ij(l))))
    - sum((l,bb)$(fr_bus(l,bb) and to_bus(l,b)), VM(b)*VM(bb)*((gij(l)/tij(l))*cos(VA(b) - VA(bb) + phi_ij(l)) + (bij(l)/tij(l))*sin(VA(b) - VA(bb) + phi_ij(l))));

reactive_balance_err2(b)$(not bus_ref(b))..
    (sum(g$gen_bus(g,b), baseMVA*QG(g)) - bus(b,'Qd')) / baseMVA - ERR2(b)
    =E=
    (sum(l$fr_bus(l,b), -((bij(l) + bsh_ij(l)/2)/sqr(tij(l))) + bij_m(l)) + sum(l$to_bus(l,b), -(bij(l)+bsh_ij(l)/2)) - bs(b) - CS(b))*sqr(VM(b))
    + sum((l,bb)$(fr_bus(l,b) and to_bus(l,bb)), VM(b)*VM(bb)*((bij(l)/tij(l))*cos(VA(b) - VA(bb) - phi_ij(l)) - (gij(l)/tij(l))*sin(VA(b) - VA(bb) - phi_ij(l))))
    + sum((l,bb)$(fr_bus(l,bb) and to_bus(l,b)), VM(b)*VM(bb)*((bij(l)/tij(l))*cos(VA(b) - VA(bb) + phi_ij(l)) - (gij(l)/tij(l))*sin(VA(b) - VA(bb) + phi_ij(l))));

vm_q_plus_mpec_ineq(b)$bus_reg(b)..
    sum(g$gen_bus(g,b), (QG(g) - gen(g,'Qmin'))*VM_Q_PLUS(b)) =L= 0;

vm_q_minus_mpec_ineq(b)$bus_reg(b)..
    sum(g$gen_bus(g,b), (gen(g,'Qmax') - QG(g))*VM_Q_MINUS(b)) =L= 0;

vm_setpoint_mpec(b)$bus_reg(b)..
    VM(b) - vm_setpoint(b) - VM_Q_PLUS(b) + VM_Q_MINUS(b) =E= 0;

defobj_mpec_ineq..
    OBJ_MPEC_INEQ =E= sum(b$(not bus_ref(b)), sqr(ERR1(b)) + sqr(ERR2(b)));

model m_mpec_ineq / real_balance_err1, reactive_balance_err2,
                    vm_q_plus_mpec_ineq, vm_q_minus_mpec_ineq, vm_setpoint_mpec,
                    defobj_mpec_ineq /;

* ----- MPEC reformulation using Fischer-Burmeister function -----
vm_q_plus_mpec_fb(b)$(bus_reg(b))..
    sum(g$gen_bus(g,b), (QG(g) - gen(g,'Qmin')) + VM_Q_PLUS(b) -
        sqrt(sqr(QG(g) - gen(g,'Qmin')) + sqr(VM_Q_PLUS(b)) +2*xi)) =E= 0;

vm_q_minus_mpec_fb(b)$(bus_reg(b))..
    sum(g$gen_bus(g,b), (gen(g,'Qmax') - QG(g)) + VM_Q_MINUS(b) -
        sqrt(sqr(gen(g,'Qmax') - QG(g)) + sqr(VM_Q_MINUS(b)) +2*xi)) =E= 0;

defobj_mpec_fb..
    OBJ_MPEC_FB
    =E= alpha_q*(card(b) / card(bus_reg))*sum(b$bus_reg(b), sum(g$gen_bus(g,b),
                                      sqr((2*QG(g) - gen(g,'Qmax') - gen(g,'Qmin'))/(gen(g,'Qmax') - gen(g,'Qmin')))))
        + alpha_r*(card(b) / card(bus_reg))*sum(b$bus_reg(b), sqr((VM(b) - vm_setpoint(b)) / delta_v))
        + alpha_u*(card(b) / (card(b) - card(bus_reg)))*sum(b$(not bus_reg(b)), sqr((VM(b) - 1.0) / delta_v))
        + alpha_va*(card(b) / (card(b) - 1))*sum(b$(not bus_ref(b)), sqr(VA(b) / pi))
        + alpha_d*(card(b) / card(l))*sum(b$(not bus_ref(b)),
                                          sum((l,bb)$(fr_bus(l,b) and to_bus(l,bb)), sqr((VA(b) - VA(bb) - phi_ij(l))/pi)))
        + alpha_y*(card(b) / card(bus_reg))*sum(b$bus_reg(b), sqr(VM_Q_PLUS(b) / delta_v))
        + alpha_z*(card(b) / card(bus_reg))*sum(b$bus_reg(b), sqr(VM_Q_MINUS(b) / delta_v));

model m_mpec_fb / real_balance, reactive_balance,
                  vm_q_plus_mpec_fb, vm_q_minus_mpec_fb, vm_setpoint_mpec,
                  defobj_mpec_fb /;

vm_setpoint(b) = min(bus(b,'Vmax'), max(bus(b,'Vmin'), bus(b,'Vm')));
vm_lo(b) = bus(b,'Vmin');
vm_up(b) = bus(b,'Vmax');

VM.L(b) = vm_setpoint(b);
VA.L(b) = bus(b,'Va');

cs_setpoint(b) = min(bus(b,'Csmax'), max(bus(b,'Csmin'), bus(b,'Cs')));
cs_lo(b) = bus(b,'Csmin');
cs_up(b) = bus(b,'Csmax');

CS.L(b) = cs_setpoint(b);
CS.FX(b) = CS.L(b);

QG.L(g) = min(gen(g,'Qmax'), max(gen(g,'Qmin'), gen(g,'Qg')));
QG.LO(g) = gen(g,'Qmin');
QG.UP(g) = gen(g,'Qmax');

* Although we do not include slack bus in our model, we fix
* variable values associated with it to properly define complementarity.
PG.FX(g) = min(gen(g,'Pmax'), max(gen(g,'Pmin'), gen(g,'Pg')));
QG.FX(g)$(gen_ref(g)) = QG.L(g);

VM.FX(b)$(bus_ref(b)) = VM.L(b);
VA.FX(b)$(bus_ref(b)) = VA.L(b);

* If feas is specified, try to find a feasible solution and
* fix PG to the feasible value.
if (%feas% > 0,
    VM.LO(b) = bus(b,'Vmin');
    VM.UP(b) = bus(b,'Vmax');
    CS.LO(b) = bus(b,'Csmin');
    CS.UP(b) = bus(b,'Csmax');

    if (%relax_pg% > 0,
        PG.LO(g) = gen(g,'Pmin');
        PG.UP(g) = gen(g,'Pmax');
    );

    solve m_feas minimizing obj using nlp;
    abort$(m_feas.modelstat ne 2) "Infeasible";
    execute_unload '%case%_start', PG, VM, CS;

* Reset values.
    PG.FX(g) = PG.L(g);
    QG.L(g) = min(gen(g,'Qmax'), max(gen(g,'Qmin'), gen(g,'Qg')));
    QG.FX(g)$(gen_ref(g)) = QG.L(g);
    VM.FX(b)$(bus_ref(b)) = VM.L(b);
    VM.L(b)$(not bus_ref(b)) = vm_setpoint(b);
    VM.LO(b)$(not bus_ref(b)) = -inf;
    VM.UP(b)$(not bus_ref(b)) = inf;
    VA.L(b)$(not bus_ref(b)) = bus(b,'Va');
    cs_setpoint(b) = CS.L(b);
    CS.FX(b) = cs_setpoint(b);
);

if (%load_start% > 0,
    execute_load '%case%_start', PG, VM, CS;
    PG.FX(g) = PG.L(g);
    VM.FX(b)$(bus_ref(b)) = VM.L(b);
    VM.L(b)$(not bus_ref(b)) = vm_setpoint(b);
    VM.LO(b)$(not bus_ref(b)) = -inf;
    VM.UP(b)$(not bus_ref(b)) = inf;
    cs_setpoint(b) = CS.L(b);
    CS.FX(b) = cs_setpoint(b);
);

if (%pf% > 0,
    VM.FX(b)$(bus_ref(b) or bus_reg(b)) = VM.L(b);

    m_pf.optfile = 1;
    solve m_pf using mcp;

* Reset bounds.
    VM.LO(b)$(not bus_ref(b)) = -inf;
    VM.UP(b)$(not bus_ref(b)) = inf;
);

if (%qlim% > 0,
    m_qlim.optfile = 1;
    solve m_qlim using mcp;
);

if (%cs% > 0,
$ontext
* For case3120sp
    bus_cs(b) = no;
    bus_cs('2306') = yes;
    bus_cs('2830') = yes;
$offtext
$ontext
* For Network_06O-124_scenario17
    bus_cs(b) = no;
    bus_cs('522') = yes; bus_cs('523') = yes; bus_cs('530') = yes; bus_cs('531') = yes;
    bus_cs('532') = yes; bus_cs('533') = yes; bus_cs('536') = yes; bus_cs('537') = yes;
*    bus_cs('539') = yes; bus_cs('542') = yes; bus_cs('543') = yes; bus_cs('544') = yes;
*    bus_cs('545') = yes; bus_cs('546') = yes; bus_cs('547') = yes; bus_cs('558') = yes;
*    bus_cs('559') = yes; bus_cs('560') = yes; bus_cs('561') = yes; bus_cs('562') = yes;
*    bus_cs('563') = yes; bus_cs('567') = yes; bus_cs('568') = yes; bus_cs('569') = yes;
*    bus_cs('570') = yes; bus_cs('594') = yes; bus_cs('596') = yes; bus_cs('597') = yes;
*    bus_cs('1325') = yes; bus_cs('1326') = yes; bus_cs('1360') = yes; bus_cs('1361') = yes;
*    bus_cs('1362') = yes; bus_cs('1363') = yes; bus_cs('1364') = yes; bus_cs('1365') = yes;
$offtext
    cs_lo(b) = cs_setpoint(b) - %cs_dev%$bus_cs(b);
    cs_up(b) = cs_setpoint(b) + %cs_dev%$bus_cs(b);
    CS.LO(b) = -inf;
    CS.UP(b) = inf;
    CS.FX(b)$(not bus_cs(b)) = cs_setpoint(b);
    CS_PLUS.FX(b)$(not bus_cs(b)) = 0;
    CS_MINUS.FX(b)$(not bus_cs(b)) = 0;
    VM_CS_PLUS.FX(b)$(not bus_cs(b)) = 0;
    VM_CS_MINUS.FX(b)$(not bus_cs(b)) = 0;
);

if (%tau% > 0,
$ontext
* For case3120sp
    tij_reg(l) = no;
    tij_reg_bus(l,b) = no;
    tij_reg('215') = yes;
    tij_reg_bus('215','316') = yes;
    tij_reg('2594') = yes;
    tij_reg_bus('2594','2774') = yes;
    tij_reg('2652') = yes;
    tij_reg_bus('2652','2930') = yes;
$offtext
$ontext
* For Network_06O-124_scenario17
    tij_reg(l) = no;
    tij_reg_bus(l,b) = no;
    tij_reg('872') = yes;
    tij_reg_bus('872','518') = yes;
    tij_reg('894') = yes;
    tij_reg_bus('894','554') = yes;
    tij_reg('914') = yes;
    tij_reg_bus('914','581') = yes;
    tij_reg('915') = yes;
    tij_reg_bus('915','582') = yes;
    tij_reg('921') = yes;
    tij_reg_bus('921','590') = yes;
    tij_reg('931') = yes;
    tij_reg_bus('931','599') = yes;
$offtext
    tau_lo(l) = tij(l) - %tau_dev%$tij_reg(l);
    tau_up(l) = tij(l) + %tau_dev%$tij_reg(l);

    TAU.L(l) = tij(l);
    TAU.LO(l) = -inf;
    TAU.UP(l) = inf;
    TAU.FX(l)$(not tij_reg(l)) = tij(l);
    TAU_PLUS.FX(l)$(not tij_reg(l)) = 0;
    TAU_MINUS.FX(l)$(not tij_reg(l)) = 0;
    VM_TAU_PLUS.FX(l)$(not tij_reg(l)) = 0;
    VM_TAU_MINUS.FX(l)$(not tij_reg(l)) = 0;
);

if (%cs% > 0 or %tau% > 0,
    if (%tau% eq 0,
        m_sshunt.optfile = 1;
        solve m_sshunt using mcp;
    elseif %cs% eq 0,
        m_tau.optfile = 1;
        solve m_tau using mcp;
    else
        m_sshunt_tau.optfile = 1;
        solve m_sshunt_tau using mcp;
    );
);

if (%mpec_ineq% > 0,
    m_mpec_ineq.dictfile = 0;
    if (%mpec_ineq_noerr% > 0,
        ERR1.FX(b) = 0;
        ERR2.FX(b) = 0;
    );
    gen(g,'Qmin')$(gen(g,'Qmin') eq -inf) = -99999.0;
    gen(g,'Qmax')$(gen(g,'Qmax') eq inf) = 99999.0;
    solve m_mpec_ineq minimizing OBJ_MPEC_INEQ using nlp;
);

if (%mpec_fb% > 0,
    m_mpec_fb.dictfile = 0;
    gen(g,'Qmin')$(gen(g,'Qmin') eq -inf) = -99999.0;
    gen(g,'Qmax')$(gen(g,'Qmax') eq inf) = 99999.0;
    solve m_mpec_fb minimizing OBJ_MPEC_FB using nlp;
);

* Check if voltage bounds are violated.
scalar idx, max_pv_viol, pv_viol, pq_viol, pv_set, binding, ndegen;
scalar sp, lev, lo, up, max_vm_viol, is_fr, is_to;

max_pv_viol = 0;
pv_viol = 0;
pq_viol = 0;
pv_set = 0;
binding = 0;
ndegen = 0;
max_vm_viol = 0;
loop(b$(not bus_ref(b)),
    if (VM.L(b) < bus(b,'Vmin') or VM.L(b) > bus(b,'Vmax'),
        idx = b.val;
        sp = vm_setpoint(b); lev = VM.L(b); lo = bus(b,'Vmin'); up = bus(b,'Vmax');
        if (bus_reg(b),
            is_fr = sum(l$fr_bus(l,b), 1);
            is_to = sum(l$to_bus(l,b), 1);
            display "VM (REG) out of bounds.", idx, lo, up, lev, sp, is_fr, is_to;
            lev = sum(g$gen_bus(g,b), QG.L(g));
            lo =  sum(g$gen_bus(g,b), QG.LO(g));
            up =  sum(g$gen_bus(g,b), QG.UP(g));
            display "QG associated with this bus.", lo, up, lev;
            pv_viol = pv_viol + 1;
            max_pv_viol = max(max_pv_viol, bus(b,'Vmin') - VM.L(b));
            max_pv_viol = max(max_pv_viol, VM.L(b) - bus(b,'Vmax'));
        else
            display "VM (NON-REG) out of bounds.", idx, lo, up, lev, sp;
            pq_viol = pq_viol + 1;
        );
    );

    if (bus_reg(b),
        idx = b.val;
        if (round(VM.L(b) - vm_setpoint(b), 6),
            sp = vm_setpoint(b); lev = VM.L(b); lo = bus(b,'Vmin'); up = bus(b,'Vmax');
            display "VM (REG) out of set point.", idx, lo, up, lev, sp;
            pv_set = pv_set + 1;
            if (abs(sp - lev) > max_vm_viol,
                max_vm_viol = abs(sp - lev);
            );
        );
        lev = sum(g$gen_bus(g,b), QG.L(g));
        lo = sum(g$gen_bus(g,b), QG.LO(g));
        up = sum(g$gen_bus(g,b), QG.UP(g));
        if ((VM.L(b) eq vm_setpoint(b)) and (lev eq lo or lev eq up),
            sp = vm_setpoint(b);
            display "VM (REG) degenerate.", idx, sp;
            display "QG associated with this bus.", lo, up, lev;
            ndegen = ndegen + 1;
        );
    );
);

* Check if there is a binding reactive power.
loop(g$(not gen_ref(g)),
    if (QG.L(g) eq QG.LO(g) or QG.L(g) eq QG.UP(g),
        idx = g.val;
        lev = QG.L(g); lo = QG.LO(g); up = QG.UP(g);
        display "Q is binding.", idx, lo, up, lev;
        binding = binding + 1;
    );
);

file out / '' /;
put out;
put 'Number of regulated buses  . . . . . . . . . ', card(bus_reg) /;
put 'Number of VM bounds violations (REG) . . . . ', pv_viol /;
put 'Number of VM bounds violations (NON-REG) . . ', pq_viol /;
put 'Number of VM out of set-point. . . . . . . . ', pv_set /;
put 'Number of QG binding . . . . . . . . . . . . ', binding /;
put 'Number of degenerate cases . . . . . . . . . ', ndegen /;
put 'Maximum voltage magnitude. . . . . . . . . . ', smax(b, VM.L(b)) /;
put 'Minimum voltage magnitude. . . . . . . . . . ', smin(b, VM.L(b)) /;
out.nr = 2;
out.nd = 15;
out.nw = 0;
if (%tau% > 0,
    put 'Real power mismatch. . . . . . . . . . . . . ', smax(b$(not bus_ref(b)), real_balance_tau.L(b) - bus(b,"Pd")/baseMVA) /;
    put 'Reactive power mismatch. . . . . . . . . . . ', smax(b$(not bus_ref(b)), reactive_balance_tau.L(b) - bus(b,"Qd")/baseMVA) /;
elseif %mpec_ineq% > 0,
    put 'Real power mismatch. . . . . . . . . . . . . ', smax(b$(not bus_ref(b)), real_balance_err1.L(b) + ERR1.L(b) - bus(b,"Pd")/baseMVA) /;
    put 'Reactive power mismatch. . . . . . . . . . . ', smax(b$(not bus_ref(b)), reactive_balance_err2.L(b) + ERR2.L(b) - bus(b,"Qd")/baseMVA) /;
else
    put 'Real power mismatch. . . . . . . . . . . . . ', smax(b$(not bus_ref(b)), real_balance.L(b) - bus(b,"Pd")/baseMVA) /;
    put 'Reactive power mismatch. . . . . . . . . . . ', smax(b$(not bus_ref(b)), reactive_balance.L(b) - bus(b,"Qd")/baseMVA) /;
);
if (%qlim% > 0,
    put 'Complementarity error  . . . . . . . . . . . ', m_qlim.objval /;
elseif %cs% > 0 or %tau% > 0,
    if (%tau% eq 0,
        put 'Complementarity error  . . . . . . . . . . . ', m_sshunt.objval /;
    elseif %cs% eq 0,
        put 'Complementarity error  . . . . . . . . . . . ', m_tau.objval /;
    else
        put 'Complementarity error  . . . . . . . . . . . ', m_sshunt_tau.objval /;
    );
);

put 'Maximum deviation from setpoint (VM) . . . . ', max_vm_viol /;
put 'Maximum band violations  . . . . . . . . . . ', max_pv_viol /;
putclose;


