# vPF
## Generating data files in GDX format
```python
$ python parse_matpower_gams.py matpowerfile
```
where
* matpowerfile: the path to the MATPOWER file, such as case1354pegase.m.

```python
$ python parse_go_gams.py gofile
```
where
* gofile: the path to the scenario directory of the GO data files.

These scripts will generate data files in GDX format that subsequently are read
by our `vPF.gms` model.

## Solving power flow with voltage regulation by generator reactive power
```gams
$ gams vPF.gms --case=casename --qlim=1
```
where
* casename: the path to the case GDX file
* qlim: solve power flow with voltage regulation by generator reactive power

Optional arguments `--feas=1` and `--relax_pg=1` enable us to compute a feasible
real power and voltage reference value.
Many of the GO data are infeasible if we were to solve power flow starting from
the given initial point.
In that case, we compute a feasible point first and solve the problem from that
initial point.
For example, the following commands performs this:
```gams
$ gams vPF.gms --case=casename --qlim=1 --feas=1 --relax_pg=1 nlp=ipopt
```

## Solving power flow with voltage regulation by generator reactive power, transformer, and switched shunt devices
```gams
$ gams vPF.gms --case=casename --tau=1 --tau_dev=dev --cs=1 --cs_dev=dev
```
where
* tau: solve power flow with voltage regulation by generator reactive power and
transformers
* tau_dev: specify the allowed range of tap values.
* cs: solve power flow with voltage regulation by generator reactive power and
switched shunt devices
* cs_dev: specify the allowed range of switched shunt device values.


