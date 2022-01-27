import os
import numpy as np

mol = 'BDP.xyz'
sol = 'MECN.xyz'
nat = 31
nats = 6 
lr = 6.5
pss = 6.2
method = '--gfn 0'
solvent = 'acetonitrile'
mdparams = "$metadyn\n save=10\n kpush=0.02\n alp=1\n$end\n$md\n time=10\n step=1\n temp=298.15\n$end\n"
mdparams += "$wall\n potential=logfermi\n sphere:auto, all\n autoscale=1.1\n$end"
cleanup = "rm gfnff_*;rm xtbrestart"

#data for 1 Mole of liquid at 1 atm to 1M solution conversion
den = 0.786 # CHCl3 
mmass = 41.05 # CHCl3

selfc = den * 1000 / mmass
expls = 29 
lstdfac = 1.987204 * 298.15 * np.log(selfc/expls) /1000

generate = [cleanup, 
"python bunch_generator.py --mol %s --sol %s --cr 2.5 --lr 999 --mr 25 -pss %f" %(mol, sol, pss),
"echo \'%s\' > microdrop.inp; cat .xcontrol >> microdrop.inp; echo '$end' >> microdrop.inp" %(mdparams),
"xtb slab.xyz --input microdrop.inp --metadyn --gff > STG1.out",
"cat xtb.trj | grep energy | cut -d\' \' -f3 > mdtrjenergy.csv",
"head -1 xtb.trj > lasttrjpoint.xyz",
"tail -r xtb.trj | grep -m 1 -B 9999 energy | tail -r >> lasttrjpoint.xyz",
"python bunch_generator.py --mol lasttrjpoint.xyz --rt trim --nat %i --nats %i --lr %f" %(nat, nats, lr)]

optimize = [cleanup,
"xtb slabtrimmed.xyz --opt vtight %s > STG2.out" %(method),
cleanup,
"mv xtbopt.xyz STG2.xyz",
"xtb slabtrimmed_solvonly.xyz --opt vtight %s > STG2_S.out" %(method),
cleanup,
"mv xtbopt.xyz STG2_S.xyz",
"xtb STG2.xyz --opt vtight %s --gbsa %s > STG2_%s.out" %(method, solvent, solvent),
cleanup,
"mv xtbopt.xyz STG2_%s.xyz" %(solvent),
"xtb STG2_S.xyz --opt vtight %s --gbsa %s > STG2_%s_S.out" %(method, solvent, solvent),
cleanup,
"mv xtbopt.xyz STG2_%s_S.xyz" %(solvent),
"xtb %s --opt vtight %s > MOL.out" %(mol, method),
cleanup,
"mv xtbopt.xyz %s" %(mol)]

thermochem = ["xtb %s --bhess  %s > MOL_THERMO.out" %(mol, method),
"cat MOL_THERMO.out | grep \'TOTAL FREE\' | tr -s \' \' | cut -d\' \' -f 6 > datapoints.txt",
cleanup,
"xtb STG2.xyz --bhess  %s > STG2_THERMO.out" %(method),
"cat STG2_THERMO.out | grep \'TOTAL FREE\' | tr -s \' \' | cut -d\' \' -f 6 >> datapoints.txt",
cleanup,
"xtb STG2_S.xyz --bhess  %s > STG2_S_THERMO.out" %(method),
"cat STG2_S_THERMO.out | grep \'TOTAL FREE\' | tr -s \' \' | cut -d\' \' -f 6 >> datapoints.txt",
cleanup,
"xtb STG2_%s.xyz --bhess  %s --gbsa %s > STG2_%s_THERMO.out" %(solvent, method, solvent, solvent),
"cat STG2_%s_THERMO.out | grep \'TOTAL FREE\' | tr -s \' \' | cut -d\' \' -f 6 >> datapoints.txt" %(solvent),
cleanup,
"xtb STG2_%s_S.xyz --bhess  %s --gbsa %s > STG2_%s_S_THERMO.out" %(solvent, method, solvent, solvent),
"cat STG2_%s_S_THERMO.out | grep \'TOTAL FREE\' | tr -s \' \' | cut -d\' \' -f 6 >> datapoints.txt" %(solvent),
cleanup]

print('--- GENERATION STEP ---')
for i in generate:
    os.system(i)
print('--- OPTIMIZATION STEP ---')
for i in optimize:
    os.system(i)
print('--- THERMOCHEMISTRY STEP ---')
for i in thermochem:
    os.system(i)

path = os.path.abspath(__file__)[:-len(os.path.basename(__file__))]
gees = np.array(np.genfromtxt(str(path) + 'datapoints.txt'))
deltagee = (gees[1] - gees[0] - gees[2] + (gees[3] - gees[1]) - (gees[4] - gees[2])) * 627.5
deltageecorr = ( deltagee - 1.89 - lstdfac ) * 4.184
print('solvation free energy equals ' + str(deltageecorr) + ' kJ/mol')
