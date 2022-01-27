import os, sys, argparse
import numpy as np

path = os.path.abspath(__file__)[:-len(os.path.basename(__file__))]

# these are to be commandline parameters

parser = argparse.ArgumentParser(description='Generate molecule with explicit solvation environment')
parser.add_argument('-mol','--mol', metavar='molecule', type=str, default='molecule.xyz',
                    help='name of molecule file')
parser.add_argument('-sol','--sol', metavar='solvent', type=str, default='chcl3',
                    help='name of solvent, available values are ...')
parser.add_argument('-method','--method', metavar='--gfn x', type=string, default='--gff',
                    help='xTB method to use for calculation')
parser.add_argument('-rt','--rt', metavar='type', type=str, default='normal',
                    help='Runtype, use normal for full run of solvation free energy')
parser.add_argument('-nat','--nat', metavar='N', type=int, default=31,
                    help='number of atoms in the beginning of .xyz to be used as comparator in trimming procedure')

args = parser.parse_args()

cl_solvent = args.sol 
mol = args.mol
nat = args.nat
method = args.method
runtyp = args.rt

#    flat file with all of the parameters (tab-separated):
#    solvent,nats,sol,den,mmass,lr,pss
#    chcl3,5,CHCl3.xyz,1.49,119.38,6.5,6.2

# parse solvent parameters or die

solvent, nats, sol, den, mmass, lr, pss = np.array(str(path) + ".dropparams", delimiter = ',', unpack=True)

notfound = True
for i, j in enumerate(solvent):
    if j != cl_solvent:
        notfound = True
        continue
    else:
        notfound = False
        solvent = solvent[i]
        nats = nats[i]
        sol = sol[i]
        den = den[i]
        mmass = mmass[i]
        lr = lr[i]
        pss = pss[i]
sys.exit('solvent not found!') if notfound

# i prefer to leave these hardcoded just for now

mdparams = "$metadyn\n save=10\n kpush=0.02\n alp=1\n$end\n$md\n time=10\n step=1\n temp=298.15\n$end\n"
mdparams += "$wall\n potential=logfermi\n sphere:auto, all\n autoscale=1.1\n$end"

cleanup = "rm gfnff_*;rm xtbrestart"

### definition of the external command lists for subroutines

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
"xtb STG2_%s.xyz --bhess  %s --gbsa %s reference > STG2_%s_THERMO.out" %(solvent, method, solvent, solvent),
"cat STG2_%s_THERMO.out | grep \'TOTAL FREE\' | tr -s \' \' | cut -d\' \' -f 6 >> datapoints.txt" %(solvent),
cleanup,
"xtb STG2_%s_S.xyz --bhess  %s --gbsa %s reference > STG2_%s_S_THERMO.out" %(solvent, method, solvent, solvent),
"cat STG2_%s_S_THERMO.out | grep \'TOTAL FREE\' | tr -s \' \' | cut -d\' \' -f 6 >> datapoints.txt" %(solvent),
cleanup]

def generate(generate)
    print('--- GENERATION STEP ---')
    for i in generate:
        os.system(i)

def omtimize(optimize)
    print('--- OPTIMIZATION STEP ---')
    for i in optimize:
        os.system(i)

def geecalc(path, thermochem)
    print('--- THERMOCHEMISTRY STEP ---')
    for i in thermochem:
        os.system(i)
        gees = np.array(np.genfromtxt(str(path) + 'datapoints.txt'))
    return gees

def dgeecalc(path, gees, den, mmass):
    # amount of explicit solvent molecules, unoptimal way to get it is to parse file generated upon 
    # bunch_generator.py run each time.
    expls = np.genfromtxt(str(path) + '.samnt')
    # equilibrium concentration of 1 mole of solution at 1 atm
    selfc = den * 1000 / mmass
    # term to correct to 1 mol/L solution standard state
    lstdfac = 1.987204 * 298.15 * np.log(selfc/expls) /1000
    # calculate to thermodynamic cycle as in reference [], 627.5 is for conversion from eH to kcal/mol
    deltagee = (gees[1] - gees[0] - gees[2] + (gees[3] - gees[1]) - (gees[4] - gees[2])) * 627.5
    # adding all correction factors and converting to kJ/mol
    deltageecorr = ( deltagee - 1.89 - lstdfac ) * 4.184
    print('solvation free energy equals ' + str(deltageecorr) + ' kJ/mol')

def normal(path, generate, optimize, thermochem, den, mmass):
    generate(generate)
    optimize(optimize)
    gees = thermochem(path, thermochem)
    dgeecalc(path, gees, den, mmass)
