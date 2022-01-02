# Bunch generator

depending on the runtype:
- \-rt gen: generates bunch of solvent around molecule in cube of dimensions (ss\*pss)^3 cut off by radius -mr or by layer of thickness -lr in .xyz format
- \-rt trim: trims bunch of solvent around the molecule (specify number of atoms in solvent and number of atoms in molecule via -nats and -nat keywords) according to the -lr keyword

for full list of commandline arguments use ./bunch_generator.py --help
