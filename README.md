# Bunch generator

depending on the runtype:
- \-rt gen: generates bunch of solvent around molecule in sphere of radius -ss * -pss or in layer of thickness -mr in .xyz format
- \-rt trim: trims bunch of solvent around the molecule (specify number of atoms in solvent and number of atoms in molecule via -nats and -nat keywords) according to the same keywords

total list of cl-arguments:
  -h, --help            show this help message and exit
  -mol molecule, --mol molecule
                        name of molecule file
  -sol solvent, --sol solvent
                        name of solvent file
  -ss N, --ss N         slab size in Angstroem
  -pss N, --pss N       size of single molecule occupied space in Angstroem
  -mr N, --mr N         cutoff radius in Angstroem
  -cr N, --cr N         collide solvent-molecule radius in Angstroem
  -lr N, --lr N         layer radius in Angstroem, set huge number to generate sphere of --mr
  -rt type, --rt type   Runtype, use gen for generation or trim for trimming with -nat set to number of atoms in the beginning of .xyz as comparator
  -nat N, --nat N       number of atoms in the beginning of .xyz to be used as comparator in trimming procedure
  -nats N, --nats N     number of atoms in each solvent molecule

