## Running MrBayes ##

# Protein:

exe *nex
lset Nst=1 Rates=Invgamma nucmodel=Protein
set Usebeagle=yes Beagleresource=1 Beagledevice=GPU Beagleprecision=single
mcmc ngen=100000

sump
sumt