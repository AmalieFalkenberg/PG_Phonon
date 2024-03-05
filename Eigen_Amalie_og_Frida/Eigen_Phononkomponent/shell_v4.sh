#!/usr/bin/env bash

mcrun -c --mpi=6 phonon_eigenvector.instr -d PH3_h0_l3.95 -n 1000000 -N 41  E=0,20  Ef=25  Dlambda=0.1  h=0  l=3.95  dA3=-90 Temp=2 width=0.005  phononmode=3  E_steps_high=20  E_steps_low=20  Verbose=0

mcrun -c --mpi=6 phonon_eigenvector.instr -d PH3_h0_l4.0 -n 1000000 -N 41  E=0,20  Ef=25  Dlambda=0.1  h=0  l=4.0  dA3=-90 Temp=2 width=0.005  phononmode=3  E_steps_high=20  E_steps_low=20  Verbose=0

mcrun -c --mpi=6 phonon_eigenvector.instr -d PH3_h0_l4.05 -n 1000000 -N 41  E=0,20   Ef=25  Dlambda=0.1  h=0  l=4.05  dA3=-90 Temp=2 width=0.005  phononmode=3  E_steps_high=20  E_steps_low=20  Verbose=0



mcrun -c --mpi=6 phonon_eigenvector.instr -d PH4_h0_l3.95 -n 1000000 -N 41  E=0,20   Ef=25  Dlambda=0.1  h=0  l=3.95  dA3=-90 Temp=2 width=0.005  phononmode=4  E_steps_high=20  E_steps_low=20  Verbose=0

mcrun -c --mpi=6 phonon_eigenvector.instr -d PH4_h0_l4.0 -n 1000000 -N 41  E=0,20  Ef=25  Dlambda=0.1  h=0  l=4.0  dA3=-90 Temp=2 width=0.005  phononmode=4  E_steps_high=20  E_steps_low=20  Verbose=0

mcrun -c --mpi=6 phonon_eigenvector.instr -d PH4_h0_l4.05 -n 1000000 -N 41  E=0,20  Ef=25  Dlambda=0.1  h=0  l=4.05  dA3=-90 Temp=2 width=0.005  phononmode=4  E_steps_high=20  E_steps_low=20  Verbose=0



mcrun -c --mpi=6 phonon_eigenvector.instr -d PH5_h0_l3.95 -n 1000000 -N 41  E=0,20   Ef=25  Dlambda=0.1  h=0  l=3.95  dA3=-90 Temp=2 width=0.005  phononmode=5  E_steps_high=20  E_steps_low=20  Verbose=0

mcrun -c --mpi=6 phonon_eigenvector.instr -d PH5_h0_l4.0 -n 1000000 -N 41  E=0,20   Ef=25  Dlambda=0.1  h=0  l=4.0  dA3=-90 Temp=2 width=0.005  phononmode=5  E_steps_high=20  E_steps_low=20  Verbose=0

mcrun -c --mpi=6 phonon_eigenvector.instr -d PH5_h0_l4.05 -n 1000000 -N 41  E=0,20  Ef=25  Dlambda=0.1  h=0  l=4.05  dA3=-90 Temp=2 width=0.005  phononmode=5  E_steps_high=20  E_steps_low=20  Verbose=0
