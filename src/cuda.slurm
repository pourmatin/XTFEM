#!/bin/bash
#----------------------------------------------------
# Example SLURM job script to run MPI applications
#----------------------------------------------------
#SBATCH -J kokkos-xtfem3d     	# Job name
#SBATCH -o xtfem.o%j      	# Name of stdout output file
#SBATCH -e xtfem.o%j      	# Name of stdout output file
#SBATCH -p gpu	 		# Queue name
#SBATCH -N 4               	# Total number of nodes requested
#SBATCH -n 4              	# Total number of mpi tasks requested
#SBATCH -t 02:00:00        	# Run time (hh:mm:ss)



# Launch the MPI application using ibrun
ibrun tacc_affinity ./xtfem3d.cuda

#Note:  ibrun does not bind tasks or threads by default
#       To bind task/threads to sockets or cores, you must use
#       ibrun with tacc_affinity or use the "srun" command directly
#       with the "--cpu_bind" option.
