#!/bin/csh
#PBS -N peakmem
#PBS -A NMMM0015
#PBS -l walltime=03:00:00
#PBS -q regular
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -l select=1:ncpus=1

setenv TMPDIR /glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

### Load the module environment
module purge
##ignore that the sticky module ncarenv/... is not unloaded
setenv LMOD_TMOD_FIND_FIRST yes
module load ncarenv/23.09
module use /glade/work/epicufsrt/contrib/spack-stack/derecho/modulefiles
module load ecflow/5.8.4
module load mysql/8.0.33

module use /glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-1.6.0/envs/unified-env/install/modulefiles/Core
module load stack-gcc/12.2.0
module load stack-cray-mpich/8.1.25
module load stack-python/3.10.13

#module available
module load jedi-fv3-env jedi-mpas-env soca-env

ulimit -s unlimited
export GFORTRAN_CONVERT_UNIT='big_endian:101-200'



module load peak_memusage


peak_memusage.exe ./goes2mpas.x

