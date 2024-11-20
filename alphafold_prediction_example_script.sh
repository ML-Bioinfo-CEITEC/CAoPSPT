#PBS -N Alphafold-OOD-JOB
#PBS -q default
#PBS -l select=1:ncpus=8:ngpus=1:mem=100gb:scratch_local=10gb:gpu_cap=sm_75
#PBS -l walltime=100:00:00
#PBS -m ae

cd $PBS_O_WORKDIR

###############################################
#change paths
OUTPUTS_DIR=/storage/brno2/home/annaa/outputs/batch_01
# here you can set one file, directory or wildcard for files inside dir, e.g. '/dir/ecto/ry/ will process all files in the path, '/dir/*.fasta' will process all *.fasata files in dir    
FASTA_INPUTS=/storage/brno2/home/annaa/alphafold_batches/batch_01
###############################################

#choose one of the Alfafold prediction method
METHOD="AF_MONOMER_REDUCED_DBS"
#METHOD="AF_MONOMER_FULL_DBS"
#METHOD="AF_MULTIMER_FULL_DBS"
#METHOD="AF_MULTIMER_REDUCED_DBS"

# load methods
source ./lib.sh

# execute method
$METHOD

# copy results from scratch
mkdir -p $OUTPUTS_DIR
cp -r  $SCRATCHDIR/* $OUTPUTS_DIR/ &&  rm -r $SCRATCHDIR/*
