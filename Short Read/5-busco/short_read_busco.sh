
#!/bin/bash
#SBATCH --job-name=busco_initial_tutorial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --mail-user=jacqueline.guillemin@uconn.edu
#SBATCH -o busco_%j.out
#SBATCH -e busco_%j.err

module load busco/3.0.2b
module unload augustus
export PATH=/home/CAM/jguillemin/augustus/bin:/home/CAM/jguillemin/augustus/scripts:$PATH
export AUGUSTUS_CONFIG_PATH=$HOME/augustus/config
###if your run crashes uncomment the following:
module unload blast/2.7.1
module load blast/2.2.31

run_BUSCO.py -i /home/CAM/jguillemin/Assembly_tutorial/Assembly_tutorial/Assembly/SPAdes_out/scaffolds.fasta -l /home/CAM/jguillemin/Assembly_tutorial/Assembly_tutorial/bacteria_odb9 -o bacterial_short_read_tutorial_busco -m geno -c 1
