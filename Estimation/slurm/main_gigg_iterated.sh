#!/bin/bash -l
#SBATCH --time=0:45:00
#SBATCH --cpus-per-task=25
#SBATCH --mem=45000M
#SBATCH --output=Outfiles/gigg_saveak.out

export t_bash=0 sv_bash=1 trend_bash=1 almonrest_bash=1  group_sparse_bash=1 ortho_choice_bash=1 dat_choice_bash=1

module load matlab

srun matlab -nodisplay -r "run_main_iterated($t_bash,$sv_bash,$trend_bash,$almonrest_bash,$group_sparse_bash,$ortho_choice_bash,$dat_choice_bash) ; exit(0)"
