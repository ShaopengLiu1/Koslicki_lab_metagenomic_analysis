# interactive code to run CAMI simulation


### go to the output folder
/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210130_CAMISIM_errorfree_simulation/output

### generate regular simulation config file
bash /data/sml6467/github/Koslicki_lab_metagenomic_analysis/src/CAMISIM_set_config.sh  simulation_input_data.txt /data/sml6467/github/Koslicki_lab_metagenomic_analysis/src/git_repo/CAMISIM/defaults/personal_run/templete_mini_config.ini


### generate error-free simulation config file
bash /data/sml6467/github/Koslicki_lab_metagenomic_analysis/src/CAMISIM_set_config.sh  simulation_input_data.txt /data/sml6467/github/Koslicki_lab_metagenomic_analysis/src/git_repo/CAMISIM/defaults/personal_run/templete_mini_config.ini yes


### manually set size from 0.1 to 0.5 (500MB) and set seed to 1234567


### separately move those files into CAMI folder and run simulation
mv new*   /data/sml6467/github/Koslicki_lab_metagenomic_analysis/src/git_repo/CAMISIM/defaults/personal_run/
cd /data/sml6467/github/Koslicki_lab_metagenomic_analysis/src/git_repo/CAMISIM/
conda activate /data/sml6467/github/Koslicki_lab_metagenomic_analysis/src/conda_env/cami_env_py27
python metagenomesimulation.py ./defaults/personal_run/new_mini_config.ini


### input files are stored for future usage.
### mv the "out" folder back here, that's the results


