

# installation instuctions for setting up python environment:
cd ~/cajalcourse
bash Miniconda3-latest-Linux-x86_64.sh
conda create -n cajal python=3.7 anaconda
conda activate cajal
conda install -c conda-forge jupyterlab 
conda install -c conda-forge nilearn
pip install --user datalad 
sudo apt-get install git-annex


# set environment:
source activate cajal
# or
conda activate cajal

# course directory:
cd ~/cajalcourse

# download dataset:
datalad install -g https://github.com/OpenNeuroDatasets/ds000245.git

# preprocess dataset:
docker run -ti --rm \
  -v ~/cajalcourse/ds000245:/data:ro \
  -v ~/cajalcourse/ds000245_preproc:/out \
  -v ~/cajalcourse:/fs poldracklab/fmriprep:latest \
  --fs-no-reconall \
  --output-spaces MNI152NLin2009cAsym \
  --force-no-bbr \
  --dummy-scans 4 \
  --sloppy \
  --notrack \
  --n_cpus 2 \
  --mem_mb 16000 \
  --participant_label sub-CTL01 sub-CTL02 \
  --fs-license-file /fs/license.txt \
  /data /out participant

#  docker run -ti --rm -v /home/padawan/cajalcourse/ds000245:/data:ro -v /home/padawan/cajalcourse/ds000245_preproc:/out -v /home/padawan/cajalcourse:/fs poldracklab/fmriprep:latest --fs-no-reconall --output-spaces MNI152NLin2009cAsym --force-no-bbr --dummy-scans 4 --sloppy --notrack --n_cpus 2 --mem_mb 16000 --participant_label sub-CTL01 sub-CTL02 --fs-license-file /fs/license.txt /data /out participant


# start python notebook:
jupyter lab

# further instructions available there 


