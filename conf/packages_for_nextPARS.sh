
# Add bin to PYTHONPATH for modules made here (genome_annotations.py)
cwd=`pwd`
export PYTHONPATH="$cwd/bin:$PYTHONPATH"


## Install pip for python packages
sudo apt-get update
sudo apt-get -y install python-pip
sudo apt-get -y install default-jre


## To install the required python packages
pip install --user argparse numpy biopython datetime pysam termcolor pandas keras tensorflow dask h5py

