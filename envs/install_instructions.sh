# Disclaimer: 
# Install instructions are truncated to ensure functionality of the dependencies needed to run the BalticFishDB workflow.
# Full installation of the tools after which the environments were names may require further steps.


### align-env
conda env create -f align-env.yaml
# install MAFFT extensions
conda activate align-env
wget https://mafft.cbrc.jp/alignment/software/mafft-7.525-with-extensions-src.tgz
tar -xzvf mafft-7.525-with-extensions-src.tgz
cd mafft-7.525-with-extensions/extensions
# modify the PREFIX line in the Makefile to point to the location of the conda environment
# run 'echo $CONDA_PREFIX' to obtain the required value
make clean
make
make install

### antismash-7.1.0
conda env create -f antismash-7.1.0.yaml

### anvio-8
conda env create -f anvio-8.yaml

### meme-5.5.9
conda env create -f meme-5.5.9.yaml

### orthofinder-2.5.5
conda env create -f orthofinder-2.5.5.yaml

### paprica-env
conda env create -f paprica-env.yaml

### seqkit-2.8.2
conda env create -f seqkit-2.8.2.yaml

### taxonkit-0.20.0
conda env create -f taxonkit-0.20.0.yaml

