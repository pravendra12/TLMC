## Installation (HPC-friendly)

> These steps use **Miniconda + conda env** and build scripts (`build.sh`) for both LMC and TLMC.

### 0) Install Miniconda3 (one-time)
```bash
cd $HOME
wget -O Miniconda3-latest-Linux-x86_64.sh \
  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# restart shell (or)
source ~/.bashrc

which conda
conda --version
````

### 1) Create conda env + Python deps

```bash
conda create -n LMC python=3.11 -y
conda activate LMC

conda install -y -c conda-forge ase spglib numpy scipy

python -c "import ase, spglib, numpy, scipy; print('Python deps OK')"
```

### 2) Load compiler / build modules (edit for your cluster)

```bash
module load gcc/13.2.0
module load cmake/3.26.3
module load boost/1.88.0
module load impi/2021.5.1
module list

# optional
module spider gcc
```

### 3) Build LMC (dependency / baseline)

```bash
root="$PWD"
mkdir -p "$root/codeLMC"
cd "$root/codeLMC"

git clone https://github.com/pravendra12/LatticeMonteCarlo-eigen.git
cd LatticeMonteCarlo-eigen

# recommended
bash build.sh

# fallback manual build:
# mkdir -p cmake-build && cd cmake-build
# cmake .. -DCMAKE_BUILD_TYPE=Release \
#   -DCMAKE_C_COMPILER="$(which gcc)" \
#   -DCMAKE_CXX_COMPILER="$(which g++)"
# make -j2
```

### 4) Build TLMC

```bash
cd "$root"
mkdir -p "$root/codeTLMC"
cd "$root/codeTLMC"

git clone https://github.com/pravendra12/TLMC.git
cd TLMC

# recommended
bash build.sh

# fallback manual build:
# mkdir -p cmake-build && cd cmake-build
# cmake .. -DCMAKE_BUILD_TYPE=Release \
#   -DCMAKE_C_COMPILER="$(which gcc)" \
#   -DCMAKE_CXX_COMPILER="$(which g++)"
# make -j2
```