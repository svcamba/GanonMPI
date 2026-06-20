# 🧬 GanonMPI
Parallelization of the ganon tool for taxonomic classification with MPI and C++ native threads.

![C++](https://img.shields.io/badge/C%2B%2B-17%20%7C%2020-blue.svg)
![MPI](https://img.shields.io/badge/MPI-Message_Passing_Interface-orange.svg)
![Status](https://img.shields.io/badge/Status-Bachelor_Thesis_Completed-brightgreen.svg)

**GanonMPI** is a distributed hybrid architecture for the ultrafast taxonomic classification of metagenomic sequences (Big Data). 

This project expands the state-of-the-art classifier *ganon*, enabling it to overcome shared-memory limitations by deploying it on distributed memory infrastructures, HPC clusters, and supercomputers.

---

## 🏗️ Hybrid Architecture: MPI + C++ Threads

To mitigate bottlenecks during the metagenomic classification phase, GanonMPI employs a two-level hierarchical parallel topology:

1. **Inter-node Level (MPI):** Utilizes the *Message Passing Interface* to perform dynamic partitioning of input files (FASTQ files) and consolidate complex data structures (such as hash maps) through collective operations and binary tree reduction across the cluster's network.
2. **Intra-node Level (C++ Threads):** Locally, each MPI process spawns multiple native C++ execution threads (`std::thread`) that communicate via concurrent queues and process sequence batches in parallel, exploiting all available CPU cores on the node.

---

## 🌿 Branching Model

The development of this project is structured into several branches to isolate experimental features and performance benchmarks. *(Navigate to the specific technical branches to explore the parallel engine source code)*:

* [`main`](../../tree/main): Main branch. Hosts the general documentation.
* [`separacionArchivosMmap`](../../tree/separacionArchivosMmap): **Main technical development branch.** Contains the implementation of Memory-Mapped Files, simultaneous read pointer splitting, and the core parallel classification engine.
* [`proceso0Trabaja`](../../tree/proceso0Trabaja): First approach to the hybrid development. The root process exclusively reads the sequence files (FASTQ) and distributes them to the other processes so each can classify its own set of sequences. Finally, the root process gathers the results and writes the final output.
* [`separacionArchivosNoRam`](../../tree/separacionArchivosNoRam): Second technical development branch. Contains the first approach to the distributed read model, where each process reads a specific set of the sequence files independently.

---

## 🚀 Requirements and Deployment

To compile and execute **GanonMPI** on an **HPC cluster**, the following components are required:
* A **POSIX-compatible environment**.
* A **C++ compiler** with **C++20** support (e.g., GCC, Clang).
* An **MPI-3 implementation** (e.g., **OpenMPI** or MPICH).
* Legacy structural libraries from the original ganon tool:
    * **SeqAn3** (biological sequence processing).
    * **Cereal** (data serialization for MPI communication).
* **CMake** build system.
* **Python 3**.
* **OpenSSL** (required to download the indexes).

### Basic Installation Example

First, if you are using a system with **Lmod** environments, load the required modules using the following commands:
```bash
module load gnu12
module load openssl/1.1.1w
module load python/3.12.9
module load cmake/3.25.2
```
Clone and compile the original tool:
```bash
git clone --recurse-submodules https://github.com/pirovc/ganon.git
cd ganon
pip install .
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DVERBOSE_CONFIG=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCONDA=OFF -DLONGREADS=OFF ..
make -j 4
make install
```
Upon completion, the executable binaries will be generated inside the `build/` directory.
Next, install **Raptor**:
```bash
git clone --branch raptor-v3.0.1 --recurse-submodules https://github.com/seqan/raptor
cd raptor
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-std=c++23 -Wno-interference-size" ..
make -j 4
```
Finally, install de parallel GanonMPI version:
```bash
git clone https://github.com/svcamba/GanonMPI.git
cd GanonMPI
# Warning: Replace 'branch_name' with the specific branch you want to use
git checkout branch_name 
cd src/ganon_classify/
make clean && make
cp ganon_classify ../../build/ganon-classify
```

## Basic Database Creation Example
```bash
ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 48 --complete-genomes --db-prefix abfv_rs_cg
```

## Sequence Download Example
For instance, to download the *SRR29606654* sequence from the *ENA* website (https://www.ebi.ac.uk/), retrieve the files and decompress them using:
```bash
gunzip SRR29606654_1.fastq.gz
gunzip SRR29606654_2.fastq.gz
```

## Execution on the Pluton Cluster (Slurm)
Create a `bash` script (e.g., `job.sh`) with the following content:
```bash
#!/bin/bash
#SBATCH --job-name=GanonMPI # Job name displayed in the Slurm queue
#SBATCH --output=SRR37421749_PairedReads_BD_GRANDE_%a.out
#SBATCH -N 8 # Number of nodes
#SBATCH -n 32 # Total number of MPI processes (across all nodes)
#SBATCH -t 70:20:00 # Estimated execution time; Slurm will kill the job if exceeded
#SBATCH --mem=0 # Request all available RAM on the assigned nodes

# MPI environment variables and local installation paths
export MPI_HOME=$HOME/software/openmpi-5.0.3-gnu12
export PATH=$MPI_HOME/bin:$PATH
export LD_LIBRARY_PATH=$MPI_HOME/lib:$LD_LIBRARY_PATH
export MANPATH=$MPI_HOME/share/man:$MANPATH
export PKG_CONFIG_PATH=$MPI_HOME/lib/pkgconfig:$PKG_CONFIG_PATH

echo "--- START MPI BENCHMARK ---"
echo "Nodes: 8"
echo "Threads per MPI process: 7; Total MPI processes: 32"

# Execution
# Using /usr/bin/time -v to benchmark memory and real time
/usr/bin/time -v $MPI_HOME/bin/mpirun -x LD_LIBRARY_PATH ../build/ganon-classify \
    --paired-reads \
    ../../secuencias/SRR37421749/SRR37421749_1.fastq,\
    ../../secuencias/SRR37421749/SRR37421749_2.fastq \
    --ibf ../ganon/abfv_rs_cg.hibf \
    --tax ../ganon/abfv_rs_cg.tax \
    --rel-cutoff 0.75 --rel-filter 0.1 --fpr-query 1e-05 \
    --skip-lca --output-all --hibf \
    --output-prefix SRR37421749_BD_GRANDE \
    --threads 7
```
Finally, submit the job to the cluster's queue system:
```bash
sbatch job.sh
```
