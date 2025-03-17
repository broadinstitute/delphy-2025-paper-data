# Delphy paper data

Datasets, scripts for analyzing them, and scripts to generate figures.

Although these scripts were last successfully run in early 2025 to prepare, execute and
analyze the Delphy and BEAST2 runs, there is no guarantee that they will run perfectly in
the future and/or on your machine.  We do not intend to maintain these scripts moving
forward.  Instead, the aim of this repo is to be living documentation of all the details
of the runs and plots.

Some files in this repo (notably `*.dphy` run results file) were to big to push to GitHub;
their sizes ranged from 24 MiB to 302 MiB.  They are available upon request to
`pvarilly@broadinstitute.org`.  Suggestions for a long-term repository for these files are
also welcome.

## Usage

Create a Python virtualenv and install the associated requirements:
```
python3 -m venv delphy-env
source ./delphy-env/bin/activate
pip install -r requirements.txt
```

Make symbolic links to the helper binaries you want to use here (you can run `setup_default_links.sh` to look up binaries via `which`, but you'll most likely need to adapt these steps manually):
```
mafft
delphy
delphy_mcc
treeannotator2 (treeannotator in BEAST2 binaries)
loganalyser2 (loganalyser in BEAST2 binaries)
sapling (https://github.com/broadinstitute/sapling)
iqtree2
```

Then enter each of the dataset directory and run the numbered scripts in order.  Warning: the BEAST runs can take a very long time to complete.

In the present repo, all the files created by these runs and analyses have been uploaded.  If you want to regenerate them from scratch, you should delete every *directory* inside each of the dataset directory (e.g., in `sars-cov-2-lemieux`, delete things like `delphy_outputs_a` and `raw`, but not `00_prepare_runs.py` nor `sample_ids.csv`).

The final plots that were composed into the paper figures are in the `plots` directory of each dataset.


## Software versions used

- Delphy Version 1.0 (build 2036, commit 06a7ee4)
- MAFFT Version 7.505 (2022/Apr/10)
- BEAST2 v2.6.2
- BEAGLE commit `3a8d3e6` (Sun Mar 10 2024)
- BEAST X 10.5.0-beta5 Prerelease #1d511b10c2 (for mpox-parker-2025 comparisons)
- Sapling Version 0.1.1 (build 2, commit a0b9da1)
- IQ-TREE 2.3.6
- TreeTime 0.11.4

# Preparing the AWS machine for the benchmarks

- Launch an Ubuntu 24.04 LTS x86-64 instance of type `c7a.2xlarge` (8 vCPUs & 16GB memory) with 8GB gp3 storage
  (for the `sars-cov-2-gisaid-week-by-week` dataset, launch an instance of type `c7a.4xlarge` (16 vCPUs & 32GB memory))
  (for the `sims` dataset, launch an instance of type `c7a.24xlarge` (96 vCPUs & 192GB memory))
- Install BEAST2 (downloaded from BEAST2 releases page: [https://github.com/CompEvol/beast2/releases])
```
  scp -i "~/.ssh/2023-01-29-aws-vs.pem" BEAST.v2.6.2.Linux.tgz ubuntu@ec2-3-78-245-33.eu-central-1.compute.amazonaws.com:.
```
- SSH into the machine, e.g.
```
  ssh -i "~/.ssh/2023-01-29-aws-vs.pem" ubuntu@ec2-3-78-245-33.eu-central-1.compute.amazonaws.com
```
- Upgrade Ubuntu packages
```
    sudo apt update
    sudo apt upgrade  # instance may need to be restarted; do it (`sudo shutdown -r now`) and log back in
```
- Install latest available Java LTS release (21 as of this writing)
```
    sudo apt install openjdk-21-jdk
```
- Check Java works and print version:
```
    java -version

    > openjdk version "21.0.6" 2025-01-21
    > OpenJDK Runtime Environment (build 21.0.6+7-Ubuntu-124.04.1)
    > OpenJDK 64-Bit Server VM (build 21.0.6+7-Ubuntu-124.04.1, mixed mode, sharing)
```
- Download, unpack BEAST X and make a symbolic link to it
```
    wget https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0-beta5/BEAST_X_v10.5.0-beta5.tgz
    tar -xzvf BEAST_X_v10.5.0-beta5.tgz
    ln -s BEASTv10.5.0/bin/beast beast1
```
- Test it
```
    ./beast1 -version
    >
    >  BEAST X v10.5.0-beta5 Prerelease #1d511b10c2, 2002-2024
    >  ...
```
- Unpack BEAST2 and make a symbolic link to it
```
    tar -xvzf BEAST.v2.6.2.Linux.tgz
    ln -s BEASTv10.5.0/bin/beast beast1
```
- Test it
```
    ./beast2 -version
    > v2.6.2
```
- Build and install BEAGLE from source (following these instructions: https://github.com/beagle-dev/beagle-lib/wiki/LinuxInstallInstructions)
```
    # Don't download JDK 11, already got JDK 21 above
    sudo apt-get install cmake build-essential autoconf automake libtool git pkg-config # openjdk-11-jdk
    export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64/  # Need this for CMake to find JDK libs below
    # Also add that same line to ~/.bashrc

    git clone --depth=1 https://github.com/beagle-dev/beagle-lib.git
    cd beagle-lib
    mkdir build
    cd build
    cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF -DCMAKE_INSTALL_PREFIX:PATH=$HOME ..
    make -j install

    export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH  # So that BEAST finds BEAGLE
    # Also add that same line to ~/.bashrc
    
    make test  # Should work
    cd ../..  # Back to home
```
- Ensure that BEAST1 finds BEAGLE
```
    ./beast1 -beagle_info
    > ...
    > --- BEAGLE RESOURCES ---
    > 
    >0 : CPU (x86_64)
    > Flags: PRECISION_SINGLE PRECISION_DOUBLE COMPUTATION_SYNCH EIGEN_REAL EIGEN_COMPLEX SCALING_MANUAL SCALING_AUTO SCALING_ALWAYS SCALERS_RAW SCALERS_LOG VECTOR_SSE VECTOR_NONE THREADING_NONE PROCESSOR_CPU FRAMEWORK_CPU
```
- Ensure that BEAST2 finds BEAGLE
```
    ./beast2 -beagle_info
    > ...
    > --- BEAGLE RESOURCES ---
    > 
    >0 : CPU (x86_64)
    > Flags: PRECISION_SINGLE PRECISION_DOUBLE COMPUTATION_SYNCH EIGEN_REAL EIGEN_COMPLEX SCALING_MANUAL SCALING_AUTO SCALING_ALWAYS SCALERS_RAW SCALERS_LOG VECTOR_SSE VECTOR_NONE THREADING_NONE PROCESSOR_CPU FRAMEWORK_CPU
```
- Now upload any BEAST2.6.2 XML file and run it with the following command (`-threads -1` uses as many threads as there are CPUs, `-beagle` enforces the use of BEAGLE)
```
    cd path/containing/beast/xml
    time ~/beast2 -threads -1 -beagle <beast_input.xml>
``` 
