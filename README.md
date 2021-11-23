## ERT_ratesamp

This code samples correlated measurement perturbations that we can use to perturb well rates in the schedule file.

---

If you plan to collaborate or contribute anything to the project, use the <a href="#1b-advanced-installation">Advanced Installation</a> option.

## 1a. Basic installation

Create a directory to clone the two repositories:

```bash
git clone git@github.com:geirev/ERT_ratesamp.git
git clone git@github.com:geirev/EnKF_sampling.git
```

After cloning, the directory structure should look like:

```bash
.
├── EnKF_sampling
└── ERT_ratesamp
```

## 1b. Advanced installation

Make a personal github account unless you already have one.
Fork the two repositorys listed above.
Next clone the repositories and set upstream to the original repositories where
you need to replace <userid> with your github userid

```bash
git clone git@github.com:<userid>/ERT_ratesamp.git
pushd ERT_ratesamp
git remote add upstream https://github.com/geirev/ERT_ratesamp
#or, if you have set up git-ssh
#git remote add upstream git://github.com:geirev/ERT_ratesamp
popd

git clone git@github.com:<userid>/EnKF_sampling.git
pushd EnKF_sampling
git remote add upstream https://github.com/geirev/EnKF_sampling
#or, if you have set up git-ssh
#git remote add upstream git://github.com:geirev/EnKF_sampling
popd
```

## 2. Required Packages


```bash
sudo apt-get -y update
sudo apt-get -y install libblas-dev liblapack-dev libfftw3-dev gfortran
sudo apt-get -y install gnuplot  # Needed if you want to use the gnuplot plotting macro
```

## 3. Compile the `EnKF_sampling` library

Navigate to the `lib` folder of the `EnKF_sampling` repository:

```bash
cd EnKF_sampling/lib
```

then compile and place all the `.o` files as well as `libsampling.a` into
the `build` directory of the `ERT_ratesamp` repository using:

```bash
make BUILD=../../ERT_ratesamp/build
```

## 4. Compile `ERT_ratesamp` code

Navigate to the `src` folder of the `ERT_ratesamp` repository:

```bash
cd ERT_ratesamp/src
```

then compile and install the executable in the target directory, defaulting to
`$HOME/bin`:

```bash
make BINDIR=$HOME/bin
```

## 5. To use
run:

```bash
ertsampling
```

The code will use precompiled paths and filenames and output an ensemble of realizations of perturbed wellrates
CONTROL_0 -- CONTROL_N, a file wells.txt with wellnames, and an "obsolete" file EPERT_0 containing all the realizations.

You have to provide the precompiled paths to files
```bash
   obshistfile='/home/geve/Dropbox/Statoil/erterr/observations/obshistnew.txt' # input to ERT
```
Example file (note that this is the only accepted version/format of this file)
```bash
HISTORY_OBSERVATION WOPR:OP_1 { ERROR = 0.05 ; ERROR_MODE = RELMIN; ERROR_MIN =   10.0; };
HISTORY_OBSERVATION WOPR:OP_2 { ERROR = 0.05 ; ERROR_MODE = RELMIN; ERROR_MIN =   10.0; };
HISTORY_OBSERVATION WOPR:OP_3 { ERROR = 0.05 ; ERROR_MODE = RELMIN; ERROR_MIN =   10.0; };
HISTORY_OBSERVATION WOPR:OP_4 { ERROR = 0.05 ; ERROR_MODE = RELMIN; ERROR_MIN =   10.0; };
HISTORY_OBSERVATION WOPR:OP_5 { ERROR = 0.05 ; ERROR_MODE = RELMIN; ERROR_MIN =   10.0; };
HISTORY_OBSERVATION WGPR:OP_1 { ERROR = 0.10 ; ERROR_MODE = RELMIN; ERROR_MIN =  100.0; };
HISTORY_OBSERVATION WGPR:OP_2 { ERROR = 0.10 ; ERROR_MODE = RELMIN; ERROR_MIN =  100.0; };
HISTORY_OBSERVATION WGPR:OP_3 { ERROR = 0.10 ; ERROR_MODE = RELMIN; ERROR_MIN =  100.0; };
HISTORY_OBSERVATION WGPR:OP_4 { ERROR = 0.10 ; ERROR_MODE = RELMIN; ERROR_MIN =  100.0; };
HISTORY_OBSERVATION WGPR:OP_5 { ERROR = 0.10 ; ERROR_MODE = RELMIN; ERROR_MIN =  100.0; };
HISTORY_OBSERVATION WWPR:OP_1 { ERROR = 0.10 ; ERROR_MODE = RELMIN; ERROR_MIN =   10.0; };
HISTORY_OBSERVATION WWPR:OP_2 { ERROR = 0.10 ; ERROR_MODE = RELMIN; ERROR_MIN =   10.0; };
HISTORY_OBSERVATION WWPR:OP_3 { ERROR = 0.10 ; ERROR_MODE = RELMIN; ERROR_MIN =   10.0; };
HISTORY_OBSERVATION WWPR:OP_4 { ERROR = 0.10 ; ERROR_MODE = RELMIN; ERROR_MIN =   10.0; };
HISTORY_OBSERVATION WWPR:OP_5 { ERROR = 0.10 ; ERROR_MODE = RELMIN; ERROR_MIN =   10.0; };
```


A standard schedule file from where we extract wellnames and rates.
```bash
   schedulefile='/home/geve/Dropbox/Statoil/eclipse/include/history.sch'
```
The place to save the CONTROL files with perturbed rates and a copy of the wells.txt file
   controlpath='/home/geve/Dropbox/Statoil/eclipse/Priors/control/'

Note that the wells.txt file may need to be edited before use.



