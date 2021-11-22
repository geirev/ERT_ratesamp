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

## 4. Compile the `ERT_ratesamp` library

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





