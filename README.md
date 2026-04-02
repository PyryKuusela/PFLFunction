The following package is based on the following work: [https://arxiv.org/abs/2604.01191](https://arxiv.org/abs/2604.01191)

## Usage
Run the following commands in a terminal window to prepare the virtual environment.
1. Clone the repo: `git clone https://github.com/PyryKuusela/PFLFunction.git`
2. With Conda, create a virtual environment: `conda env create -f environment.yml`
3. Activate the environment: `conda activate pflfunction`
4. Install the package while in the root directory: `pip install .`
   (For an editable install instead, i.e. for people who wish to develop the code, one can instead run: `pip install -e .`)

In the scripts folder, make the shell files executables:
1. In a terminal window, run: `chmod +x parallel_zeta.sh`

In the scripts folder, activate the conda environment and run the script or example notebook:
1. Usage of the script: `./parallel_zeta.sh operator_name primes scaling label acc nadd`

## To-Do
The following are works-in-progress with regards to upgrading the code.
1. Implement the $$\text{ord}_p(\text{W}(\varphi^p)^{-1})$$ into the default value of $$A$$ when no values for $$\texttt{acc}$$ and $$\texttt{nadd}$$ are given.
2. Upgrade the form of the available denominators for $$\text{U}$$-matrix computations.
3. Allow for elliptic curves and one-parameter fourfolds, i.e. $$b=2$$ and $$b=5$$.
4. Implement our own polynomial code to minimize dependencies on e.g. $$\texttt{SymPy}$$.
5. Run through all of the code with a profiler.
6. Work towards multi-parameter threefolds, i.e. Hodge type $$(1,m,m,1)$$.
