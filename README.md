# SubPEx 1.0

## What is it?

Subpocket explorer (SubPEx) is a tool that uses weighted ensemble (WE) path
sampling, as implemented in WESTPA, to maximize pocket conformational search.
SubPEx's goal is to produce a diverse ensemble of protein conformations for use
in ensemble docking.

As with any WE implementation, SubPEx uses a progress coordinate to focus
computational power on sampling phase space. The available progress coordinates
are:

- composite RMSD (a linear combination of backbone and pocket heavy-atom RMSD)
- pocket heavy atoms RMSD
- backbone RMSD
- Jaccard distance of pocket volumes (jd)

We highly recommend using the composite RMSD progress coordinate. Use of other
coordinates is not officially supported.

## SubPEx use

### Downloading and installing SubPEx

The first step is to download, clone, or copy the repository.

```bash
git clone https://git.durrantlab.pitt.edu/erh91/SubPEx-Erich.git
```

The repository includes a `yaml` file to create a `conda` environment with all
the necessary dependencies. Simply run this command:

```bash
conda env create -f environment.yaml
```

To activate the new `conda` environment, run:

```bash
conda activate westpa
```

Some users may wish to create their own environments or to use an existing
WESTPA environment. If so, install the following packages so SubPEx can
calculate the progress coordinate:

- Westpa
- MDAnalysis
- NumPy
- SciPy
- scikit-learn
- yaml

Note that installing SubPEx into an existing environment is not officially
supported.

### Configuring a SubPEx run

Users should take advantage of our autobuilder (`wizard.py`) to setup their
SubPEx simulations. In some cases, however, users may wish to manually set up
their simulations by editing key SubPEx/WESTPA files. This approach is not
officially supported, but we provide the below instructions for
advanced/adventurous users.

___Link your preliminary, equilibrated simulation___

1. SubPEx assumes you have already run preliminary simulations to equilibrate
   your system. Soft link or copy your preliminary, equilibrated trajectories
   and necessary restart files to the `./reference/` directory. Rename the files
   `mol` with the appropriate extension. (Note that `./reference/` already
   contains the `namd.md.conf` and `amber.prod_npt.in` template files, which
   SubPEX uses to interface with the NAMD and AMBER MD engines, respectively.)
    - If using NAMD, soft link the `.dcd` file of the final equilibration run.
      NAMD requires other files to restart simulations as well. Be sure to soft
      link the `.xsc`, `.coor`, and `.inpcrd` files as well. Remember the
      `.prmtop` file as well.

      ```bash
      ln -s /file/path/to/simulation/my_namd_file.dcd /WEST/ROOT/reference/mol.dcd
      ln -s /file/path/to/simulation/my_namd_file.xsc /WEST/ROOT/reference/mol.xsc
      ln -s /file/path/to/simulation/my_namd_file.coor /WEST/ROOT/reference/mol.coor
      ln -s /file/path/to/simulation/my_namd_file.inpcrd /WEST/ROOT/reference/mol.inpcrd
      ln -s /file/path/to/simulation/my_namd_file.prmtop /WEST/ROOT/reference/mol.prmtop
      ```

    - If using Amber, the filetype that works with the SubPEx algorithm is
      `.nc`. You need to soft link the `.rst` file of the final equilibration
      run as well. Remember the `.prmtop` file as well.

      ```bash
      ln -s /file/path/to/trajectory/my_amber_file.nc /WEST/ROOT/reference/mol.nc
      ln -s /file/path/to/trajectory/my_amber_file.rst /WEST/ROOT/reference/mol.rst
      ln -s /file/path/to/trajectory/my_amber_file.prmtop /WEST/ROOT/reference/mol.prmtop
      ```

2. Extract the last frame of the preliminary, equilibrated trajectory as a `pdb`
   file with your preferred molecular analysis program (e.g., VMD). Soft link
   that to the `./reference/` directory as well, and name the link
   `last_frame.pdb`.

      ```bash
      ln -s /file/path/to/last/frame/my_last_frame.pdb /WEST/ROOT/reference/last_frame.pdb
      ```

___Edit the `west.cfg` file___

1. Edit the following parameters in the `west.cfg` file:
    - the directory portion of the path variables, though the basename itself
      should not change. ___NOTE: Be sure to use full (not relative) paths.___
       - `reference`: the PDB file that will be used in EVERY SINGLE
         progress-coordinate calculation (the last frame of the preliminary,
         equilibrated simulation mentioned above).
       - `selection_file`: path to a text file that will contain the pocket
         selection string (MDAnalysis selection notation). This file will be
         automatically generated in a subsequent step, but specify its future
         path here.
       - `reference_fop`: path to an `xyz` file that will contain the field of
         points needed to calculate the `jd` progress coordinate. This file is
         also useful for visualizing the selected pocket. It will be
         automatically generated in a subsequent step.
       - `west_home`: home directory of the SubPEx run. You'll most likely want
         to use the same directory that contains the `west.cfg` file itself.
       - `topology`: topology file needed for the MD simulations (likely the
         same topology file used in the preliminary, equilibrated simulations).
    - the progress coordinate (`pcoord`) to use.
       - `composite`: composite RMSD (recommended)
       - `prmsd`: pocket heavy atoms RMSD (not officially supported)
       - `bb`: backbone RMSD (not officially supported)
       - `jd`: Jaccard distance (not officially supported)
    - the auxiliary data (`auxdata`) to calculate and save.
       - `composite`: composite RMSD
       - `prmsd`: pocket heavy atoms RMSD*
       - `pvol`: pocket volume (requires `jd` too)
       - `bb`: backbone RMSD
       - `rog`: radius of gyration of the pocket (requires `jd` too)
       - `jd`: Jaccard distance
    - make sure that the WESTPA progress coordinate and auxdata match the SubPEx
      ones (these sections are both found in the `west.cfg` file).
       - The WESTPA progress coordinate is specified at `west -> data ->
         datasets`, `subpex -> pcoord`, and in `adaptive_binning/adaptive.py`
       - The WESTPA auxiliary data is at `west -> executable -> datasets`
       - The SubPEx progress coordinate is at `subpex -> pcoord`
       - The SubPEx auxiliary data is at `subpex -> auxdata`

___Define the pocket to sample___

1. You must define the location of the binding pocket you wish to sample. Find
   the coordinates of the pocket center and radius using the extracted last
   frame.
    - Visual inspection is often useful at this step. You might first create a
      PDB file with a CA dummy atom. Load that together with the extracted last
      frame of the previous step into your preferred visualization software
      (ChimeraX, PyMol, VMD, etc.). Then manually move the dummy atom to the
      pocket center and measure its location. Similarly, use the dummy atom to
      determine the radius from that center required to encompass the pocket of
      interest.
2. Return to the `west.cfg` file and edit the following parameters:
   - `center`: the pocket center
   - `radius`: the pocket radius
   - `resolution`: the distance between adjacent pocket-filling grid points
     (especially important if using the `jd` progress coordinate)
3. Run `python westpa_scripts/get_reference_fop.py west.cfg`. This script will
   generate the files specified by the `selection_file` and `reference_fop`
   parameters in the `west.cfg` file.
4. Visually inspect the pocket field of points (fop) and/or the selection string
   (MDAnalysis selection syntax).
   - Ensure that the points in the fop (`reference_fop`) file entirely fill the
     pocket of interest.
   - Ensure that the residues (`selection_file`) truly line the pocket of
     interest.
   - Note that the popular molecular visualization program VMD can load `xyz`
     files and select residues.
5. After visual inspection, adjust the `west.cfg` file (`center`, `radius`, and
     `resolution` parameters) and re-run the
     `westpa_scripts/get_reference_fop.py` script. Continue to recalculate the
     pocket as needed to fine-tune your pocket.

___Setup the progress coordinate calculations___

1. Update the variables in the `adaptive_binning/adaptive.py` file to indicate
   the number of walkers per bin, the bins' minimum and maximum values, etc.
    - This file controls the adaptive binning scheme that SubPEx uses.
    - A detailed description of each variable is given in the file itself.
2. Change `westpa_scripts/get_pcoord.sh`.
    - This script runs when calculating initial progress coordinates for new
      initial states (istates).
    - A the beginning of the file, modify the line `export ENGINE="NAMD"` to
      match your MD engine (`NAMD` or `AMBER`, in capital letters).
3. Modify the `westpa_scripts/runseg.sh` file.
    - This file runs each WESTPA/SubPEx walker (segment). It creates the needed
      directory, runs the walker simulation, and calculates the progress
      coordinate.
    - At the beginning of the file, modify the line `export ENGINE="NAMD"` to
      match your MD engine (`NAMD` or `AMBER`, in capital letters).

___Setup the environment___

1. Revise the `env.sh` file.
   - The file itself contains further instructions as comments.
   - Among other things, be sure to set the environmental variables required to
     run the NAMD or AMBER executables, as well as the appropriate WORKMANAGER.
   - Setting the appropriate variables may be complicated if using a
     supercomputing center. You may need to consult with an IT administrator.
2. Modify the appropriate MD configuration file in `./reference/` directory
   (`./reference/amber.prod_npt.in` if using AMBER, `./reference/namd.md.conf`
   if using NAMD).
   - Make sure the number of frames saved per simulation equals `pcoordlength`
     minus one (`pcoordlength` is defined in the `adaptive_binning/adaptive.py`
     file). For example:
     - If using AMBER, modify `./reference/amber.prod_npt.in` to make sure
       `nstlim` / `ntwx` = `pcoordlength` - 1.
     - If using NAMD, modify `./reference/namd.md.conf` to make sure `run` /
       `dcdfreq` = `pcoordlength` - 1.
3. Activate the WESTPA conda environment and source the init.sh file.
4. Execute the `. init.sh` file. Note that this will delete any data from
   previous SubPEx runs.

```bash
conda activate westpa
```

___Running SubPEx___

1. To run SubPEx, execute the `./run.sh` file from the command line.
   - You can also run SubPEx on a supercomputing cluster. See the
     `./aux_scripts/run.slurm.sh` for an example submission script for the slurm
     job scheduler. Note that you will likely need to modify the submission
     script for your specific cluster. Please check with your IT administrator
     to troubleshoot any cluster-specific problems.
2. If errors occur during execution, check the `./job_logs` directory. (If there
   is no `./job_logs` directory, that alone will cause WESTPA/SubPEx to fail.)

__Notes__: WESTPA simulations are not easy to set up. You are likely to
encounter errors on your first attempt. To debug, check the `west.log` file and
the output in the `job_logs` directory.

## Important scripts and files and what they do

- __env.sh__ sets up some environmental variables.
- __init.sh__ initializes the run. Creates the basis and initial states and
  calculates progress coordinates for those, using the `bstate.py` script.
- __gen_istate.sh__ makes `istates` directory.
- __get_pcoord.sh__ copies and links necessary files. Then calls on `bstate.py`
  script. This script gets called by `init.sh`.
- __run.sh__ starts the run.
- __runseg.sh__ WESTPA runs this script for each trajectory segment. The script
  has three jobs:
    1. Link the necessary files for MD simulations
    2. Run the MD simulation
    3. Calculate the progress coordinate using the `pcoord.py` script.
- __west.cfg__ the file containing the run's configuration.
- __west.h5__ contains all the results of the WE run.
- __get_reference_fop.py__ calculates the initial field of points for JD
  progress coordinate and creates a selection string for MDAnalysis.
- __pcoord_istate.py__ calculates the progress coordinate for the initial
  states.
- __pcoord.py__ calculates the progress coordinate for the production run.
