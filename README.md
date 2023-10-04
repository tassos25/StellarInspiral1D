# StellarInspiral1D

This is an extension for MESA to simulate an in-spiraling object in a star. Hence, you need to have an installed version of MESA (we recommend to use the MESA release r23.05.1).

Its main application is meant for simulating a common envelope phase where the envelope of one star engulfs its companion. This companion can be any kind of object. While for the in-spiral it will be treated as a simple object.

In its current version are three steps:

1. [Pre-evolution](#step-1:-pre-evolution)
2. [Relaxation](#step-2:-relaxation)
3. [In-spiral](#step-3:-in-spiral)

Those steps are meant to be run consecutively. To do so there are shell scripts which do this:

* rn_HMS-CO.sh
* rn_HMS-HMS.sh
* rn_HMS-sinlge.sh

They differ only by using a different first step.

To run a full simulation you simple run one of those scripts, e.g. run `./rn_HMS-CO.sh` in a shell/terminal opened in the main directory of this extension.

Those will move to the different subdirectories and call the scripts in there. Each subdirectory for steps contain the common MESA scripts `clean`, `mk`, `rn`, `re`, and a make directory with the `makefile`. The make files are adopted to  compile the needed code for each step. Beside this, there is another script in each subdirectory for steps `cpy` it takes care of moving the finial models to a common place, see [models section](#models). The files starting with the key word `inlist` are described in the [inlists section](#inlists).

## Step 1: Pre-evolution

Currently there are three options for the evolution until the onset of the engulfment. Those are based on the framework use in POSYDON. Each case has its own subdirectory starting with the key word `step1_`.

### HMS-CO

This starts with a zero-age main sequence hydrogen rich star and a compact object as companion. The companion is treated as a point mass in a binary with the stellar component. Hence, it uses the binary module of MESA to generate the stellar model at the onset of the engulfment.

### HMS-HMS

This starts with two zero-age main sequence hydrogen rich stars in a binary. Here, both components are simulates as stars until a condition for an unstable engulfment of both stellar components is reached.

### HMS-single

This starts with a single zero-age main sequence hydrogen rich star and evolves it as a single star. This is aimed for use case, where the in-spiraling component has a negligible effect on the stellar evolution until the engulfment.

### other cases

If there are other use cases, where a different progenitor evolution is needed a new step 1 can be introduced.

## Step 2: Relaxation

This step is aimed to make a smooth transition to the in-spiraling phase. It places the in-spiraling object into the star and allow the stellar profile to adjust. So it evolves the star a bit, while keeping the new object at a given position.

This phase acts in the subdirectory `step2_RelaxModel` and used only the star module of MESA together with the in-spiral extension.

## Step 3: In-spiral

Here we enable the drag to act on the in-spiraling object. This drag changes the orbit of the in-spiraling object and injects energy in the stellar envelope and imparts a torque.

This phase acts in the subdirectory `step3_inspiral` and used only the star module of MESA together with the in-spiral extension.

## models

There is a subdirectory `models`. In there are the models sorted, which are needed at the beginning of each step. Hence, it initially contains a ZAMS file only. When running after each step the finial model of this step gets copied in here to serve as the input for the next step. An already existing model is saved by extending `.old` to it. Hence, there is a saved copy of the last run (only).

## inlists

There is a subdirectory `inlists`. It contains the common inlists:

* `inlist_star_common`: this is the common in-list loaded for all star modules. It contains all the parameters, which get (potentially) changed from the MESA defaults. It although contains all the stuff, which might get changed for individual steps only. This should serve as an overview, hence it shows the set value, a short description, a statement, what the MESA default is and a comment, why it got changed.
* `inlist_binary_common`: this is similar to the previous one, but for the binary module of MESA.
* `inlist_x_ctrls_POSYDON`: this is a shared in-list for all the extra value needed for the steps which are based on POSYDON modifications in the `run_binary_extras` or `run_star_extras`, see [common_code section](#POSYDON)
* `inlist_x_ctrls_CE`: this is a shared in-list for the common envelope part. It contains the extra values needed for this code part, see [common_code section](#inspiral).
* `history_columns_common.list`: this is a shard list for columns writen to the star's history files.
* `profile_columns_common.list`: this is a shard list for columns writen to the star's profile files.

Each step contains main in-list(s) to take care of the stacking of in-lists, those are called `inlist` (, `inlist1`, `inlist2`). Additionally, there are in-lists specific for each step and MESA module are therefore called `inlist_#module_step?` (`#module` is `binary` or `star`; `?` is `1`, `2`, or `3`). Here you can find the step specific changes. In the case there is more than one star, there are last level in-lists `inlist_star?` (`?` is `1`, `2`) which only do the naming of files and directories to differentiate the two stars.

## common_code

There is a subdirectory `common_code`. It contains different bunches of code and data.

### inspiral

In this subdirectory is the code needed for the in-spiraling steps:

* `CE_run_star_extras`: this file contains the common functions usually in the `run_star_extras`, which are common for all the in-spiraling steps
* `CE_adjust_mdot`: tbw
* `CE_after_struct_burn_mix`: tbw
* `CE_before_struct_burn_mix`: tbw
* `CE_energy`: tbw
* `CE_orbit`: tbw
* `CE_timestep`: tbw
* `CE_torque`: tbw

### POSYDON

* `POSYDON_run_star_extras`: this file contains the common functions usually in the `run_star_extras`, which are common for all the steps based on POSYDON
* `POSYDON_run_binary_extras`: this file contains the common functions usually in the `run_binary_extras`, which are common for all the steps based on POSYDON

### ionization

This is the required code from the ionization branch of MESA which is not part of the MESA revision r23.05.1, while needed for this extension.

### ionization_data

This is the required data needed for the ionization code. It is copied from MESA revision r9793, which contained this the ionization code and this data set in `${MESA_DIR}/data`.

## FAQ

tbw
