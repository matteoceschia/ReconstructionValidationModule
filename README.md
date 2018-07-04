# ValidationModule readme


Cheryl Patrick (UCL)
Last updated July 4th 2018

ValidationModule  is a Falaise pipeline module to process a selection of basic validation data and output a ROOT ntuple file, to be used for standard plots and metrics to validate data quality and basic reconstruction. The code is based on the SensitivityModule, but outputs completely different data. It should serve as an example for further validation modules; this one looks at reconstructed quantities like the number of reconstructed hits in an event or the reconstructed energy in a calorimeter.

## Files:

- ValidationModule.cpp
- ValidationModule.h
- CMakeLists.txt
- ValidationModuleExample.conf.in


## Description

Add to an flreconstruct pipeline to generate a ROOT ntuple file with some pertinent branches. To build it, do

``` console
$ ls
CMakeLists.txt                   SensitivityModule.h
README.md                        SensitivityModuleExample.conf.in
SensitivityModule.cpp
$ mkdir build
$ cd build
$ cmake -DCMAKE_PREFIX_PATH=<pathtoyourfalaiseinstall> ..
...
$ make
...
... If you are developing the module, you can test it by doing ...
$ make test
```

Note: if you get a QT5 error, you may need to specify the QT5 path when you run the cmake line, as given by `brew --prefix qt5-base`. For example, you can run:
``` console
$ cmake -DCMAKE_PREFIX_PATH="$(brew --prefix qt5-base);$(brew --prefix)" ..
``` 

The build will create the `libValidationModule` shared library plus the example `flreconstruct` pipeline
script `ValidationModuleExample.conf`. Assuming that you have an `input.brio` file that contains
the `SD`, `CD`, `TCD`, `TTD` and `PTD` banks from the full reconstruction pipeline of `flreconstruct`
(up to and including gamma clustering), this can be run as:

``` console
... Assume we run in the build dir ...
$ ls
CMakeCache.txt                SensitivityModuleExample.conf
CMakeFiles                    cmake_install.cmake
Makefile
...
$ flreconstruct -i /path/to/input.brio -p ValidationModuleExample.conf
...
$ ls
CMakeCache.txt                ValidationModuleExample.conf
CMakeFiles                    cmake_install.cmake
Makefile                      Validation.root
```

The output file will by default be called `Validation.root` so don’t run it multiple times concurrently in the same directory
or you will overwrite the previous file! Use the falaise flreconstruct pipeline instructions to see how to integrate this module in your pipeline.

There is now the option to configure the output filename in the module configuration file.
The final two lines of the configuration file must read:

[name="processing" type="ValidationModule"]
filename_out : string[1] = "my_filename.root"


## Branches in the tuple -

**h_calorimeter_hit_count** : Total number of reconstructed calorimeter hits

**h_calo_hits_over_threshold** : Total number of reconstructed calorimeter hits above the 50keV trigger threshold

**h_cluster_count** : Total number of reconstructed clusters

**h_track_count** : Total number of reconstructed tracks (not including gamma 'tracks')

**h_negative_track_count** : Tracks with negative charge (from curvature, assuming coming from foil)

**h_positive_track_count** : Tracks with positive charge (from curvature, assuming coming from foil)

**h_associated_track_count** : Number of tracks with associated calo hit

**h_geiger_hit_count** : Total number of tracker hits

**v_all_track_hit_counts** : Vector of the number of tracker hits in each individual track

**h_total_calorimeter_energy** : Summed energy in all calorimeters (MeV)

**h_calo_energy_over_threshold** : Summed energy in all calorimeters (MeV) of calorimeter hits above the 50keV trigger threshold

**h_unassociated_calorimeter_energy** : Summed energy in calorimeters not associated to tracks (considered to be gammas) (MeV)

**h_unassociated_energy_over_threshold** : Summed energy in calorimeters for hits over 50keV that are not associated to tracks (considered to be gammas) (MeV)

**h_associated_calorimeter_energy** : Summed calo energy associated to tracks (MeV)

**h_associated_energy_over_threshold** : Summed calo energy for hits over 50keV that are associated to tracks (MeV)

**h_calo_hit_time_separation** Time in ns between first and last calorimeter hits

**t_cell_hit_count** Vector of tracker cells that have a geiger hit. Encoded using the EncodeLocation function

**tm_average_drift_radius.t_cell_hit_count** Vector of drift radii for each Geiger hit in mm. The order of the hits corresponds to the order of cell locations in t_cell_hit_count

**c_calorimeter_hit_map** Vector with all calorimeter hit locations, encoded using EncodeLocation

**cm_average_calorimeter_energy.c_calorimeter_hit_map** Vector of energies of calorimeter hits in MeV. The order of the hits corresponds to the order of hit locations in c_calorimeter_hit_map
