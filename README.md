# SensitivityModule readme


Cheryl Patrick (UCL)
Last updated March 6, 2018

ValidationModule  is a Falaise pipeline module to process a selection of basic validation data and output a ROOT ntuple file, to be used for standard plots and metrics to validate data quality and basic reconstruction. The code is based on the SensitivityModule, but outputs completely different data


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


## Output tuple structure -

**calorimeter_hit_count** : Total number of calorimeter hits

**cluster_count** : Total number of clusters

**track_count** : Total number of tracks (not including gamma 'tracks')

**negative_track_count** : Tracks with negative charge (from curvature, assuming coming from foil)

**positive_track_count** : Tracks with positive charge (from curvature, assuming coming from foil)

**associated_track_count** : Number of tracks with associated calo hit

**geiger_hit_count** : Total number of tracker hits

**all_track_hit_counts** : Vector of the number of tracker hits in each individual track

**total_calorimeter_energy** : Summed energy in all calorimeters (MeV)

**unassociated_calorimeter_energy** : Summed energy in calorimeters not associated to tracks (considered to be gammas) (MeV)

**associated_calorimeter_energy** : Summed calo energy associated to tracks (MeV)

**calo_hit_time_separation** Time in ns between first and last calorimeter hits

