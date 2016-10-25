## CRIV_app installation: 

Criv is based on the Python language, version 2.7. To install Python 2.7 : [https://www.python.org/](https://www.python.org/)

Download CRIV_app as a zip archive or clone it with git.

Criv depends on Gmsh (mesh creation) and Sutra (flow model). You need to place the two executables in the ./bin directory of the CRIV_app directory depending on your Operating System (Microsoft, macOS, Linux).

The executables and source codes can be found at following locations:
**SUTRA**: 
[http://water.usgs.gov/nrp/gwsoftware/sutra.html] (http://water.usgs.gov/nrp/gwsoftware/sutra.html)

**Gmsh**: 
[http://geuz.org/gmsh/#Download](http://geuz.org/gmsh/#Download)

## Using CRIV

You first need to set your input parameters in the ``./model/user_param.txt`` text file following : 

```
cell_thickness[L] 30 
cell_width[L] 10 
river_width[L] 5 
river_depth[L] 0.5
bank_angle[°] 45 
riverbed_thickness[L] 1 
anisotropy[-] 0.6 
aquifer_hydraulic_conductivity[L/T] 5e-04 
riverbed_hydraulic_conductivity[L/T] 5e-04 

```

[[https://github.com/rivtools/Criv/fig/model_shem.png|alt=octocat]]

To launch the app you need to open a python terminal in the ``CrivApp/model/`` directory

Now you can use the different functions of CRIV_app by copy/pasting the commands from ``CrivApp/execut.py`` in the Python terminal. 

First you need to past:
```python
import prep_file
param = prep_file.read_input_file()
```
**Plot model geometry**: 
past :
```python
prep_file.build_model(param)
```
then, you can take a look to ``CrivApp/output/model.png`` so as to check the validity of the created geometry.

**Calculate Xfar**: 
The grid cell of the horizontal model where a Cauchy-type boundary condition is applied should be large enough to include all the converging / diverging flows associated with the stream. For these reasons, the distance away from the stream where groundwater flow is horizontal, Xfar, must be estimated as fallowing. You first need to set the maximum expected groundwater - surface water head gradient in the ``./model/user_param.txt`` text file :

```
maximum_gw-sw_head_gradient[%] 2

```

Then, past in the Python terminal :

```python
prep_file.compute_Xfar(param)
```

then, you can observed the horizontal component of flow from rivers center in the  ``CrivApp/output/model.png`` so as to adapt your cell width.

**Compute CRIV**: 

After preliminary steps you can coompute CRIV by pasting in the Python terminal :

```python
prep_file.compute_CRIV(param)
```

regression line, CRIV and R2 value can then be obtain in the ``CrivApp/output/`` folder.

**CRIV distribution from parameter distribution**:

The probabilistic CRIV distribution can also be obtained by random sampling from the prior statistical distributions set in ``./model/user_param.txt`` :

```
SD_cell_thickness[L] 1 
SD_cell_width[L] 10
SD_river_width[L] 1
SD_river_depth[L] 1  
SD_bank_angle[°] 1  
SD_riverbed_thickness[L] 3
SD_anisotropy[-] 2 
SD_aquifer_hydraulic_conductivity[L/T] 5
SD_reverbed_hydraulic_conductivity[L/T] 5 

```

by pasting : 

```python
prep_file.CRIV_distrib(param)
```
and 

```python
prep_file.CRIV_dist_plot(param)
```

to plot the probabilistic CRIV distribution in the ``CrivApp/output/`` folder.
