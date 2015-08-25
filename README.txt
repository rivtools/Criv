this app work with Python language.


To install Python 2.7 : https://www.python.org/

CRIV_app installation: 

Download CRIV_app as zip archive in your computer

You need to install two programm in the ./bin directory of the CRIV_app directory depending on your OS:

SUTRA : http://water.usgs.gov/nrp/gwsoftware/sutra.html

Gmsh : http://geuz.org/gmsh/#Download

Launch the app:

You need to inplement your parameters in the ./model/user_param.txt file

To launch the app you need to open a python terminal in the CrivApp/model/ directory

Now you can launch the different tools of CRIV_app by past differents line of the CrivApp/execut.py in the python terminal. 

#Plot model geometry
prep_file.build_model(param)

---> to show the geometry of your problem (model dimensions) 
     in the CrivApp/output/model.png

#Compute CRIV
prep_file.compute_CRIV(param)

---> to compute the CRIV

#Plot regression line
prep_file.plot_CRIV()

---> to see the regression line associate to the correlation between stream-aquifer discharge and head difference between stream and aquifer
     in the  CrivApp/output/regrassion_line.png

#Get CRIV and R2 from regression line
prep_file.write_CRIV(param)

---> obtain files with values of CRIV and regression coeffecient of the regression line
     in the CrivApp/output/CRIV_value.txt &  CrivApp/output/R2_value.txt

#calculate Xfar
prep_file.compute_Xfar(param)

---> Obtain the plot of groundwater horizontal fluxes componants vs distance from the stream to determin Xfar
     in the  CrivApp/output/plot_xfar.png
