import prep_file

param = prep_file.read_input_file()

#Plot model geometry
prep_file.build_model(param)

#Compute CRIV
prep_file.compute_CRIV(param)

#Plot regression line
prep_file.plot_CRIV()

#Get CRIV and R2 from regression line
prep_file.write_CRIV(param)

#calculate Xfar
prep_file.compute_Xfar(param)

