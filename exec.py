# import Criv pre-processing script
import prep_file

# -- Load parameters from user_param.txt
param = prep_file.read_input_file()

# -- Plot model geometry
# output saved in CrivApp/output/model.png
prep_file.build_model(param)

# -- Calculate Xfar
# output saved in CrivApp/output/plot_xfar.png
prep_file.compute_Xfar(param)

# -- Compute CRIV
prep_file.compute_CRIV(param)

# -- Plot regression line
# output saved in CrivApp/output/regression_line.png
prep_file.plot_CRIV()

# -- Get CRIV and R2 from regression line
# output saved in  CrivApp/output/CRIV_value.txt and CrivApp/output/R2_value.txt
prep_file.write_CRIV(param)

# -- Get CRIV distribution from parameter distribution
prep_file.CRIV_distrib(param)

#plot CRIV and parameter distribution
prep_file.CRIV_dist_plot(param)

