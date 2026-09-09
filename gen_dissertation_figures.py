# %% Generate dissertation figures
import gen_dissertation_figures as figures

figures.configure(SHOW_PLOTS=True) # uncomment to show plots 
# figures.configure(WAIT_FOR_PLOTS=True) # uncomment to not wait for plots to be closed before continuing

# Comment out completed figures while working on the remaining ones.
# figures.color_palette.generate()
# figures.losses.generate()
# figures.error_heatmaps.generate()
# figures.control_simulation.generate()
figures.pid_shape_control_simulation.generate()
# figures.mpc_shape_control_simulation.generate()
# figures.fmag_control_simulation.generate()
# figures.tcv_geometry.generate()
# figures.net_example1.generate()
# figures.net_example2.generate()
# figures.net_example3.generate()
# figures.comparison.generate()

figures.pid_shape_control_simulation_v2.generate()
# figures.mpc_shape_control_simulation_v2.generate()
# figures.fmag_control_simulation_v2.generate()


# figures.tcv_geometry.generate()

import matplotlib.pyplot as plt
plt.show()
