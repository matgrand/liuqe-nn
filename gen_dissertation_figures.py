# %% Generate dissertation figures
import gen_dissertation_figures as figures

# figures.configure(SHOW_PLOTS=True) # uncomment to show plots
# figures.configure(WAIT_FOR_PLOTS=True) # uncomment to not wait for plots to be closed before continuing

# Comment out completed figures while working on the remaining ones.

# figures.error_heatmaps.generate()
# figures.color_palette.generate()

figures.tcv_geometry.generate()
figures.swish_act.generate()
figures.sched_losses.generate()
figures.mag_meas.generate()
figures.raw_ds_example.generate()
figures.final_ds_example.generate()
figures.importance_sampling.generate()
figures.net_example.generate()
figures.comparison.generate()

figures.pid_shape_control_simulation.generate()
figures.mpc_shape_control_simulation.generate()
figures.fmag_control_simulation.generate()

figures.compile_inspection_pdf()

import matplotlib.pyplot as plt
plt.show()
