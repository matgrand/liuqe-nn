# Figure prompts

## tcv_geometry
Change figure tcv geometry
To contextualize where it will be in the narration read /home/mg/repos/phd_thesis/ch2.md.
Figure idea: the goal of this figure is to give an idea of what liuqe is doing, its native grid, an example of flux map with the LCFS, see an example of ctrl points location to understand what the overall goal of the surrogate is, also show the location of the poloidal field coils, magnetic measurements Bm and flux loops
Modifications: (1) change current titles/descriptions + (2) add 1 more panel:
1. current (a) -> Liuqe grid. curr (b) Poloidal flux \psi and LCFS. curr (c) Shape control example points. Also change the description below to something coeherent with the narration, like do not mention the PID. something like this is like this is liuqe, ctrl pts, mag meas and coils positions
2. add a new panel (-> 4 horiz panels) with the mag sensors (find the correct name) Bm, flux loops Ff, and Pfc positions Ia. Read the guides to understand how to do it, I suspect that some data might be missing and I need to retrieve from lac. If that's the case, prepare a .mat script that save the necessary data to a mat file, that I will run in the lac server, I will then move the data in the same place in the directory u saved it and you will be able to load it from the tcv_geometry.py script. is it clear? (additional note: the guides might not be enough to understand how to plot those posisitions, if that's the case the meq/liuqe repo is at /home/mg/repos/meq (also /home/mg/repos/hybrid) might e useful)
Other mods:
- drop the colorbar and the assciated Poloidal flux label
- for the flux map, instead of showing all the samples, only show points inside the vessel (see how gen_dissertation_figures/net_example.py does it)
- do not change data sources from the current panels, just pretend it's from liuqe. Complete the top level description accordingly. It should match the narrative


## TCV shapes
The goal is to recreate figure gen_dissertation_figures/fig_useful_material/tcv_shape_flexibility_draft.png, we need to recreate it exactly, the only difference will be the colortyle and layout (2x6 -> 3x4).
Additional notes:
- read the guides
- I want you to create 2 scripts:
    1. a matlab script that I will run myself in the tcv server that will load and save the data in mat files.
    2. a python script that creates the final figure using the mat files, using the dissertation format, be inspired by gen_dissertation_figures/net_example.py for the layout (3x4 tough)
- I want the same shots and instants.
- call it tcv_shape_flexibility_draft.py
- the shot numbers and instants are (in order): [(9849,0.5), (6010,0.8), (5650,0.8), (8856,0.29), (11962,1.0), (11368,0.65), (18548,1.5), (6442,0.5), (10159,0.01), (19373,0.42), (11928,0.67), (36151,0.457)]
- in the m file, try to get the already saved eq first, re-run the eq as a backup. And do not interrupt the script if one of the eq doesnt load (there are some very old shots, Im not sure it will be very easy to load them)