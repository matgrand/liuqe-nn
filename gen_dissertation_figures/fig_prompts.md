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