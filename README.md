# HPLSA
High-Pressure Linear Stability Analysis (HPLSA)

HPLSA was designed to serve as a flexible tool to develop linear stability methods for ideal-gas and real-gas thermodynamic frameworks, including wrappers for CoolProp and RefProp libraries. The code is equipped with several comments for readability. 
Additionally, the repository includes a description of the code, instructions and a guide for users.
Therefore, the linear stability modal and non-modal analysis can be fully reproduced by the interested reader.
Although the code is suitable for Poiseuille flows, it can be easily adapted to other wall-bounded cases, such as Couette flow, by adjusting the initial and boundary conditions.
Furthermore, the operator is built for temporal eigenproblems (prescribing streamwise and spanwise wavenumbers and solving the eigenproblem for the angular frequency and growth rate of the perturbation), but it is prepared to be expanded to solve spatial problems also, which is typical in the case of external flows.
This solver requires the previous installation of the high-pressure compressible flow solver (HPCFS) available at https://github.com/marc-bernades/HPCFS, which embeds some functionalities and thermodynamic models necessary for the complete usage of HPLSA.

Main files contained within the respository:

1. Main executable scripts

Main_LST: modal-analysis linear stability

Main_LST_TG: transient growth (non-modal) linear stability

Main_KineticEnergyBudget: performs 2D kinetic energy budget and contains the functionalities to derive the plots to (i) compare different cases, (ii) assess the energy budget and (iii) quantify compressibility effects and vorticity transport contributors.

2. Main postprocess

Main_Postprocess: Compile the modal-analysis results
- Perturbation profiles, stability diagrams, spectra, eigenvectors

Main_Posprocess_TG: Compile the non-modal-analysis results
- Perturbation profiles, optimum patterns, input/output responses, growth rate maps

2. Functions
- LST_functions: baseflow, LST opertor build, LST code with boundary conditions for Poiseuille flow and transient growth
- Plot_functions: baseflow, energy budget, growth rates, optimal perturbation, spectrum, stability diagram and transient growth (called via Main_Postprocess scripts)
- Utilities: bulks, numeric-based and thermodynamic-based Jacobians, central derivatives, normalizations, scalings, discretization
