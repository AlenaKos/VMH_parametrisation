# VMH_parametrisation
Air Parcel Model with Ice Physics: Code for Homogeneous Ice Nucleation Study.
This repository contains the full codebase used in the paper "A New Parameterization for Homogeneous Ice Nucleation Driven by Highly Variable Dynamical Forcings."

Includes:

-The core model simulates atmospheric air parcels, incorporating ice microphysics to study homogeneous ice nucleation under different dynamical forcings based on the double-moment bulk scheme [Baumgartner, M., and P. Spichtinger, 2019], [Dolaptchiev, S. I., Spichtinger, P., Baumgartner, M., & Achatz, U. (2023)]. 

-script for reading and filtering the data form ICON-MS-GWaM model. Production of the initial conditions for gravity waves associated preturbations from the global climate model coupled to the GW parametrisation MS-GWaM [Bölöni, G., Y.-H. Kim, S. Borchert, and U. Achatz, 2021].  

-script for calculation of the fitting coefficients based on the ensamle calculations

-script for producing the plots for the draft 


References
Baumgartner, M., and P. Spichtinger, 2019: Homogeneous nucleation from an asymptotic point of view. Theor. Comput. Fluid Dyn., 33, 83–106, https://doi.org/10.1007/s00162-019-00484-0.
Dolaptchiev, S. I., Spichtinger, P., Baumgartner, M., & Achatz, U. 2023. Interactions between gravity waves (GW) and cirrus clouds: Asymptotic modeling of wave-induced ice nucleation. Journal of the Atmospheric Sciences, 80(12), 2861-2879, https://doi.org/10.1175/JAS-D-22-0234.1
Bölöni, G., Y.-H. Kim, S. Borchert, and U. Achatz, 2021: Toward transient subgrid-scale gravity wave representation in atmospheric models. Part I: Propagation model including nondissipative wave–mean-flow interactions. J. Atmos. Sci., 78, 1317–1338, https://doi.org/10.1175/JAS-D-20-0065.1.
