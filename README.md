# comsol4 multipole calculations

This repository houses the code suitable for extracting multipole components from the electromagnetic responds of periodic structures under normal-incidence illumination, when simulated in COMSOL4. This code specifically was used to obtain mutlipole results for:

[Krishnamoorthy et. al. "Infrared dielectric metamaterials from high refractive index chalcogenides", Nature Comm. 11, 1692 (2020)](https://www.nature.com/articles/s41467-020-15444-0)

The code is not necessarily pretty, but it has been vetted repeatedly to:

* Work with a range of COMSOL geometries, via the COMSOLs Matlab API
* Compute the multipoles correctly
* Store them correctly and load them correctly
* Reproduce electromagnetic fields radiated by these multipoles, correctly.

## Files

* `Mult_Extract.m` - the top-level file that would be run to launch the extraction of the multipoles given COMSOL file. Includes a loop over the wavelengths, to get the spectrum
* `Get_Bi2Te3_VS_loop_z_offset_good_model.m` - launched by `Mult_Extract.m` to prepare the COMSOL geometry, which will then be used to compute multipole integrals
* `CSL4v4_ArbMed_ComputeAllMultipoleData_Callback_July19.m` - is called by `Mult_Extract.m` to actually run all the multipole integrals and store them
* `CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17.m` - called by `CSL4v4_ArbMed_ComputeAllMultipoleData_Callback_July19.m` to parse the simplistic integrand expressions into longer ones that would be accepted by COMSOL

* `Analyze_Mult.m` - loads the multipole calculations, and plots them
* `CSL4v4_RecoverMultipolesFromDataMat_Sep17.m` - called by `Analyze_Mult.m`. Does the actual loading. Simple code, valuable due to vetting
* `CSL4v4_ArbMed_GetAnalyticMuPolFields_Sep17.m` - once multupole values are loded, use their values to compute the electric field that would be emitted by and array of such multipoles
