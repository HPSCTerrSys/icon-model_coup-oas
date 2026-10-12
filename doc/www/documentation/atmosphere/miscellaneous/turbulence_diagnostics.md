```{eval-rst}
:orphan:
```

(ref_atmosphere_turbdiag)=
# Turbulence Diagnostics

To enable a more through evaluation of turbulent processes in ICON, the option exists to calculate and write **additional turbulence parameters** via the namelist parameter `ldiagnose_gstke` in the `io_nml`. When this option is enabled, two additional diagnostic 3D parameters are defined and available for output. In addition to the (subgrid-scale, parameterised) turbulent kinetic energy (SGS TKE, `tke`), which is a prognostic variable for `inwp_turb=1`(the Raschendorfer Turbdiff scheme), the eddy dissipation rate (`edr`) and turbulent length scale (`tur_len_scale`) can be added to the output namelist. This writes the instantaneous value of these parameters.

The diagnostic additionally calculates a **grid-scale contribution to the TKE** from the temporal fluctuation of the resolved winds (`gs_tke`) over a pre-specified time period `gstke_interval`, also set in the `io_nml`. A temporal average of the SGS TKE and the EDR over the same specified period is also calculated and available as output variables with the names `sgs_tke` and `avg_edr`.

The sum of `gs_tke` and `sgs_tke` is then a more accurate approximation of the actual TKE of the atmosphere. As grid resolution increases, it becomes more relevant to include the contribution from the grid-scale, resolved TKE.
