```{eval-rst}
:orphan:
```

(ref_atm_nwpmisc)=
# Further Options for ICON with NWP physics

(ref_sstsic_ext)=
## Sea-Surface Temperature and Sea-Ice Fraction

The sea-surface temperature (SST) has to be provided to the atmosphere model. The physics package has a sea-ice model
that provides prognostic ice surface temperatures and is able to melt sea ice but not to create new ice from water.
Thus, both the initial SST and ice fraction have to be provided externally, with optional updates during the
simulation. Information about SST and sea-ice fraction can be obtained from different sources, which are selected by
the {term}`sstice_mode` parameter:

1. SST and sea-ice fraction are read from the analysis. The SST is kept constant. This mode also applies to coupled
   atmo/ocean simulations.
2. SST and sea-ice fraction are read from the analysis. The SST is updated by climatological increments (from the
   extpar file) on a daily basis.
3. SST and sea-ice fraction are updated daily, based on climatological monthly means.
4. SST and sea-ice fraction are updated daily, based on actual monthly means.
6. SST and sea-ice fraction are updated with a user-defined interval.

In modes 3, 4 and 6, SST and ice fraction are read from NetCDF files specified in {term}`sst_td_filename` and
{term}`ci_td_filename`. Modes 3 and 4 interpret the given fields as monthly means. They are associated with the
middle of the month and ICON interpolates linearly between them. The filename parameters support keywords `<path>`,
`<gridfile>`, `<year>`, and `<month>`. These get replaced by the `master_nml:model_base_dir`, the current grid
filename, or the current simulation year or month, respectively. In climatology mode, the year is replaced by the
string `CLIM`, and the twelve months form a repeating cycle.

Mode 6 provides more fine-grained control: the SST and ice fraction fields assume the value stored in the file at
the corresponding time stamp, with linear interpolation in between. The durations between time stamps can be
nonuniform. This allows emulation of mode 4 with arbitrary averaging periods by setting the time stamp to the
middle of the averaging interval. The SST and ice-cover files are a sequence with an interval given by
{term}`sst_file_interval`, anchored on midnight of Jan 1st of the current simulation year. The filename patterns
support all keywords available for extpar files (like domain number `DOM<idom>` or grid resolution `R<nroot0>B<jlev>`)
and the time-stamp components `<year>`, `<month>`, `<day>`, `<hh>`, `<mm>`, and `<ss>`. The NetCDF variables have to be
named `SST` and `SIC`, respectively.


## Glossary of Namelist Parameters

_Operational NWP setting marked by {material-regular}`settings;1em;pst-color-secondary`_

:::{glossary}
sstice_mode
  (`&lnd_nml`) Sea surface temperature / sea ice cover origin ([description see above](ref_sstsic_ext)) 2:{material-regular}`settings;1em;pst-color-secondary`

sst_td_filename
  (`lnd_nml`) Filename for time dependent sea surface temperature data ([description see above](ref_sstsic_ext))

ci_td_filename
  (`&lnd_nml`) Filename for time dependent sea ice cover data ([description see above](ref_sstsic_ext))

sst_file_interval
  (`&lnd_nml`) Interval for `sst_td_filename` and `ci_td_filename` ([description see above](ref_sstsic_ext))
:::
