```{eval-rst}
:orphan:
```

(ref_output_coupling)=
# Output Coupling

The **output coupling** feature provides a generic interface for connecting ICON to external output components using the YAC (Yet Another Coupler) library. This allows users to develop custom output components that run independently of ICON, enabling flexible post-processing, analysis, and visualization workflows.

## Key Features

- **Generic Interface**: Works with any YAC-compatible output component
- **Flexible Grid Interpolation**: YAC provides various interpolation methods
- **Temporal Interpolation**: YAC handles temporal coupling automatically
- **Metadata**: Transfers CF conventions and GRIB2 metadata defined in ICON to the output component
- **3D Field Support**: Handles multi-level atmospheric and oceanic data with collection sizes
- **Fractional Masking**: Optional validity masks for handling ocean data

## Configuration

### Enabling Output Coupling

ICON must be compiled with YAC support enabled:

```bash
./configure --enable-coupling
```

Add to your ICON namelist (atmosphere or ocean):

```fortran
&coupling_mode_nml
  coupled_to_output = .TRUE.
/
```

### Running in MPI MPMD Mode

The output coupling requires running ICON and the output component together in **MPI MPMD (Multiple Program Multiple Data) mode**. This allows both executables to run simultaneously and communicate via YAC.

#### Basic MPMD Execution

The general syntax for running ICON with an output component is:

```bash
mpirun -n <N_ICON> <icon_executable> : \
       -n <N_OUTPUT> <output_component> <output_args>
```

Where:
- `<N_ICON>`: Number of MPI ranks for ICON
- `<icon_executable>`: Path to ICON executable (e.g., `icon`)
- `<N_OUTPUT>`: Number of MPI ranks for the output component (typically 1)
- `<output_component>`: Path to output component executable (e.g., `python simple_output.py`)
- `<output_args>`: Arguments for the output component

#### Example: ICON with simple_output.py

```bash
mpirun -n 4 ./icon : \
       -n 1 python run/checksuite.clm/simple_output.py output.nc \
                  atm,icon_grid,temp \
                  atm,icon_grid,u \
                  atm,icon_grid,v \
                  --bounds -180 180 -90 90 \
                  --res 360 180 \
                  --output-interval PT1H
```

This runs ICON on 4 MPI ranks and the output component on 1 rank.

#### MPI Implementation Specifics

Different MPI implementations have slightly different MPMD syntax:

**OpenMPI:**
```bash
mpirun -n 4 ./icon : -n 1 python simple_output.py [args]
```

**SLURM with srun:**
```bash
srun --multi-prog mpmd.conf
```

Where `mpmd.conf` contains:
```
0-3   ./icon
4     python simple_output.py output.nc atm,icon_grid,temp [...]
```

### Component Naming Convention

:::{admonition} Important: Component Naming
:class: admonition-icontheme

{material-regular}`warning;2em;pst-color-secondary` When output coupling is enabled, ICON automatically registers an additional YAC component with an `_output` suffix appended to the main component name (defined in the `master.nml`). The component without the suffix is used for the physical coupling. For example:

- If the main ICON atmosphere component is named `atm`, the output component will be `atm_output`
- If the main ICON ocean component is named `oce`, the output component will be `oce_output`

**All fields are exposed on this `_output` component**, not on the main component. Your output component must reference the correct component name when defining couples and accessing fields.
:::

For example, in `simple_output.py`:
```bash
python simple_output.py output.nc \
  atm_output,icon_grid,temp \
  atm_output,icon_grid,u \
  atm_output,icon_grid,v
```

Note the use of `atm_output` instead of `atm`.

:::{admonition} Accessing Diagnostitc Variables
:class: admonition-icontheme

{material-regular}`warning;2em;pst-color-secondary` Some diagnositic variables are only evaluated if an `output_nml` exists that uses this variable.
If you want to access these variables over the output_coupling, make sure that such an `output_nml` exists in your configuration, otherwise the fields would provide invalid data.
See for example [](ref_output_coupling_vertical_grid_type).

:::

### Field Selection Criteria

Variables in the ICON variable list (`add_var`) are automatically exposed if they meet ALL of the following criteria:

1. **Output flag set**: `info%loutput = .TRUE.`
2. **Not a container**: Container variables (tracer collections) are excluded
3. **Correct data type**: Only `REAL_T` (floating-point) variables
4. **Valid horizontal grid**: Either
   - `GRID_UNSTRUCTURED_CELL` (cell-centered data)
   - `GRID_UNSTRUCTURED_VERT` (vertex-centered data)
5. **Correct dimensions**:
   - First dimension must be `nproma` (blocking factor)
   - For 2D variables: 2nd dimension must be `≥ nblks`
   - For 3D variables: 3rd dimension must be `≥ nblks`

(ref_output_coupling_vertical_grid_type)=
### Vertical Grid Types

The coupling supports multiple vertical coordinate systems, identified by prefixes:

- **(no prefix)**: Model levels (`level_type_ml`)
- `pl::`: Pressure levels (`level_type_pl`)
- `hl::`: Height levels (`level_type_hl`)
- `il::`: Isentropic levels (`level_type_il`)

Example: A temperature field on pressure levels would be exposed as `pl::temp`.

:::{admonition} Tip: Vertical Interpolation
:class: admonition-icontheme

{material-outlined}`info;2em;pst-color-primary` To expose model level variables on pressure levels (or other vertical coordinates), you can use a dummy `output_nml` entry that performs the vertical interpolation. For example:

```fortran
&output_nml
 output_filename  = "pressure_level_for_output_coupling"
 filename_format  = "<output_filename>_<levtype_l>"
 filetype         = 999                       ! prohibits any file output
 remap            = 0
 output_start     = "${start_date}"           ! output_start = output_end
 output_end       = "${end_date}"             ! --> write once only irrespective of
 output_interval  = "${atmTimeStep}"          !     the output interval and
 file_interval    = "P999Y"    !     the file interval
 pl_varlist       = 'temp'
 p_levels         =  20000, 50000, 85000
/
```

This configuration triggers ICON's vertical interpolation machinery, making the pressure-level fields available in the variable list with the `pl::` prefix for output coupling, without actually writing any output files.
:::

:::{admonition} Warning: Diagnostic Variable Timing
:class: admonition-icontheme

{material-regular}`warning;2em;pst-color-secondary` Make sure you are using the correct `output_start`, `output_end` and `output_interval` settings. The diagnostic variables are only updated based on these values.
:::

:::{admonition} Collection Selection
:class: admonition-icontheme

{material-outlined}`info;2em;pst-color-primary` If you dont want to access all vertical levels of a field you can set the "collection_selection" configuration in the YAC coupling. See [YAC documentation](https://dkrz-sw.gitlab-pages.dkrz.de/yac/dd/dfa/yaml_file.html).
:::

## Metadata

The coupling automatically transfers metadata to output components in JSON format. For more details on metadata handling in YAC, see the [YAC metadata documentation](https://dkrz-sw.gitlab-pages.dkrz.de/yac/da/d1b/usage.html#phase_def_meta).

### CF Conventions
- `standard_name`: CF standard name
- `units`: Physical units
- `short_name`: Abbreviated name
- `long_name`: Descriptive name

### GRIB2 Information
- `discipline`: GRIB2 discipline code
- `category`: Parameter category
- `number`: Parameter number
- `bits`: Number of bits for packing
- `gridtype`: Grid type identifier
- `subgridtype`: Sub-grid type identifier
- `additional_keys`: Any extra GRIB2 key-value pairs

### Variable Groups
Variables are tagged with their group memberships for easier filtering in output components.

## Fractional Masking

Three standard mask fields are automatically created:

1. **`valid_mask_sfc`**: Surface (2D) land-sea mask
2. **`valid_mask`**: Full model levels (nlev)
3. **`valid_mask_half`**: Half levels (nlev+1)

### Mask Values

- `1.0`: valid cells
- `0.0`: invalid cells
- For vertex-centered data: A vertex is valid if ANY adjacent cell is valid

:::{admonition} Ocean Component Note
:class: admonition-icontheme

{material-outlined}`info;2em;pst-color-primary` In the ocean component, only valid cells/vertices carry valid data.
:::

## Performance Considerations

### Timers

The coupling tracks performance with these timers:

- `timer_coupling_output`: Total coupling time
- `timer_coupling_output_buf_prep`: Buffer preparation
- `timer_coupling_output_1stput`: First put operation
- `timer_coupling_output_put`: Subsequent put operations

Check these in ICON's timing output to assess coupling overhead.

## Example Output Component

ICON includes `run/checksuite.clm/simple_output.py` as a basic reference implementation:

```bash
python simple_output.py output.nc \
  atm,icon_grid,temp \
  atm,icon_grid,u \
  atm,icon_grid,v \
  --bounds -180 180 -90 90 \
  --res 360 180 \
  --output-interval PT1H
```

This interpolates ICON fields to a regular lat-lon grid and writes NetCDF output.

**Command-line Arguments:**

- `output.nc`: Output filename
- Positional: `component,grid,field` tuples of source fields (component name, grid name, field name)
- `--bounds`: Lon/lat bounds (min_lon, max_lon, min_lat, max_lat)
- `--res`: Resolution (n_lon, n_lat)
- `--output-interval`: ISO 8601 duration string

## Common Issues and Solutions

:::{admonition} Troubleshooting: Variables Not Appearing
:class: admonition-icontheme

If variables are not appearing in the coupling:

1. Verify that the `add_var` call was made by ICON
2. Set `msg_level >= 15` to see detailed messages explaining why variables in the var list are omitted
3. Check for correct component name (e.g. `atm` vs. `atmo`)
4. Check if all diagnostic variables are in a output_nml, to ensure that they are updated.
:::

## Custom Output Components

To create your own output component:

1. Initialize YAC and define your component
2. Define grids and points where you want to receive data
3. Define the target fields
4. Call `yac_fdef_couple()` to establish connections (optional, couplings can also be defined in a yaml file)
5. Call `yac_fenddef()` to finalize setup
6. In your time loop, call `yac_fget()` to receive data from ICON

See the [YAC documentation](https://gitlab.dkrz.de/dkrz-sw/yac) for details.

## Limitations

- **Single patch only**: Currently only `patch_id = 1` is supported

## References

- YAC Documentation: [https://dkrz-sw.gitlab-pages.dkrz.de/yac/](https://dkrz-sw.gitlab-pages.dkrz.de/yac/)
- Hiopy: [https://nils.gitlab-pages.dkrz.de/coyote/](https://gitlab.dkrz.de/dkrz-sw/yac)
