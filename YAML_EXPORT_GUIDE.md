# YAML Export for Oscillation Fitter

This code provides functionality to convert `SystConfig` objects into YAML format suitable for the oscillation fitter.

## Features

- Convert all systematics from a FHiCL configuration to YAML format
- Customize parameters for each systematic
- Filter and export specific systematics
- Two usage patterns: method on `SystConfig` or standalone function

## Quick Start

```python
from systs import SystConfig

# Load configuration
config = SystConfig('systs.fcl')

# Export all systematics to YAML
yaml_str = config.to_oscillation_yaml(output_file='oscillation_config.yaml')
```

## API Reference

### `SystConfig.to_oscillation_yaml()`

Convert the entire SystConfig to YAML format for the oscillation fitter.

**Parameters:**
- `output_file` (str, optional): Path to save the YAML file. If None, returns the YAML string.
- `sample_names` (list[str], optional): List of sample names. Default: `["ND_*", "FD_*", "ATM"]`
- `error` (float, optional): Error value for the systematic. Default: `1.0`
- `flat_prior` (bool, optional): Whether to use a flat prior. Default: `False`
- `param_bounds` (list[float], optional): Parameter bounds [min, max]. Default: `[-9999.0, 9999.0]`
- `generated_value` (float, optional): Generated parameter value. Default: `0.0`
- `prefit_value` (float, optional): Pre-fit parameter value. Default: `0.0`
- `spline_mode` (list[int], optional): Mode for spline. Default: `[0]`
- `spline_type` (str, optional): Type of spline. Default: `"PerEvent"`
- `mcmc_step_scale` (float, optional): MCMC step scale. Default: `0.4`
- `syst_type` (str, optional): Type of systematic. Default: `"Spline"`
- `param_group` (str, optional): Parameter group. Default: `"Xsec"`

**Returns:**
- `str`: YAML formatted string

### `create_oscillation_yaml_from_systematics()`

Create a YAML configuration for the oscillation fitter from a list of `Systematic` objects.

This standalone function provides the same functionality as the method but can be used with filtered or custom lists of systematics.

**Parameters:** Same as `to_oscillation_yaml()` plus:
- `systematics` (list[Systematic]): List of Systematic objects to convert

## Usage Examples

### Example 1: Export All Systematics

```python
from systs import SystConfig

config = SystConfig('systs.fcl')
yaml_str = config.to_oscillation_yaml(
    output_file='all_systematics.yaml',
    sample_names=["ND_*", "FD_*", "ATM"],
    error=1.0,
    mcmc_step_scale=0.4
)
```

### Example 2: Export Specific Systematics

```python
from systs import SystConfig, create_oscillation_yaml_from_systematics

config = SystConfig('systs.fcl')
all_systs = config.get_all_systs()

# Filter for MEC systematics only
mec_systs = [s for s in all_systs if 'MEC' in s.name.upper()]

yaml_str = create_oscillation_yaml_from_systematics(
    mec_systs,
    output_file='mec_systematics.yaml',
    param_group='MEC'
)
```

### Example 3: Custom Parameters for Different Groups

```python
from systs import SystConfig, create_oscillation_yaml_from_systematics

config = SystConfig('systs.fcl')
all_systs = config.get_all_systs()

# Export cross-section systematics with custom bounds
xsec_systs = [s for s in all_systs if 'xsec' in s.name.lower()]
yaml_str = create_oscillation_yaml_from_systematics(
    xsec_systs,
    output_file='xsec_systematics.yaml',
    param_bounds=[-5.0, 5.0],
    mcmc_step_scale=0.3,
    param_group='Xsec'
)

# Export flux systematics with different parameters
flux_systs = [s for s in all_systs if 'flux' in s.name.lower()]
yaml_str = create_oscillation_yaml_from_systematics(
    flux_systs,
    output_file='flux_systematics.yaml',
    param_bounds=[-3.0, 3.0],
    mcmc_step_scale=0.5,
    param_group='Flux'
)
```

### Example 4: Corrections vs Systematics

The code automatically handles corrections differently from regular systematics. For corrections:
- `Generated` and `PreFitValue` are set to the correction's central value (`syst.cv`)
- For regular systematics, they use the provided defaults

```python
config = SystConfig('systs.fcl')
all_systs = config.get_all_systs()

# Separate corrections from systematics
corrections = [s for s in all_systs if s.isCorrection]
systematics = [s for s in all_systs if not s.isCorrection]

print(f"Found {len(corrections)} corrections and {len(systematics)} systematics")
```

## Output Format

The generated YAML follows this structure:

```yaml
- Systematic:
    SampleNames: ["ND_*", "FD_*", "ATM"]
    Error: 1.0
    FlatPrior: false
    Names:
      FancyName: systematic_name
      ParameterName: systematic_name
    ParameterBounds:
    - -9999.0
    - 9999.0
    ParameterValues:
      Generated: 0.0
      PreFitValue: 0.0
    SplineInformation:
      Mode: [0]
      SplineName: systematic_name
      Type: PerEvent
    StepScale:
      MCMC: 0.4
    Type: Spline
    ParameterGroup: Xsec
```

## Testing

Run the test script to see examples:

```bash
python test_yaml_export.py
```

This will generate several example YAML files demonstrating different use cases.

## Requirements

- `yaml` (PyYAML): For YAML serialization
- All existing dependencies in `systs.py`

Install with:
```bash
pip install pyyaml
```
