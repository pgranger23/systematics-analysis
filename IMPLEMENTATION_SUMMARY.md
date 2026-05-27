# YAML Export Implementation Summary

## What Was Implemented

I've added functionality to convert `SystConfig` objects into YAML format compatible with the oscillation fitter.

## Files Modified

### 1. **systs.py** (Modified)
Added two main functions:

#### `SystConfig.to_oscillation_yaml()` method
- Converts all systematics from a SystConfig object to YAML
- Located at line ~355
- Returns YAML string and optionally saves to file

#### `create_oscillation_yaml_from_systematics()` standalone function  
- Converts a list of Systematic objects to YAML
- Located at line ~819
- More flexible - works with filtered lists

Both functions support customization of:
- Sample names (default: `["ND_*", "FD_*", "ATM"]`)
- Error values
- Parameter bounds
- MCMC step scales
- Spline configuration
- Parameter groups

### 2. **test_yaml_export.py** (New)
Complete test script demonstrating:
- Exporting all systematics
- Filtering specific systematics (e.g., MEC only)
- Custom parameter configurations

### 3. **example_yaml_simple.py** (New)
Minimal standalone example that:
- Shows the YAML structure clearly
- Works without needing FCL files
- Demonstrates both regular systematics and corrections

### 4. **YAML_EXPORT_GUIDE.md** (New)
Comprehensive documentation including:
- API reference
- Usage examples
- Output format description
- Installation requirements

## Usage Examples

### Basic Usage
```python
from systs import SystConfig

# Load configuration
config = SystConfig('systs.fcl')

# Export all systematics
yaml_str = config.to_oscillation_yaml(output_file='oscillation_config.yaml')
```

### Filter Specific Systematics
```python
from systs import SystConfig, create_oscillation_yaml_from_systematics

config = SystConfig('systs.fcl')
all_systs = config.get_all_systs()

# Export only MEC systematics
mec_systs = [s for s in all_systs if 'MEC' in s.name.upper()]
yaml_str = create_oscillation_yaml_from_systematics(
    mec_systs,
    output_file='mec_systematics.yaml',
    param_group='MEC'
)
```

### Custom Parameters
```python
yaml_str = config.to_oscillation_yaml(
    output_file='custom_config.yaml',
    sample_names=["ND_numu", "FD_nue"],
    error=0.5,
    param_bounds=[-3.0, 3.0],
    mcmc_step_scale=0.2
)
```

## Output Format

The generated YAML matches the required format:

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

## Special Features

1. **Automatic Correction Handling**: Corrections automatically use their central value (`syst.cv`) for Generated and PreFitValue
2. **Flexible Filtering**: Easy to filter systematics by name, type, or any other property
3. **Batch Export**: Can export all systematics or specific subsets
4. **Type Safety**: Uses proper Optional type hints for Python type checking

## Testing

Run the example to verify:
```bash
python3 example_yaml_simple.py
```

This will generate `example_oscillation_config.yaml` demonstrating the output format.

## Dependencies

Requires PyYAML:
```bash
pip install pyyaml
```

All other dependencies are already in systs.py.
