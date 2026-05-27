#!/usr/bin/env python3
"""
Simple example demonstrating the YAML export functionality.
This creates a minimal example without needing to load an actual FCL file.
"""

import yaml
from typing import Optional

# Minimal Systematic class for demonstration
class Systematic:
    def __init__(self, name: str, is_correction: bool = False, cv: float = 0.0):
        self.name = name
        self.isCorrection = is_correction
        self.cv = cv
        self.id = hash(name) % 1000  # Simple ID
        self.paramVariations = [-1, 0, 1] if not is_correction else [1]

def create_single_systematic_yaml(systematic: Systematic,
                                 sample_names: Optional[list[str]] = None,
                                 error: float = 1.0,
                                 flat_prior: bool = False,
                                 param_bounds: Optional[list[float]] = None,
                                 generated_value: float = 0.0,
                                 prefit_value: float = 0.0,
                                 spline_mode: Optional[list[int]] = None,
                                 spline_type: str = "PerEvent",
                                 mcmc_step_scale: float = 0.4,
                                 syst_type: str = "Spline",
                                 param_group: str = "Xsec") -> dict:
    """Create a YAML dict for a single systematic."""
    
    # Set defaults
    if sample_names is None:
        sample_names = ["ND_*", "FD_*", "ATM"]
    if param_bounds is None:
        param_bounds = [-9999.0, 9999.0]
    if spline_mode is None:
        spline_mode = [0]
    
    return {
        'Systematic': {
            'SampleNames': sample_names,
            'Error': error,
            'FlatPrior': flat_prior,
            'Names': {
                'FancyName': systematic.name,
                'ParameterName': systematic.name
            },
            'ParameterBounds': param_bounds,
            'ParameterValues': {
                'Generated': generated_value if not systematic.isCorrection else systematic.cv,
                'PreFitValue': prefit_value if not systematic.isCorrection else systematic.cv
            },
            'SplineInformation': {
                'Mode': spline_mode,
                'SplineName': systematic.name,
                'Type': spline_type
            },
            'StepScale': {
                'MCMC': mcmc_step_scale
            },
            'Type': syst_type,
            'ParameterGroup': param_group
        }
    }


def main():
    print("=" * 70)
    print("YAML Export Example for Oscillation Fitter")
    print("=" * 70)
    print()
    
    # Create some example systematics
    systematics = [
        Systematic("MaCCQE", is_correction=False),
        Systematic("MaRES", is_correction=False),
        Systematic("MEC_C", is_correction=True, cv=1.0),
        Systematic("FSI_pi_abs", is_correction=False),
    ]
    
    print(f"Creating YAML for {len(systematics)} systematics:")
    for s in systematics:
        correction_label = " (correction)" if s.isCorrection else ""
        print(f"  - {s.name}{correction_label}")
    print()
    
    # Convert to YAML
    yaml_list = []
    for syst in systematics:
        yaml_dict = create_single_systematic_yaml(
            syst,
            sample_names=["ND_*", "FD_*", "ATM"],
            error=1.0,
            mcmc_step_scale=0.4,
            param_group="Xsec"
        )
        yaml_list.append(yaml_dict)
    
    # Convert to YAML string
    yaml_str = yaml.dump(yaml_list, default_flow_style=False, sort_keys=False)
    
    # Print the YAML
    print("Generated YAML:")
    print("-" * 70)
    print(yaml_str)
    print("-" * 70)
    
    # Save to file
    output_file = "example_oscillation_config.yaml"
    with open(output_file, 'w') as f:
        f.write(yaml_str)
    print(f"\nSaved to: {output_file}")
    print()
    
    # Show one systematic in detail
    print("=" * 70)
    print("Example: Single Systematic Structure")
    print("=" * 70)
    print("\nFor systematic 'MaCCQE':")
    print(yaml.dump([yaml_list[0]], default_flow_style=False, sort_keys=False))


if __name__ == "__main__":
    main()
