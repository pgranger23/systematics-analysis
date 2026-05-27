#!/usr/bin/env python3
"""
Test script to demonstrate YAML export functionality for oscillation fitter.
"""

from systs import SystConfig, create_oscillation_yaml_from_systematics

def main():
    # Load the systematics configuration
    print("Loading systematics configuration from systs.fcl...")
    config = SystConfig('systs.fcl')
    
    # Get all systematics
    all_systs = config.get_all_systs()
    print(f"Found {len(all_systs)} systematics")
    print()
    
    # Example 1: Export all systematics using the built-in method
    print("=" * 70)
    print("Example 1: Export all systematics using SystConfig.to_oscillation_yaml()")
    print("=" * 70)
    yaml_str = config.to_oscillation_yaml(
        output_file='oscillation_all_systematics.yaml',
        sample_names=["ND_*", "FD_*", "ATM"],
        error=1.0,
        mcmc_step_scale=0.4
    )
    print("First systematic in YAML format:")
    print("-" * 70)
    # Print just the first systematic as an example
    lines = yaml_str.split('\n')
    for i, line in enumerate(lines):
        print(line)
        if i > 20:  # Print first ~20 lines
            print("... (truncated)")
            break
    print()
    
    # Example 2: Export only specific systematics (e.g., those containing 'MEC')
    print("=" * 70)
    print("Example 2: Export only MEC systematics")
    print("=" * 70)
    mec_systs = [s for s in all_systs if 'MEC' in s.name.upper()]
    if mec_systs:
        print(f"Found {len(mec_systs)} MEC systematics:")
        for s in mec_systs:
            print(f"  - {s.name}")
        
        yaml_str = create_oscillation_yaml_from_systematics(
            mec_systs,
            output_file='oscillation_mec_systematics.yaml',
            param_group='MEC'
        )
        print()
    else:
        print("No MEC systematics found")
        print()
    
    # Example 3: Export with custom parameters
    print("=" * 70)
    print("Example 3: Export first 3 systematics with custom parameters")
    print("=" * 70)
    if len(all_systs) >= 3:
        first_three = all_systs[:3]
        yaml_str = create_oscillation_yaml_from_systematics(
            first_three,
            output_file='oscillation_custom.yaml',
            sample_names=["ND_numu", "FD_nue"],
            error=0.5,
            flat_prior=True,
            param_bounds=[-3.0, 3.0],
            mcmc_step_scale=0.2,
            param_group='CustomGroup'
        )
        print("Exported systematics:")
        for s in first_three:
            print(f"  - {s.name}")
    
    print()
    print("=" * 70)
    print("All files generated successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
