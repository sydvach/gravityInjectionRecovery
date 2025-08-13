import numpy as np
import pandas as pd

def make_injection_params(separations_mas, contrasts):
    # Prepare an empty list for all parameter sets
    rows = []
    
    # Loop over separation and magnitude, and sample random phases
    for sep in separations_mas:
        for contrast in contrasts:
            # Random position angles in degrees (0–360)
            phases = np.random.uniform(0, 2*np.pi, 5) if sep > 0 else [0]*5
            for phase in phases:
                delta_ra = sep * np.cos(phase)
                delta_dec = sep * np.sin(phase)
                rows.append({
                    "contrast": contrast,
                    "sep_mas": round(sep, 3),
                    "phase_deg": round(np.degrees(phase), 3),
                    "deltaRA_mas": round(delta_ra, 4),
                    "deltaDec_mas": round(delta_dec, 4)
                })
    
    # Convert to DataFrame
    param_table = pd.DataFrame(rows)
    return param_table

if __name__ == '__main__':

separations_mas = np.linspace(0, 65, 5)   # fibre FWHM is 65 mas(?)
contrasts = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
param_table = make_injection_params(separations_mas, contrasts)
param_table.to_csv("injection_parameters.csv", index=False)
