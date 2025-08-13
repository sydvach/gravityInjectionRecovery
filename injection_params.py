import numpy as np
import pandas as pd

def make_injection_params(separations_mas, contrasts):
    rows = []
    for sep in separations_mas:
        for contrast in contrasts:
            # Random position angles in degrees (0–360)
            phases = np.random.uniform(0, 2*np.pi, 5) if sep > 0 else [0]*5
            for phase in phases:
                delta_ra = sep * np.cos(phase)
                delta_dec = sep * np.sin(phase)
                rows.append({
                    "Contrast": contrast,
                    "Sep_mas": round(sep, 3),
                    "Phase_deg": round(np.degrees(phase), 3),
                    "deltaRA_mas": round(delta_ra, 4),
                    "deltaDec_mas": round(delta_dec, 4)
                })
    
    
    param_table = pd.DataFrame(rows)
    param_table['folderName'] = np.arange(len(param_table))
    return param_table

separations_mas = np.linspace(0, 40, 5)   # fibre FWHM is 65 mas(?)
contrasts = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
param_table = make_injection_params(separations_mas, contrasts)
param_table.to_csv("injection_parameters.csv", index=False)
