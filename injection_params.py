import numpy as np
import pandas as pd

def make_injection_params(separations_mas, contrasts):
    rows = []
    for i in range(len(separations_mas)):
        for contrast in contrasts:
            # Random position angles in degrees (0–360)
            phases = np.random.uniform(0, 2*np.pi, 5)
            for phase in phases:
                try:
                    sep = np.random.uniform(separations_mas[i], separations_mas[i+1])
                except:
                    sep = np.random.uniform(40,45)
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

if __name__ == "__main__":
    separations_mas = np.linspace(0, 40, 5)   # fibre FWHM is 65 mas(?)
    contrasts = np.array([1e-4, 1e-5, 1e-6])
    param_table = make_injection_params(separations_mas, contrasts)
    param_table.to_csv("/Users/svach/gravityInjectionRecovery/injection_parameters.csv", index=False)
