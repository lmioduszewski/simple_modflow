import numpy as np


def calculate_calibration_statistics(observed, simulated):
    """
    Calculate common calibration statistics.

    Parameters:
    observed (np.array): Array of observed values.
    simulated (np.array): Array of simulated values.

    Returns:
    dict: Dictionary of calibration statistics.
    """
    if len(observed) != len(simulated):
        raise ValueError("Observed and simulated arrays must have the same length.")

    # Mean Error (ME)
    me = np.mean(simulated - observed)

    # Mean Absolute Error (MAE)
    mae = np.mean(np.abs(simulated - observed))

    # Root Mean Square Error (RMSE)
    rmse = np.sqrt(np.mean((simulated - observed) ** 2))

    # Normalized Root Mean Square Error (NRMSE)
    observed_range = np.max(observed) - np.min(observed)
    nrmse = rmse / observed_range

    # Nash-Sutcliffe Efficiency (NSE)
    nse = 1 - (np.sum((simulated - observed) ** 2) / np.sum((observed - np.mean(observed)) ** 2))

    # Coefficient of Determination (R^2)
    r_squared = np.corrcoef(observed, simulated)[0, 1] ** 2

    stats = {
        'ME': me,
        'MAE': mae,
        'RMSE': rmse,
        'NRMSE': nrmse,
        'NSE': nse,
        'R^2': r_squared
    }

    return stats