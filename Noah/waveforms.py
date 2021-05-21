from enum import Enum

import numpy as np

def get_prescaler(sample_rate_MSps):
    """
    Calculate the prescaler for the given sample rate (in MS/s).
    
    See section 1.2.3 of the SD1 3.x Software for M320xA / M330xA Arbitrary Waveform Generators User's Guide
    
    The effective sampling rate is calculated as:
    
        1 GS/s              , prescaler == 0
        200 MS/s            , prescaler == 1
        100 MS/s / prescaler, prescaler > 1
    """
    
    if sample_rate_MSps == 1000:
        return 0
    elif sample_rate_MSps == 200:
        return 1
    elif sample_rate_MSps <= 50:
        if 100 % sample_rate_MSps == 0:
            return 100 // sample_rate_MSps
    
    raise Exception(f"Invalid sample rate: {sample_rate_MSps} MS/s. Sample rate must be 1 GS/s, 200 MS/s or 100 MS/s divided by an integer")

def get_num_points(time_ns, sample_rate_MSps):
    time_step_ns = 1000 / sample_rate_MSps
    return int(time_ns / time_step_ns)
    
def create_sine(period_ns, repetition, sample_rate_MSps):
    n_pts = get_num_points(period_ns * repetition, sample_rate_MSps)
    phi = np.linspace(0, np.pi*2*repetition, n_pts)
    w = np.sin(phi)
    
    if len(w) < 2000:
        raise Exception("not enough data")
        
    return w



# def get_time_step(prescaler):
#     """
#     Calculate the time step between points in nanoseconds based on the prescaler
#     """
    
#     if prescaler == 0:
#         return 1 # ns
#     elif prescaler == 1:
#         return 5 # ns
#     elif (prescaler > 1 && prescaler < 4096):
#         return prescaler * 10 # ns
    
#     raise Exception(f"Invalid prescaler: {prescaler}. Prescaler must be between 0 and 4095")

# def get_time_step(sample_rate_Msps):
#     """
#     Calculate the time step between points in nanoseconds
#     """
#     return 1000 / sample_rate_Msps # ns
    