# Relevant Documentation: SD1 3.x Software for M320xA / M330xA Arbitrary Waveform Generators User's Guide

from enum import Enum

import numpy as np

# Path must be setup before this import
import keysightSD1

# See section 1.2.2 of the documentation
MIN_CYCLIC_QUEUE_DURATION_NS = 2000

def get_prescaler(sample_rate_MSps):
    """
    Calculate the prescaler for the given sample rate (in MS/s).
    
    See section 1.2.3 of the documentation
    
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

def check_duration(duration_ns):
    if len(duration_ns) < MIN_CYCLIC_QUEUE_DURATION_NS:
        raise Exception(f"Not enough data: {duration_ns} ns. Minimum {MIN_CYCLIC_QUEUE_DURATION_NS} ns required.")
    
def create_sine(period_ns, repetition, sample_rate_MSps):
    n_pts = get_num_points(period_ns * repetition, sample_rate_MSps)
    
    phi = np.linspace(0, np.pi*2*repetition, n_pts)
    w = np.sin(phi)
    
    return w

def create_gaussian(a, b, c, duration_ns, sample_rate_MSps):
    n_pts = get_num_points(duration_ns, sample_rate_MSps)
    
    x = np.linspace(-(duration_ns // 2), (duration_ns // 2), n_pts)
    w = a * np.exp(-(x - b)**2 / (2 * c**2))
    
    return w
    
def queue_oneshot(awg, channel, wave, sample_rate_MSps, delay=0, cycles=1, trigger=keysightSD1.SD_TriggerModes.AUTOTRIG, trigger_pxi_line=1, trigger_mode=keysightSD1.SD_TriggerBehaviors.TRIGGER_FALL):
    awg.set_channel_offset(0.0, channel)
    awg.set_channel_amplitude(0.5, channel)

    awg.set_channel_wave_shape(keysightSD1.SD_Waveshapes.AOU_AWG, channel)
    awg.awg_queue_config(channel, keysightSD1.SD_QueueMode.ONE_SHOT)
    
    if trigger == keysightSD1.SD_TriggerModes.EXTTRIG:
        awg.awg_config_external_trigger(channel, 4000 + trigger_pxi_line, trigger_mode) # See Table 36 of the documentation for why 4000 is added
    
    wave_awg = awg.upload_waveform(wave)
    
    prescaler    = get_prescaler(sample_rate_MSps)
    awg.awg_queue_waveform(channel, wave_awg, trigger, delay, cycles, prescaler)
    
    awg.awg_start(channel)

def trigger_fall(awg, trigger_pxi_line):
    awg.set_pxi_trigger(0, trigger_pxi_line) # 0 = ON
    awg.set_pxi_trigger(1, trigger_pxi_line) # 1 = OFF