# Relevant Documentation: SD1 3.x Software for M310xA / M330xA Digitizers User's Guide

import numpy as np

# Path must be setup before this import
import keysightSD1

DIG_SAMPLE_RATE_MSPS = 500 # MS/s

def setup_dig_pxi(
    dig,
    channel,
    duration_ns,
    full_scale,
    prescaler=0,
    impedance=keysightSD1.AIN_Impedance.AIN_IMPEDANCE_HZ,
    coupling=keysightSD1.AIN_Coupling.AIN_COUPLING_DC,
    trigger_pxi_line=1,
    trigger_behavior=keysightSD1.SD_TriggerBehaviors.TRIGGER_FALL):
    
    n_points = int(DIG_SAMPLE_RATE_MSPS * duration_ns // 1000)
    
    # Stop and flush DAQ
    dig.daq_stop(channel)
    dig.daq_flush(channel)
        
    # Setup input channel
    dig.set_prescaler(prescaler, channel)
    dig.set_full_scale(full_scale, channel)
    dig.set_impedance(impedance, channel)
    dig.set_coupling(coupling, channel)
    
    # Setup DAQ
    dig.set_points_per_cycle(n_points, channel)
    dig.set_n_cycles(1, channel)
    dig.set_daq_trigger_delay(0, channel)
    dig.set_daq_trigger_mode(keysightSD1.SD_TriggerModes.HWDIGTRIG, channel)
    
    # Setup DAQ digital trigger
    dig.set_digital_trigger_source(trigger_pxi_line + 4000, channel)
    dig.set_digital_trigger_behaviour(trigger_behavior, channel)
    
    # Setup DAQ read
    dig.set_n_points(n_points, channel)
    dig.set_timeout(5000, channel)
    
    # Start DAQ
    dig.daq_start(channel)
