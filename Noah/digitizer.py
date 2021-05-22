import numpy as np

# Path must be setup before this import
import keysightSD1

def setup_dig(dig, channel, full_scale, impedance=keysightSD1.AIN_Impedance.AIN_IMPEDANCE_HZ, coupling=keysightSD1.AIN_Coupling.AIN_COUPLING_DC, trigger_pxi_line=1, trigger_mode=keysightSD1.SD_AIN_TriggerMode.FALLING_EDGE):
    dig.daq_stop(channel)
    dig.daq_flush(channel)
    
    # TODO:
    # - full_scale
    # - prescaler
    # - trigger_mode
    
    
    dig.set_full_scale(full_scale, channel)
    dig.set_impedance(impedance, channel)
    dig.set_coupling(coupling, channel)
    
    dig.set_prescaler(prescaler, channel)
    # TODO not sure about triggers
    dig.set_digital_trigger_source(trigger_pxi_line, channel)
    dig.set_digital_trigger_mode(trigger_mode, channel)
    
    moduleIn.channelPrescalerConfig(ch, DIG_PRESCALER)
    moduleIn.DAQdigitalTriggerConfig(ch, pxi1, trigger_mode)
    moduleIn.DAQconfig(ch, TOT_POINTS_IN, 1, DELAY_IN, ext_trigger)