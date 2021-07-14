import matplotlib.pyplot as plt
import numpy as np
import time
from tqdm import tqdm

"""
High level classes and functions for common processes for
controlling the Keysight FieldFpx RF Analyzer N9914A VNA

Note that the code has only been tested on this particular VNA model
"""
import pyvisa

rm = pyvisa.ResourceManager()


def list_resources():
    """
    List all resources connected to the device

    :return: tuple of strings
    """
    return rm.list_resources()


class InvalidFileTypeError(Exception):
    pass


class NA:
    """
    Operate device in Network Analyzer (NA) mode
    """
    def __init__(self, resource: str, timeout: int = 10_000):
        """
        :param resource: str, VISA resource. see `func:list_resources`
        :param timeout: int, operations timeout in ms, defaults to 10000
        """
        self.na = rm.open_resource(resource)
        self.na.timeout = timeout
        self.reset()

    def reset(self, force: bool = False):
        """
        Reset the device

        :param force: *RST if True, defaults to False
        """
        if force:
            self.na.write("*RST")
        self.na.write("*CLS")
        self.na.write("SYST:PRES")
        # Set Mode to NA
        self.na.write("INST:SEL 'NA'")

    def single_sweep(self):
        self.na.write('INIT:IMM;*OPC?')

    def trigger_mode(self, mode: bool):
        """
        Set trigger mode to hold for trigger synchronization

        :param mode: True for ON, False for OFF
        """
        mode = "ON" if mode else "OFF"
        self.na.write(f"INIT:CONT {mode};*OPC?")  # OFF for hold and ON for continue

    def output_mode(self, mode: "S21" or "S11"):
        self.na.write(f'CALC:PAR:DEF {mode}')

    def trigger(self):
        """
        Trigger a single measurement and wait
        """
        self.na.write("INIT: IMMediate;*OPC?")

    def measurement(
            self,
            output_mode: "S21" or "S11",
            num_points: int,
            start_freq: float,
            end_freq: float,
            power: float,
            bandwidth: float,
            average_pts: int,
            single_sweep: bool = True,
    ):
        """
        Trigger a single measurement with the specified configuration,
        and wait until it is complete

        :param output_mode: "S21" or "S11"
        :param num_points: int, number of data points to acquire
        :param start_freq: float, starting frequency
        :param end_freq: float, ending frequency
        :param power: float, magnitude of signal power in mDb
        :param bandwidth: float, measurement bandwidth
        :param average_pts: int, number of measurements to average over at each point
        :param single_sweep: bool, single measurement if True, else continuous. Defaults to True.
        """
        self.output_mode(output_mode)
        self.trigger_mode(single_sweep)
        self.na.write("SENS:SWE:POIN " + str(num_points))
        self.na.write("SENS:FREQ:START " + str(start_freq))
        self.na.write("SENS:FREQ:STOP " + str(end_freq))
        self.na.write("SOUR:POW -" + str(power))
        self.na.write("BWID " + str(bandwidth))
        self.na.write("AVER:COUN " + str(average_pts))
        self.trigger()

    def save_data(self, fname: str, fmt: "csv" or "s1p" or "png") -> None:
        """

        :param fname: local filename without extension
        :param fmt: "csv" or "s1p" or "png")
        :return: None
        """
        commands = {
            "png": "MMEM:STOR:IMAG",
            "csv": "MMEM:STOR:FDAT",
            "s1p": "MMEM:STOR:SNP"
        }
        if fmt not in ("csv", "s1p", "png"):
            raise InvalidFileTypeError
        fname = fname+"."+fmt
        command = commands[fmt]
        self.na.write(f'MMEM:DEL "temp.{fmt}"')
        self.na.write(f'{command} "temp.{fmt}"')

        # BINBLOCK transfer
        data = self.na.query_binary_values(
            f'MMEM:DATA? "temp.{fmt}"',
            datatype="B",
            container=bytearray
        )

        with open(fname, "wb") as f:
            f.write(data)

        self.na.write(f'MMEM:DEL "temp.{fmt}"') # delete from VNA to save space

    def get_handle(self) -> pyvisa.Resource:
        """
        Return the pyVisa resource
        :return: pyvisa.Resource
        """
        return self.na

print(list_resources())
na = NA('ASRL/dev/ttyACM0::INSTR') 

# Parameters
powers = [25]

output_mode= "S21"
num_points = 200
start_freq = 4e9
end_freq = 6e9
bandwidth = 10e3
average_pts = 10

# for each iteration: measure -> save data
for i, power in tqdm(enumerate(powers)):
    na.measurement(
        output_mode,
        num_points,
        start_freq,
        end_freq,
        power,
        bandwidth,
        average_pts,
    )
    
    na.save_data(f"Power_Sweep/power_{i}", "png")