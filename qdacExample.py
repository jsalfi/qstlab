# qdacExample.py
# Copyright QDevil ApS, 2018, 2019
VERSION = "1.21"

import qdac
import random
import math

with qdac.qdac('/dev/cu.usbserial-FTA634UE') as q:
    print("QDAC Serial number: %s" % q.getSerialNumber())
    print("Number of channels: %d" % q.getNumberOfChannels())

    print("-----------------------------------------------")
    print("Setting Channel 1 voltage range to 10 V")
    result = q.setVoltageRange(channel=1, theRange=10)
    print("Result: %s" % result)
    print("Setting Channel 1 DC voltage to 1.23456 V")
    result = q.setDCVoltage(channel=1, volts=1.23456)
    print("Result: %s" % result)
    voltage1 = q.getDCVoltage(1)
    print("Voltage output on channel 1 is %f V" % voltage1)
    current1 = q.getCurrentReading(1)
    print("Current on channel 1 is %e A" % current1)

    print("-----------------------------------------------")
    print("Defines triangle function generator for generator 1 and starts it on channel 2")
    result = q.defineFunctionGenerator(generator=qdac.Generator.generator1, waveform=qdac.Waveform.triangle, period=100, dutycycle=50)
    print("Result: %s" % result)
    q.setChannelOutput(channel=2, generator=qdac.Generator.generator1, amplitude=1.0, offset=0.0)

    print("-----------------------------------------------")
    print("Defines pulse train and starts it on channel 3")
    result = q.definePulsetrain(lowDuration=100, highDuration=2000, lowVolts=0.01, highVolts=0.03)
    print("Result: %s" % result)
    q.setChannelOutput(channel=3, generator=qdac.Generator.pulsetrain, amplitude=1.0, offset=0.0)

    print("-----------------------------------------------")
    print("Sets Channel 4 current range to +/- 1uA, before changing the Voltage range")
    result = q.setCurrentRange(channel=4, theRange=1.0e-6)
    print("Result: %s" % result)
    print("Sets Channel 4 voltage range to +/- 1 Volt and sets the DC voltage to 0.123456 V")
    result = q.setVoltageRange(channel=4, theRange=1.0)
    print("Result: %s" % result)
    result = q.setDCVoltage(channel=4, volts=0.123456)
    print("Result: %s" % result)

    print("-----------------------------------------------")
    print("Defines arbitrary waveform and start it on channel 5")
    n = 800
    samples = samples = [-5*(0.1*i/n)*math.sin(math.pi*i/(n/20)) for i in range(n)]
    result = q.defineAWG(samples)
    print("Result: %s" % result)
    q.setChannelOutput(channel=5, generator=qdac.Generator.AWG, amplitude=1.0, offset=0.0)

    print("-----------------------------------------------")
    print("Sets sync output channel 1 to trigger on start of AWG waveform")
    result = q.setSyncOutput(1, qdac.Generator.AWG)
    print("Result: %s" % result)

    print("-----------------------------------------------")
    print("Soft-waits for beginning of pulsetrain")
    if q.waitForSync(qdac.Generator.pulsetrain, timeout=15):
        print("Got a sync message - beginning of pulsetrain")
    else:
        print("WaitForSync timeout")
    input("Press enter to set all channels to DC, 1V range and 0.123456V output")
    for ch in range(1, q.getNumberOfChannels() + 1):
        q.setChannelOutput(ch, qdac.Generator.DC)
        q.setCurrentRange(channel=ch, theRange=1.0e-6)   # Current range HAS to be 1uA befoe switching voltage range to 1V
        q.setVoltageRange(channel=ch, theRange=1.0)
        q.setDCVoltage(channel=ch, volts=.123456)

