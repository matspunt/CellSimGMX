Introduction
============

``simpleble`` is a high-level OO Python package which aims to provide an easy and intuitive way of interacting with nearby Bluetooth Low Energy (BLE) devices (GATT servers). In essence, this package is an extension of the ``bluepy`` package created by Ian Harvey (see `here <https://github.com/IanHarvey/bluepy/>`_)

The aim here was to define a single object which would allow users to perform the various operations performed by the ``bluepy.btle.Peripheral``, ``bluepy.btle.Scanner``, ``bluepy.btle.Service`` and ``bluepy.btle.Characteristic`` classes of ``bluepy``, from one central place. This functionality is facilitated by the ``simpleble.SimpleBleClient`` and ``simpleble.SimpleBleDevice`` classes, where the latter is an extention/subclass of ``bluepy.btle.Peripheral``, combined with properties of ``bluepy.btle.ScanEntry``.

The current implementation has been developed in Python 3 and tested on a Raspberry Pi Zero W, running Raspbian 9 (stretch), but should work with Python 2.7+ (maybe with minor modifications in terms of printing and error handling) and most Debian based OSs.

Motivation
**********

As a newbie experimenter/hobbyist in the field of IoT using BLE communications, I found it pretty hard to identify a Python package which would enable one to use a Raspberry Pi (Zero W inthis case) to swiftly scan, connect to and read/write from/to a nearby BLE device (GATT server).

This package is intended to provide a quick, as well as (hopefully) easy to undestand, way of getting a simple BLE GATT client up and running, for all those out there, who, like myself, are hands-on learners and are eager to get their hands dirty from early on.


.. code-block:: python

    """This example demonstrates a simple BLE client that scans for devices,
    connects to a device (GATT server) of choice and continuously reads a characteristic on that device.

    The GATT Server in this example runs on an ESP32 with Arduino. For the
    exact script used for this example see `here <https://github.com/nkolban/ESP32_BLE_Arduino/blob/6bad7b42a96f0aa493323ef4821a8efb0e8815f2/examples/BLE_notify/BLE_notify.ino/>`_
    """

    from bluepy.btle import *
    from simpleble import SimpleBleClient, SimpleBleDevice
