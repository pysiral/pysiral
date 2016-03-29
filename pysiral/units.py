# -*- coding: utf-8 -*-

from construct import Adapter

import numpy as np


class OneHundredth(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)/100)


class TenThousands(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*10000)


class OneHundredthDecibel(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)/100)


class OneTenthMicroDeg(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*1e-7)


class MicroDeg(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*1e-6)


class TenMicroDeg(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*1e-5)


class Centimeter(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*1e-2)


class MilliMeter(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*1e-3)


class Micrometer(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*1e-6)


class PicoSecond(Adapter):
    def _decode(self, obj, context):
        return float(obj)*1e-12


class MicroWatts(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*1e-6)


class MicroRadians(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*1e-6)


class OneTenthMicroRadians(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*1e-6)


class OneTenths(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*1e-1)


class OneThousands(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*1e-3)


class TenPascal(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)*1e-1)


class Per256(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)/256.)


class Per2048(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)/2048.)


class Per8096(Adapter):
    def _decode(self, obj, context):
        return np.float32(float(obj)/8096.)
