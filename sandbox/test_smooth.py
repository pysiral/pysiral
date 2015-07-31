# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


def test_smooth():

    x = np.linspace(0, 1, 500, endpoint=False)
    # y = signal.square(2*np.pi*5*x)
    y = np.zeros(500)
    y[100:500:50] = 100

    plt.figure(facecolor="white")
    plt.plot(x, y)
    plt.plot(x, smooth(y, 3))
    #plt.ylim(-1.5, 1.5)
    plt.show()


def smooth(x, window):
    return np.convolve(x, np.ones(window)/window, mode='same')

if __name__ == "__main__":
    test_smooth()
