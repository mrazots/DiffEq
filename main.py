import sys
import matplotlib.pyplot as plt
import math


def func(x, y):
    return (1-2*y) * math.exp(x) + y*y + math.exp(2*x)


def isnan(num):
    if num != num or num == float('Inf'):
        return True


def construction(xx0, yy0, hh, xx):

    plt.ylim(-120, 1000)

    e = math.exp(1)

    x0 = xx0
    y0 = yy0
    h = hh
    x = xx

    #IVP

    xdots = []
    ydots = []

    x0 = xx0
    y0 = yy0

    xdots1 = [x0]
    ydots1 = [y0]

    c = 1/(2-math.exp(-5))-5
    i = x0 + h
    while i <= x:
        xdots1.append(i)
        ydots1.append(1/(c - i) + math.exp(i))
        i = i + h

    plt.figure(1)
    plt.plot(xdots1, ydots1, label = "IVP")
    plt.legend()


    #Euler’s method

    errorx = []
    errory = []

    i = x0
    while i <= x:
        errorx.append(i)
        i += h

    f = func(x0, y0)
    hf = h * f
    xdots.append(x0)
    ydots.append(y0)

    i = x0 + h
    while i < x:
        x1 = i
        y1 = y0 + hf
        if isnan(y1):
            y1 = sys.float_info.max
        xdots.append(x1)
        ydots.append(y1)

        f = func(x1, y1)
        hf = h * f
        y0 = y1
        i = i + h

    plt.figure(1)
    plt.plot(xdots, ydots, label = "Euler’s method")
    plt.legend()


    #Euler’s method error

    i = 0
    while i < len(ydots):
        errory.append(ydots[i] - ydots1[i])
        i += 1

    plt.figure(2)
    plt.plot(errorx, errory, label = "Euler's method errors")
    plt.legend()


    #Improved Euler’s method

    x0 = xx0
    y0 = yy0
    xdots = []
    ydots = []

    xdots.append(x0)
    ydots.append(y0)

    i = x0
    j = 0
    while i < x:
        arg1 = xdots[j] + h/2
        f = func(xdots[j], ydots[j])
        arg2 = ydots[j] + h * f/2
        f = func(arg1, arg2)
        hf = h * f
        y = ydots[j] + hf
        if isnan(y):
            y = sys.float_info.max
        xdots.append(i + h)
        ydots.append(y)

        i = i + h
        j = j + 1

    plt.figure(1)
    plt.plot(xdots, ydots, label = "Improved Euler’s method")
    plt.legend()


    #Improved Euler’s method error

    errory = []
    i = 0
    while i < len(ydots) - 1:
        errory.append(ydots[i] - ydots1[i])
        i += 1

    plt.figure(2)
    plt.plot(errorx, errory, label = "Improved Euler’s method errors")
    plt.legend()


    #Runge-Kutta method

    x0 = xx0
    y0 = yy0

    xdots = []
    ydots = []

    xdots.append(x0)
    ydots.append(y0)

    i = x0
    j = 0
    while i < x:
        k1 = func(xdots[j], ydots[j])
        k2 = func(xdots[j] + h/2, ydots[j] + h*k1/2)
        k3 = func(xdots[j] + h/2, ydots[j] + h*k2/2)
        k4 = func(xdots[j] + h, ydots[j] + h*k3)
        del_y = h * (k1 + 2*k2 + 2*k3 + k4) / 6
        y = ydots[j] + del_y
        if isnan(y):
            y = sys.float_info.max
        xdots.append(i+h)
        ydots.append(y)

        i = i + h
        j = j + 1

    plt.figure(1)
    plt.plot(xdots, ydots, label = "Runge-Kutta method")
    plt.legend()


    #Runge-Kutta method errors

    errory = []
    i = 0
    while i < len(ydots) - 1:
        errory.append(ydots[i] - ydots1[i])
        i += 1


    plt.figure(2)
    plt.plot(errorx, errory, label = "Runge-Kutta method errors")
    plt.legend()


    plt.ylim(-600, 1000)
    plt.show()


construction(-5, 2, 0.02, 0)