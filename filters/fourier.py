import numpy as np
import matplotlib.pyplot as plt
import PyQt6
_ = PyQt6

W = 1
W_f = 3
P = 14

#"""
a = 0.53836
def f(x): # hamming window
    return a + (1 - a) * np.cos(np.pi * x)
#"""

n = 2 ** P + 1
dx = (2 * W) / (n - 1)
x = np.linspace(-W, W, n)
y = [f(np.abs(_)) for _ in x]
k = np.linspace(-W_f, W_f, n)

F = np.abs([np.sum(y * np.exp(-2j * np.pi * f * x) * dx) for f in k])

_, axes = plt.subplots(2, 1)

formatter = lambda x, _: "" if x == 0 else f"{x:.12g}"
def setplot(ax):
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)

setplot(axes[0])
axes[0].plot(x, y, color="red")
axes[0].set_title("f(x)")

setplot(axes[1])
axes[1].plot(k, F, color="blue")
axes[1].set_title("F(k)")
axes[1].set_xlim(-1, 1)

plt.show()
