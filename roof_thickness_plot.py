
import numpy as np
import matplotlib.pyplot as plt


def calc(idx=0):
    b = [1.107,0.99227,1.181]
    tr = [286.5,286.5,382]
    td = [36.5,18.3,33.51]
    ri = [56.4,31.8,48.62]

    x = np.arange(0,b[idx],0.001)

    result = tr[idx] * 1 / np.sqrt(1 - np.power(td[idx]/ri[idx] * np.sin(x), 2))
    return x, result


if __name__ == "__main__":

    fig, ax = plt.subplots()
    for i in range(3):
        x, res = calc(idx=i)
        ax.plot(x, res, label=i)
    ax.legend()
    plt.show()