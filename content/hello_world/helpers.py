import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def check_first():
    a = np.loadtxt("./out/a.txt")
    b = np.loadtxt("./out/b.txt")
    c = np.loadtxt("./out/c.txt")
    if np.all(c == a + b):
        print("Success (`c == a + b`)!")
    else:
        print("At some point there is `c != a + b`...")


def plot_c():
    T = np.loadtxt("out/meshT.txt")
    x, y = np.loadtxt("out/meshP.txt", unpack=True)
    c = np.loadtxt("out/c.txt")

    T_ref = np.loadtxt("refs/contours/meshT.txt")
    x_ref, y_ref = np.loadtxt("refs/contours/meshP.txt", unpack=True)
    c_ref = np.loadtxt("refs/contours/c.txt")

    fig = plt.figure(figsize=(10, 4), dpi=100)

    ax1 = fig.add_subplot(121)
    ax1.axis("equal")
    ax1.set_title("Computed")
    tp1 = ax1.tripcolor(x, y, triangles=T - 1, facecolors=c)
    plt.colorbar(tp1, ax=ax1)

    ax2 = fig.add_subplot(122)
    ax2.axis("equal")
    ax2.set_title("Expected")
    tp2 = ax2.tripcolor(x_ref, y_ref, triangles=T_ref - 1, facecolors=c_ref)
    plt.colorbar(tp2, ax=ax2)

    plt.show()


w = h = 77


def image_to_data(filename):
    global w, h
    im = plt.imread(filename)
    w, h, c = im.shape
    im = im.reshape((w * h * c,))
    s = f"{w}\n{h}\n" + ("\n".join(map(str, im)))
    with open("picture.txt", "w") as f:
        f.write(s)


def make_color_map(r, g, b):
    c = np.vstack((r, g, b)).T / 255
    c = ListedColormap(c)

    v = np.arange(r.shape[0])

    return c, v


def load_image(path):
    T = np.loadtxt(path + "/meshT.txt")
    x, y = np.loadtxt(path + "/meshP.txt", unpack=True)
    # need to correct scales
    x = x / np.max(x) * w
    y = (1 - y / np.max(y)) * h

    r = np.loadtxt(path + "/r.txt")
    g = np.loadtxt(path + "/g.txt")
    b = np.loadtxt(path + "/b.txt")

    c, v = make_color_map(r, g, b)

    return T, x, y, v, c


def plot_saturation_curve():

    x = np.linspace(0, 1, 50)
    plt.plot(x, x, ".", label="Identity")

    def curve(x, p, q):
        return ((1 - x) ** p * (x / q) + x ** p * ((x - 1) / q + 1)) / (
            x ** p + (1 - x) ** p
        )

    def plot_curve(x, p, q):
        y = np.array([curve(i, p, q) for i in x])
        plt.plot(x, y, label=f"p = {p}, q = {q}")

    plot_curve(x, 4, 3)
    plot_curve(x, 2, 2)

    plt.legend()
    plt.show()


def plot_image(part):
    T, x, y, v, c = load_image("./out/")
    T_orig, x_orig, y_orig, v_orig, c_orig = load_image("./refs/logo_triangles/")
    T_ref, x_ref, y_ref, v_ref, c_ref = load_image("./refs/" + part)

    fig = plt.figure(figsize=(12, 4), dpi=100)

    ax1 = fig.add_subplot(131)
    ax1.axis("equal")
    ax1.set_title("Computed")
    ax1.tripcolor(x, y, triangles=T - 1, facecolors=v, cmap=c)

    ax2 = fig.add_subplot(132)
    ax2.axis("equal")
    ax2.set_title("Reference")
    ax2.tripcolor(x_ref, y_ref, triangles=T_ref - 1, facecolors=v_ref, cmap=c_ref)

    ax3 = fig.add_subplot(133)
    ax3.axis("equal")
    ax3.set_title("Original")
    ax3.tripcolor(x_orig, y_orig, triangles=T_orig - 1, facecolors=v_orig, cmap=c_orig)

    plt.show()


def plot_single_image():
    T, x, y, v, c = load_image("./out/")

    fig = plt.figure(figsize=(4, 4), dpi=100)

    ax1 = fig.add_subplot()
    ax1.axis("equal")
    ax1.set_title("Computed")
    ax1.tripcolor(x, y, triangles=T - 1, facecolors=v, cmap=c)

    plt.show()
