import chemist


def h_nucleus(x, y, z):
    return chemist.Nucleus('H', 1, 1836.15, x, y, z)


def h2_nuclei():
    h0 = h_nucleus(0.0, 0.0, 0.0)
    h1 = h_nucleus(0.0, 0.0, 1.3984)

    return chemist.Nuclei(h0, h1)
