from numpy import log10, sqrt, abs


def f_haaland(Re, D, e):
    return (1 / (-1.8 * log10(((e / D) / 3.7) ** 1.11 + 6.9 / Re))) ** 2


def f_colebrook(Re, D, e, prec=1e-8):
    f_aux = f_haaland(Re, D, e)  # Primeira aproximação para f

    def f_colebrook_step(f, Re, D, e):
        return (1 / (-2 * log10(e / (3.7 * D) + 2.51 / (Re * sqrt(f))))) ** 2

    f_next = f_colebrook_step(f_aux, Re, D, e)

    while abs(f_next - f_aux) > prec:
        f_aux = f_next
        f_next = f_colebrook_step(f_aux, Re, D, e)

    return f_next


def reynolds(rho, V, D, mu):
    return rho * V * D / mu
