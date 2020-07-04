import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def simple_plot(name, plot_option=False):
    data = np.genfromtxt(name)

    plt.title("Cv vs T from Dissipation Theorem")
    L = 5
    T, E, E2 = data[1:, 0], data[1:, 1], data[1:, 2]
    cv = 1 / (T ** 2 * L * L) * (E2 - E ** 2)
    plt.plot(T, cv, '--', label=str(L))

    plt.xlabel(r"$T$")
    plt.ylabel(r"$c_v$")
    plt.savefig("comparisonL5.png")

    if(plot_option):
        plt.legend()
        plt.show()


def simple_plot2(name):
    data = np.genfromtxt(name)

    L = 2.
    T, M = data[5:, 0], data[5:, 1]
    print(M)
    M = M / (L * L)
    print(M)
    plt.plot(T, M, '--', label="6")

    plt.show()


def plot_cv(systems):
    """
    This routine plot the heat capacities as a function of temperature for
      each system from the Energies moments obtained.
    """
    # Intended lables known beforehand
    labels = [2, 4, 6, 8, 16, 32, 64]

    plt.title("Cv vs T from Dissipation Theorem")
    # systems = systems[:5]
    # Going over avery system, calculate cv and plot
    for i in range(len(systems)):
        data = np.genfromtxt(systems[i])
        T, E, E2 = data[1:, 0], data[1:, 1], data[1:, 2]
        cv = 1 / (T ** 2 * labels[i] * labels[i]) * (E2 - E ** 2)
        plt.plot(T, cv, '--', label=labels[i])

    # General labels
    plt.xlabel(r"$T$")
    plt.ylabel(r"$c_v$")
    plt.xlim(T.min(), T.max())
    plt.legend()
    plt.savefig("periodic2.png")
    plt.show()


def plot_E(systems):
    """
    This routine plot the heat capacities as a function of temperature for
      each system from the Energies moments obtained.
    """
    # Intended lables known beforehand
    labels = [2, 4, 6, 8, 16, 32, 64]

    plt.title("E vs T")
    # systems = systems[:5]
    # Going over every system, calculate cv and plot
    for i in range(len(systems)):
        data = np.genfromtxt(systems[i])
        T, E = data[1:, 0], data[1:, 1]
        plt.plot(T, E / (labels[i] * labels[i]), '--', label=labels[i])

    # General labels
    plt.xlabel(r"$T$")
    plt.ylabel(r"$E$")
    plt.xlim(T.min(), T.max())
    plt.legend()
    plt.savefig("free2.png")
    plt.show()


def plot_2cv(systems):
    """
    This routine plot the heat capacities as a function of temperature for
      each system from the energies obtained and compare if this one cv were
      obtained by a derivate process respect T in the E
    """

    # Mostrly repeating plot_cv routine; only change finite differences.
    # systems2 = [L2_system, L4_system, L8_system, L16_system, L32_system,
    #             L64_system]
    labels = [2, 4, 8, 16, 32, 64]

    fig, ax = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(15, 10))
    plt.suptitle("Dissipation Theorem vs Finite Differences", size="18")
    for i in range(len(systems)):
        data = np.genfromtxt(systems[i])
        T, E, E2 = data[1:, 0], data[1:, 1], data[1:, 2]
        cv = 1 / (T ** 2 * labels[i] * labels[i]) * (E2 - E ** 2)

        # Appling finite differences for the derivate
        hN = (T[1:] - T[:-1]) * labels[i] * labels[i]  # Not uniform mapping
        cvd = (E[1:] - E[:-1]) / hN

        plt.subplot(int('23' + str(i + 1)))
        plt.title("L = " + str(labels[i]))
        plt.plot(T[:-1], cv[:-1], '--', label="Dis. T.")
        plt.plot(T[:-1], cvd, '--', label="Fin. D.")
        plt.legend()

    fig.text(0.5, 0.04, r"$T$", ha='center', size="16")
    fig.text(0.04, 0.5, r"$c_v$", va='center', rotation='vertical', size="16")
    plt.savefig("free3.png")
    plt.show()


def peaks_anaysis(systems, plot_option=False):
    """
    This routine evaluates the nature of the peaks that arise as L increase.
    """

    def Tc(L, b, a, nu):
        return b + a * L ** (-1 / nu)

    def lineal(lL, la, s):
        return la + s * lL

    T_max = np.ones(len(systems))
    Ls = np.array([2, 4, 6, 8, 16, 32, 64])
    Ls = Ls[:len(systems)]

    for i in range(len(systems)):
        data = np.genfromtxt(systems[i])
        T, E, E2 = data[1:, 0], data[1:, 1], data[1:, 2]
        cv = 1 / (T ** 2 * Ls[i] * Ls[i]) * (E2 - E ** 2)
        con_max = (cv == cv.max())
        T_max[i] = T[con_max]
        print(T[con_max])

    # per_errors = np.array([0.02, 0.05, 0.05, 0.05, 0.05, 0.05, 0.02])
    popt, pcov = curve_fit(Tc, Ls, T_max, sigma=T_max * 0.005,
                           absolute_sigma=True)

    if(plot_option):
        plt.title("Tc vc L")
        plt.plot(Ls[:len(systems)], T_max, '.')
        Ls_con = np.linspace(Ls.min(), Ls.max())
        plt.plot(Ls_con, Tc(Ls_con, *popt), '-')
        plt.ylabel("Tc")
        plt.xlabel("L")
        plt.savefig("free5.png")
        plt.show()

    print("Parameters: Tc %.2f, nu; %.2f" % (popt[0], popt[2]))
    print("Covariance Matrix:", pcov)
    print(popt[2], pcov[2][2] ** 0.5)
    print(popt[0], pcov[0][0] ** 0.5)

    print("---------------------------------------------------------\n\n")

    return popt[2], pcov[2][2] ** 0.5


def Critical_exponents_cv(systems, systemsM, nu, nu_err):
    cv_max = np.ones(len(systems))
    M_max = np.ones(len(systems))
    Ls = np.array([2, 4, 6, 8, 16, 32, 64])
    Ls = Ls[:len(systems)]
    T_max = np.ones(len(systems))

    def Power_law(L, a, s):
        return abs(a) * L ** s

    def lineal(L, la, s):
        return la + s * L

    for i in range(len(systems)):
        # Calculate Cv and retrieve max
        data = np.genfromtxt(systems[i])
        T, E, E2 = data[1:, 0], data[1:, 1], data[1:, 2]
        cv = 1 / (T ** 2 * Ls[i] * Ls[i]) * (E2 - E ** 2)
        i_max = np.argmax(cv_max)
        cv_max[i] = cv[i_max]
        T_max[i] = T[i_max]

        # Calculate M and retrive max
        dataM = np.genfromtxt(systemsM[i])
        TM, M = dataM[1:, 0], dataM[1:, 1]
        h = TM[1] - TM[0]
        con_max = (TM < T[i_max] + h / 2) * (TM > T[i_max] - h / 2)
        M_max[i] = M[con_max] / (Ls[i] * Ls[i])

    # popt, pcov = curve_fit(Power_law, Ls / T_max, cv_max)
    # plt.plot(Ls, cv_max, '.')
    # Ls_con = np.linspace(Ls.min(), Ls.max())
    # plt.plot(Ls_con, Power_law(Ls_con, *popt), '-')

    # s, s_err = popt[1], pcov[1][1] ** 0.5
    # alpha_err = ((s * nu_err) ** 2 + (nu * s_err)**2)**0.5

    # print(nu)
    # print("---------------------------------------------------------")
    # print("alpha: %.2f +/- %.2f" % (s * nu, alpha_err))
    # print("---------------------------------------------------------\n")
    # plt.show()

    # popt, pcov = curve_fit(Power_law, Ls / T_max, M_max)
    # plt.plot(Ls[:-1], M_max[:-1], '.')
    # Ls_con = np.linspace(Ls.min(), Ls.max())
    # plt.plot(Ls_con, Power_law(Ls_con, *popt), '-')

    # print(nu)
    # s, s_err = popt[1], pcov[1][1] ** 0.5
    # beta_err = ((s * nu_err) ** 2 + (nu * s_err)**2)**0.5
    # print("---------------------------------------------------------")
    # print("beta: %.3f +/- %.3f" % (-s * nu, beta_err))
    # print("---------------------------------------------------------\n")
    # plt.show()

    popt, pcov = curve_fit(lineal, np.log(Ls[1:] / T_max[1:]), np.log(cv_max[1:]))

    plt.title("ln(cv) vs ln(L)")
    print("values:", popt[1] * nu, pcov[0][0]**0.5 * nu + nu_err * popt[1])
    plt.plot(np.log(Ls[:]), np.log(cv_max[:]), '.')
    Ls_con = np.linspace(Ls.min(), Ls.max())
    plt.plot(np.log(Ls_con), lineal(np.log(Ls_con), *popt), '-')
    plt.xlabel("ln(L)")
    plt.ylabel("ln(cv)")
    plt.savefig("Free6.png")
    plt.show()

    plt.title("ln(M) vs ln(L)")
    popt, pcov = curve_fit(lineal, np.log(Ls[1:]), np.log(M_max[1:]))
    print("values2:  ", popt[1] * nu, pcov[0][0]**0.5 * nu + nu_err * popt[1])
    plt.plot(np.log(Ls[:] / T_max[:]), np.log(M_max[:]), '.')
    Ls_con = np.linspace(Ls.min(), Ls.max())
    plt.plot(np.log(Ls_con), lineal(np.log(Ls_con), *popt), '-')
    plt.xlabel("ln(L)")
    plt.ylabel("ln(M)")
    plt.savefig("Free7.png")
    plt.show()


def plot_M(systems):
    """
    This routine plot the heat capacities as a function of temperature for
      each system from the Energies moments obtained.
    """
    # Intended lables known beforehand
    labels = [2, 4, 6, 8, 16, 32, 64]

    plt.title("M vs T")
    # systems = systems[:5]
    # Going over avery system, calculate cv and plot
    for i in range(len(systems)):
        data = np.genfromtxt(systems[i])
        T, M = data[12:, 0], data[12:, 1]
        plt.plot(T, M / (labels[i] * labels[i]), '--', label=labels[i])

    # General labels
    plt.xlabel(r"$T$")
    plt.ylabel(r"$M$")
    plt.legend()
    plt.savefig("free4.png")
    plt.show()


################################################
# Heat Capacity Analisys #######################
################################################
# directory = "Periodic_data/"
# name = directory + "EmvsTn_steps80000000L4n_points50.txt"
# name = "Free_data/EmvsTn_steps1000000L16n_points50.txt"
# name = "Free_data/EmvsTn_steps4000000L8n_points50.txt"
# name = "Periodic_data/EmvsTn_steps1000001L2n_points50.txt"
# simple_plot(name)


################################################
# Energy plot ##################################
################################################
# Making list
# systems = np.genfromtxt("periodic_names640Mp50.txt", dtype="str")
# plot_E(systems[:])
# systems = np.genfromtxt("free_names160Mp50.txt", dtype="str")
# plot_E(systems[:])


#################################
# Graphs of Cv ##################
#################################
# Making list of periodic systems
# systems = np.genfromtxt("periodic_names640Mp50.txt", dtype="str")
# plot_cv(systems[:])
# systems = np.genfromtxt("free_names160Mp50.txt", dtype="str")
# plot_cv(systems[:])


#################################
# Graphs of Cv with finite dif. #
################################
# systems = np.genfromtxt("periodic_names640Mp50.txt", dtype="str")
# systems = np.concatenate((systems[:2], systems[3:]))
# plot_2cv(systems)
# systems = np.genfromtxt("free_names160Mp50.txt", dtype="str")
# systems = np.concatenate((systems[:2], systems[3:]))
# plot_2cv(systems)


#########################################################
# Analyzing critic temperatures  and Critical exponents #
#########################################################
# systems = np.genfromtxt("periodic_names640Mp50.txt", dtype="str")
# systems2 = np.genfromtxt("periodic_namesM160Mp50.txt", dtype="str")
# nu, nu_err = peaks_anaysis(systems[:])
# Critical_exponents_cv(systems, systems2, nu, nu_err)

# systems = np.genfromtxt("free_names160Mp50.txt", dtype="str")
# systems2 = np.genfromtxt("free_namesM128Mp70.txt", dtype="str")
# nu, nu_err = peaks_anaysis(systems)
# Critical_exponents_cv(systems, systems2, nu, nu_err)


################################################
# Magnetization Analisys #######################
# ################################################
# name1 = "Free_data/MmvsTn_steps2000000L2n_points70.txt"
# name2 = "Free_data/MmvsTn_steps4000000L4n_points70.txt"
# name3 = "Free_data/MmvsTn_steps6000000L6n_points70.txt"
# name4 = "Free_data/MmvsTn_steps16000000L8n_points70.txt"
# name5 = "Free_data/MmvsTn_steps64000000L16n_points70.txt"
# name6 = "Free_data/MmvsTn_steps128000000L32n_points70.txt"
# print(name1)
# print(name2)
# print(name3)
# print(name4)  # from 5
# print(name5)  # from 5
# print(name6)

# name7 = "Free_data/MmvsTn_steps480000000L64n_points70.txt"
# print(name7)
# simple_plot2(name7)


################################################
# Magnetization plots ##########################
################################################
# systems = np.genfromtxt("periodic_namesM160Mp50.txt", dtype="str")
# plot_M(systems)
# systems = np.genfromtxt("free_namesM128Mp70.txt", dtype="str")
# plot_M(systems)


# Same system different seeds
# name = "Periodic_data/EmvsTn_steps10000000L5n_points50.txt"
# simple_plot(name)
# name = "Periodic_data/EmvsTn_steps10000001L5n_points50.txt"
# simple_plot(name, True)


################################
# Analisis of Markov processes #
################################


plt.title("Magnetization vs steps")
M = np.array([9, 9, 9, 9, 9, 11, 13, 15, 15, 17, 19, 19, 19, 19, 21, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25])
M2 = np.array([-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 5, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 9, 9, 9, 11, 11, 11, 11, 13, 13, 13, 13, 15, 15, 15, 15, 15, 15, 15, 15, 15, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 19, 19, 19, 21, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25])
plt.plot(M)
plt.plot(M2)
plt.xlabel("steps")
plt.savefig("Magnetization.png")
plt.show()
