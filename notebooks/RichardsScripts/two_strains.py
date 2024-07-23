import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def dXdt(X,t,alpha, gamma, delta, Kab, Kba):
    Ia,Ib,Ra,Rb,Rab = X
    R0 = 1-Ra-Rb-Rab

    dIa = alpha*(R0+Rb)*Ia - delta*Ia
    dIb = alpha*(R0+Ra)*Ib - delta*Ib
    dRa = alpha*(R0*(1-Kba)*Ia - Ra*Ib) - gamma*Ra
    dRb = alpha*(R0*(1-Kab)*Ib - Rb*Ia) - gamma*Rb
    dRab = alpha*(R0*Kab*Ib + R0*Kba*Ia + Ra*Ib + Rb*Ia) - gamma*Rab

    return np.array([dIa, dIb, dRa, dRb, dRab])


if __name__=="__main__":
    Kab = 0.8
    Kba = 0.65
    delta = 1.0
    gamma = 1/20
    alpha = 3.0

    X0 = np.array([1e-2, 0, 0.05, 0, 0.5])
    T = np.linspace(0,100,300)

    Ieq = gamma/delta*(1-delta/alpha)*np.dot(np.linalg.inv([[1, Kab], [Kba, 1]]), [1,1])

    sol = odeint(dXdt, X0, T, args=(alpha, gamma, delta, Kab, Kba))

    plt.figure()
    plt.plot(T, sol[:,0])
    plt.plot(T, sol[:,1])

    plt.figure()
    plt.plot(T, sol[:,2]+sol[:,4])
    plt.plot(T, sol[:,3]+sol[:,4])
    plt.plot(T, sol[:,4])
    plt.axhline(1-delta/alpha)

    X1 = sol[-1]
    X1[1]=1e-6
    sol = odeint(dXdt, X1, T, args=(alpha, gamma, delta, Kab, Kba))

    plt.figure()
    plt.plot(T, sol[:,0])
    plt.plot(T, sol[:,1])
    plt.axhline(Ieq[0])
    plt.axhline(Ieq[1])

    plt.figure()
    plt.plot(T, sol[:,2]+sol[:,4])
    plt.plot(T, sol[:,3]+sol[:,4])
    plt.plot(T, sol[:,4])
    plt.axhline(1-delta/alpha)

    plt.figure()
    plt.plot(T, sol[:,1]/(sol[:,0] + sol[:,1]))
    plt.axhline((1-Kba)/(2-Kab-Kba))

    K=31
    vals = np.linspace(0,1,K)
    fraction = np.zeros((K,K))
    fraction_theory = np.zeros((K,K))
    T = np.linspace(0,1000,300)
    for i,Kab in enumerate(vals):
        for j,Kba in enumerate(vals):
            X0 = np.array([1e-2, 1e-2, 0.05, 0.05, 0.5])
            sol = odeint(dXdt, X1, T, args=(alpha, gamma, delta, Kab, Kba))
            fraction[i,j] = sol[-1,1]/(sol[-1,0] + sol[-1,1])
            fraction_theory[i,j] = (1-Kba)/(2-Kab-Kba)

    plt.figure()
    for i, v in enumerate(vals):
        plt.plot(1/(1+(1-vals)/(1-v)), fraction[i], 'o')
        plt.plot(1/(1+(1-vals)/(1-v)), fraction_theory[i], ls='--')

    plt.plot([0,1], [0.75, 0.25])
    x = np.linspace(0,1,101)
    plt.plot(x, 1-x-1.5*x*(x-1)*(x-0.5))
    plt.ylabel('frequency')
    plt.xlabel('(1-f)/(2-f-b)')
