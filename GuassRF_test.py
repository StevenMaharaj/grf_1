import matplotlib.pyplot as plt
import numpy as np
from GuassRF.simulation import KL

N = 100 # order of the KL expansion
M = 200 # M+1 quadrature points
def Bm(t,s):
    return np.minimum(t,s)
a, b = 0., 1. # domain of simulation

kl_instance = KL()
X,phi,L = kl_instance.one_d_nys(N,M,a,b,Bm)


print(L[:10])
print([ 1./((k+0.5)**2*np.pi**2) for k in range(10)])
# plot eigenvalues: pi/L = (k-0.5)**2 for BM
L_ex = [(k+0.5)**2 for k in range(10)]
L_app = 1./(L[:10]*np.pi**2)
plt.plot(L_ex, label = "exact eigenvalues")
plt.plot(L_app,'x', label = "numerical eigenvalues")
plt.legend()
plt.ylabel(r' $\frac{1}{\lambda_k\pi^2}$')
plt.title(' Eigenvalues')
plt.savefig("BM_EV_eg.pdf")
plt.close()
t= np.linspace(a,b,M+1) # t-grid
print(np.linalg.norm(phi[:,4],2))
e = phi[:,4]
plt.plot(t, np.sqrt(2)*np.sin(4.5*np.pi*t),'x',label= "Exact")
plt.plot(t, e, label = "Numerical")
plt.title("Eigenfunction, k = {}".format(5))
plt.legend()
plt.savefig("BM_EF_eg.pdf")