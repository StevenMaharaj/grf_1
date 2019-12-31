import matplotlib.pyplot as plt
import numpy as np
from GuassRF.simulation import one_d_nys
import argparse



parser = argparse.ArgumentParser()
parser.add_argument('-t', choices=['1d','2d'], type=str, required=True)

args = parser.parse_args()

print(args.t)

# if args._get_args()=='f':
#     print('the f')
# else:
#     print(' nt the f')
# N = 100 # order of the KL expansion
# M = 200 # M+1 quadrature points
# def Bm(t,s):
#     return np.minimum(t,s)
# a, b = 0., 1. # domain of simulation


# X,phi,L = one_d_nys(N,M,a,b,Bm)


# # print(L[:10])
# # print([ 1./((k+0.5)**2*np.pi**2) for k in range(10)])
# # plot eigenvalues: pi/L = (k-0.5)**2 for BM
# L_ex = [(k+0.5)**2 for k in range(10)]
# L_app = 1./(L[:10]*np.pi**2)
# plt.plot(L_ex, label = "exact eigenvalues")
# plt.plot(L_app,'x', label = "numerical eigenvalues")
# plt.legend()
# plt.ylabel(r' $\frac{1}{\lambda_k\pi^2}$')
# plt.title(' Eigenvalues')
# plt.savefig("BM_EV_eg.pdf")
# plt.close()
# t= np.linspace(a,b,M+1) # t-grid
# # print(np.linalg.norm(phi[:,4],2))
# e = phi[:,4]
# # plt.plot(t, np.sqrt(2)*np.sin(4.5*np.pi*t),'x',label= "Exact")
# # plt.plot(t, e, label = "Numerical")
# # plt.title("Eigenfunction, k = {}".format(5))
# # plt.legend()
# # plt.savefig("BM_EF_eg.jpg")
# plt.plot(t,X)

# plt.title("Standard Brownian Motion")
# # plt.legend()
# plt.savefig("SBM2_eg.jpg")