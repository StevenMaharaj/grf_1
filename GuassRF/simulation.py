import numpy as np
class KL:
    def __init__(self):
        pass

    def one_d_nys(self,N,M,a,b,Cov,quad = "EOLE"):
        """
        Karhunen-Loeve in 1-Dimension using Nystrom method.
        -----
        input
        -----
        N: Order of the Karhunen-Loeve expansion.
        M: number of quadrature intervals . N <=M
        a,b: domain of simulation, X_t for t in [a,b]
        Cov: The covariance function, a bivariate function
        quad: Quadrature used."EOLE" for the EOLE method. I tried Gauss-Legendre
        before and there was an issue with inaccurate simulation at the end
        points of the simulation domain
        -----
        output
        -----
        X: a 1-D array of the random field
        phi: a 2-D arrray whose columns are the eigenfunctions
        L: an 1-D array of the eigenvalues.
        """
        if N > M:
            raise ValueError('Order of expansion N should be less than quadrature\
        points used')
        if quad == "EOLE": # EOLE method
            x = np.linspace(a,b,M+1) # EOLE uniform grid.
            W = (1./M)*(b-a)*np.eye(M+1) #EOLE weight matrix
            x1,x2 = np.meshgrid(x,x)
            C = Cov(x1,x2) # covariance matrix
            B = np.dot(np.dot(np.sqrt(W),C),np.sqrt(W)) #symmetric B matrix.
            L,y = np.linalg.eig(B) # eigenvalues and -vectors of B
            X = np.zeros(M+1)
            W_inv = np.sqrt((float(M)/(b-a)))*np.eye(M+1) # weights matrix.
            phi = np.dot(W_inv,y) # original eigenvector problem.
            Z = np.random.randn(M+1)
            for i in range(N):
                X += Z[i]*np.sqrt(L[i])*phi[:,i]
            return X, phi, L
        else:
            raise ValueError('We only have EOLE quadrature for now.')

#
# Execution guard
#
if __name__ == '__main__':
	pass
