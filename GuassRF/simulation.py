import numpy as np

def one_d_nys(N,M,a,b,Cov,quad = "EOLE"):
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



def two_d_nys(N,n,m,lims,Cov,quad = "EOLE"):
    """
    Solver using the Nystrom method for finding the Karhunen-Loeve expansion	-----
    input
    -----
    N: The order of the Karhunen-Loeve expansion.
    n,m: n and m are the number of grid points along x and y direction respectively. 
    lims: lims=[a,b,c,d] simulation domain is [a,b] x [c,d]
    Cov: the covariance function should. Should be given as c(x,y), x and y bivariate vectors.
    quad: The quadrature method used. EOLE will be the only implemented for now 
    """
    a,b,c,d = lims # extract domain limits
    A = (b-a)*(d-c) # Omega area of the rectangular domain.
    x, y = np.linspace(a,b,n), np.linspace(a,b,m) 
    W =(A/((n-1)*(m-1)))*np.eye(n*m)
    xx = np.hstack([np.repeat(x,m).reshape(n*m,1),np.tile(y,n).reshape(n*m,1)])
    xxx = np.hstack([np.repeat(xx,n*m,axis=0),np.tile(xx,[n*m,1])])
    C = Cov(xxx[:,0:2],xxx[:,2:]).reshape(n*m,n*m) #Covariance matrix, check this.
    B = np.dot(np.dot(np.sqrt(W),C),np.sqrt(W)) # symmetric pos def B
    L,y = np.linalg.eig(B) # eigeevalues and vectors of B.
    W_inv = np.sqrt(float((n-1)*(m-1))/A)*np.eye(n*m) # invert W matrix.
    phi = np.dot(W_inv,y)
    Z = np.random.randn(N) # array of standard normal RVs.
    X = np.zeros(n*m)
    for i in range(N):
        X += np.sqrt(L[i].real)*Z[i]*phi[:,i].real
    return X.reshape(n,m), np.array(phi),L # just return eigensuite for now
#
# Execution guard
#
if __name__ == '__main__':
	pass
