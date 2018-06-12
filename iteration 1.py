import itertools
# Phi is an array-valued function of rho
# that gives the dimensionless G potetial phi.
# rho is the dimensionless density distribution.
# rho_max (with dimension) is the maximum density from the first guess rho distribution.
# A,B=two fixed (r,mu) coordinates. B is on the rotation axis.
def iterate(rho,rho_max,R_e,A,B):
    R_max=R_e
    for i in range(2):
        phi=np.zeros(KDIV*NDIV).reshape((NDIV,KDIV))
        for v in range(NDIV):
            for w in range(KDIV):
                phi[v,w]=PHI(v,w)
        phi=phi/(rho_max*G*R_e**2)
        q=rho_max/b #b=1.964e6 g cm^-3
        h_0=2*(phi[A]-phi[B])/square(R_max/R_e)
        F=-phi+(phi[A]-phi[B])
        F_max=max(F)
        C=(F_max-F[A]*sqrt(1+(q)**(2./3.)))/(sqrt(1+(q)**(2./3.))-1)
        H=F+C
        H_max=max(H)
        if i==0:
            h_n,C_n,H_n=(h_0,C,H_max)
        else:
            h_m,C_m,H_m=(h_0,C,H_max) # m=n+1
            H_ratio,h_ratio,C_ratio=(abs(1-(H_n/H_m)),abs(1-(h_n/h_m)),abs(1-(C_n/C_m)))
        for t in list(itertools.product(NDIV,KDIV)): # t is an element of Cartesian product
            if H[t]>=0:
                rho[t]=(max((1+(q)**(2./3.))*square(H[t]/H_max)-1,0)**1.5)/q
            else:
                rho[t]=0.
        rho_max=rho_max*max(rho)
        rho=rho/max(rho)
        R_e=sqrt(8.*a*sqrt(1+(rho_max/b)**(2./3.))/(b*G*rho_max*H_max)) #a=6.00e22 dynes cm^-2
    while H_ratio>=0.0001 and h_ratio>=0.0001 and C_ratio>=0.0001:
    	h_n,C_n,H_n=(h_m,C_m,H_m)
    	phi=np.zeros(KDIV*NDIV).reshape((NDIV,KDIV))
        for i in range(NDIV):
            for j in range(KDIV):
                phi[i,j]=PHI(i,j)
        phi=phi/(rho_max*G*R_e**2)
        q=rho_max/b
        h_0=2*(phi[A]-phi[B])/square(R_max/R_e)
        F=-phi+(phi[A]-phi[B])
        F_max=max(F)
        C=(F_max-F[A]*sqrt(1+(q)**(2./3.)))/(sqrt(1+(q)**(2./3.))-1)
        H=F+C
        H_max=max(H)
        h_m,C_m,H_m=(h_0,C,H_max)
            H_ratio,h_ratio,C_ratio=(abs(1-(H_n/H_m)),abs(1-(h_n/h_m)),abs(1-(C_n/C_m)))
        for t in list(itertools.product(NDIV,KDIV)):
            if H[t]>=0:
                rho[t]=(max((1+(q)**(2./3.))*square(H[t]/H_max)-1,0)**1.5)/q
            else:
                rho[t]=0.
        rho_max=rho_max*max(rho)
        rho=rho/max(rho)
        R_e=sqrt(8.*a*sqrt(1+(rho_max/b)**(2./3.))/(b*G*rho_max*H_max))

    
