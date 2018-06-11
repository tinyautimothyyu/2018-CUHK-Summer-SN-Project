import itertools
# Phi is an array-valued function of rho 
# that gives the dimensionless G potetial phi. 
# rho is the dimensionless density distribution.
# rho_max (with dimension) is the maximum density from the first guess rho distribution. 
# Given R_p, R_e,and A,B=two fixed (i,j) coordinates 
def iterate(rho,rho_max,R_e):
    for i in range(2):
        phi=Phi(rho,rho_max,R_e)
        q=rho_max/b #b=1.964e6 g cm^-3
        h_0=2*(phi[A]-phi[B])/square(X[A]/R_e) #X is the meshgrid of X coordinate (with dimension length)
        F=-phi+(phi[A]-phi[B]) 
        F_max=max(F[:int(A[0]),:int(B[1])]) #restrict F_max in a rectangle enclosing the wd 
        C=(F_max-F[A]*sqrt(1+(q)**(2./3.)))/(sqrt(1+(q)**(2./3.))-1)
        H=F+C
        H_max=max(H[:int(A[0]),:int(B[1])]) #or just max(H); I dont know H values outside the wd
        if i==0:
            h_n,C_n,H_n=(h_0,C,H_max) 
        else:
            h_m,C_m,H_m=(h_0,C,H_max) # m=n+1
            H_ratio,h_ratio,C_ratio=(abs(1-(H_n/H_m)),abs(1-(h_n/h_m)),abs(1-(C_n/C_m)))
        for t in list(itertools.product(range(shape(rho)[0]),range(shape(rho)[1]))): # t is an element of Cartesian product
            if H[t]>=0:
                rho[t]=(max((1+(q)**(2./3.))*square(H[t]/H_max)-1,0)**1.5)/q
            else:
                rho[t]=0.
        rho_max=rho_max*max(rho)
        rho=rho/max(rho)
        R_e=sqrt(8.*a*sqrt(1+(rho_max/b)**(2./3.))/(b*G*rho_max*H_max)) #a=6.00e22 dynes cm^-2
    while H_ratio>=0.0001 and h_ratio>=0.0001 and C_ratio>=0.0001:
    	h_n,C_n,H_n=(h_m,C_m,H_m)
    	phi=Phi(rho,rho_max,R_e)
        q=rho_max/b 
        h_0=2*(phi[A]-phi[B])/square(X[A]/R_e)
        F=-phi+(phi[A]-phi[B])
        F_max=max(F[:int(A[0]),:int(B[1])]) 
        C=(F_max-F[A]*sqrt(1+(q)**(2./3.)))/(sqrt(1+(q)**(2./3.))-1)
        H=F+C
        H_max=max(H[:int(A[0]),:int(B[1])]) 
        h_m,C_m,H_m=(h_0,C,H_max)
            H_ratio,h_ratio,C_ratio=(abs(1-(H_n/H_m)),abs(1-(h_n/h_m)),abs(1-(C_n/C_m)))
        for t in list(itertools.product(range(shape(rho)[0]),range(shape(rho)[1]))): 
            if H[t]>=0:
                rho[t]=(max((1+(q)**(2./3.))*square(H[t]/H_max)-1,0)**1.5)/q
            else:
                rho[t]=0.
        rho_max=rho_max*max(rho)
        rho=rho/max(rho)  
        R_e=sqrt(8.*a*sqrt(1+(rho_max/b)**(2./3.))/(b*G*rho_max*H_max))
    	
    
