import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
get_col = lambda col: (line.split(' ')[col-1] for line in open('Star_WENO_CentralDensity_1.dat')) #
T=[] #Time
D=[] #Central Density
for i in get_col(14):
    T.append(float(i))

for i in get_col(27):
    D.append(float(i))

# Number of samplepoints
N = len(T)
# sample spacing
t=T[-1]/N
C=D[0] #Initial central rho
D=np.array(D)-C
T=np.array(T)

yf = scipy.fftpack.fft(D)
xf = np.linspace(0.0, 1.0/(2.0*t), N/2)

Y=2.0/N * np.abs(yf[:N//2])
fig, ax = plt.subplots()
ax.plot(T,D)
#ax.plot(xf, Y)
plt.show()
mx=[] #all maxima 
for i in range(len(Y)):
    if Y[i]==max(Y):
        x1=xf[i]
        mx.append([xf[i],Y[i]])
    elif i<len(Y)-1 and i>0:
        if Y[i]>Y[i+1] and Y[i]>Y[i-1]:
            mx.append([xf[i],Y[i]])

print('The initial density is '+str(C))
print('The frequency is '+str(x1))