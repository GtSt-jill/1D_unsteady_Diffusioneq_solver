import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# import os
# print(os.getcwd())

n_p=1001 # number of nodes
n_e=n_p-1 # number of elements

frame_number = 100

fig = plt.figure()
data=np.loadtxt("u_dynamics_ad=0.50_Di=0.00.dat")
x=data[:,0]
u=data[:,1]

ims = []
for i in range(frame_number):
    im = plt.plot(x[n_p*i:n_p*(i+1)],u[n_p*i:n_p*(i+1)],"r")
    ims.append(im)

ani = animation.ArtistAnimation(fig, ims, interval=frame_number)

plt.show()

# ani.save('animation_tmp.gif', writer='imagemagick')