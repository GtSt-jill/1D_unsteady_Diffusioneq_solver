import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# import os
# print(os.getcwd())

n_p = 5001 # number of nodes
n_e = n_p-1 # number of elements

frame_number = 100

fig = plt.figure()
data = np.loadtxt("u_dynamics_ad=0.50_Di=0.00.dat")
x = data[:,0]
u = data[:,1]

ims = []
for i in range(frame_number):
    im = plt.plot(x[n_p*i:n_p*(i+1)],u[n_p*i:n_p*(i+1)], color = "r")
    ims.append(im)

ani = animation.ArtistAnimation(fig, ims, interval=20)
# plt.plot(x[0:n_p],u[0:n_p], "r")
# plt.savefig("initial condition.pdf")
plt.show()

# ani.save('animation_0615_ver3.gif', writer='imagemagick')