import numpy as np
phi = np.pi/4
theta = 0.0
psi = 0.0

cphi_2 = np.cos(phi/2)
sphi_2 = np.sin(phi/2)
ctheta_2 = np.cos(theta/2)
stheta_2 = np.sin(theta/2)
cpsi_2 = np.cos(psi/2)
spsi_2 = np.sin(psi/2)

e0 = cphi_2*ctheta_2*cpsi_2 + sphi_2*stheta_2*spsi_2
ex = sphi_2*ctheta_2*cpsi_2 - cphi_2*stheta_2*spsi_2
ey = cphi_2*stheta_2*cpsi_2 + sphi_2*ctheta_2*spsi_2
ez = cphi_2*ctheta_2*spsi_2 - sphi_2*stheta_2*cpsi_2

ca = np.cos(10*np.pi/180)
sa = np.sin(10*np.pi/180)
cb = np.cos(15*np.pi/180)
sb = np.sin(15*np.pi/180)

L = 2000
D = 200
S = 20

Fx = -D*ca*cb -S*ca*sb + L*sa
print(Fx)
Fy = S*cb -D*sb
print(Fy)
Fz = -D*sa*cb -S*sa*sb - L*ca
print(Fz)
print(np.cos(np.pi/2))

print(np.rad2deg(np.arctan(np.tan(45*np.pi/180)*np.cos(45*np.pi/180))))
