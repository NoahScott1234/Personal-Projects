
import matplotlib.pyplot as plt
import math
G = 6.67e-11
M = 5.972e24
x = 6371000.0 + float(input('Enter Height Above Earth: '))
y = 0
v0 = math.sqrt(G*M/x)
vx = 0
vy = v0
dt = 0.01
t = 0
r = (x**2 + y**2)**0.5
T = 2*math.pi*math.sqrt(r**3/(G*M))
x_val = []
y_val = []
r_val = []
t_val = []
r_earth = 6371000.0
while t < T: 
    ax= -G*M*x/r**3
    ay = -G*M*y/r**3
    vx = vx + ax * dt
    vy = vy + ay * dt
    x = x + vx * dt
    y = y + vy * dt
    r = (x**2 + y**2)**0.5
    x_val.append(x)
    y_val.append(y)
    r_val.append(r)
    t_val.append(t)
    t += dt
print(f'Time Taken to Complete Orbit: {t:.2f} seconds')
print(f'Velocity of object: {v0:.2f} m/s')

fig, axs = plt.subplots(2,1)

axs[0].plot(x_val, y_val)
axs[0].set_title("Orbit")

axs[1].plot(t_val, r_val)
axs[1].set_title("Radius vs Time")

plt.show()



