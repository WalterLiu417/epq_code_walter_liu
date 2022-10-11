import math
import matplotlib.pyplot as plt

# Initial setup and constants
h = 0
vel = 0.00001
a = 0
r = 25
p_h = 0.08375
p_a1 = 1.25787
p_a2 = 0.999864
v_a1 = 0.0000181385
v_a2 = -3.3757*10**-10
dt = 0.01
t = 0
T = [0]
H = [0]
G = 6.67*10**-11
M_E = 5.97*10**24
vol = (4/3)*math.pi*r**3
r_E = 6378000
m_payload = 10000

while True: # Loop runs until the velocity < 0 which causes the log function to throw an error
    # Add half of the change in velocity first
    vel += 0.5*a*dt
    h += vel*dt
    p_air = p_a1 * p_a2 ** h
    v_air = v_a1 + v_a2 * h
    reynolds = (p_air * vel * r) / v_air
    try:
        drag_coeff = 4.808*10**-2*(math.log(reynolds)**2) - 1.406 *  math.log(reynolds) + 10.490
    except ValueError:
        break
    drag_force = 0.5*drag_coeff*p_air*(2*math.pi*r**2)*(vel**2)
    g = G*M_E/((r_E + h)**2)
    upthrust = p_air * vol * g
    weight = p_h * vol * g + m_payload * g
    # Newton's second law
    a = (upthrust - weight - drag_force) / (p_h * vol)
    # Second half of velocity
    vel += 0.5*a*dt
    h += vel*dt
    t += dt
    T.append(t)
    H.append(h)
    print(t, h, vel)

# Plot a height against time graph
plt.plot(T, H)
plt.show()

