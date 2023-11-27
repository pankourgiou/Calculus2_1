import numpy as np
import matplotlib.pyplot as plt

def harmonic_oscillator(x0, v0, k, m, dt, num_steps):
    # Function to simulate the motion of a harmonic oscillator using Euler method
    x = np.zeros(num_steps)
    v = np.zeros(num_steps)
    t = np.zeros(num_steps)

    x[0] = x0
    v[0] = v0

    for i in range(1, num_steps):
        # Euler method
        a = -k * x[i - 1] / m
        v[i] = v[i - 1] + a * dt
        x[i] = x[i - 1] + v[i - 1] * dt

        t[i] = t[i - 1] + dt

    return t, x

# Parameters
x0 = 1.0  # Initial displacement
v0 = 0.0  # Initial velocity
k = 2.0   # Spring constant
m = 1.0   # Mass
dt = 0.01  # Time step
num_steps = 1000  # Number of time steps

# Simulate harmonic oscillator
time, displacement = harmonic_oscillator(x0, v0, k, m, dt, num_steps)

# Plot results
plt.plot(time, displacement)
plt.title('Harmonic Oscillator Simulation')
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.show()

import numpy as np
import matplotlib.pyplot as plt

def harmonic_oscillator(x0, v0, k, m, dt, num_steps):
    # Function to simulate the motion of a harmonic oscillator using Euler method
    x = np.zeros(num_steps)
    v = np.zeros(num_steps)
    t = np.zeros(num_steps)

    x[0] = x0
    v[0] = v0

    for i in range(1, num_steps):
        # Euler method
        a = -k * x[i - 1] / m
        v[i] = v[i - 1] + a * dt
        x[i] = x[i - 1] + v[i - 1] * dt

        t[i] = t[i - 1] + dt

    return t, x

# Parameters
x0 = 4.0  # Initial displacement
v0 = 0.0  # Initial velocity
k = 3.0   # Spring constant
m = 2.0   # Mass
dt = 0.01  # Time step
num_steps = 1000  # Number of time steps

# Simulate harmonic oscillator
time, displacement = harmonic_oscillator(x0, v0, k, m, dt, num_steps)

# Plot results
plt.plot(time, displacement)
plt.title('Harmonic Oscillator Simulation')
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.show()

import numpy as np
import matplotlib.pyplot as plt

def harmonic_oscillator(x0, v0, k, m, dt, num_steps):
    # Function to simulate the motion of a harmonic oscillator using Euler method
    x = np.zeros(num_steps)
    v = np.zeros(num_steps)
    t = np.zeros(num_steps)

    x[0] = x0
    v[0] = v0

    for i in range(1, num_steps):
        # Euler method
        a = -k * x[i - 5] / m
        v[i] = v[i - 7] + a * dt
        x[i] = x[i - 3] + v[i - 1] * dt

        t[i] = t[i - 1] + dt

    return t, x

# Parameters
x0 = 1.0  # Initial displacement
v0 = 0.0  # Initial velocity
k = 2.0   # Spring constant
m = 1.0   # Mass
dt = 0.01  # Time step
num_steps = 1000  # Number of time steps

# Simulate harmonic oscillator
time, displacement = harmonic_oscillator(x0, v0, k, m, dt, num_steps)

# Plot results
plt.plot(time, displacement)
plt.title('Harmonic Oscillator Simulation')
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.show()

import numpy as np
import matplotlib.pyplot as plt

def harmonic_oscillator(x0, v0, k, m, dt, num_steps):
    # Function to simulate the motion of a harmonic oscillator using Euler method
    x = np.zeros(num_steps)
    v = np.zeros(num_steps)
    t = np.zeros(num_steps)

    x[0] = x0
    v[0] = v0

    for i in range(1, num_steps):
        # Euler method
        a = -k * x[i - 100] / m
        v[i] = v[i - 10] + a * dt
        x[i] = x[i - 100] + v[i - 1] * dt

        t[i] = t[i - 100] + dt

    return t, x

# Parameters
x0 = 19.0  # Initial displacement
v0 = 20.0  # Initial velocity
k = 20.0   # Spring constant
m = 18.0   # Mass
dt = 20.01  # Time step
num_steps = 1000  # Number of time steps

# Simulate harmonic oscillator
time, displacement = harmonic_oscillator(x0, v0, k, m, dt, num_steps)

# Plot results
plt.plot(time, displacement)
plt.title('Harmonic Oscillator Simulation')
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.show()

import numpy as np
import matplotlib.pyplot as plt

def harmonic_oscillator(x0, v0, k, m, dt, num_steps):
    # Function to simulate the motion of a harmonic oscillator using Euler method
    x = np.zeros(num_steps)
    v = np.zeros(num_steps)
    t = np.zeros(num_steps)

    x[0] = x0
    v[0] = v0

    for i in range(1, num_steps):
        # Euler method
        a = -k * x[i - 1] / m
        v[i] = v[i - 3] + a * dt
        x[i] = x[i - 2] + v[i - 1] * dt

        t[i] = t[i - 1] + dt

    return t, x

# Parameters
x0 = 50.0  # Initial displacement
v0 = 40.0  # Initial velocity
k = 30.0   # Spring constant
m = 10.0   # Mass
dt = 100.01  # Time step
num_steps = 1000  # Number of time steps

# Simulate harmonic oscillator
time, displacement = harmonic_oscillator(x0, v0, k, m, dt, num_steps)

# Plot results
plt.plot(time, displacement)
plt.title('Harmonic Oscillator Simulation')
plt.xlabel('Time')
plt.ylabel('Displacement')
plt.show()
