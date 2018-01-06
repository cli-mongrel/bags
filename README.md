BAGS - btrettel’s air gun simulation
---

BAGS is a simulator script written by btrettel and is being maintained by yours truly.

Right now, the simulation is configred simply be editing the script with values:
```# inputs that describe the gun

gas = 'air'
F   = 2.5     # flow coefficient (typically called C_v)
d   = 13.4e-3 # barrel diameter, m
L   = 0.3048  # barrel length, m
m   = 1.3e-3  # mass of projectile, kg
dP  = 5.e5    # gauge pressure, Pa
P_f = 14.e3   # minimum pressure to move projectile, Pa
t_o = 5.e-3   # valve opening time, s
V_c = 3.3e-5  # gas chamber volume, m^3
V_d = 8.2e-6  # dead volume, m^3

# configuration

P_ref = 101325.0 # reference/standard atmospheric pressure in Pa
R_bar = 8314.472 # universal gas constant in J/kmol*K
g     = 9.80665  # gravitational acceleration in m/s^2
T_ref = 293.15   # reference temperature in K
dt    = 1.e-6    # time interval in seconds

# atmospheric conditions (pressure and temperature) available from
# http://weather.gov in the USA and a barometer there and elsewhere; see
# http://en.wikipedia.org/wiki/Inch_of_mercury to convert barometer measurements
# use T_atm = T_ref for the reference temperature (293.15 K / 20.15 C / 68 F)
# use P_atm = P_ref for the standard atmospheric pressure
# use T_c = T_atm and T_b = T_atm unless gas in gun is hot or cold

T_atm = 293.0 # atmospheric temperature in K
T_c   = T_atm # initial gas chamber temperature, K
T_b   = T_atm # initial barrel gas temperature, K
P_atm = P_ref # atmospheric temperature in Pa
```

After doing this, running the script gives the following output:
```
mon@expedit ~/bags $ python3 bags.py 
V_m = 102.0 m/s (muzzle velocity)
eta = 34.2% (energy efficiency)
t_m = 6.6 ms (dwell time)
T_b = -41.8 C (final barrel gas temperature)
a   = 304.9 m/s (speed of sound)
Ma  = 0.335 (Mach number)
```

That's it for now. Thanks to Ben for the cool work. We discussed some ideas of how to expand the script in the future, which I plan to implement.
