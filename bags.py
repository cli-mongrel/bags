# BAGS version 1.0.4

# Original work by Ben Trettel (http://trettel.org/)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# BAGS assumes a perfect gas (constant specific heats plus ideal gas) for
# simplicity. The errors associated with a constant specific heat are
# likely small in comparison to the errors associated with the valve
# model, for example. Also, BAGS uses a "zone" approximation for the gas
# in the barrel and gas chamber. The pressure, mass density, temperature,
# etc., are treated as uniformly distributed in space in each zone in
# this approximation. As a consequence, BAGS is inaccurate above Mach 0.4
# or so. Accuracy at moderate or higher Mach numbers requires at least
# taking account of inhomogeneities and as a consequence requires a much
# more computationally expensive approach. The simplest would be a 1D
# Lagrangian blob approximation.

from math import *
import csv

# inputs that describe the gun

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

def gasprop(prop,gas):
   # reads gas properties from gasdata.csv
   # syntax is gasprop('desired property','gas name')
   # ex. gasprop('MW','N2')

   reader = csv.reader(open("gasdata.csv", "r"), delimiter = ',')

   data = ['C_p', 'C_v', 'MW']

   item = data.index(prop) + 1

   for row in reader:
      if row[0][0] != '#':
         if row[0] == gas:
            return float(row[item])

def P_r_crit(Gamma):
   # returns the critical pressure ratio, above which flow velocity is limited to the sonic velocity
   # Gamma is the ratio of specific heats of the gas
   # See http://en.wikipedia.org/wiki/Choked_flow for discussion

   return (2. / (Gamma + 1.))**(Gamma / (Gamma - 1.))

def a(Gamma, MW, T, R_bar):
   # returns the speed of sound in m/s
   # Gamma is the ratio of specific heats of the gas
   # MW is the molecular weight of the gas in kg/kmol
   # T is the temperature of the gas in K
   # R_bar is the gas constant in J/kmol*K
   # Is the simple speed of sound equation for a perfect gas
   # See http://en.wikipedia.org/wiki/Speed_of_sound#Speed_in_ideal_gases_and_in_air

   return sqrt(Gamma * R_bar * T / MW)

def m_dot(P_in, P_out, T_in, Gamma, t, t_o, F, s, MW, r):
   # returns the mass flow through the valve in kg/s
   # P_in is the pressure of the gas stream flowing in in Pa
   # P_out is the pressure of the gas stream flowing out in Pa
   # T_in is the temperature of the gas stream flowing in in K
   # k is the ratio of specific heats of the gas
   # t is the current time in s
   # t_o is the opening time of the valve in s
   # F is the valve flow coefficient
   # s is the specific gravity of the gas
   # MW is the molecular weight of the gas in kg/kmol
   # r is the critical pressure ratio

   # Determine if the flow is sonic (which occurs if the pressure ratio is above the critical pressure ratio) and from there calculate the mass flow rate.
   # I took these equations from a book and do not know what assumptions they make.
   # Similar equations can be found at http://www.engineeringtoolbox.com/flow-coefficients-d_277.html
   # In the future I plan to make this part more accurate by using a better valve model and a force balance on the piston.

   if r * P_in <= P_out:
      flow = 8.189e-7 * F * sqrt((P_in**2. - P_out**2.) / (T_in * s))
   else:
      flow = (5.612e-7 * F * P_in) / sqrt(s * T_in)

   # if the valve is still opening, scale the flow linearly against the fraction the valve is open
   if t < t_o:
      flow = flow * (t / t_o)

   return (P_ref * flow * MW) / (R_bar * T_ref)

def V_b(x, A, V_d):
   # returns the barrel volume at projectile position x (m) in m^3
   # A is the cross-sectional area of the barrel in m^2
   # V_d is the dead volume between the valve and the barrel in m^3

   return A * x + V_d

def x_2(P_b, T_b, x_dot, P_f, Gamma, M, R_bar, P_atm, T_atm, A, m):
   # returns the current projectile acceleration in m/s^2
   # P_b is the current barrel pressure in Pa
   # T_b is the current barrel temperature in K
   # x_dot is last iteration's velocity in m/s
   # P_f is the minimum pressure to move the projectile in Pa
   # k is the ratio of specific heats of the gas
   # M is the molecular weight of the gas in kg/kmol
   # R_bar is the gas constant in J/kmol*K
   # P_atm is the atmospheric pressure in Pa
   # T_atm is the atmopsheric temperature in K
   # A is the cross-sectional area of the barrel in m^2
   # m is the projectile mass in kg
   # Uses basic Newtonian mechanics (i.e. F = m * a)

   # calculate the total pressure and friction forces

   force = A * (P_b - P_atm)
   friction = A * P_f

   # use some logic to determine if the projectile will accelerate and if so, how the force of pressure and friction com_bine

   if x_dot == 0.:
      # if the projectile is static and if the total friction force is greater than the total pressure force, then the projectile will not accelerate because the static friction force has not been exceeded

      if abs(force) > abs(friction):
         return (force - friction) / m
      else:
         return 0.
   else:
      # if the velocity is negative, the direction of the friction force is opposite of the normal direction

      if x_dot >= 0.:
         return (force - friction) / m
      else:
         return (force + friction) / m

def T_update(m_dot_in, m_cv, T_flow, T_cv, work, dt, C_v):
   # returns the current barrel or chamber gas temperature in K depending on the direction of flow neglecting the kinetic energy of the flow
   # m_dot_in is the mass flow rate into the control volume in kg/s
   # m_cv is the total mass of the control volume gas at time t - deltat in kg/s
   # T_flow is the temperature of the flowing gas in K
   # T_cv is the temperature of the control volume gas at time t - deltat in K
   # work is the work done during the last time interval in J
   # dt is the time step size in s
   # C_v is the specific heat capacity at constant volume of the gas
   # Derived from the Reynolds transport theorem
   # See http://en.wikipedia.org/wiki/Reynolds_transport_theorem

   return (m_dot_in * T_flow * dt + m_cv * T_cv) / (m_dot_in * dt + m_cv) - work / (C_v * (m_dot_in * dt + m_cv))

def PE(P_start, P_min, V_start, Gamma):
   # returns the potential energy of a gas in J
   # P_start is the maximum (starting) pressure of the gas in Pa
   # P_min is the pressure the gas will be expanded to in Pa
   # V_start is the starting volume of the gas in m^3
   # Gamma is the ratio of specific heats of the gas
   # Derived from polytropic process relationships
   # See http://en.wikipedia.org/wiki/Polytropic_process

   return (P_min * (P_start / P_min)**(1. / Gamma) * V_start - P_start * V_start) / (1. - Gamma)

def P_ideal_gas(m, T, V, R_bar, M):
   # returns the pressure of a gas according to the ideal gas law and the conditions given
   # m is the gas mass in kg
   # T is the gas temperature in K
   # V is the volume in m^3
   # R_bar is the gas constant in J/kmol*K
   # M is the molecular weight of the gas in kg/kmol
   # Derived from the mass form of the ideal gas law
   # See http://en.wikipedia.org/wiki/Ideal_gas_law

   return (m * R_bar * T) / (V * M)

# get gas properties

C_p   = gasprop('C_p', gas)
C_v   = gasprop('C_v', gas)
Gamma = round(C_p / C_v, 3)
MW    = gasprop('MW', gas)
MWair = gasprop('MW', 'air')

# calculate some shortcuts

s = MW / MWair        # calculate the specific gravity of the gas from known values
A = (pi / 4.) * d**2. # calculate the cross-sectional area of the barrel
r = P_r_crit(Gamma)   # calculate the critical pressure ratio

# initial conditions before looping

t     = 0.
x     = 0.
work  = 0.
x_dot = 0.
P_b   = P_atm
P_c   = dP + P_atm
m_b   = (P_b * V_b(x, A, V_d) * MW) / (R_bar * T_b)
m_c   = (P_c * V_c * MW) / (R_bar * T_c)

# run the loop

while x < L and t < 0.1:
   if P_b <= P_c:
      # flow into the barrel
      m_dot_in = m_dot(P_c, P_b, T_c, Gamma, t, t_o, F, s, MW, r)
      T_b = T_update(m_dot_in, m_b, T_c, T_b, work, dt, C_v)
      m_b = m_b + m_dot_in * dt
      m_c = m_c - m_dot_in * dt
   else:
      # flow into the gas chamber
      m_dot_in = m_dot(P_b, P_c, T_b, Gamma, t, t_o, F, s, MW, r)
      T_c = T_update(m_dot_in, m_c, T_b, T_c, 0, h, C_v)
      T_b = T_update(-m_dot_in, m_b, T_c, T_b, work, dt, C_v)
      m_b = m_b - m_dot_in * dt
      m_c = m_c + m_dot_in * dt

   # calculate the new chamber and barrel pressures

   P_b = P_ideal_gas(m_b, T_b, V_b(x, A, V_d), R_bar, MW)
   P_c = P_ideal_gas(m_c, T_c, V_c, R_bar, MW)

   # save the current velocity for the work calculation

   x_dot_prev = x_dot

   # Compute the acceleration, then use the forward Euler method to compute the new velocity

   x_dotdot = x_2(P_b, T_b, x_dot, P_f, Gamma, MW, R_bar, P_atm, T_atm, A, m)
   x_dot = x_dot + dt * x_dotdot

   # the increase in kinetic energy between the last step and now is the work done on the projectile
   # TODO: Account for friction

   work = m * (x_dot**2. - x_dot_prev**2.) / 2.

   # move the projectile down the barrel

   x = x + dt * x_dot

   # increment the time counter

   t = t + dt

print("V_m =", round(x_dot, 1), "m/s (muzzle velocity)")
print("eta =", str(100. * round(0.5 * m * x_dot**2. / PE(P_atm + dP, P_atm, V_c, Gamma), 3))+"% (energy efficiency)")
print("t_m =", round(t * 1.e3, 2), "ms (dwell time)")
print("T_b =", round(T_b - 273.15, 1), "C (final barrel gas temperature)")
print("a   =",   round(a(Gamma, MW, T_b, R_bar), 1), "m/s (speed of sound)")
print("Ma  =",  round(x_dot / a(Gamma, MW, T_b, R_bar), 3), "(Mach number)")
