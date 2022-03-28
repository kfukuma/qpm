import numpy as np
import quaternion

sind   = lambda x : np.sin(np.deg2rad(x))
cosd   = lambda x : np.cos(np.deg2rad(x))
tand   = lambda x : np.tan(np.deg2rad(x))
asind  = lambda y : np.rad2deg(np.arcsin(y))
acosd  = lambda y : np.rad2deg(np.arccos(y))
atan2d = lambda y, x : np.rad2deg(np.arctan2(y, x))


def qisb(dec, inc, stk, dip):  # transformation of a block sample
    pvec = np.array([cosd(dec)*cosd(inc), sind(dec)*cosd(inc), sind(inc)])
    p = np.quaternion(*pvec)

    # transform -dip around x-axis
    theta1 = np.deg2rad(-dip)
    u1 = np.quaternion(0,1,0,0)        # x-axis [1,0,0]
    q1 = np.exp(0.5*theta1*u1)

    # transform -strike around z'(Z)-axis           
    theta2 = np.deg2rad(-stk)
    u2 = np.quaternion(0,0,0,1)        # z'(Z)-axis [0,0,1]
    q2 = np.exp(0.5*theta2*u2)

    q = q1*q2
    pprime = q.inverse()*p*q

    return atan2d(pprime.imag[1], pprime.imag[0]), asind(pprime.imag[2])


def qisd(dec, inc, azm, plg):   # transformation of a drill core
    pvec = np.array([cosd(dec)*cosd(inc), sind(dec)*cosd(inc), sind(inc)])
    p = np.quaternion(*pvec)

    # transform plg around y-axis
    theta1 = np.deg2rad(plg)
    u1 = np.quaternion(0,0,1,0)        # y-axis [0,1,0]
    q1 = np.exp(0.5*theta1*u1)

    # transform -azm around z'(Z)-axis            
    theta2 = np.deg2rad(-azm)
    u2 = np.quaternion(0,0,0,1)        # z'(Z)-axis [0,0,1]
    q2 = np.exp(0.5*theta2*u2)

    q = q1*q2
    pprime = q.inverse()*p*q

    return atan2d(pprime.imag[1], pprime.imag[0]), asind(pprime.imag[2])


def qisp(dec, inc):   # transformation of a piston core
    pvec = np.array([cosd(dec)*cosd(inc), sind(dec)*cosd(inc), sind(inc)])
    p = np.quaternion(*pvec)

    # transform -90 deg around y-axis
    theta = np.deg2rad(-90.)
    u = np.quaternion(0,0,1,0)        # y-axis [0,1,0]
    q = np.exp(0.5*theta*u)

    pprime = q.inverse()*p*q

    return atan2d(pprime.imag[1], pprime.imag[0]), asind(pprime.imag[2])


def qtilt(dec, inc, stkf, dipf):   # tilt-correction
    pvec = np.array([cosd(dec)*cosd(inc), sind(dec)*cosd(inc), sind(inc)])
    p = np.quaternion(*pvec)

    # rotate -dip around strike direction
    theta = np.deg2rad(-dipf)
    u = np.quaternion(0, cosd(stkf), sind(stkf), 0)
    q = np.exp(0.5*theta*u)

    pprime = q*p*q.inverse()

    return atan2d(pprime.imag[1], pprime.imag[0]), asind(pprime.imag[2])


def qvgp(dec, inc, lat, lon):  # VGP calculation
    clat = 90. - lat
    mclat = np.arctan2(2., tand(inc))

    pvec = np.array([cosd(lon)*sind(clat), sind(lon)*sind(clat), cosd(clat)])
    p = np.quaternion(*pvec)

    # rotate magnetic colatitude around n,
    #   n: normal to the meridian plane passing site
    zc = np.exp(1j*np.deg2rad(lon))*(-1j)
    n = np.quaternion(0., np.real(zc), np.imag(zc), 0.)
    theta1 = mclat
    q1 = np.exp(0.5*theta1*n)

    # rotate -dec around site                   
    theta2 = np.deg2rad(-dec)
    q2 = np.exp(0.5*theta2*p)

    q = q2*q1
    pprime = q*p*q.inverse()

    return  90.-acosd(pprime.imag[2]), atan2d(pprime.imag[1], pprime.imag[0])


'''
# an example script to execute 
import qpm

dec = 346.8                    # qisb
inc =  68.5
stk = 216.7
dip =  76.
print("dec_i, inc_i = ", qpm.qisb(dec, inc, stk, dip))
# dec_i, inc_i = 147.8, 8.3

dec =   8.1                    # qisd
inc =  45.2
azm = 347.
plg =  27.
print("dec_i, inc_i = ", qpm.qisd(dec, inc, azm, plg))
# dec_i, inc_i = 5.3, 71.6

dec = 346.8                    # qisp
inc =  68.5
print("dec_i, inc_i = ", qpm.qisp(dec, inc))
# dec_i, inc_i = -5.1, -20.9

dec =   5.3                    # qtilt
inc =  71.6
stkf = 135.
dipf =  21.
print("dec_t, inc_t = ", qpm.qtilt(dec, inc, stkf, dipf))
# dec_t, inc_t = -74.3, 76.6

dec = 154.0                   # qvgp
inc = -58.0                    
lat =  45.5                   # in deg, +, -: north, south
lon = -73.                    # in deg, +, -: east, west or 0-360
print("lat_p, lon_p = ", qpm.qvgp(dec, inc, lat, lon))
# lat_p, lon_p =  -69.6, 6.6
'''