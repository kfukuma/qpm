# version 1.2
import numpy as np
import quaternion

sind     = lambda th   : np.sin(np.deg2rad(th))
cosd     = lambda th   : np.cos(np.deg2rad(th))
tand     = lambda th   : np.tan(np.deg2rad(th))
asind    = lambda x    : np.rad2deg(np.arcsin(x))
acosd    = lambda x    : np.rad2deg(np.arccos(x))
atand    = lambda x    : np.rad2deg(np.arctan(x))
atan2d   = lambda y, x : np.rad2deg(np.arctan2(y, x))
di2xyz   = lambda d,i  : np.array([cosd(d)*cosd(i), sind(d)*cosd(i), sind(i)])
clph2xyz = lambda cl,ph: np.array([cosd(ph)*sind(cl), sind(ph)*sind(cl), cosd(cl)])
xyz2di   = lambda x,y,z: np.array([atan2d(y,x), asind(z)])
xyz2lph  = lambda x,y,z: np.array([90.-acosd(z), atan2d(y,x)])


def qisb(dec, inc, stk, dip):  # transformation of a block sample
    pvec = di2xyz(dec, inc)
    p = np.quaternion(*pvec)

    # transform -dip around x-axis
    theta1 = np.deg2rad(-dip)
    u1 = np.quaternion(0,1,0,0)        # x-axis [1,0,0]
    q1 = np.exp(0.5*theta1*u1)

    # transform -strike around z-axis           
    theta2 = np.deg2rad(-stk)
    u2 = np.quaternion(0,0,0,1)        # z-axis [0,0,1]
    q2 = np.exp(0.5*theta2*u2)

    q = q1*q2
    pprime = q.inverse()*p*q

    return xyz2di(pprime.x,pprime.y,pprime.z)


def qisd(dec, inc, azm, plg):   # transformation of a drill core
    pvec = di2xyz(dec, inc)
    p = np.quaternion(*pvec)

    # transform -(90-plg) around y-axis
    theta1 = np.deg2rad(-(90.-plg))
    u1 = np.quaternion(0,0,1,0)        # y-axis [0,1,0]
    q1 = np.exp(0.5*theta1*u1)

    # transform -azm around z-axis            
    theta2 = np.deg2rad(-azm)
    u2 = np.quaternion(0,0,0,1)        # z-axis [0,0,1]
    q2 = np.exp(0.5*theta2*u2)

    q = q1*q2
    pprime = q.inverse()*p*q

    return xyz2di(pprime.x,pprime.y,pprime.z)


def qisp(dec, inc):   # transformation of a piston core
    pvec = di2xyz(dec, inc)
    p = np.quaternion(*pvec)

    # transform -90 deg around y-axis
    theta = np.deg2rad(-90.)
    u = np.quaternion(0,0,1,0)        # y-axis [0,1,0]
    q = np.exp(0.5*theta*u)

    pprime = q.inverse()*p*q

    return xyz2di(pprime.x,pprime.y,pprime.z)


def qtilt(dec, inc, stkf, dipf):   # tilt-correction
    pvec = di2xyz(dec, inc)
    p = np.quaternion(*pvec)

    # rotate -dip around strike direction
    theta = np.deg2rad(-dipf)
    u = np.quaternion(0, cosd(stkf), sind(stkf), 0)
    q = np.exp(0.5*theta*u)

    pprime = q*p*q.inverse()

    return xyz2di(pprime.x,pprime.y,pprime.z)


def qvgp(dec, inc, lat, lon):  # VGP determination
    clat = 90. - lat
    mclat = atan2d(2., tand(inc))

    pvec = clph2xyz(clat,lon)
    p = np.quaternion(*pvec)

    ppvec = clph2xyz(clat-mclat,lon)
    pp = np.quaternion(*ppvec)

    # rotate -dec around site                   
    theta = np.deg2rad(-dec)
    q = np.exp(0.5*theta*p)

    pprime = q*pp*q.inverse()

    return  xyz2lph(pprime.x, pprime.y, pprime.z)


def qvgpx(slat, slon, plat, plon):  # expected direction from VGP
    cslat = 90. - slat
    cplat = 90. - plat

    psvec = clph2xyz(cslat,slon)
    ps = np.quaternion(*psvec)       # position vector of site
    ppvec = clph2xyz(cplat,plon)
    pp = np.quaternion(*ppvec)       # position vector of VGP
    northp = np.quaternion(0,0,0,1)  # north pole

    pspp = ps*pp
    psnp = ps*northp

    mclat = np.arccos(-pspp.w)         # dot product of psvec and ppvec
    incx = atand(2./np.tan(mclat))

    psppc = np.quaternion(*pspp.imag)  # cross product of psvec and ppvec
    psnpc = np.quaternion(*psnp.imag)  # cross product of psvec and north pole 
    npsppc = psppc.normalized()*psnpc.normalized()

    return  np.array([acosd(-npsppc.w), atand(2./np.tan(mclat))])


'''
# an example script to execute 
import qpm

dec = 346.8                    # qisb
inc =  68.5
stk = 216.7
dip =  76.
print("dec_i, inc_i = ", qpm.qisb(dec, inc, stk, dip))
# dec_i, inc_i = 147.8, 8.3

dec =  332.0                   # qisd
inc =  46.0
azm =  25.
plg =  53.
print("dec_i, inc_i = ", qpm.qisd(dec, inc, azm, plg))
# dec_i, inc_i = 5.5, 11.9

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

slat =  55.                    # qvgpx in deg, +, -: north, south
slon =  13.                    #       in deg, +, -: east, west or 0-360
plat =  77.3                   #       in deg, +, -: north, south
plon =  154.7                  #       in deg, +, -: east, west or 0-360
print("dec_x, inc_x = ", qpm.qvgpx(slat, slon, plat, plon))
# dec_x inc_x = 11.0, 63.0

'''