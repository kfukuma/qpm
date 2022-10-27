# version 1.1
import numpy as np
import quaternion

sind   = lambda x : np.sin(np.deg2rad(x))
cosd   = lambda x : np.cos(np.deg2rad(x))
tand   = lambda x : np.tan(np.deg2rad(x))
asind  = lambda y : np.rad2deg(np.arcsin(y))
acosd  = lambda y : np.rad2deg(np.arccos(y))
atand  = lambda y : np.rad2deg(np.arctan(y))
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

    return atan2d(pprime.y, pprime.x), asind(pprime.z)


def qisd(dec, inc, azm, plg):   # transformation of a drill core
    pvec = np.array([cosd(dec)*cosd(inc), sind(dec)*cosd(inc), sind(inc)])
    p = np.quaternion(*pvec)

    # transform -(90-plg) around y-axis
    theta1 = np.deg2rad(-(90.-plg))
    u1 = np.quaternion(0,0,1,0)        # y-axis [0,1,0]
    q1 = np.exp(0.5*theta1*u1)

    # transform -azm around z'-axis            
    theta2 = np.deg2rad(-azm)
    u2 = np.quaternion(0,0,0,1)        # z'(Z)-axis [0,0,1]
    q2 = np.exp(0.5*theta2*u2)

    q = q1*q2
    pprime = q.inverse()*p*q

    return atan2d(pprime.y, pprime.x), asind(pprime.z)


def qisp(dec, inc):   # transformation of a piston core
    pvec = np.array([cosd(dec)*cosd(inc), sind(dec)*cosd(inc), sind(inc)])
    p = np.quaternion(*pvec)

    # transform -90 deg around y-axis
    theta = np.deg2rad(-90.)
    u = np.quaternion(0,0,1,0)        # y-axis [0,1,0]
    q = np.exp(0.5*theta*u)

    pprime = q.inverse()*p*q

    return atan2d(pprime.y, pprime.x), asind(pprime.z)


def qtilt(dec, inc, stkf, dipf):   # tilt-correction
    pvec = np.array([cosd(dec)*cosd(inc), sind(dec)*cosd(inc), sind(inc)])
    p = np.quaternion(*pvec)

    # rotate -dip around strike direction
    theta = np.deg2rad(-dipf)
    u = np.quaternion(0, cosd(stkf), sind(stkf), 0)
    q = np.exp(0.5*theta*u)

    pprime = q*p*q.inverse()

    return atan2d(pprime.y, pprime.x), asind(pprime.z)


def qsli(dec1, inc1, dec2, inc2, t):  # spherical linear interpolation
    pvec1 = np.array([cosd(dec1)*cosd(inc1), sind(dec1)*cosd(inc1),  sind(inc1)])
    p1 = np.quaternion(*pvec1)
    pvec2 = np.array([cosd(dec2)*cosd(inc2), sind(dec2)*cosd(inc2),  sind(inc2)])
    p2 = np.quaternion(*pvec2)

    pprime = np.power(p2*p1.inverse(), t)*p1

    return  atan2d(pprime.y, pprime.x), asind(pprime.z)


def qvgp(dec, inc, lat, lon):  # VGP calculation
    clat = 90. - lat
    mclat = atan2d(2., tand(inc))

    pvec = np.array([cosd(lon)*sind(clat), sind(lon)*sind(clat),  cosd(clat)])
    p = np.quaternion(*pvec)

    # slerp of site and north pole by magnetic colatitude 
    northp = np.quaternion(0,0,0,1)         # north pole
    t = mclat/clat
    ptmp = np.power(northp*p.inverse(), t)*p

    # rotate -dec around site                   
    theta = np.deg2rad(-dec)
    q = np.exp(0.5*theta*p)

    pprime = q*ptmp*q.inverse()

    return  90.-acosd(pprime.z), atan2d(pprime.y, pprime.x)


def qvgpx(slat, slon, plat, plon):  # expected direction from VGP
    cslat = 90. - slat
    cplat = 90. - plat

    psvec = np.array([cosd(slon)*sind(cslat), sind(slon)*sind(cslat),  cosd(cslat)])
    ps = np.quaternion(*psvec)       # position vector of site
    ppvec = np.array([cosd(plon)*sind(cplat), sind(plon)*sind(cplat),  cosd(cplat)])
    pp = np.quaternion(*ppvec)       # position vector of VGP
    northp = np.quaternion(0,0,0,1)  # north pole

    pspp = ps*pp
    psnp = ps*northp

    mclat = np.arccos(-pspp.w)         # dot product of psvec and ppvec
    incx = atand(2./np.tan(mclat))

    psppc = np.quaternion(*pspp.imag)  # cross product of psvec and ppvec
    psnpc = np.quaternion(*psnp.imag)  # cross product of psvec and north pole 
    npsppc = psppc.normalized()*psnpc.normalized()
    decx = acosd(-npsppc.w)

    return  acosd(-npsppc.w), atand(2./np.tan(mclat))


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

dec1 =   5.3                    # qsli
inc1 =  71.6
dec2 =   5.3 + 180.0
inc2 = -81.6     
t =   0.5                    # 0 =< t =< 1
print("dec_l, inc_l = ", qpm.qsli(dec1, inc1, dec2, inc2, t))
# dec_l, inc_l = 5.3, -13.4

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