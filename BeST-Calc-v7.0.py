#this code is for surface energy calculation of Butler equation
#Considering Excess Molar Volume (XA*VA+XB*VB(Vi=a+b*T+c*T**2) or XA*VA+XB*VB+XA*XB*Vhi(Vi=V0*exp(VA)))
#Considering the Size Dependent Surface Tension for pure Component (Xiong's model(Phys. Chem. Chem. Phys., 2011, 13, 10648–10651), or Nanda's model(Phys. Lett. A, 2012, 376(19), 1647-1649), or Jiang'model(J. Phys. Chem. B 2004, 108, 18, 5617–5619))
#Date: 20210524
#Date: 20210625
#update: add the user-defined composition step; modify the x2val to 0.01; add the size effect of the surface tension for pure elements (contribution to the Gibbs energy of pure element);
#Date: 20210702
#update: add the exception of None type; add the surface Gibbs energy;
#Date: 20210819
#update: adjust the code of dichotomy;
#Date: 20211115
#upate: pste_r[j] should not consider the excess term effect on the molar volume and surface tension (especially for the molar volume, since there is not excess term in surfaca tension in the present database), but it does not influence these which don't have the excess term.
#Date: 20220214
#update: add the judgement of surface energy (not negative value);
#Date: 20220815
#update: The surface Gibbs energy, SGB[j], should multiply 10**9 since the unit of radius is nanometer.

import numpy as np
import csv

#necessary input variables
R = 8.31451
Nav = 6.022*10**(23)
#which phase to calculate
ph = input("Phase: ")
print('\n')
#elements
el1 = input("The first element in the solution model: ")
print('\n')
el2 = input("The second element in the solution model: ")
print('\n')
#shape factor
c = float(eval(input("Shape factor for Gsurf=2*C*σ*V/r: ")))
print('\n')
#packing fraction
pf = float(eval(input('The 2-D interfacial packing fraction of the interface, f(Ai=f*(Nav**1/3)*(Vi**2/3), f usually is 1.091): ')))
print('\n')
#ratio of coordination number
ra = float(eval(input('The ratio of the coordination number in the surface to that in the bulk, βmix(Gsurf,ex=βmix*Gbulk,ex, βmix usually is 0.85 for liquid and 0.84 for solid phase): ')))
print('\n')
#Interaction parameters
a = [None]*5
b = [None]*5
#The maximum order of interation parameters
od = int(eval(input("The maximum order interaction parameters(do not excess 4): ")))
print('\n')
if od <= 4:
    for k in range(0, od + 1, 1):
        print('The ' + str(k) + ' order interaction parameters, a+b*T: ')
        a[k]=float(eval(input('a' + str(k) + '= ')))
        b[k]=float(eval(input('b' + str(k) + '= ')))
        print('\n')
    for t in range(k+1, 5, 1):
        a[t] = 0
        b[t] = 0
else:
    print('Incorrect order for interation parameters. \nThe maximum order should be among 0 to 4')
    print('\n')
    input("Press <enter> to close program")
    exit()

#size dependent model for surface tension of pure component
print('This code has three models for Size dependent Surface Tension of pure component, σ=a+b*T:')
print("Jiang'model(J. Phys. Chem. B 2004, 108, 18, 5617–5619): σ(r)=σ*(1-2*hi/r), hi: atomic diameter of element i, r: particle radius")
print("Xiong's model(Phys. Chem. Chem. Phys., 2011, 13, 10648–10651): σ(r)=σ*(1-0.725*hi/r), hi: atomic diameter of element i, r: particle radius")
print("Nanda's model(Phys. Lett. A, 2012, 376(19), 1647-1649): σ(r)=σ*(1-2*c*βi/r), βi: constant of bulk element i, c: shape factor, unity for spherical nanoparticles, r: particle radius")
print('The surface tension for pure ' + el1 + ': ')
pste1a = float(eval(input('The constant term, a = ')))
pste1t1 = float(eval(input('The first order term on Temperature, b = ')))
print('\n')
print('The surface tension for pure ' + el2 + ': ')
pste2a = float(eval(input('The constant term, a = ')))
pste2t1 = float(eval(input('The first order term on Temperature, b = ')))
print('\n')
#This flag determines to consider the size dependent surface tension for alloy or pure component or not
sd = input('Considering the Size Dependent Surface Tension[y/n]: ')
if sd in ['y', 'Y']:
    sds = input('Size dependent Surface Tension for[A: pure element/B: binary alloy]: ')
    sdm = input('Which model[Jiang/Xiong/Nanda]: ')
    if sdm == 'Jiang':
        hel1 = float(eval(input("Atomic diameter (nm) of " + el1 + ": ")))
        print('\n')
        hel2 = float(eval(input("Atomic diameter (nm) of " + el2 + ": ")))
        print('\n')
    elif sdm == 'Xiong':
        hel1 = float(eval(input("Atomic diameter (nm) of " + el1 + ": ")))
        print('\n')
        hel2 = float(eval(input("Atomic diameter (nm) of " + el2 + ": ")))
        print('\n')
    elif sdm == 'Nanda':
        sf = float(eval(input('The factor for different shape particle of β, c(unity for spherical nanoparticle): ')))
        hel1 = float(eval(input("β of " + el1 + ": ")))
        print('\n')
        hel2 = float(eval(input("β of " + el2 + ": ")))
        print('\n')
    else:
        print('Invalid input')
        print('\n')
        input("Press <enter> to close program")
        exit()
elif sd in ['n', 'N']:
    print('\n')
else:
    print('Invalid input')
    print('\n')
    input("Press <enter> to close program")
    exit()

#This flag determines to take the excess molar volume into account or not
print('This code has two models to describe the molar volume:')
print('Negelect the excess molar volume: VAB=XA*VA+XB*VB, Vi=a+b*T+c*T**2')
print('Considering the excess molar volume: VAB=XA*VA+XB*VB+XA*XB*[Vij0+(XA-XB)*Vij1+...], VA=V0*exp(VA), Vij=V0*exp(VA)')
em = input('Considering the Excess Molar Volume[y/n]:')
if em in ['n', 'N']:
    print('The molar volume for pure ' + el1 + ': ')
    mve1a = float(eval(input('The constant term, a = ')))
    mve1t1 = float(eval(input('The first order term on Temperature, b = ')))
    mve1t2 = float(eval(input('The second order term on Temperature, c = ')))
    print('\n')
    print('The molar volume for pure ' + el2 + ': ')
    mve2a = float(eval(input('The constant term, a = ')))
    mve2t1 = float(eval(input('The first order term on Temperature, b = ')))
    mve2t2 = float(eval(input('The second order term on Temperature, c = ')))
    print('\n')
elif em in ['y', 'Y']:
    print('The molar volume for pure ' + el1 + ': ')
    mve1V0 = float(eval(input('V0(constant, m3/mol): ')))
    print('VA(a+b*T+c*T^2+d*T^3+e*T^4+f*T^-1): ')
    mve1VAa = float(eval(input('a= ')))
    mve1VAb = float(eval(input('b= ')))
    mve1VAc = float(eval(input('c= ')))
    mve1VAd = float(eval(input('d= ')))
    mve1VAe = float(eval(input('e= ')))
    mve1VAf = float(eval(input('f= ')))
    print('\n')
    print('The molar volume for pure ' + el2 + ': ')
    mve2V0 = float(eval(input('V0(constant, m3/mol): ')))
    print('VA(a+b*T+c*T^2+d*T^3+e*T^4+f*T^-1): ')
    mve2VAa = float(eval(input('a= ')))
    mve2VAb = float(eval(input('b= ')))
    mve2VAc = float(eval(input('c= ')))
    mve2VAd = float(eval(input('d= ')))
    mve2VAe = float(eval(input('e= ')))
    mve2VAf = float(eval(input('f= ')))
    print('\n')
    ov = int(eval(input("The maximum order interaction parameters for excess molar volume(do not excess 4): ")))
    print('\n')
    vijV0 = [None]*5
    vijVAa = [None]*5
    vijVAb = [None]*5
    vijVAc = [None]*5
    vijVAd = [None]*5
    vijVAe = [None]*5
    vijVAf = [None]*5
    if ov <= 4:
        for k in range(0, ov + 1, 1):
            print('The ' + str(k) + ' order interaction parameters for excess molar volume')        
            vijV0[k]=float(eval(input('V0(constant, m3/mol)= ')))
            print('VA(a+b*T+c*T^2+d*T^3+e*T^4+f*T^-1): ')
            vijVAa[k]=float(eval(input('VAa= ')))
            vijVAb[k]=float(eval(input('VAb= ')))
            vijVAc[k]=float(eval(input('VAc= ')))
            vijVAd[k]=float(eval(input('VAd= ')))
            vijVAe[k]=float(eval(input('VAe= ')))
            vijVAf[k]=float(eval(input('VAf= ')))
            print('\n')
        for t in range(k+1, 5, 1):
            vijV0[t] = 0
            vijVAa[t] = 0
            vijVAb[t] = 0
            vijVAc[t] = 0
            vijVAd[t] = 0
            vijVAe[t] = 0
            vijVAf[t] = 0
    else:
        print('Incorrect order for interation parameters of excess molar volume. \nThe maximum order should be among 0 to 4')
        print('\n')
        input("Press <enter> to close program")
        exit()
else:
    print('Incorrect order for interation parameters of excess molar volume. \nThe maximum order should be among 0 to 4')
    print('\n')
    input("Press <enter> to close program")
    exit()

#The calculation radius range
Rmin = float(eval(input("The minimum radius, nm: ")))
print('\n')
Rmax = float(eval(input("The maximum radius, nm: ")))
print('\n')
Rsp = float(eval(input("The calculation step for radius, nm: ")))
print('\n')

#The calculation temperature range
Tmin = float(eval(input("The minimum Temperature, K: ")))
print('\n')
Tmax = float(eval(input("The maximum Temperature, K: ")))
print('\n')
Tsp = float(eval(input("The calculation step for Temperature, K: ")))
print('\n')

#The calculation composition range
Xmin = float(eval(input("The minimum x("+ el2 + ") (Not be 0): ")))
print('\n')
Xmax = float(eval(input("The maximum x("+ el2 + ") (Not be 1): ")))
print('\n')
Xsp = float(eval(input("The calculation step for composition: ")))
print('\n')

print("**********Let's start calculating!**********")
print("**********It may take some minutes!**********")
print("**********Please wait!**********")
print('\n')

#define some functions
#define function for partial molar excess Gibbs energy of element A
def PMEGE1(xb, T, L0val, L1val, L2val, L3val, L4val):
    PGE1 = xb**2*(L0val+(3-4*xb)*L1val+L2val*(5-16*xb+12*xb**2)+L3val*(7-36*xb+60*xb**2-32*xb**3)+L4val*(9-64*xb+168*xb**2-192*xb**3+80*xb**4))
    return PGE1
#define function for partial molar excess Gibbs energy of element B
def PMEGE2(xb, T, L0val, L1val, L2val, L3val, L4val):
    PGE2 = (1-xb)**2*(L0val+L1val*(1-4*xb)+L2val*(1-8*xb+12*xb**2)+L3val*(1-12*xb+36*xb**2-32*xb**3)+L4val*(1-16*xb+72*xb**2-128*xb**3+80*xb**4))
    return PGE2
#define the function of molar surface area
def PMIA(f, MV):
    pmia = f*Nav**(1/3)*MV**(2/3)
    return pmia
#define the function of alloy surface tension from element A
def DABA(da, miaa, T, sxa, bxa, bpga, spga):
    daba = da + (R*T/miaa)*np.log(sxa/bxa)+(spga-bpga)/miaa
    return daba
#define the function of alloy surface tension from element B
def DABB(db, miab, T, sxb, bxb, bpgb, spgb):
    dabb = db + (R*T/miab)*np.log(sxb/bxb)+(spgb-bpgb)/miab
    return dabb
#size dependent surface tension for pure element in Jiang's model
def SIGMA_JIANG(r,hi,sigma):
    sigma_r = sigma*(1-2*hi/r)
    return sigma_r
#size dependent surface tension for pure element in Xiong's model
def SIGMA_XIONG(r,hi,sigma):
    sigma_r = sigma*(1-0.725*hi/r)
    return sigma_r
#size dependent surface tension for pure element in Nanda's model
def SIGMA_NANDA(r,c,bi,sigma):
    sigma_r = sigma*(1-2*c*bi/r)
    return sigma_r
#define the funciton for Excess molar volume
#define the function of molar volume for pure component, VI
def VI(T, V0, VAa, VAb, VAc, VAd, VAe, VAf):
    VI = V0*np.exp(VAa + VAb*T + VAc*T**2 + VAd*T**3 + VAe*T**4 + VAf*T**(-1))
    return VI
#Interaction parameter for excess molar volume
def VIJ(T, V0, VAa, VAb, VAc, VAd, VAe, VAf):
    VIJ = V0*np.exp(VAa + VAb*T + VAc*T**2 + VAd*T**3 + VAe*T**4 + VAf*T**(-1))
    return VIJ
#define the partial molar volume, Vhi=Vh+(1-xi)*aVh/axi, element A
def PMVE1(xb, T, V1, Vij0, Vij1, Vij2, Vij3, Vij4):
    PVE1 = V1+xb**2*(Vij0+(3-4*xb)*Vij1+Vij2*(5-16*xb+12*xb**2)+Vij3*(7-36*xb+60*xb**2-32*xb**3)+Vij4*(9-64*xb+168*xb**2-192*xb**3+80*xb**4))
    return PVE1
#define the partial molar volume, Vhi=Vh+(1-xi)*aVh/axi, element B
def PMVE2(xb, T, V2, Vij0, Vij1, Vij2, Vij3, Vij4):
    PVE2 = V2+(1-xb)**2*(Vij0+Vij1*(1-4*xb)+Vij2*(1-8*xb+12*xb**2)+Vij3*(1-12*xb+36*xb**2-32*xb**3)+Vij4*(1-16*xb+72*xb**2-128*xb**3+80*xb**4))
    return PVE2

n = 1E10
if sd in ['y', 'Y']:
    f = open(el1.upper() + el2.upper() + '-' + ph.upper() + '-SurfaceEnergy-ButlerEquation' + '-SizeDependent-' + sd.upper() + '-' + sds.upper() + '-' + sdm.upper() + '-ExcessMolarVolume-' + em.upper() + '.csv', 'w', newline='', encoding='utf-8')
else:
    f = open(el1.upper() + el2.upper() + '-' + ph.upper() + '-SurfaceEnergy-ButlerEquation' + '-SizeDependent-N-ExcessMolarVolume-' + em.upper() + '.csv', 'w', newline='', encoding='utf-8')
csv_writer = csv.writer(f)
csv_writer.writerow(['System', el1 + el2])
csv_writer.writerow(['Phase', ph])
r = Rmin
while r <= Rmax:
    x2val = Xmin
    if sd in ['y', 'Y']:
        if sds in ['a', 'A']:
            if sdm == 'Jiang':
                pste1a_r = SIGMA_JIANG(r,hel1,pste1a)
                pste1t1_r = SIGMA_JIANG(r,hel1,pste1t1)
                pste2a_r = SIGMA_JIANG(r,hel2,pste2a)
                pste2t1_r = SIGMA_JIANG(r,hel2,pste2t1)
            elif sdm == 'Xiong':
                pste1a_r = SIGMA_XIONG(r,hel1,pste1a)
                pste1t1_r = SIGMA_XIONG(r,hel1,pste1t1)
                pste2a_r = SIGMA_XIONG(r,hel2,pste2a)
                pste2t1_r = SIGMA_XIONG(r,hel2,pste2t1)
            else:
                pste1a_r = SIGMA_NANDA(r, sf, hel1,pste1a)
                pste1t1_r = SIGMA_NANDA(r, sf, hel1,pste1t1)
                pste2a_r = SIGMA_NANDA(r, sf, hel2,pste2a)
                pste2t1_r = SIGMA_NANDA(r, sf, hel2,pste2t1)
        else:
            pass
    else:
        pass
    while x2val <= Xmax:
        x1val = 1-x2val
        ds0 = 1
        ds1 = x1val-x2val
        ds2 = ds1**2
        ds3 = ds1**3
        ds4 = ds1**4
        Tval = Tmin 
        T = []
        pge1 = []
        pge2 = []
        i = 0
        if sd in ['y', 'Y']:
            if sds in ['b', 'B']:
                hal = x1val*hel1+x2val*hel2 ##h or β for alloy
            else:
                pass
        else:
            pass
        while Tval <= Tmax:
            l0val = a[0] + b[0]*Tval
            l1val = a[1] + b[1]*Tval
            l2val = a[2] + b[2]*Tval
            l3val = a[3] + b[3]*Tval
            l4val = a[4] + b[4]*Tval
            T.append(Tval)
            pge1.append(PMEGE1(x2val, Tval, l0val, l1val, l2val, l3val, l4val)) #partial excess Gibbs energy
            pge2.append(PMEGE2(x2val, Tval, l0val, l1val, l2val, l3val, l4val))
            Tval = round(Tval + Tsp, 6)
            i = i + 1
        spge1 = [] #surface partial excess Gibbs energy for element 1
        spge2 = [] #surface partial excess Gibbs energy for element 2
        pste1_r = [None]*i #surface tension of pure element 1
        pste1_rr = [None]*i # surface tension of pure element 1 including the size effect
        pste2_r = [None]*i #surface tension of pure element 2
        pste2_rr = [None]*i # surface tension of pure element 2 including the size effect
        pste_r = [None]*i #the sencond term in Y
        mve1 = [None]*i #molar volume of pure element 1
        mve2 = [None]*i #molar volume of pure element 2
        pmve1 = [None]*i #partial molar volume of pure element 1
        pmve2 = [None]*i #partial molar volume of pure element 2
        mvph = [None]*i #Molar volume for phase
        stph = [None]*i #surface energy contribution of alloy, a part of the first term in Y
        YB = [None]*i #Y for multiple linear regression
        SGB = [None]*i #the surface Gibbs energy
        #MLR = [None]*i #for multiple linear regression, not used        
        miae1 = [None]*i #molar surface area for element 1
        miae2 = [None]*i #molar surface area for element 2
        sx1 = [None]*i #surface mole fraction of element 1
        sx2 = [None]*i #surface mole fraction of element 2
        st = [None]*i #surface tension of alloy or size dependent surface tension of binary alloy
        for t in range(0,i,1):
            spge1val = ra * pge1[t]
            spge1.append(spge1val)
            spge2val = ra * pge2[t]
            spge2.append(spge2val)
        for j in range(0,i,1):
            if sd in ['y', 'Y']:
                if sds in ['a', 'A']:
                    pste1_r[j] = pste1a_r + pste1t1_r*T[j]
                    pste2_r[j] = pste2a_r + pste2t1_r*T[j]
                else:
                    pste1_r[j] = pste1a + pste1t1*T[j]
                    pste2_r[j] = pste2a + pste2t1*T[j]
            else:
                pste1_r[j] = pste1a + pste1t1*T[j]
                pste2_r[j] = pste2a + pste2t1*T[j]
            if em in ['n', 'N']:
                mve1[j] = mve1a + mve1t1*T[j] + mve1t2*T[j]**2
                mve2[j] = mve2a + mve2t1*T[j] + mve2t2*T[j]**2
                pmve1[j] = mve1[j] 
                pmve2[j] = mve2[j]
                mvph[j] = mve1[j]*x1val+mve2[j]*x2val #Molar volume for phase=X1V1+X2V2
            else:
                mve1[j] = round(VI(T[j], mve1V0, mve1VAa, mve1VAb, mve1VAc, mve1VAd, mve1VAe, mve1VAf), 12)
                mve2[j] = round(VI(T[j], mve2V0, mve2VAa, mve2VAb, mve2VAc, mve2VAd, mve2VAe, mve2VAf), 12)
                vij = [] #Interaction parameters for excess molar volume
                for k in range(0, 5, 1):
                    vij.append(round(VIJ(T[j], vijV0[k], vijVAa[k], vijVAb[k], vijVAc[k], vijVAd[k], vijVAe[k], vijVAf[k]), 12))
                pmve1[j] = round(PMVE1(x2val, T[j], mve1[j], vij[0], vij[1], vij[2], vij[3], vij[4]), 12)
                pmve2[j] = round(PMVE2(x2val, T[j], mve2[j], vij[0], vij[1], vij[2], vij[3], vij[4]), 12)
                mvph[j] = round(x1val*mve1[j]+x2val*mve2[j]+x1val*x2val*(ds0*vij[0]+ds1*vij[1]+ds2*vij[2]+ds3*vij[3]+ds4*vij[4]), 12)
            stph[j] = mvph[j]*2*c/(x1val*x2val) #surface energy contribution for all pure components
            if sd in ['y', 'Y']:
                if sds in ['b', 'B']:
                    if sdm == 'Jiang':
                        pste1_rr[j] = SIGMA_JIANG(r, hel1, pste1_r[j])
                        pste2_rr[j] = SIGMA_JIANG(r, hel2, pste2_r[j])
                    elif sdm == 'Xiong':
                        pste1_rr[j] = SIGMA_XIONG(r, hel1, pste1_r[j])
                        pste2_rr[j] = SIGMA_XIONG(r, hel2, pste2_r[j])
                    else:
                        pste1_rr[j] = SIGMA_NANDA(r, sf, hel1, pste1_r[j])
                        pste2_rr[j] = SIGMA_NANDA(r, sf, hel2, pste2_r[j])
                else:
                    pste1_rr[j] = pste1_r[j]
                    pste2_rr[j] = pste2_r[j]
            else:
                pste1_rr[j] = pste1_r[j]
                pste2_rr[j] = pste2_r[j]
            pste_r[j] = 2*c*(pste1_rr[j]*x1val*mve1[j]+pste2_rr[j]*x2val*mve2[j])/(x1val*x2val) #the sencond term in Y including the size effect
            miae1[j] = PMIA(pf, pmve1[j])
            miae2[j] = PMIA(pf, pmve2[j])    
            if x2val < 0.000001:
                sx2[j] = 0
                sx1[j] = 1
                st[j] = pste1_r[j]
                if sd in ['y', 'Y']:
                    if sds in ['b', 'B']:
                        if sdm == 'Jiang':
                            st[j] = SIGMA_JIANG(r, hal, st[j])
                        elif sdm == 'Xiong':
                            st[j] = SIGMA_XIONG(r, hal, st[j])
                        else:
                            st[j] = SIGMA_NANDA(r, sf, hal, st[j])
                    else:
                        pass
                else:
                    pass
            elif x2val > 0.999999:
                sx1[j] = 0
                sx2[j] = 1
                st[j] = pste2_r[j]
                if sd in ['y', 'Y']:
                    if sds in ['b', 'B']:
                        if sdm == 'Jiang':
                            st[j] = SIGMA_JIANG(r, hal, st[j])
                        elif sdm == 'Xiong':
                            st[j] = SIGMA_XIONG(r, hal, st[j])
                        else:
                            st[j] = SIGMA_NANDA(r, sf, hal, st[j])
                    else:
                        pass
                else:
                    pass
            else:
                st1val0 = []
                st2val0 = []
                sx1val = 10**(-4)
                sx1val0 = []
                dstval0 = []
                while sx1val <= 1:
                    sx2val = 1 - sx1val
                    sx1val0.append(sx1val)
                    st1val = DABA(pste1_r[j], miae1[j], T[j], sx1val, x1val, pge1[j], spge1[j])
                    st1val0.append(st1val)
                    st2val = DABB(pste2_r[j], miae2[j], T[j], sx2val, x2val, pge2[j], spge2[j])
                    st2val0.append(st2val)
                    dstval0.append(st1val - st2val)
                    sx1val = sx1val + 10**(-4)
                M = len(sx1val0)
                boolflag = False
                for k in range(1, M, 1):
                    if (dstval0[0]*dstval0[k])<0:
                        sx10 = sx1val0[0]
                        dst0 = dstval0[0]
                        sx11 = sx1val0[k]
                        dst1 = dstval0[k]                            
                        for t in range(0, int(n), 1):
                            sx1h = (sx10+sx11)/2
                            sx2h = 1 - sx1h
                            dsth = DABA(pste1_r[j], miae1[j], T[j], sx1h, x1val, pge1[j], spge1[j])-DABB(pste2_r[j], miae2[j], T[j], sx2h, x2val, pge2[j], spge2[j])
                            if dsth*dst0<0:
                                sx11 = sx1h
                                dst1 = dsth
                                if abs(sx10-sx11)<10**(-6):
                                    boolflag = True
                                    break
                            else:
                                sx10 = sx1h
                                dst0 = dsth
                                if abs(sx10-sx11)<10**(-6):
                                    boolflag = True
                                    break
                        sx1[j] = (sx10+sx11)/2
                        sx2[j] = 1 - sx1[j]
                        st[j] = DABA(pste1_r[j], miae1[j], T[j], sx1[j], x1val, pge1[j], spge1[j])                        
                        if sd in ['y', 'Y']:
                            if sds in ['b', 'B']:
                                if sdm == 'Jiang':
                                    st[j] = SIGMA_JIANG(r, hal, st[j])
                                elif sdm == 'Xiong':
                                    st[j] = SIGMA_XIONG(r, hal, st[j])
                                else:
                                    st[j] = SIGMA_NANDA(r, sf, hal, st[j])
                            else:
                                pass
                        else:
                            pass
                    if boolflag == True:
                        break
            if st[j] is not None:
                if st[j] < 0:
                    st[j] = None
                else:
                    YB[j] = stph[j]*st[j]-pste_r[j]
                    SGB[j] = 2*c*st[j]*mvph[j]/r*10**9
            else:
                pass
        csv_writer.writerow(['Nanoparticle radius (nm)', r])
        csv_writer.writerow(['x('+ el2 + ')', x2val])
        csv_writer.writerow(['T (K)', 'G(' + el1 + ', E, bulk) (J/mol)', 'G(' + el2 + ', E, bulk) (J/mol)', 'G(' + el1 + ', E, surface) (J/mol)', 'G(' + el2 + ', E, surface) (J/mol)', 'x(' + el1 + ', surface)', 'x(' + el2 + ', surface)', 'Partial molar volume of ' + el1 + ' in ' + ph + ' (m3/mol)', 'Partial molar volume of ' + el2 + ' in ' + ph + ' (m3/mol)', 'Molar volume for ' + ph + ' (m3/mol)', 'Surface tension of ' + el1 + ' (N/m)', 'The size effect Surface tension of ' + el1 + ' (N/m)', 'Surface tension of ' + el2 + ' (N/m)', 'The size effect Surface tension of ' + el2 + ' (N/m)', 'Surface tension (N/m)', 'The surface Gibbs energy (J/mol)', 'Y (for multiple linear regression)', '(x1-x2)**0', '(x1-x2)**1', '(x1-x2)**2', '(x1-x2)**3', '(x1-x2)**4'])
        for j in range(0,i,1):
            csv_writer.writerow([T[j], pge1[j], pge2[j], spge1[j], spge2[j], sx1[j], sx2[j], pmve1[j], pmve2[j], mvph[j], pste1_r[j], pste1_rr[j], pste2_r[j], pste2_rr[j], st[j], SGB[j], YB[j], ds0, ds1, ds2, ds3, ds4])
        #f.close()
        x2val = round(x2val+Xsp, 6)
    r = round(r+Rsp, 6)
f.close()
print("**********The calculation is done!**********")
input("Press <enter> to close program")
exit()
