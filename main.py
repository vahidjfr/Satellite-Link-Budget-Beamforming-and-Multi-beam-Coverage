import numpy as np
import math
from scipy.constants import c
from scipy.constants import pi
from scipy import special
from scipy.integrate import quad


## Distribution of cloud attenuation along slant path based on global digital maps
def cloud_loss(f,e,typ):

    k=0.0819*f
    if typ=='h':
        m=0.5
    elif typ=='m':
        m=0.05
    elif typ=='l':
        m=0.001
    l=(k*m* 0.45 )/math.sin(math.radians(e))
    l= math.exp(l/10)
    return l

## main code structure
## ask about parameters
h_sat = int(input("Determine altitude of SAT: in KM"))
elevation = int(input("Determine elevation angle: in degree"))

## confirm input parameters

while True:
    answ = input("Do you want to use beamforming [y/n]?")
    if answ == "n":
        f = 1
        break
    elif answ == "y":
        tet0 = float(input("tetasteering lobe of antenna:"))
        phi0 = float(input("Phi of main steering lobe of antenna:"))
        phi = float(input("Determine azimuth of antenna:"))
        answ2 = input("what type of beamforming do you use? [Linear=li, planar=pl, circular=ci]")

    ## info about frequency and antenna
    while True:
        print("General parameters about link")
        fc = input("Carrier frequency:")
        fc = int(fc)
        #det_band(fc)
        landa = c/fc
        print("Lambda is:", landa)
        answ = input ("are you sure about general parameters []y/n ")
        if answ == "y":
            del answ
            break
        elif answ == "n":
            del answ
            continue
        else:
            del answ
            break

## info about position of satellite
    while True:
        h_sat = int(input("Determine altitude of Sat: in Km"))
        elevation = int(input("Determine elevation angle: [in degree] "))
        #di = get_slant_range(elevation, h_sat)
        answ = input(" Are sure about Satellite position parameters:[y/n")
        if answ == "y":
            del answ
            break
        elif answ == "n":
            del answ
            continue
        else:
            del answ
            break


## General information of link
    while True:

        print("General parameters about Link")
        fc ,fcu = input("Carrier frequency with its prefix:['g' = giga,'m'= meg a,'k' =kilo](ex: 40 g) ").split()
        fc= int(fc)
        if fcu =='g' or fcu == 'G':
            fc = fc * 1e9
            print(" Your frequency is:", fc)
        elif fcu == 'm':
            fc = fc * 1e6
            print (" Your frequency is:", fc)
        elif fcu == 'k':
            fc = fc * 1e3
            print(" Your frequency is:", fc)
        else:
            print(" Your prefix is not defined, or You don't put any prefix ", det_band(fc))
        landa=c/fc
        print ("Lambda is: ",landa)
        answ=input ("Are you sure about General parameters[y/n] ")
        if answ == 'y':
            del answ
            break
        elif answ == 'n':
            del answ
            continue
        else:
            del answ
            break



### Information about the position of the satellite
while True:

    h_sat = float(input("Determine altitude of SAT: in Km"))
    elevation = float(input("Determine me elevation angle:[in Degree] "))
    di = get_slant_range(elevation,h_sat)[0]
    print("Slant Range is: ", di )
    print("earth centeral angle is: ",get_slant_range(elevation,h_sat)[1])
    print(" Nadir is: ",get_slant_range(elevation,h_sat)[2])
    IAA=2.55604187 * (1-math.acos(math.radians(get_slant_range(elevation,h_s at)[2])))
    print("Instaneuse Acess Area in KM :",IAA)

    answ=input("Are you sure about Satelite position parameters:[y/n] ")

    if answ == 'y':
        del answ
        break
    elif answ == 'n':
        del answ
        continue
    else:
        del answ
    break



print("***Operation is started***")
answ = input("Downlink or Uplink:[dl/ul]")

if answ == "ul":
    ### Uplink parameters

while True:

    print("Enter UPLINK Parameters:")
    print("________")
    print("***Satellite Parameters***")
    print("________")
    print("**** Satellite Antenna Parameters*****")

    teta_tx = float(input(" Satellite antenna direction angle: in degree"))
    da = float(input("Determine reflector diametere in meter:"))
    r = da / 2
    tet_3db = 70 * (landa / da)
    eff = float(input("Input efficiency 0f antena [0.01-1]:"))
    gt_max = eff * (70 * pi / tet_3db)
    ltet = w2db(12 * ((teta_tx / tet_3db) ** 2))
    ltet = w2db(ltet)
    print("Maximum gain is:", gt_max)
    gtx = ntn_patern *(teta_tx, landa, r) * gt_max

    print("Tx Antenna gain in dB is:", gtx)
    answ = input("Are you sure about Tx antenna parameters:[y/n]")
    if answ == 'y':
        del answ
        break
    elif answ == 'n':
        del answ
        continue
    else:
        del answ
        break



def ntn_patern(teta,landa,r):
    k=2*pi*landa

    if teta==0 :
        patt=1
        print("pattern is:",patt)

    else:
        patt=4*pow(abs((special.jv(1, k*r*math.sin(math.radians(teta))))/k*r*math.sin(math.radians(teta))),2)
        print("pattern is:", patt)



##Depointing margin
ltet=w2db(12*((teta_tx/tet_3db)**2))
ltet=w2db(ltet)



### page 19

while True:

        answ = input("Do you want to use beamforming[y/n]?")
        if answ == "n":
            f = 1
            break

        elif answ == "y":
            tet0 = float(input("tetasteering lobe of antenna"))
            phi0 = float(input(" Phi of main steering lobe of antenna:"))
            phi = float(input(" Determine azimuth of antenna:"))
            answ2 = str(input("what type of beamforming use?[Linear=li,Plenary=pl,Circular = ci]"))

        if answ2 == "li":
            del answ2
            k = int(input("Number of elements"))
            da = k * pi * math.sin(teta_tx)
            psi = 2 * pi * da * ((math.sin(teta_tx) - math.sin(tet0)) / landa)
            f = abs((1 / k) * (math.sin(k * psi / 2) / math.sin(psi / 2)))
            print("Linear beam forming array factor is: ", f)

        elif answ2 == "pl":
            del answ2
            k, l = input("how many elements you have in length and width? x[5 6]").split()
            k = int(k)
            l = int(l)
            du = k * pi * math.sin(teta_tx)
            dv = l * pi * math.sin(teta_tx)
            psiu = (2 * pi * du * ((math.sin(teta_tx) * math.cos(phi)) - (math.sin(tet0) * math.cos(phi0))) / landa)
            psiv = (2 * pi * du * ((math.sin(teta_tx) * math.sin(phi)) - (math.sin(tet0) * math.sin(phi0))) / landa)
            f = abs((1 / (k * l)) * (math.sin(k * psiu / 2) * math.sin(l * psiv / 2)) / (math.sin(psiu / 2) * math.sin(psiv / 2)))
            print("Planer beamforming array factor is", f)

        elif answ2 == "ci":
            del answ2
            f = 0
            phik = 0
            R = float(input("Determine radius in m:"))

            for i in range(1, k + 1):
                phik = 360 * (i / k)
                fk = math.exp(complex(0, (2 * pi / landa) * R * ((math.sin(teta_tx) * math.cos(phi - phik)) - (math.sin(tet0) * math.cos(phi0 - phik)))))
                f = fk + f
            print(f)

        else:
            print("if yes put 'y' and press enter if your answer is no put 'n'and press enter")

        answ = input("Are you sure about Beamforming antenna parameters: [y / n]")

        if answ == 'y':
            del answ
            break

        elif answ == 'n':
            del answ
            continue
        else:
            del answ
            break
    gtx = gtx *f
    gtx=w2db(gtx)
    print("TX antenna gain in dB is: " , gtx)




## power parameters of transmitter

while True:
    pt , ptu = input("Earth station Transmitter output power with its unit:[dB='db',watt='w'] ex: 9 w ").split()
    pt = float(pt)
    if ptu == 'w':
        pt = w2db(pt)
    bl = float(input(" Earth station back off loss in dB: "))
    fel = float(input("Earth station branching and feeder loss in dB: "))
    answ = input("Are you sure about power of antenna parameters:[y/n]")

    if answ == 'y':
        del answ
        break
    elif answ == 'n':
        del answ
        continue
    else:
        del answ
        break




## atmospheric loss (gas, cloud, rain)
while True:
        answ=input("Do you want to consider Atomospheric loss?[y/n]")
        if answ == 'n':
            del answ
            latm = 0
            break
        elif answ=='y':
            del answ

            while True:
                answ = input("Do you consider cloud attenuation?[y/n]")
                if answ == 'n':
                    del answ
                    lcloud = 0
                elif answ == 'y':
                    del answ
                    typ= input("Determine density of fog:[low='l',medium = 'm',high='h']")

                    lcloud = cloud_loss(fc,elevation,typ)
                answ = input("Do you consider rain attenuation?[y/n]")
                if answ == 'n':
                    del answ
                    lrain = 0
                    break
                elif answ == 'y':
                    lrain=rain_loss(fc,elevation)
                    del answ
                    break

        latm= 0.3 + w2db(math.sqrt(pow(lcloud,2)+pow(lrain,2)))
        print("Atomospheric loss is dB: ",latm)
        answ=input("are yoy sure about Atomospheric loss paramerer?[y/n] ")
        if answ == 'n':
            continue
            del answ
        elif answ == 'y':
            del answ
            break
        else:
            break




### Receiver parameters calculation loop

    while True:
        print(" Antenna Parameters for earth estation")
        teta_rx = float(input(" Earth station transmit antenna direction angle in degree:"))
        darx = float(input("Determine reflector diametere in meter:"))
        rrx=da/2
        tet_3db_rx=70*(landa/darx)
        eff_rx = float(input("Input efficiency 0f antena [0.01-1]:"))
        gt_max_rx= eff_rx*(70*pi/tet_3db_rx)
        print("Maximum gain in RX is:",gt_max_rx)
        grx=ntn_patern(teta_rx,landa,rrx)*gt_max

        while True:
            answ=input("Do you want to use beamforming[y/n]?")
            if answ=="n":
                f=1
                break
            elif answ=="y":
                tet0=float(input("tetafteering lobe of antenna"))
                phi0=float(input(" Phi of main steering lobe of antenna:"))
                phi=float(input(" Determine azimuth of antenna:"))
                answ2=input("what type of beamforming use?[Linear=li,Planery=pl,Circular=ci]")

                if answ2=="li":
                    del answ2
                    k=int(input("Number of elemnets"))
                    da=k*pi*math.sin(teta_tx)
                    psi=2*pi*da*((math.sin(teta_tx)- math.sin(tet0))/landa)
                    f=abs((1/k)*(math.sin(k*psi/2)/math.sin(psi/2)))
                    print("Linear beam forming array factor is: ",f)

                elif answ2=="pl":
                    del answ2
                    k,l =input("how many elements you have in length and width?ex[5 6]").split()
                    k=int(k)
                    l=int(l)
                    du=k*pi*math.sin(teta_tx)
                    dv=l*pi*math.sin(teta_tx)
                    psiu=(2*pi*du*((math.sin(teta_tx)*math.cos(phi) - math.sin(tet0)*math.cos(phi0)))/landa)
                    psiv=(2*pi*du*((math.sin(teta_tx)*math.sin(phi))-(math.sin(tet0)*math.sin(phi0)))/landa)
                    f=abs((1/(k*l))*(math.sin(k*psiu/2) * math.sin(l*psiv/2))/(math.sin(psiu/2) * math.sin(psiv/2)))

                    print("Planer beamforming array factor is",f)

                elif answ2=="ci":
                    del answ2
                    f=0
                    phik = 0
                    R=float(input("Determine radius in m:"))
                    for i in range(1,k+1):
                    phik = 360 * (i / k)
                    fk = math.exp(complex(0, (2 * pi / landa) * R((math.sin(teta_tx) * math.cos(phi - phik)) - (math.sin(tet0) * math.cos(phi0 - phik)))))
                    f = fk + f
                print("Circular beamforming array factor is", f)

            grx = grx * f
            grx = w2db(grx)
            print("TX antenna gain in dB is: ", grx)

        else:
        print("if yes put 'y' and press enter if your answer is no put'n' and pressenter")

        answ = input("Are you sure about Beamforming antenna parameters:[y/n]")

        if answ == 'y':
            del answ
            break

        elif answ == 'n':
            del answ
            continue
        else:
            del answ
            break

    print("RX gain is : ", grx)
    print("****Earth Station calculation*****")

        sfl = float(input("Earth station branching and feeder loss:"))
        answ = input("are yoy sure about RX antenna parameter?[y/n] ")
        if answ == 'n':
            del answ
            continue
        elif answ == 'y':
            del answ
            break
        else:
            break




#### noise parameters

while True:
    print("Noise Parameters for earth estation")
    ta , tau = input("Determine antenna tempretur with its unit[kelvin='k',Celsius='c'],(ex; 270 k): ")
    ta = float(ta)
    if tau == 'c' or tau == 'C':
        ta=c2k(ta)

    tf , tfu = input("Determine antenna tempretur with its unit[kelvin='k',Celsius='c'],(ex; 270 k): ")
    tf = float(tf)
    if tfu == 'c' or tfu == 'C':
        tf=c2k(tf)

    te , teu = input("Determine Equivalent with its unit[kelvin='k',Celsius='c'],(ex; 270 k): ")
    te = float(te)
    if teu == 'c' or teu == 'C':
        te=c2k(te)
    te = int(input("Determine feeder temprature:"))

    n = (ta/sfl)+tf*(1-1/sfl)+te
    n=w2db(n)
    print("Noise power is in dB: ",n)

    answ=input("are yoy sure about Noise antenna paramerer?[y/n] ")
    if answ == 'n':
        continue
        del answ
    elif answ == 'y':
        del answ
        break
    else:
        break



#### Multibeam (compare performance of beam with and without multibeam)

while True:
    answ=input(" Do you want use multi-beam analyse?[y/n] ")
    if answ == 'n' :
        del answ
        r_edge=math.sin(math.radians(teta_tx))
        print("Beam radius is: ",r_edge)

        ABS=math.sqrt(3) * math.radians(teta_tx)
        print("Adjacent Beam Spacing is:",ABS)
        print("Gain at SAT in dB is without any change: ",gtx)
        print("EIRP without beamforming is:",eirp)
        break

    elif answ=='y' :
        nb=int(input("Determine number of beam in multi-beam system"))
        tet_edge = teta_tx/nb
        r_edge=math.sin(math.radians(tet_edge))
        print("Beam radius is: ",r_edge)

        ABS = math.sqrt(3) * math.radians(teta_tx)
        print("Adjacent Beam Spacing is: ",ABS)
        print("Gain at SAT in dB is withou beam forming is: ",gtx)

        gtx_edge=ntn_patern(tet_edge,landa,r)*gt_max
        print("Gain of satelite after beamforming: ", gtx_edge)
        print("EIRP without beamforming is: ",eirp)

        eirp_edge = pt + gtx_edge - l
        print("EIRP with beamforming is: ",eirp_edge)
        print("*****Link budget Performances*****")

        prxn0_edge = eirp-l+grx-n+228.2
        print("Total power over noise ratio is with Multibeam in dB: ",prxn0_edge)
        break

    answ = input("Are you sure about parameters related to beamforming?[y/n]")

    if answ == 'y':
        del answ
        break
    elif answ == 'n':
        del answ
        continue
    else:
        break



#### Coverage (calculate parameters about max and min coverage probability)

while True:
        answ = input("Do you live in Urban or Rural ? [u/r][urban='u',rural='r'")
        if answ == 'r':
            alpha = 2.5
            del answ

         elif answ =='u':
            alpha = 3.25
            del answ

        Pl= pt - l +alpha*10*math.log10(r_edge)
        print("Path loss probability is:",Pl)

        if det_band(fc) == 's':
            sigma_max = 15.6
            sigma_min = 1.2

        else:
            sigma_max = 17.1
            sigma_min = 2.4

        st_max= sigma_max + 10*math.log10(r_edge)
        st_min= sigma_min + 10*math.log10(r_edge)
        integrand=0.5* special.erfc(-st_max)
        integrand1=0.5* special.erfc(-st_min)
        Pcov1= quad(integrand, 0, r_edge)
        print("Max coverage probability:" , Pcov1)
        Pcov2= quad(integrand1, 0, r_edge)
        print("Max coverage probability:" , Pcov2)

            answ = input("Are you sure about parameters related to Coverage?[y/n]")
            if answ == 'y':
                del answ
                break
            elif answ == 'n':
                del answ
                continue
            else:
                break















