#!/usr/bin/python

# ============ cubex-cl.py ====================
# Version of CubeX for command line operation
# Use command python ./cubex-cl.py <filename>
# <filename> should refer to an input file containing tab-delimited columns of data
# where each column is one of nine diplotypes.
# =============================================

import cgi,os,time,sys,math

datafilename = sys.argv[1]

def cubex(inputarray):
    result = {}
    n1111 = float(inputarray[0])
    n1112 = float(inputarray[1])
    n1122 = float(inputarray[2])
    n1211 = float(inputarray[3])
    n1212 = float(inputarray[4])
    n1222 = float(inputarray[5])
    n2211 = float(inputarray[6])
    n2212 = float(inputarray[7])
    n2222 = float(inputarray[8])
    n = (n1111 + n1112 + n1122 + n1211 + n1212 + n1222 + n2211 + n2212 + n2222)
    result["n"] = n
    p = ((n1111+n1112+n1122)*2.0+(n1211+n1212+n1222))/(2.0 * n)
    q = ((n1111+n1211+n2211)*2.0+(n1112+n1212+n2212))/(2.0 * n)
    n11 = (2.0*n1111 + n1112 + n1211) 
    n12 = (2.0*n1122 + n1112 + n1222)
    n21 = (2.0*n2211 + n2212 + n1211)
    n22 = (2.0*n2222 + n2212 + n1222) 
    a0 = -n11*p*q
    a1 = -n11*(1.0 - 2.0*p - 2.0*q) - n1212*(1.0 - p - q) + 2.0*n*p*q
    a2 = 2.0*n*(1.0 - 2.0*p - 2.0*q) - 2.0*n11 - n1212
    a3 = 4.0 * n
    minhap = n11 / (2.0 * float(n))
    maxhap = (n11 + n1212) / (2.0 * float(n))
    result["minhap"] = minhap
    result["maxhap"] = maxhap
    if p < 1.0 and p > 0.0:
        result["hwchisnp1"] = ((((n1111 + n1112 + n1122) - ((p ** 2)*n)))**2)/((p ** 2)*n)+\
                    ((((n1211 + n1212 + n1222) - ((2 * p * (1.0-p))*n)))**2)/((2 * p * (1.0-p))*n)+\
                    ((((n2211 + n2212 + n2222) - (((1-p) ** 2)*n)))**2)/(((1-p) ** 2)*n)
    else:
        result["hwchisnp1"] = 0.0
    if q < 1.0 and q > 0.0:
        result["hwchisnp2"] = ((((n1111 + n1211 + n2211) - ((q ** 2)*n)))**2)/((q ** 2)*n)+\
                    ((((n1112 + n1212 + n2212) - ((2 * q * (1.0-q))*n)))**2)/((2 * q * (1.0-q))*n)+\
                    ((((n1122 + n1222 + n2222) - (((1-q) ** 2)*n)))**2)/(((1-q) ** 2)*n)
    else:
        result["hwchisnp2"] = 0.0
    
    a = a3
    b = a2
    c = a1
    dee = a0
    
    xN = -b/(3.0*a)
    d2 = (math.pow(b,2)-3.0*a*c)/(9*math.pow(a,2))
    yN = a * math.pow(xN,3) + b * math.pow(xN,2) + c * xN + dee
    yN2 = math.pow(yN,2)

    h2 = 4 * math.pow(a,2) * math.pow(d2,3)
    result["realnoproblem"] = 0
    if abs(yN2-h2) <= 0.00001:
        result["realnoproblem"] = 1

    if yN2 > h2:
        # option 1
        number1 = 0.0
        number2 = 0.0
        if (1.0/(2.0*a)*(-yN + math.pow((yN2 - h2),0.5))) < 0:
            number1 = -math.pow(-(1.0/(2.0*a)*(-yN + math.pow((yN2 - h2),0.5))),1.0/3.0)
        else: number1 = math.pow((1.0/(2.0*a)*(-yN + math.pow((yN2 - h2),0.5))),1.0/3.0)
        if (1.0/(2.0*a)*(-yN - math.pow((yN2 - h2),0.5))) < 0:
            number2 = -math.pow(-(1.0/(2.0*a)*(-yN - math.pow((yN2 - h2),0.5))),1.0/3.0)
        else: number2 = math.pow((1.0/(2.0*a)*(-yN - math.pow((yN2 - h2),0.5))),1.0/3.0)
        alpha = xN + number1 + number2
        result["inputdata"] = str(inputarray)
        result["alpha"] = alpha
        result["beta"] = "Not a real root"
        result["gamma"] = "Not a real root"
        result["p"] = p
        result["q"] = q
        result["pq"] = str(p * q)
        if result["alpha"] >= minhap - 0.00001 and result["alpha"] <= maxhap + 0.00001:
            result["alphaposs"] = 1
            f11 = result["alpha"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["alphaDprime"] = round(Dprime,3)
            result["alpharsquared"] = round(rsquared,4)
            result["alphaf11"] = f11
            result["alphaf12"] = f12
            result["alphaf21"] = f21
            result["alphaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["alphaposs"] = 0
            else:
                result["alpha-e1111"] = n * f11**2
                result["alpha-e1112"] = 2 * n * f11 * f12
                result["alpha-e1122"] = n * f12**2
                result["alpha-e1211"] = 2 * n * f11 * f21
                result["alpha-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["alpha-e1222"] = 2 * n * f12 * f22
                result["alpha-e2211"] = n * f21**2
                result["alpha-e2212"] = 2 * n * f21 * f22
                result["alpha-e2222"] = n * f22**2
                if result["alpha-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["alpha-e1111"])**2/result["alpha-e1111"])
                else: chisq1111 = 0.0
                if result["alpha-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["alpha-e1112"])**2/result["alpha-e1112"])
                else: chisq1112 = 0.0
                if result["alpha-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["alpha-e1122"])**2/result["alpha-e1122"])
                else: chisq1122 = 0.0
                if result["alpha-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["alpha-e1211"])**2/result["alpha-e1211"])
                else: chisq1211 = 0.0
                if result["alpha-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["alpha-e1212"])**2/result["alpha-e1212"])
                else: chisq1212 = 0.0
                if result["alpha-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["alpha-e1222"])**2/result["alpha-e1222"])
                else: chisq1222 = 0.0
                if result["alpha-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["alpha-e2211"])**2/result["alpha-e2211"])
                else: chisq2211 = 0.0
                if result["alpha-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["alpha-e2212"])**2/result["alpha-e2212"])
                else: chisq2212 = 0.0
                if result["alpha-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["alpha-e2222"])**2/result["alpha-e2222"])
                else: chisq2222 = 0.0
                result["alpha-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else:
            result["alphaposs"] = 0
        result["betaposs"] = 0
        result["gammaposs"] = 0

    elif yN2 == h2:
        # option 2
        base = yN/2.0*a
        delta = math.copysign(math.pow(abs(base),(1.0/3.0)), base) ## changed from math.pow((yN/2.0*a),(1.0/3.0)) to avoid failing when (yN/2.0*a) is negative
        result["inputdata"] = str(inputarray)
        result["alpha"] = xN + delta
        result["beta"] = xN + delta
        result["gamma"] = xN - 2.0*delta
        result["p"] = str(p)
        result["q"] = str(q)
        result["pq"] = str(p * q)
        if result["alpha"] >= minhap - 0.00001 and result["alpha"] <= maxhap + 0.00001:
            result["alphaposs"] = 1
            f11 = result["alpha"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["alphaDprime"] = round(Dprime,3)
            result["alpharsquared"] = round(rsquared,4)
            result["alphaf11"] = f11
            result["alphaf12"] = f12
            result["alphaf21"] = f21
            result["alphaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["alphaposs"] = 0
            else:
                result["alpha-e1111"] = n * f11**2
                result["alpha-e1112"] = 2 * n * f11 * f12
                result["alpha-e1122"] = n * f12**2
                result["alpha-e1211"] = 2 * n * f11 * f21
                result["alpha-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["alpha-e1222"] = 2 * n * f12 * f22
                result["alpha-e2211"] = n * f21**2
                result["alpha-e2212"] = 2 * n * f21 * f22
                result["alpha-e2222"] = n * f22**2
                if result["alpha-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["alpha-e1111"])**2/result["alpha-e1111"])
                else: chisq1111 = 0.0
                if result["alpha-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["alpha-e1112"])**2/result["alpha-e1112"])
                else: chisq1112 = 0.0
                if result["alpha-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["alpha-e1122"])**2/result["alpha-e1122"])
                else: chisq1122 = 0.0
                if result["alpha-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["alpha-e1211"])**2/result["alpha-e1211"])
                else: chisq1211 = 0.0
                if result["alpha-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["alpha-e1212"])**2/result["alpha-e1212"])
                else: chisq1212 = 0.0
                if result["alpha-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["alpha-e1222"])**2/result["alpha-e1222"])
                else: chisq1222 = 0.0
                if result["alpha-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["alpha-e2211"])**2/result["alpha-e2211"])
                else: chisq2211 = 0.0
                if result["alpha-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["alpha-e2212"])**2/result["alpha-e2212"])
                else: chisq2212 = 0.0
                if result["alpha-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["alpha-e2222"])**2/result["alpha-e2222"])
                else: chisq2222 = 0.0
                result["alpha-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else:
            result["alphaposs"] = 0
        if result["beta"] >= minhap - 0.00001 and result["beta"] <= maxhap + 0.00001:
            result["betaposs"] = 1
            f11 = result["beta"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["betaDprime"] = round(Dprime,3)
            result["betarsquared"] = round(rsquared,4)
            result["betaf11"] = f11
            result["betaf12"] = f12
            result["betaf21"] = f21
            result["betaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["betaposs"] = 0
            else:
                result["beta-e1111"] = n * f11**2
                result["beta-e1112"] = 2 * n * f11 * f12
                result["beta-e1122"] = n * f12**2
                result["beta-e1211"] = 2 * n * f11 * f21
                result["beta-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["beta-e1222"] = 2 * n * f12 * f22
                result["beta-e2211"] = n * f21**2
                result["beta-e2212"] = 2 * n * f21 * f22
                result["beta-e2222"] = n * f22**2
                if result["beta-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["beta-e1111"])**2/result["beta-e1111"])
                else: chisq1111 = 0.0
                if result["beta-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["beta-e1112"])**2/result["beta-e1112"])
                else: chisq1112 = 0.0
                if result["beta-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["beta-e1122"])**2/result["beta-e1122"])
                else: chisq1122 = 0.0
                if result["beta-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["beta-e1211"])**2/result["beta-e1211"])
                else: chisq1211 = 0.0
                if result["beta-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["beta-e1212"])**2/result["beta-e1212"])
                else: chisq1212 = 0.0
                if result["beta-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["beta-e1222"])**2/result["beta-e1222"])
                else: chisq1222 = 0.0
                if result["beta-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["beta-e2211"])**2/result["beta-e2211"])
                else: chisq2211 = 0.0
                if result["beta-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["beta-e2212"])**2/result["beta-e2212"])
                else: chisq2212 = 0.0
                if result["beta-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["beta-e2222"])**2/result["beta-e2222"])
                else: chisq2222 = 0.0
                result["beta-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else:
            result["betaposs"] = 0
        if result["gamma"] >= minhap - 0.00001 and result["gamma"] <= maxhap + 0.00001:
            result["gammaposs"] = 1
            f11 = result["gamma"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["gammaDprime"] = round(Dprime,3)
            result["gammarsquared"] = round(rsquared,4)
            result["gammaf11"] = f11
            result["gammaf12"] = f12
            result["gammaf21"] = f21
            result["gammaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["gammaposs"] = 0
            else:
                result["gamma-e1111"] = n * f11**2
                result["gamma-e1112"] = 2 * n * f11 * f12
                result["gamma-e1122"] = n * f12**2
                result["gamma-e1211"] = 2 * n * f11 * f21
                result["gamma-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["gamma-e1222"] = 2 * n * f12 * f22
                result["gamma-e2211"] = n * f21**2
                result["gamma-e2212"] = 2 * n * f21 * f22
                result["gamma-e2222"] = n * f22**2
                if result["gamma-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["gamma-e1111"])**2/result["gamma-e1111"])
                else: chisq1111 = 0.0
                if result["gamma-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["gamma-e1112"])**2/result["gamma-e1112"])
                else: chisq1112 = 0.0
                if result["gamma-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["gamma-e1122"])**2/result["gamma-e1122"])
                else: chisq1122 = 0.0
                if result["gamma-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["gamma-e1211"])**2/result["gamma-e1211"])
                else: chisq1211 = 0.0
                if result["gamma-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["gamma-e1212"])**2/result["gamma-e1212"])
                else: chisq1212 = 0.0
                if result["gamma-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["gamma-e1222"])**2/result["gamma-e1222"])
                else: chisq1222 = 0.0
                if result["gamma-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["gamma-e2211"])**2/result["gamma-e2211"])
                else: chisq2211 = 0.0
                if result["gamma-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["gamma-e2212"])**2/result["gamma-e2212"])
                else: chisq2212 = 0.0
                if result["gamma-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["gamma-e2222"])**2/result["gamma-e2222"])
                else: chisq2222 = 0.0
                result["gamma-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else:
            result["gammaposs"] = 0
        
    elif yN2 < h2:
        #option 3
        h = math.pow(h2, 0.5)
        theta = ((math.acos(-yN/h))/3.0)
        #delta = math.pow((yN/2.0*a),(1.0/3.0)) # is it correct to reuse this?
        delta = math.pow(d2,0.5)
        result["inputdata"] = str(inputarray)
        result["alpha"] = xN + 2.0 * delta * math.cos(theta)
        result["beta"] = xN + 2.0 * delta * math.cos(2.0 * math.pi/3.0 + theta)
        result["gamma"] = xN + 2.0 * delta * math.cos(4.0 * math.pi/3.0 + theta)
        result["p"] = p
        result["q"] = q
        result["pq"] = p * q
        if result["alpha"] >= minhap - 0.00001 and result["alpha"] <= maxhap + 0.00001:
            result["alphaposs"] = 1
            f11 = result["alpha"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["alphaDprime"] = round(Dprime,3)
            result["alpharsquared"] = round(rsquared,4)
            result["alphaf11"] = f11
            result["alphaf12"] = f12
            result["alphaf21"] = f21
            result["alphaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["alphaposs"] = 0
            else:
                result["alpha-e1111"] = n * f11**2
                result["alpha-e1112"] = 2 * n * f11 * f12
                result["alpha-e1122"] = n * f12**2
                result["alpha-e1211"] = 2 * n * f11 * f21
                result["alpha-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["alpha-e1222"] = 2 * n * f12 * f22
                result["alpha-e2211"] = n * f21**2
                result["alpha-e2212"] = 2 * n * f21 * f22
                result["alpha-e2222"] = n * f22**2
                if result["alpha-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["alpha-e1111"])**2/result["alpha-e1111"])
                else: chisq1111 = 0.0
                if result["alpha-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["alpha-e1112"])**2/result["alpha-e1112"])
                else: chisq1112 = 0.0
                if result["alpha-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["alpha-e1122"])**2/result["alpha-e1122"])
                else: chisq1122 = 0.0
                if result["alpha-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["alpha-e1211"])**2/result["alpha-e1211"])
                else: chisq1211 = 0.0
                if result["alpha-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["alpha-e1212"])**2/result["alpha-e1212"])
                else: chisq1212 = 0.0
                if result["alpha-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["alpha-e1222"])**2/result["alpha-e1222"])
                else: chisq1222 = 0.0
                if result["alpha-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["alpha-e2211"])**2/result["alpha-e2211"])
                else: chisq2211 = 0.0
                if result["alpha-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["alpha-e2212"])**2/result["alpha-e2212"])
                else: chisq2212 = 0.0
                if result["alpha-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["alpha-e2222"])**2/result["alpha-e2222"])
                else: chisq2222 = 0.0
                result["alpha-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else: result["alphaposs"] = 0
        if result["beta"] >= minhap - 0.00001 and result["beta"] <= maxhap + 0.00001:
            result["betaposs"] = 1
            f11 = result["beta"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["betaDprime"] = round(Dprime,3)
            result["betarsquared"] = round(rsquared,4)
            result["betaf11"] = f11
            result["betaf12"] = f12
            result["betaf21"] = f21
            result["betaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["betaposs"] = 0
            else:
                result["beta-e1111"] = n * f11**2
                result["beta-e1112"] = 2 * n * f11 * f12
                result["beta-e1122"] = n * f12**2
                result["beta-e1211"] = 2 * n * f11 * f21
                result["beta-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22
                result["beta-e1222"] = 2 * n * f12 * f22
                result["beta-e2211"] = n * f21**2
                result["beta-e2212"] = 2 * n * f21 * f22
                result["beta-e2222"] = n * f22**2
                if result["beta-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["beta-e1111"])**2/result["beta-e1111"])
                else: chisq1111 = 0.0
                if result["beta-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["beta-e1112"])**2/result["beta-e1112"])
                else: chisq1112 = 0.0
                if result["beta-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["beta-e1122"])**2/result["beta-e1122"])
                else: chisq1122 = 0.0
                if result["beta-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["beta-e1211"])**2/result["beta-e1211"])
                else: chisq1211 = 0.0
                if result["beta-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["beta-e1212"])**2/result["beta-e1212"])
                else: chisq1212 = 0.0
                if result["beta-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["beta-e1222"])**2/result["beta-e1222"])
                else: chisq1222 = 0.0
                if result["beta-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["beta-e2211"])**2/result["beta-e2211"])
                else: chisq2211 = 0.0
                if result["beta-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["beta-e2212"])**2/result["beta-e2212"])
                else: chisq2212 = 0.0
                if result["beta-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["beta-e2222"])**2/result["beta-e2222"])
                else: chisq2222 = 0.0
                result["beta-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else: result["betaposs"] = 0
        if result["gamma"] >= minhap - 0.00001 and result["gamma"] <= maxhap + 0.00001:
            result["gammaposs"] = 1
            f11 = result["gamma"]
            f12 = p - f11
            f21 = q - f11
            f22 = 1 - (f11 + f12 + f21)
            D = (f11 * f22) - (f12 * f21)
            if D >= 0.0:
                Dmax = min(p * (1.0-q), q * (1.0-p))
            else:
                Dmax = min(p*q,(1-p)*(1-q))
            Dprime = D / Dmax
            rsquared = (D ** 2) / (p * (1-p) * q * (1-q))
            result["gammaDprime"] = round(Dprime,3)
            result["gammarsquared"] = round(rsquared,4)
            result["gammaf11"] = f11
            result["gammaf12"] = f12
            result["gammaf21"] = f21
            result["gammaf22"] = f22
            if min(f11,f12,f21,f22) < -0.0000001:
                result["gammaposs"] = 0
            else:
                result["gamma-e1111"] = n * f11**2 
                result["gamma-e1112"] = 2 * n * f11 * f12
                result["gamma-e1122"] = n * f12**2
                result["gamma-e1211"] = 2 * n * f11 * f21
                result["gamma-e1212"] = 2 * n * f12 * f21 + 2 * n * f11 * f22 
                result["gamma-e1222"] = 2 * n * f12 * f22
                result["gamma-e2211"] = n * f21**2
                result["gamma-e2212"] = 2 * n * f21 * f22
                result["gamma-e2222"] = n * f22**2
                if result["gamma-e1111"] > 0.0:
                    chisq1111 = ((n1111 - result["gamma-e1111"])**2/result["gamma-e1111"])
                else: chisq1111 = 0.0
                if result["gamma-e1112"] > 0.0:
                    chisq1112 = ((n1112 - result["gamma-e1112"])**2/result["gamma-e1112"])
                else: chisq1112 = 0.0
                if result["gamma-e1122"] > 0.0:
                    chisq1122 = ((n1122 - result["gamma-e1122"])**2/result["gamma-e1122"])
                else: chisq1122 = 0.0
                if result["gamma-e1211"] > 0.0:
                    chisq1211 = ((n1211 - result["gamma-e1211"])**2/result["gamma-e1211"])
                else: chisq1211 = 0.0
                if result["gamma-e1212"] > 0.0:
                    chisq1212 = ((n1212 - result["gamma-e1212"])**2/result["gamma-e1212"])
                else: chisq1212 = 0.0
                if result["gamma-e1222"] > 0.0:
                    chisq1222 = ((n1222 - result["gamma-e1222"])**2/result["gamma-e1222"])
                else: chisq1222 = 0.0
                if result["gamma-e2211"] > 0.0:
                    chisq2211 = ((n2211 - result["gamma-e2211"])**2/result["gamma-e2211"])
                else: chisq2211 = 0.0
                if result["gamma-e2212"] > 0.0:
                    chisq2212 = ((n2212 - result["gamma-e2212"])**2/result["gamma-e2212"])
                else: chisq2212 = 0.0
                if result["gamma-e2222"] > 0.0:
                    chisq2222 = ((n2222 - result["gamma-e2222"])**2/result["gamma-e2222"])
                else: chisq2222 = 0.0
                result["gamma-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222)
        else: result["gammaposs"] = 0
        
    else: result["Error"] = "No answer"
    # Return the result as a dictionary array
    return result

datafile = open(datafilename, 'r')
print """<?xml version="1.0" encoding="UTF-8"?>\n<dataset filename=\"""" + datafilename + """\">"""
for lines in datafile:
    lines = lines.split()
    if len(lines) == 9:
        result = cubex([lines[0],lines[1],lines[2],lines[3],lines[4],lines[5],lines[6],lines[7],lines[8]])
        resultkeys = result.keys()
        resultkeys.sort()
        print """    <result>\n        <inputdata>""" + result["inputdata"] + """</inputdata>"""
        for item in resultkeys:
            if item != "inputdata":
                print """        <""" + str(item) + """>""" + str(result[item]) + """</""" + str(item) + """>"""
        print """    </result>"""
    else:
        print """    <result>Input data in incorrect format</result>"""
        
print """</dataset>"""
datafile.close()
