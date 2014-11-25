"""
Next step: Compute standard deviation of the estimates computed.
"""

from scipy.optimize import minimize
from math import sqrt
from random import random
from random import seed
import numpy
import utils
reload(utils)
from utils import *

seed()

xyzAllActual = ()
for i in range(NUM_SENSORS):
    xyzAllActual = xyzAllActual + ((random() * F, random() * F, random() * (-H/40)), )

epsilon_index = 1

xEstSamples = [[] for _ in range(NUM_SENSORS)]
yEstSamples = [[] for _ in range(NUM_SENSORS)]
zEstSamples = [[] for _ in range(NUM_SENSORS)]
errors = []

for ctr in range(1000):
    #print "********************************************"
    #print "EPSILON_INDEX = ", epsilon_index
    #print "********************************************"
    EPSILON_S_REAL = EPSILON_SOIL[epsilon_index]["real"]
    EPSILON_S_IMG = EPSILON_SOIL[epsilon_index]["img"]


    ALPHA_SOIL = OMEGA * sqrt( (MU_S * EPSILON_S_REAL/2) * (sqrt(1 + (EPSILON_S_IMG/EPSILON_S_REAL) ** 2) - 1))
    RHO = ((sqrt(MU_A / EPSILON_AIR) - sqrt(MU_S / EPSILON_S_REAL)) / \
        (sqrt(MU_A / EPSILON_AIR) + sqrt(MU_S / EPSILON_S_REAL))) ** 2

    TAU = 1 - RHO

    speed_soil = ( \
        sqrt( \
        (MU_S * EPSILON_S_REAL / 2 ) * (sqrt(1 + (EPSILON_S_IMG/EPSILON_S_REAL) ** 2) + 1)\
        ) \
    ) \
        ** (-1) * (10 ** (-7))


    AnchorsXYZ = {"above": [(0., 0., H), (F, 0., H), (0., F, H), (F, F, H)], \
        "below": [(0., 0., -H/4), (F, 0., -H/4), (0., F, -H/4), (F, F, -H/4)]}
    #AnchorsXYZ = {"above": [(0., 0., H), (F, 0., H), (0., F, H), (F, F, H)]}

    obsTAnchors = {"above": (), "below": ()}

    #Anchor to sensor observations
    for xyzOneSensorActual in xyzAllActual:
        observed = {}
        i = 0
        for anchor in AnchorsXYZ["below"]:
            distance = sqrt((anchor[0] - xyzOneSensorActual[0]) ** 2 +\
                    (anchor[1] - xyzOneSensorActual[1]) ** 2 +\
                    (anchor[2] - xyzOneSensorActual[2]) ** 2)
            power = p_average_soil2soil(distance, ALPHA_SOIL)
            if (power > -110.0):
                mean = distance / speed_soil + random()/40
                sigma = sigma_soil2Soil(distance, ALPHA_SOIL)
                ob = numpy.random.normal(loc = mean, scale = sigma, size = 1)[0]
                observed[i] = numpy.random.uniform(ob-DELTA_T, ob+DELTA_T)
            i = i + 1
        obsTAnchors["below"] = obsTAnchors["below"] + (observed, )

        observed = {}
        i = 0
        for anchor in AnchorsXYZ["above"]:
                distanceAir = sqrt((anchor[0] - xyzOneSensorActual[0]) ** 2 +\
                    (anchor[1] - xyzOneSensorActual[1]) ** 2 +\
                    (anchor[2]) ** 2)
                distanceSoil = abs(xyzOneSensorActual[2])
                signalPowerdB = p_average_air2soil(distanceAir, distanceSoil, ALPHA_SOIL, TAU)
                if signalPowerdB > -110.0:
                    mean = distanceAir / speed + distanceSoil / speed_soil
                    sigma = sigma_air2Soil(distanceAir, distanceSoil, ALPHA_SOIL, TAU)
                    ob = numpy.random.normal(loc = mean, scale = sigma, size = 1)[0]
                    observed[i] = numpy.random.uniform(ob-DELTA_T, ob+DELTA_T)
                i = i + 1
                #observed = observed + ((distanceAir + distanceSoil) / speed + random()/400, )
        obsTAnchors["above"] = obsTAnchors["above"] + (observed, )

    observedTime3DSS = {}

    # Sensor to sensor observations
    for i in range(len(xyzAllActual)):
        for j in range(i+1, len(xyzAllActual)):
            xyzi = xyzAllActual[i]
            xyzj = xyzAllActual[j]
            dist = sqrt((xyzi[0] - xyzj[0]) ** 2 + (xyzi[1] - xyzj[1]) ** 2 +(xyzi[2] - xyzj[2]) ** 2)
            signalPowerdB = p_average_soil2soil(dist, ALPHA_SOIL)

            if signalPowerdB > -110.0:
                mean = dist / speed_soil
                sigma = sigma_soil2Soil(dist, ALPHA_SOIL)
                observed = numpy.random.normal(loc = mean, scale = sigma, size = 1)[0]
                observedTime3DSS[str(i) + '-' + str(j)] = numpy.random.uniform(observed-DELTA_T, observed+DELTA_T)
    # Estimation
    xyzFirstEst = ()
    for i in range(NUM_SENSORS):
        xyzInitialEstimate = (random() * F, random() * F, random() * (-H/40))
        res = minimize(timeOfArrivalMatcherAnchorToOneSensor, xyzInitialEstimate, args=(i, obsTAnchors, AnchorsXYZ, observedTime3DSS, speed_soil), method='Nelder-Mead', \
        options = {"maxiter":1e6})

       	e = res.x
       	if (e[2] > 0 ):
       	    e[2] = -e[2]
        #print res.success
        xyzFirstEst = xyzFirstEst  + ((e[0], e[1], e[2]),)

        #print xyzAllActual[i]
        #print e
        #print "*****"
    #print xyzAllActual
    #print "xyzFirstEst = ", xyzFirstEst
    res = minimize(timeOfArrialMatcher3DX, xyzFirstEst, args=(obsTAnchors, AnchorsXYZ, observedTime3DSS, speed_soil), method='Nelder-Mead', options = {"maxiter":1e6})

    e = [(res.x[i * 3: (i+1) * 3]) for i in range(NUM_SENSORS)]

    #for i in range(len(xyzAllActual)):
    #    print "||||||"
    #    print xyzAllActual[i]
    #    print e[i]

    params_text = r"\noindent$\epsilon_s'$ = " + '{0:.3g}'.format(EPSILON_S_REAL/EPSILON_0) + r"\\\\" \
    + r"$\epsilon_s''$ = " + '{0:.3g}'.format(EPSILON_S_IMG/EPSILON_0) + r"\\\\" \
    + r"$\alpha^{(s)}$ = " + '{0:.3g}'.format(ALPHA_SOIL) + r" N/m\\\\"

    xActual = [xyz[0] for xyz in xyzAllActual]
    yActual = [xyz[1] for xyz in xyzAllActual]
    zActual = [xyz[2] for xyz in xyzAllActual]

    xFirstEst = [xyz[0] for xyz in xyzFirstEst]
    yFirstEst = [xyz[1] for xyz in xyzFirstEst]
    zFirstEst = [xyz[2] for xyz in xyzFirstEst]

    xEst = [xyz[0] for xyz in e]
    yEst = [xyz[1] for xyz in e]
    zEst = [xyz[2] for xyz in e]

    err = 0.0
    for i in range(len(xActual)):
        d = (xEst[i] - xActual[i]) ** 2 + \
            (yEst[i] - yActual[i]) ** 2 + \
            (zEst[i] - zActual[i]) ** 2
        err = err + sqrt(d)

    #print "DELTA_T = ", DELTA_T, "error = ", err
    #print "T_S = ", T_S, "error = ", err
    print ctr, ". epsilon_index = ", epsilon_index, "error = ", err
    errors.append(err)
    for i in range(NUM_SENSORS):
        xEstSamples[i].append(e[i][0])
        yEstSamples[i].append(e[i][1])
        zEstSamples[i].append(e[i][2])

    #plot([xyz[0] for xyz in xyzAllActual], [xyz[1] for xyz in xyzAllActual], \
    #        [xyz[0] for xyz in xyzFirstEst], [xyz[1] for xyz in xyzFirstEst], \
    #    [xyz[0] for xyz in e], [xyz[1] for xyz in e], EPSILON_S_REAL, EPSILON_S_IMG, params_text, "X", "Y", [-1, F + 1 , -1, F + 1])
    #plot([xyz[0] for xyz in xyzAllActual], [xyz[2] for xyz in xyzAllActual], \
    #    [xyz[0] for xyz in xyzFirstEst], [xyz[2] for xyz in xyzFirstEst], \
    #    [xyz[0] for xyz in e], [xyz[2] for xyz in e], EPSILON_S_REAL, EPSILON_S_IMG, params_text, "X", "Z", [-1, F + 1, -H/40, 0])

print "X STD = ", numpy.mean([numpy.std(x) for x in xEstSamples])
print "Y STD = ", numpy.mean([numpy.std(x) for x in yEstSamples])
print "Z STD = ", numpy.mean([numpy.std(x) for x in zEstSamples])
print "error = ", numpy.mean(errors)
