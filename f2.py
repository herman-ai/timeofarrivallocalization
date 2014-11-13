"""
TODO: 1. Remove ToA observations from other sensor nodes that are too far away.
2. Work on an iterative solution
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

for epsilon_index in range(1):
#for epsilon_index in range(len(EPSILON_S_REAL_ARRAY)):
    print "********************************************"
    print "EPSILON_INDEX = ", epsilon_index
    print "********************************************"
    EPSILON_S_REAL = EPSILON_SOIL[epsilon_index]["real"]
    EPSILON_S_IMG = EPSILON_SOIL[epsilon_index]["img"]


    ALPHA_SOIL = OMEGA * sqrt( (MU_S * EPSILON_S_REAL/2) * (sqrt(1 + (EPSILON_S_IMG/EPSILON_S_REAL) ** 2) - 1))
    RHO = ((sqrt(MU_A / EPSILON_AIR) - sqrt(MU_S / EPSILON_S_REAL)) / \
        (sqrt(MU_A / EPSILON_AIR) + sqrt(MU_S / EPSILON_S_REAL))) ** 2

    TAU = 1 - RHO

    xyzAllActual = ()
    for i in range(NUM_SENSORS):
        xyzAllActual = xyzAllActual + ((random() * F, random() * F, random() * (-H/40)), )


    AnchorsXYZ = {"above": [(0., 0., H), (F, 0., H), (0., F, H), (F, F, H)], \
        "below": [(0., 0., -H/4), (F, 0., -H/4), (0., F, -H/4), (F, F, -H/4)]}
    #AnchorsXYZ = {"above": [(0., 0., H), (F, 0., H), (0., F, H), (F, F, H)]}

    obsTAnchors = {"above": (), "below": ()}

    for xyzOneSensorActual in xyzAllActual:
        observed = ()
        for anchor in AnchorsXYZ["below"]:
            distance = sqrt((anchor[0] - xyzOneSensorActual[0]) ** 2 +\
                    (anchor[1] - xyzOneSensorActual[1]) ** 2 +\
                    (anchor[2] - xyzOneSensorActual[2]) ** 2)
            observed = observed + (distance / speed + random()/40, )
        obsTAnchors["below"] = obsTAnchors["below"] + (observed,)

        observed = ()
        for anchor in AnchorsXYZ["above"]:
                distanceAir = sqrt((anchor[0] - xyzOneSensorActual[0]) ** 2 + (anchor[1] - xyzOneSensorActual[1]) ** 2 + (anchor[2]) ** 2)
                distanceSoil = abs(xyzOneSensorActual[2])
                mean = (distanceAir + distanceSoil) / speed
                sigma = sigma_air2Soil(distanceAir, distanceSoil, ALPHA_SOIL, TAU)
                #print "mean = ", (distanceAir + distanceSoil) / speed, ", sigma = ", sigma
                observed = observed + (numpy.random.normal(loc = mean, scale = sigma, size = 1)[0], )
                #observed = observed + ((distanceAir + distanceSoil) / speed + random()/400, )
        obsTAnchors["above"] = obsTAnchors["above"] + (observed,)

    observedTime3DSS = {}

    for i in range(len(xyzAllActual)):
        for j in range(i+1, len(xyzAllActual)):
            xyzi = xyzAllActual[i]
            xyzj = xyzAllActual[j]
            dist = sqrt((xyzi[0] - xyzj[0]) ** 2 + (xyzi[1] - xyzj[1]) ** 2 +(xyzi[2] - xyzj[2]) ** 2)
            # TODO: continue the loop if the two nodes are too far apart
            signalPowerdB = p_average_soil2soil(dist, ALPHA_SOIL)
            print "signal Power = ", signalPowerdB

            observed = dist / speed + random()/400
            mean = dist / speed
            sigma = sigma_soil2Soil(dist, ALPHA_SOIL)
            print "mean = ", mean, ", sigma=", sigma
            observedTime3DSS[str(i) + '-' + str(j)] = observed
    xyzFirstEst = ()
    for i in range(NUM_SENSORS):
        xyzInitialEstimate = (random() * F, random() * F, random() * (-H/40))
        res = minimize(timeOfArrivalMatcherAnchorToOneSensor, xyzInitialEstimate, args=(i, obsTAnchors, AnchorsXYZ, observedTime3DSS), method='Nelder-Mead', \
        options = {"maxiter":1e6})

    	e = res.x
    	if (e[2] > 0 ):
    	    e[2] = -e[2]
        print res.success
        xyzFirstEst = xyzFirstEst  + ((e[0], e[1], e[2]),)

        #print xyzAllActual[i]
        #print e
        print "*****"
    print xyzAllActual
    print xyzFirstEst
    res = minimize(timeOfArrialMatcher3DX, xyzFirstEst, args=(obsTAnchors, AnchorsXYZ, observedTime3DSS), method='Nelder-Mead', options = {"maxiter":1e6})

    e = [(res.x[i * 3: (i+1) * 3]) for i in range(NUM_SENSORS)]

    for i in range(len(xyzAllActual)):
        print "||||||"
        print xyzAllActual[i]
        print e[i]

    params_text = r"\noindent$\epsilon_s'$ = " + '{0:.3g}'.format(EPSILON_S_REAL/EPSILON_0) + r"\\\\" \
    + r"$\epsilon_s''$ = " + '{0:.3g}'.format(EPSILON_S_IMG/EPSILON_0) + r"\\\\" \
    + r"$\alpha^{(s)}$ = " + '{0:.3g}'.format(ALPHA_SOIL) + r" N/m\\\\"

    plot([xyz[0] for xyz in xyzAllActual], [xyz[1] for xyz in xyzAllActual], \
        [xyz[0] for xyz in e], [xyz[1] for xyz in e], EPSILON_S_REAL, EPSILON_S_IMG, params_text, "X", "Y", [-1, F + 1 , -1, F + 1])
    plot([xyz[0] for xyz in xyzAllActual], [xyz[2] for xyz in xyzAllActual], \
        [xyz[0] for xyz in e], [xyz[2] for xyz in e], EPSILON_S_REAL, EPSILON_S_IMG, params_text, "X", "Z", [-1, F + 1, -H/40, 0])
