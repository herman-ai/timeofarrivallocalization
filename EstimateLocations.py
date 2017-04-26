
from scipy.optimize import minimize
from math import sqrt
# import numpy as np
from random import seed
import csv
import utils
# reload(utils)
from utils import *

print("Number of sensors = {}".format(NUM_SENSORS))

actual_s_locations = np.random.random_sample(size=(NUM_SENSORS, 3)) * (F, F, -H/H)  # X, Y, Z coordinates

print("actual_s_locations.shape = {}".format(actual_s_locations.shape))

AnchorsXYZ = [(0., 0., H),
              (F, 0., H),
              (0., F, H),
              (F, F, H),
              (0., 0., -H/4),
              (F, 0., -H/4),
              (0., F, -H/4),
              (F, F, -H/4)]

with open("data/actual_sensor_locations.csv", "w") as f:
    csvwriter = csv.writer(f, delimiter="\t")
    csvwriter.writerow(["#X", "Y", "Z"])
    for xyz in actual_s_locations:
        csvwriter.writerow(xyz)

with open("data/anchor_locations.csv", "w") as f:
    csvwriter = csv.writer(f, delimiter="\t")
    csvwriter.writerow(["#X", "Y", "Z"])
    for xyz in AnchorsXYZ:
        csvwriter.writerow(xyz)



#
# f_x_std = open("x_standard_dev.csv", "w")
# f_y_std = open("y_standard_dev.csv", "w")
# f_z_std = open("z_standard_dev.csv", "w")
# f_err = open("error.csv", "w")
#
# #for epsilon_index in range(4):
epsilon_index = 0

# DELTA_T is a parameter
# DELTA_T_ARR = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12]   #2 nano seconds
#
DELTA_T = 0.02   # For Sensitivity analysis
# for DELTA_T in DELTA_T_ARR:
xEstSamples = [[] for _ in range(NUM_SENSORS)]
yEstSamples = [[] for _ in range(NUM_SENSORS)]
zEstSamples = [[] for _ in range(NUM_SENSORS)]
errors = []

EPSILON_S_REAL = np.real(EPSILON_SOIL[epsilon_index])
EPSILON_S_IMG = np.imag([epsilon_index])

ALPHA_SOIL = OMEGA * \
             sqrt((MU_S * EPSILON_S_REAL/2) *
                  (sqrt(1 + (EPSILON_S_IMG/EPSILON_S_REAL) ** 2) - 1))

RHO = ((sqrt(MU_A / EPSILON_AIR) - sqrt(MU_S / EPSILON_S_REAL)) / \
    (sqrt(MU_A / EPSILON_AIR) + sqrt(MU_S / EPSILON_S_REAL))) ** 2

TAU = 1 - RHO

speed_soil = ( \
    sqrt( \
    (MU_S * EPSILON_S_REAL / 2 ) * (sqrt(1 + (EPSILON_S_IMG/EPSILON_S_REAL) ** 2) + 1)\
    ) ) ** (-1) * (10 ** (-7))

obsTAnchors = {"above": (), "below": ()}

# Anchor to sensor observations
toa_observed_from_anchors = []
for xyzOneSensorActual in actual_s_locations:
    toa_from_anchors = []
    for anchor in AnchorsXYZ[:4]:     #For anchors above soil
        # assert anchor.shape == xyzOneSensorActual.shape
        distance = np.sqrt(np.sum((anchor-xyzOneSensorActual)**2))

        # power = p_average_soil2soil(distance, ALPHA_SOIL)   # Received signal power for soil anchors
        distanceAir = sqrt((anchor[0] - xyzOneSensorActual[0]) ** 2 + \
                           (anchor[1] - xyzOneSensorActual[1]) ** 2 + \
                           (anchor[2]) ** 2)
        distanceSoil = abs(xyzOneSensorActual[2])
        power = p_average_air2soil(distanceAir, distanceSoil, ALPHA_SOIL, TAU)

        if (power > -110.0):
            # Only use the time of arrival measurement if the power is meaningful
            # mean = distance / speed_soil # TODO Introduce better randomness in the calculation
            # sigma = sigma_soil2Soil(distance, ALPHA_SOIL)
            mean = distanceAir / speed + distanceSoil / speed_soil
            sigma = sigma_air2Soil(distanceAir, distanceSoil, ALPHA_SOIL, TAU)
            print("Sampling t with mean = {} and sigma = {}".format(mean, sigma))
            ob = np.random.normal(loc = mean, scale = sigma, size = 1)[0]
            toa_from_anchors.append(ob)
        else:
            toa_from_anchors.append(np.nan)
            # We will randomize twice for sensitivity analysis
            # observed[i] = np.random.uniform(ob-DELTA_T, ob+DELTA_T)   #
    toa_observed_from_anchors.append(toa_from_anchors)
    print("toa_from_anchors = {}".format(toa_from_anchors))
toa_observed_from_anchors = np.asarray(toa_observed_from_anchors)
print("toa_observed_from_anchors.shape = {}".format(toa_observed_from_anchors.shape))

# TODO from here
#         observedTime3DSS = {}
#
#         # Sensor to sensor observations
#         #f = open("SensorNeighbors_" + str(epsilon_index) + "_" + str(NUM_SENSORS) + ".csv", "w")
#         for i in range(len(actual_s_locations)):
#             #f.write(str(i))
#             for j in range(len(actual_s_locations)):
#                 if (i == j):
#                     continue
#                 xyzi = actual_s_locations[i]
#                 xyzj = actual_s_locations[j]
#                 dist = sqrt((xyzi[0] - xyzj[0]) ** 2 + (xyzi[1] - xyzj[1]) ** 2 +(xyzi[2] - xyzj[2]) ** 2)
#                 signalPowerdB = p_average_soil2soil(dist, ALPHA_SOIL)
#
#                 if signalPowerdB > -110.0:
#                     #f.write("," + str(j))
#                     mean = dist / speed_soil
#                     sigma = sigma_soil2Soil(dist, ALPHA_SOIL)
#                     observed = numpy.random.normal(loc = mean, scale = sigma, size = 1)[0]
#                     observedTime3DSS[str(i) + '-' + str(j)] = numpy.random.uniform(observed-DELTA_T, observed+DELTA_T)
#             #f.write("\n")
#         #f.close()
#         # Estimation
#         xyzFirstEst = ()
x0 = np.random.random_sample(size=(NUM_SENSORS, 3)) * (F, F, -H/H)  # X, Y, Z coordinates
result = minimize(timeOfArrivalMatcherAnchorsToSensors,
                  x0,
                  args=(toa_observed_from_anchors, AnchorsXYZ, speed_soil),
                  method="Nelder-Mead",
                  options={"maxiter": 1e6})
        # for i in range(NUM_SENSORS):
        #     xyzInitialEstimate = (np.random.random() * F, np.random.random() * F, np.random.random() * (-H/H))
        #     res = minimize(timeOfArrivalMatcherAnchorToOneSensor,
        #                       xyzInitialEstimate,
        #                       args=(i, obsTAnchors, AnchorsXYZ, speed_soil),
        #                       method='Nelder-Mead', \
        #         options = {"maxiter":1e6})
#
#             e = res.x
#             if (e[2] > 0 ):    # Sensor height cannot be greater than ground level
#                 e[2] = -e[2]
#             #print res.success
#             xyzFirstEst = xyzFirstEst  + ((e[0], e[1], e[2]),)
#
#             #print actual_s_locations[i]
#             #print e
#             #print "*****"
#         #print actual_s_locations
#         #print "xyzFirstEst = ", xyzFirstEst
#         res = minimize(timeOfArrialMatcher3DX, xyzFirstEst, args=(obsTAnchors, AnchorsXYZ, observedTime3DSS, speed_soil), method='Nelder-Mead', options = {"maxiter":1e6})
#
#         e = [(res.x[i * 3: (i+1) * 3]) for i in range(NUM_SENSORS)]
#
#         #for i in range(len(actual_s_locations)):
#         #    print "||||||"
#         #    print actual_s_locations[i]
#         #    print e[i]
#
#         params_text = r"\noindent$\epsilon_s'$ = " + '{0:.3g}'.format(EPSILON_S_REAL/EPSILON_0) + r"\\\\" \
#         + r"$\epsilon_s''$ = " + '{0:.3g}'.format(EPSILON_S_IMG/EPSILON_0) + r"\\\\" \
#         + r"$\alpha^{(s)}$ = " + '{0:.3g}'.format(ALPHA_SOIL) + r" N/m\\\\"
#
#         xActual = [xyz[0] for xyz in actual_s_locations]
#         yActual = [xyz[1] for xyz in actual_s_locations]
#         zActual = [xyz[2] for xyz in actual_s_locations]
#
#         xFirstEst = [xyz[0] for xyz in xyzFirstEst]
#         yFirstEst = [xyz[1] for xyz in xyzFirstEst]
#         zFirstEst = [xyz[2] for xyz in xyzFirstEst]
#
#         xEst = [xyz[0] for xyz in e]
#         yEst = [xyz[1] for xyz in e]
#         zEst = [xyz[2] for xyz in e]
#
#         err = 0.0
#         for i in range(len(xActual)):
#             d = (xEst[i] - xActual[i]) ** 2 + \
#                 (yEst[i] - yActual[i]) ** 2 + \
#                 (zEst[i] - zActual[i]) ** 2
#             err = err + sqrt(d)
#         #if err > 10.0:
#         #    continue
#
#         #print "DELTA_T = ", DELTA_T, "error = ", err
#         #print "T_S = ", T_S, "error = ", err
#         #print ctr, ". epsilon_index = ", epsilon_index, "error = ", err
#         errors.append(err)
#         for i in range(NUM_SENSORS):
#             xEstSamples[i].append(e[i][0])
#             yEstSamples[i].append(e[i][1])
#             zEstSamples[i].append(e[i][2])
#         #if (ctr == 0):
#         #    f = open("estimates" + str(epsilon_index) + "_" + str(NUM_SENSORS) +".csv", "w")
#         #    f.write("X,Y,Z,X (est1),Y(est1),Z(est1),X(est),Y(est),Z(est)\n")
#         #    for i in range(NUM_SENSORS):
#         #        f.write(str(actual_s_locations[i][0]) + "," + str(actual_s_locations[i][1]) + "," + str(actual_s_locations[i][2]))
#         #        f.write("," + str(xyzFirstEst[i][0]) + "," + str(xyzFirstEst[i][1]) + "," + str(xyzFirstEst[i][2]))
#         #        f.write("," + str(e[i][0]) + "," + str(e[i][1]) + "," + str(e[i][2]) + "\n")
#         #    f.close()
#
#         #plot([xyz[0] for xyz in actual_s_locations], [xyz[1] for xyz in actual_s_locations], \
#         #        [xyz[0] for xyz in xyzFirstEst], [xyz[1] for xyz in xyzFirstEst], \
#         #    [xyz[0] for xyz in e], [xyz[1] for xyz in e], EPSILON_S_REAL, EPSILON_S_IMG, params_text, "X", "Y", [-1, F + 1 , -1, F + 1])
#         #plot([xyz[0] for xyz in actual_s_locations], [xyz[2] for xyz in actual_s_locations], \
#         #    [xyz[0] for xyz in xyzFirstEst], [xyz[2] for xyz in xyzFirstEst], \
#         #    [xyz[0] for xyz in e], [xyz[2] for xyz in e], EPSILON_S_REAL, EPSILON_S_IMG, params_text, "X", "Z", [-1, F + 1, -H/40, 0])
#
#     #print "X STD = ", numpy.mean([numpy.std(x) for x in xEstSamples])
#     #print "Y STD = ", numpy.mean([numpy.std(x) for x in yEstSamples])
#     #print "Z STD = ", numpy.mean([numpy.std(x) for x in zEstSamples])
#     print "delta_t = ", DELTA_T, "error = ", numpy.mean(errors) / NUM_SENSORS
#
#     f_x_std.write(str(epsilon_index) + "," +str(numpy.mean([numpy.std(x) for x in xEstSamples])) + "\n")
#     f_y_std.write(str(epsilon_index) + "," +str(numpy.mean([numpy.std(x) for x in yEstSamples])) + "\n")
#     f_z_std.write(str(epsilon_index) + "," +str(numpy.mean([numpy.std(x) for x in zEstSamples])) + "\n")
#     f_err.write(str(DELTA_T) + "," + str(numpy.mean(errors) / NUM_SENSORS) + "\n")
#
#     print "max xEst STD = ", numpy.max([numpy.std(x) for x in xEstSamples])
#     print "max yEst STD = ", numpy.max([numpy.std(x) for x in yEstSamples])
#     print "max zEst STD = ", numpy.max([numpy.std(x) for x in zEstSamples])
#     print "max error = ", numpy.max(errors) / NUM_SENSORS
# #
# # f_x_std.close()
# # f_y_std.close()
# # f_z_std.close()
# # f_err.close()
