import numpy as np
from math import pi
from scipy.optimize import minimize
from math import log10


H = 200
F = 100
NUM_SENSORS = 5
speed = 3. * 10**8
MU_0 = 4 * pi * 10 ** (-7)
MU_A = MU_0
MU_S = 1.0084 * MU_0
LAMBDA = 0.7
PA_0 = 38 # 38 dBm
PS_0 = 10  # 19 dBm
N_0 = -110.0  #dBm
T_S = 10.0 ** (-3)

# Scaling constants for optimization stability
scale = 100.
scale_lst_sq_obj_fn = 10**14

FREQUENCY = 833 * 10 ** 6
BETA_SQ = FREQUENCY ** 2
OMEGA = 2 * pi * FREQUENCY

EPSILON_0 = 8.85 * 10 ** (-12)
EPSILON_SOIL = np.asarray((2.361970519728 + 0.096670930496j,
                            5.08697889259241 + 0.4413468113884j,
                            16.4109595611802 + 2.36876126519685j,
                            24.4855102012741 + 3.85704851601056j)) * EPSILON_0
epsilon_index = 0

EPSILON_S_REAL = np.real(EPSILON_SOIL[epsilon_index])
EPSILON_S_IMG = np.imag([epsilon_index])
EPSILON_AIR = EPSILON_0

ALPHA_SOIL = OMEGA * \
             np.sqrt((MU_S * EPSILON_S_REAL/2) *
                  (np.sqrt(1 + (EPSILON_S_IMG/EPSILON_S_REAL) ** 2) - 1))

RHO = ((np.sqrt(MU_A / EPSILON_AIR) - np.sqrt(MU_S / EPSILON_S_REAL)) / \
    (np.sqrt(MU_A / EPSILON_AIR) + np.sqrt(MU_S / EPSILON_S_REAL))) ** 2

TAU = 1 - RHO

speed_soil = ( \
    np.sqrt( \
    (MU_S * EPSILON_S_REAL / 2 ) * (np.sqrt(1 + (EPSILON_S_IMG/EPSILON_S_REAL) ** 2) + 1)\
    ) ) ** (-1)

def p_average_air2soil(d_a, d_s, ALPHA_SOIL, TAU):
  pathLossAir = 20 * log10(4 * pi * d_a / LAMBDA)
  pathLossSoil = ALPHA_SOIL * d_s
  return PA_0 - pathLossAir - pathLossSoil

def sigma_air2Soil(d_a, d_s, ALPHA_SOIL, TAU):
  specularPowerdB = p_average_air2soil(d_a, d_s, ALPHA_SOIL, TAU)
  specular_power = 10 ** ((specularPowerdB - 30) / 10.0)
  TOTAL_NOISE = 10 ** ((N_0 - 30.0) / 10.0)
  sigma = np.sqrt(8 * (pi ** 2) * specular_power * T_S * BETA_SQ / TOTAL_NOISE) ** -1
  return sigma


# Least squares objective function
def toa_squared_error(xyz, toa_observed_from_anchors, anchor_locations, speed_soil):

    xyz=xyz.reshape(-1,3)*(1,1,1/scale)
    sq_error = 0.
    for i in range(xyz.shape[-2]):
        distance_air = np.sqrt(np.sum((xyz[i,:2]-anchor_locations[:,:2])**2, axis=1)+anchor_locations[:,2]**2)
        distance_soil = abs(xyz[i,2])
        estimated_toa = distance_air / speed + distance_soil / speed_soil
        sq_error += np.sum((toa_observed_from_anchors[i,:]-estimated_toa)**2)*scale_lst_sq_obj_fn
    return sq_error

if __name__ == "__main__":
    AnchorsXYZ = np.asarray([[0., 0., H],
                  [F, 0., H],
                  [0., F, H],
                  [F, F, H],
                  # (0., 0., -H / 4),
                  # (F, 0., -H / 4),
                  # (0., F, -H / 4),
                  # (F, F, -H / 4)
                  ])
    print("speed = {}, in soil = {}".format(speed, speed_soil))
    print(AnchorsXYZ)
    actual_s_locations = np.random.random_sample(size=(NUM_SENSORS, 3)) * (F, F, -H / H)  # X, Y, Z coordinates
    print("Actual sensor location = {}".format(actual_s_locations))
    xyz0 = np.random.random_sample(size=(NUM_SENSORS, 3)) * (F, F, -H / H)  # X, Y, Z coordinates
    print("Guessed sensor location = {}".format(xyz0))

    toa_observed_from_anchors = []
    for xyzOneSensorActual in actual_s_locations:
        toa_from_anchors = []
        for anchor in AnchorsXYZ[:4]:  # For anchors above soil
            # assert anchor.shape == xyzOneSensorActual.shape
            distance = np.sqrt(np.sum((anchor - xyzOneSensorActual) ** 2))

            # power = p_average_soil2soil(distance, ALPHA_SOIL)   # Received signal power for soil anchors
            distanceAir = np.sqrt((anchor[0] - xyzOneSensorActual[0]) ** 2 + \
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
                ob = np.random.normal(loc=mean, scale=sigma, size=1)[0]
                toa_from_anchors.append(ob)
            else:
                toa_from_anchors.append(np.nan)
                # We will randomize twice for sensitivity analysis
                # observed[i] = np.random.uniform(ob-DELTA_T, ob+DELTA_T)   #
        toa_observed_from_anchors.append(toa_from_anchors)

    toa_observed_from_anchors = np.asarray(toa_observed_from_anchors)
    print("toa_observed_from_anchors = {}".format(toa_observed_from_anchors))


    neg_likelihood = toa_squared_error(actual_s_locations * (1, 1, scale),
                                       toa_observed_from_anchors,
                                       AnchorsXYZ,
                                       speed_soil)

    print("neg_likelihood at minima = {}".format(neg_likelihood))
    bnds = (((None, None),)*2+((-100,0),))*NUM_SENSORS
    result = minimize(toa_squared_error,
                      xyz0 * (1,1,scale),
                      args=(toa_observed_from_anchors, AnchorsXYZ, speed_soil),
                      # method="Nelder-Mead",
                      method="L-BFGS-B",
                      bounds=bnds,
                      options={"maxiter":1e6})
    print(result.x.reshape(-1,3)*(1,1,1/scale))
    print(actual_s_locations)
    print(result.success)

