import numpy as np
from math import pi
from scipy.optimize import minimize
NUM_SENSORS = 5

scale_lst_sq_obj_fn = 10**23

H = 200
F = 100


T_S = 10.0 ** (-3)
PS_0 = 10  # 19 dBm
MU_0 = 4 * pi * 10 ** (-7)
MU_A = MU_0
MU_S = 1.0084 * MU_0
LAMBDA = 0.7
N_0 = -110.0  #dBm
FREQUENCY = 833 * 10 ** 6
BETA_SQ = FREQUENCY ** 2
z_scale = 100.

EPSILON_0 = 8.85 * 10 ** (-12)
EPSILON_SOIL = np.asarray((2.361970519728 + 0.096670930496j,
                            5.08697889259241 + 0.4413468113884j,
                            16.4109595611802 + 2.36876126519685j,
                            24.4855102012741 + 3.85704851601056j)) * EPSILON_0

epsilon_index = 0
EPSILON_S_REAL = np.real(EPSILON_SOIL[epsilon_index])
EPSILON_S_IMG = np.imag([epsilon_index])
EPSILON_AIR = EPSILON_0

speed_soil = ( \
    np.sqrt( \
    (MU_S * EPSILON_S_REAL / 2 ) * (np.sqrt(1 + (EPSILON_S_IMG/EPSILON_S_REAL) ** 2) + 1)\
    ) ) ** (-1)

OMEGA = 2 * pi * FREQUENCY

ALPHA_SOIL = OMEGA * \
             np.sqrt((MU_S * EPSILON_S_REAL/2) *
                  (np.sqrt(1 + (EPSILON_S_IMG/EPSILON_S_REAL) ** 2) - 1))

# Reflection coefficient
RHO = ((np.sqrt(MU_A / EPSILON_AIR) - np.sqrt(MU_S / EPSILON_S_REAL)) / \
    (np.sqrt(MU_A / EPSILON_AIR) + np.sqrt(MU_S / EPSILON_S_REAL))) ** 2


def p_average_soil2soil(d, ALPHA_SOIL):
    pathLoss = 20 * np.log10(4 * pi * d / LAMBDA)
    return PS_0 - d * ALPHA_SOIL -pathLoss #in dBm

def sigma_soil2Soil(d, ALPHA_SOIL):
  specularPowerdB = p_average_soil2soil(d, ALPHA_SOIL)
  specular_power = 10 ** ((specularPowerdB - 30) / 10.0)
  TOTAL_NOISE = 10 ** ((N_0 - 30.0) / 10.0)
  sigma = (8 * (pi ** 2) * specular_power * T_S * BETA_SQ / TOTAL_NOISE) ** -1
  return sigma


def toa_squared_error(xyz, toa_observed_from_other_sensors):
    xyz = xyz.reshape(-1,3)*(1,1,1/z_scale)
    sq_error = 0.
    for i in range(xyz.shape[-2]):
        for j in range(xyz.shape[-2]):
            if not toa_observed_from_other_sensors[i,j] > 0:
                continue
            dist = np.sqrt(np.sum((xyz[i] - xyz[j]) ** 2))
            mean_toa =dist/speed_soil
            sq_error += ((toa_observed_from_other_sensors[i,j]-mean_toa)**2)*scale_lst_sq_obj_fn
            # sigma_toa = sigma_soil2Soil(dist, ALPHA_SOIL)
    return sq_error

def p_average_soil2soil_reflected(d, ALPHA_SOIL):
    pathLoss = 20 * np.log10(4 * pi * d / LAMBDA)
    return PS_0 - d * ALPHA_SOIL - pathLoss - RHO #in dBm

def sigma_soil2Soil_reflected(d, ALPHA_SOIL):
  specularPowerdB = p_average_soil2soil_reflected(d, ALPHA_SOIL)
  specular_power = 10 ** ((specularPowerdB - 30) / 10.0)
  TOTAL_NOISE = 10 ** ((N_0 - 30.0) / 10.0)
  sigma = (8 * (pi ** 2) * specular_power * T_S * BETA_SQ / TOTAL_NOISE) ** -1
  return sigma


def toa_neg_log_likelihood_soil2soil(xyz, toa_observed_from_other_sensors_dir,
                                     toa_observed_from_other_sensors_reflected):
    xyz = xyz.reshape(-1, 3) * (1, 1, 1 / z_scale)
    neg_log_likelihood = 0.
    #for direct path
    for i in range(xyz.shape[-2]):
        for j in range(xyz.shape[-2]):
            if i == j:
                continue
            if toa_observed_from_other_sensors_dir[i, j] > 0 and toa_observed_from_other_sensors_reflected[i,j] > 0:
                dist = np.sqrt(np.sum((xyz[i] - xyz[j]) ** 2))
                mean_toa = dist / speed_soil
                sigma2_toa = sigma_soil2Soil(dist, ALPHA_SOIL)**2
                # neg_log_likelihood += ((toa_observed_from_other_sensors_dir[i, j] - mean_toa) ** 2) /sigma2_toa
                # neg_log_likelihood += 0.5*np.log(sigma2_toa)

                d1 = np.sqrt((xyz[i, 2]) ** 2 + \
                             ((xyz[i, 0] - xyz[j, 0]) ** 2 + \
                              (xyz[i, 1] - xyz[j, 1]) ** 2) *
                             ((xyz[i, 2] / (xyz[i, 2] + xyz[j, 2])) ** 2))

                d2 = np.sqrt((xyz[j, 2]) ** 2 + \
                             ((xyz[i, 0] - xyz[j, 0]) ** 2 + \
                              (xyz[i, 1] - xyz[j, 1]) ** 2) *
                             ((xyz[j, 2] / (xyz[i, 2] + xyz[j, 2])) ** 2))

                dist = d1 + d2
                mean_toa_ref = dist / speed_soil
                sigma2_toa_ref = sigma_soil2Soil_reflected(dist, ALPHA_SOIL) ** 2
                # neg_log_likelihood += ((toa_observed_from_other_sensors_reflected[i, j] - mean_toa_ref) ** 2) / sigma2_toa_ref
                # neg_log_likelihood += 0.5 * np.log(sigma2_toa_ref)

                t = np.asarray([toa_observed_from_other_sensors_dir[i, j], toa_observed_from_other_sensors_reflected[i, j]]).reshape(2,1)
                t_bar = np.asarray([mean_toa, mean_toa_ref])
                sigma = np.asarray([[1/sigma2_toa[0],0], [0,1/sigma2_toa_ref[0]]])

                a = sigma.dot(t-t_bar)
                b = (t-t_bar).T
                # neg_log_likelihood += ((t-t_bar).dot(sigma)).dot((t-t_bar).T)
                neg_log_likelihood += b.dot(a)
                neg_log_likelihood -= 0.5*np.log(np.linalg.det(sigma))
            else:
                raise Exception("test")

            # if toa_observed_from_other_sensors_reflected[i,j] > 0:
                # print("i={},j={}".format(i,j))


    return neg_log_likelihood


#direct path
if __name__=="__main__":

    actual_s_locations = np.loadtxt("data/act_sensor_loc.csv")
    est_s_locations_stage1 = np.loadtxt("data/est_stage1.csv")

    print("estimated locations (stage1) = \n{}".format(est_s_locations_stage1))
    print("actual locations = \n{}".format(actual_s_locations))

    toa_observed_from_other_sensors_dir = np.zeros(shape=(NUM_SENSORS, NUM_SENSORS))

    for i in range(actual_s_locations.shape[-2]):
        for j in range(actual_s_locations.shape[-2]):
            if i == j:
                continue
            dist = np.sqrt(np.sum((actual_s_locations[i]-actual_s_locations[j])**2))

            power = p_average_soil2soil(dist, ALPHA_SOIL)

            # print("power = {} dB".format(power))

            if power > N_0:
                mean_toa = dist / speed_soil
                sigma_toa = sigma_soil2Soil(dist, ALPHA_SOIL=ALPHA_SOIL)
                ob = np.random.normal(loc=mean_toa, scale=sigma_toa, size=1)[0]
                toa_observed_from_other_sensors_dir[i,j] = ob
            else:
                print("****************power alert***********")

    print("observed toa from other sensors (dir)=\n{}".format(toa_observed_from_other_sensors_dir))

    toa_observed_from_other_sensors_reflected = np.zeros(shape=(NUM_SENSORS, NUM_SENSORS))
    for i in range(actual_s_locations.shape[-2]):
        for j in range(actual_s_locations.shape[-2]):
            if i == j:
                continue

            d1 = np.sqrt((actual_s_locations[i,2])**2 + \
                         ((actual_s_locations[i,0]-actual_s_locations[j,0])**2 + \
            (actual_s_locations[i, 1] - actual_s_locations[j, 1]) ** 2) *
                         ((actual_s_locations[i,2] / (actual_s_locations[i,2]+actual_s_locations[j,2]))**2))

            d2 = np.sqrt((actual_s_locations[j,2])**2 + \
                         ((actual_s_locations[i,0]-actual_s_locations[j,0])**2 + \
            (actual_s_locations[i, 1] - actual_s_locations[j, 1]) ** 2) *
                         ((actual_s_locations[j,2] / (actual_s_locations[i,2]+actual_s_locations[j,2]))**2))

            dist = d1+d2

            power = p_average_soil2soil_reflected(dist, ALPHA_SOIL)

            # print("power = {} dB".format(power))

            if power > N_0:
                mean_toa = dist / speed_soil
                sigma_toa = sigma_soil2Soil_reflected(dist, ALPHA_SOIL=ALPHA_SOIL)
                ob = np.random.normal(loc=mean_toa, scale=sigma_toa, size=1)[0]
                toa_observed_from_other_sensors_reflected[i,j] = ob
            else:
                print("**************2.power alert***********")
    #
    print("observed toa from other sensors =\n{}".format(toa_observed_from_other_sensors_reflected))

    # xyz0 = np.random.random_sample(size=(NUM_SENSORS, 3)) * (F, F, -H / H)  # X, Y, Z coordinates
    xyz0 = est_s_locations_stage1


    sq_error_at_actual = toa_neg_log_likelihood_soil2soil(actual_s_locations*(1,1,z_scale),
                                                          toa_observed_from_other_sensors_dir,
                                                          toa_observed_from_other_sensors_reflected)
    print("squared error at actual = {}".format(sq_error_at_actual))


    sq_error_at_seed = toa_neg_log_likelihood_soil2soil(xyz0*(1,1,z_scale),
                                                           toa_observed_from_other_sensors_dir,
                                                        toa_observed_from_other_sensors_reflected)
    print("squared error at seed = {}".format(sq_error_at_seed))

    bnds = (((None, None),)*2+((-100,0),))*NUM_SENSORS

    result = minimize(toa_neg_log_likelihood_soil2soil,
                      xyz0*(1,1,z_scale),
                      args=(toa_observed_from_other_sensors_dir,
                            toa_observed_from_other_sensors_reflected),
                      method="L-BFGS-B",
                      bounds=bnds,
                      options={"maxiter":1e6})
    print(result.success)
    est_s_locations_stage2a = result.x.reshape(-1, 3) * (1, 1, 1 / z_scale)
    print("actual =\n{}".format(actual_s_locations))
    # print("stage 1=\n{}".format(est_s_locations_stage1))
    print("stage 2 =\n{}".format(est_s_locations_stage2a))
    np.savetxt("data/est_stage2a.csv", est_s_locations_stage2a)


