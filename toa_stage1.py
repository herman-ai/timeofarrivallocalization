import numpy as np
from math import pi
from scipy.optimize import minimize
from math import log10

speed = 3. * 10**8
NUM_RUNS = 1
MU_0 = 4 * pi * 10 ** (-7)
EPSILON_0 = 8.85 * 10 ** (-12)

LAMBDA = 0.7
PA_0 = 38 # 38 dBm
PS_0 = 10  # 19 dBm
N_0 = -110.0  #dBm
T_S = 10.0 ** (-3)
z_scale = 100.

# Scaling constants for optimization stability
scale_lst_sq_obj_fn = 10**14

class ToaLocalizer:

    EPSILON_SOIL = np.asarray((2.361970519728 + 0.096670930496j,
                                5.08697889259241 + 0.4413468113884j,
                                16.4109595611802 + 2.36876126519685j,
                                24.4855102012741 + 3.85704851601056j)) * EPSILON_0
    
    def __init__(self, epsilon_index, f):
        e_s_real = np.real(ToaLocalizer.EPSILON_SOIL[epsilon_index])
        e_s_img = np.imag([epsilon_index])
        e_a = EPSILON_0
        mu_a = MU_0
        mu_s = 1.0084 * MU_0
        omega = 2 * pi * f
        self.beta_sq = f ** 2
        self.alpha_s = omega * \
                       np.sqrt((mu_s * e_s_real/2) *
                          (np.sqrt(1 + (e_s_img/e_s_real) ** 2) - 1))

        self.rho = ((np.sqrt(mu_a / e_a) - np.sqrt(mu_s / e_s_real)) / \
            (np.sqrt(mu_a / e_a) + np.sqrt(mu_s / e_s_real))) ** 2   # Reflection coefficient

        self.tau = 1 - self.rho   # Transmission coefficient

        self.speed_soil = ( \
            np.sqrt( \
            (mu_s * e_s_real / 2 ) * (np.sqrt(1 + (e_s_img/e_s_real) ** 2) + 1)\
            )) ** (-1)

    def p_average_air2soil(self, d_a, d_s, ALPHA_SOIL=None):
      pathLossAir = 20 * np.log10(4 * pi * d_a * self.tau / LAMBDA)
      if ALPHA_SOIL is None:
        ALPHA_SOIL = self.alpha_s
      pathLossSoil = ALPHA_SOIL * d_s
      return PA_0 - pathLossAir - pathLossSoil

    def sigma_air2Soil(self, d_a, d_s, ALPHA_SOIL=None):
      if ALPHA_SOIL is None:
        ALPHA_SOIL = self.alpha_s
      specularPowerdB = self.p_average_air2soil(d_a, d_s, ALPHA_SOIL=ALPHA_SOIL)
      specular_power = 10 ** ((specularPowerdB - 30) / 10.0)
      TOTAL_NOISE = 10 ** ((N_0 - 30.0) / 10.0)
      sigma = np.sqrt(8 * (pi ** 2) * specular_power * T_S * self.beta_sq / TOTAL_NOISE) ** -1
      return sigma


    # Least squares objective function
    def toa_squared_error_air2soil(self, xyz, toa_observed_from_anchors, anchor_locations, speed_soil):

        xyz=xyz.reshape(-1,3)*(1, 1, 1 / z_scale)
        sq_error = 0.
        for i in range(xyz.shape[-2]):
            distance_air = np.sqrt(np.sum((xyz[i,:2]-anchor_locations[:,:2])**2, axis=1)+anchor_locations[:,2]**2)
            distance_soil = abs(xyz[i,2])
            mean_toa = distance_air / speed + distance_soil / speed_soil
            sq_error += np.sum((toa_observed_from_anchors[i,:]-mean_toa)**2)*scale_lst_sq_obj_fn
        return sq_error

    def toa_neg_log_likelihood_air2soil(self, xyz, toa_observed_from_anchors, anchor_locations, speed_soil, ALPHA_SOIL=None):
        if ALPHA_SOIL is None:
            ALPHA_SOIL = self.alpha_s
        xyz = xyz.reshape(-1, 3) * (1, 1, 1 / z_scale)
        neg_log_likelihood = 0
        for i in range(xyz.shape[-2]):
            distance_air = np.sqrt(
                np.sum((xyz[i, :2] - anchor_locations[:, :2]) ** 2, axis=1) + anchor_locations[:, 2] ** 2)
            distance_soil = abs(xyz[i, 2])

            mean_toa = distance_air / speed + distance_soil / speed_soil
            sigma2_toa = self.sigma_air2Soil(distance_air, distance_soil, ALPHA_SOIL)**2

            neg_log_likelihood += np.sum(((toa_observed_from_anchors[i,:]-mean_toa)**2)/sigma2_toa)
            neg_log_likelihood += np.sum(np.log(sigma2_toa))
        return neg_log_likelihood

if __name__ == "__main__":
    H = 200
    F = 100
    NUM_SENSORS = 5
    num_samples = 0
    xy_error = 0
    z_error = 0
    np.random.seed(0)
    localizer = ToaLocalizer(epsilon_index=0, f=833 * 10 ** 6)    # epsilon index = 0
    speed_soil = localizer.speed_soil
    for i in range(NUM_RUNS):
        AnchorsXYZ = np.asarray([[0., 0., H],
                      [F, 0., H],
                      [0., F, H],
                      [F, F, H]
                      ])
        actual_s_locations = np.random.random_sample(size=(NUM_SENSORS, 3)) * (F, F, -H / H)  # X, Y, Z coordinates
        # print("actual_s_locations.shape = {}".format(actual_s_locations.shape))
        # print("Actual sensor location = {}".format(actual_s_locations))
        xyz0 = np.random.random_sample(size=(NUM_SENSORS, 3)) * (F, F, -H / H)  # X, Y, Z coordinates
        # print("Guessed sensor location = {}".format(xyz0))

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
                power = localizer.p_average_air2soil(distanceAir, distanceSoil)
                # print("power = {} dB".format(power))
                if (power > -110.0):
                    mean = distanceAir / speed + distanceSoil / speed_soil
                    sigma = localizer.sigma_air2Soil(distanceAir, distanceSoil)
                    # print("Sampling t with mean = {} and sigma = {}".format(mean, sigma))
                    ob = np.random.normal(loc=mean, scale=sigma, size=1)[0]
                    toa_from_anchors.append(ob)
                else:
                    toa_from_anchors.append(np.nan)
                    # We will randomize twice for sensitivity analysis
                    # observed[i] = np.random.uniform(ob-DELTA_T, ob+DELTA_T)   #
            toa_observed_from_anchors.append(toa_from_anchors)

        toa_observed_from_anchors = np.asarray(toa_observed_from_anchors)

        # print("toa_observed_from_anchors = {}".format(toa_observed_from_anchors))


        # neg_likelihood = localizer.toa_neg_log_likelihood_air2soil(actual_s_locations * (1, 1, z_scale),
        #                                                  toa_observed_from_anchors,
        #                                                  AnchorsXYZ,
        #                                                  speed_soil)
        #
        # # print("neg_likelihood at minima = {}".format(neg_likelihood))
        # neg_likelihood = localizer.toa_neg_log_likelihood_air2soil(xyz0 * (1, 1, z_scale),
        #                                                  toa_observed_from_anchors,
        #                                                  AnchorsXYZ,
        #                                                  speed_soil)

        # print("neg_likelihood at other = {}".format(neg_likelihood))


        bnds = (((None, None),)*2+((-100,0),))*NUM_SENSORS
        result = minimize(localizer.toa_neg_log_likelihood_air2soil,
                          xyz0 * (1, 1, z_scale),
                          args=(toa_observed_from_anchors, AnchorsXYZ, speed_soil),
                          # method="Nelder-Mead",
                          method="L-BFGS-B",
                          bounds=bnds,
                          options={"maxiter":1e6})
        estimated_locations = result.x.reshape(-1,3) * (1, 1, 1 / z_scale)
        #
        # np.savetxt("data/est_stage1.csv", estimated_locations)
        # np.savetxt("data/act_sensor_loc.csv", actual_s_locations)


        if result.success == True:
            print("num_samples = {}".format(num_samples))
            print("estimated locations = \n{}".format(estimated_locations))
            print("actual locations = \n{}".format(actual_s_locations))
            num_samples += 1
            xy_error += np.sum(np.sqrt(np.sum((estimated_locations[:,:2]-actual_s_locations[:,:2])**2, axis=1)))
            z_error  += np.sum(np.abs(estimated_locations[:,2]-actual_s_locations[:,2]))
            # print("z_error = {}".format(np.sum(np.abs(estimated_locations[:,2]-actual_s_locations[:,2]))))

    print("xy error = {}".format(xy_error/(num_samples*NUM_SENSORS)))
    print("z error = {}".format(z_error/(num_samples*NUM_SENSORS)))

