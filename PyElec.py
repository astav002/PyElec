

import numpy as np
import pandas as pd
import json
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from matplotlib import cm
import argparse
import os


class PyElec():
    def __init__(self):
        self.version = "2022.0.0.0"
        self.authoring = "Adam Stavola, Thomas Jefferson National Accelerator Facility"
        self.reference = "Based on work at JLAB: RCD-TBD-96 _001 (Stapleton) and  RCD-RPN-97 _001 (Degtiarenko)"

        print("******************************************")
        print("PyElect {}".format(self.version))
        print(self.authoring)
        print(self.reference)
        print("******************************************")
        # set program verbosity - 3 prints everything
        self.verbose = 3
        self.data_dir = os.path.join(os.getcwd(), "Data")

        # setup global defaults
        self.N_A = 6.022e23
        self.frac_scatter = 0.5
        self.low_e_limit = 500
        self.high_e_limit = 12000

        self.integral_low_energy = 100

        self.lambda_a = 30
        self.lambda_b = 55
        self.lambda_c = 120

        # hall weidth dimensions
        self.hall_a_w = 56.4
        self.hall_c_w = 48.62
        self.hall_b_w = 31.8

    def SheildingCoefficients(self, energy_mev):
        """
        Return the shielding coefficients based on energy of the beam
        """
        coeff_C = 0.000225
        coeff_B = 0.00225
        coeff_A = 0.009545

        energy_gev = energy_mev / 1000
        if (energy_gev  >= 1):
            pass
        elif (energy_gev < 0.5):
            coeff_B = coeff_B * 1.6 * energy_gev ** 1.5
        else:
            coeff_B = coeff_B * (0.566 + 0.434 * (energy_gev - 0.5) / 0.5)
            coeff_C = coeff_C * (-0.1668221 + 0.8889406 * energy_gev + 0.2638326 * energy_gev ** 2)      
            
        self.shld_coeff_A = coeff_A
        self.shld_coeff_B = coeff_B
        self.shld_coeff_C = coeff_C


    def RoofThickness(self, x, tr, td, ri):
        """
        This function computes the roof thickness as function of the angle x
        """
        # Public Function f6(x_) As Double
        # 'For comments look ELEC5b.FORTRAN code
        #     t = tr / Sqr(1 - (td * Sin(x_) / ri) ^ 2)
        #     f6 = a1t * Exp(-t / lamat) * Sin(x_)
        #     f6 = f6 + b1t * Exp(-t / lambt) * Sin(x_)
        #     f6 = f6 + c1t * Exp(-t / lamct) * Sin(x_)
        # End Function

        a = self.shld_coeff_A
        b = self.shld_coeff_B
        c = self.shld_coeff_C
        lambda_a = self.lambda_a
        lambda_b = self.lambda_b
        lambda_c = self.lambda_c

        # check to ensure not imaginary computation
        arg = 1 - (td*np.sin(x)/ri)**2
        arg = np.clip(arg, 1e-12, None)
        thickness = tr / np.sqrt(arg)        

        roof_thickness = np.sin(x) * (a * np.exp(-thickness /lambda_a) + 
                            b*np.exp(- thickness / lambda_b) + 
                            c * np.exp(- thickness / lambda_c))

        return roof_thickness



    def RadiationLength(self, z, a):
        """
        Calculate radiation length

        Not sure the purpose of this since we have it in lookup table, 
        but we can calculate the radiation length as the average of the following 
        
        x01 = 716.405 * at_ * (1 + 0.12 * ((z_ / 82) ^ 2)) / (z_ * (z_ + 1) * Log(183 * 1 / (z_ ^ 0.3333)))
        x02 = 716.405 * at_ / (z_ ^ 2 * Log(183 / (z_ ^ 0.3333)) + z_ * Log(1440 / (z_ ^ 0.6666)))
        f1 = (x01 + x02) / 2
        """

        x_o_1 = 716.405 * a * (1+0.12 * np.power(z/82, 2)) / (z*(z+1)*np.log(183 * 1 / np.power(z,0.3333)))
        x_o_2 = 716.405 * a / (np.power(z,2) * np.log(183 / np.power(z, 0.3333)) + z * np.log(1440 / np.power(z,  0.6666)))

        return (x_o_1 + x_o_2)/2

    def ScatterProduction(self,z, a, t, X_t_o, position, E_o_MeV,  current_uA, crit_radius_cm, crit_distance_cm, f_scatter=0.5):
        """
        The target specific values are passed in as arrays - verified with lead target 1000 mg/cm2 12/6/2021 AJS
        z - array of  z values
        a = array of A values
        t - array of target thicknesses in g/cm2
        X_t_o - array of thicknesses in  X_o
        position - array of the target locations in cm
        E_o_MeV - beam energy in MeV
        current_uA - beam current in uA
        crit_radius_cm - critical window radius in cm
        crit_distance_cm - critical window distance to pivot in cm

        Excel implementation
        f2 = 0.157 * ts_ * cur_ / (e0_ * thet_ ^ 2 + 0.157 * ts_ / e0_)
        f2 = f2 + (e0_ * cur_ - f2) * Exp(-thet_ ^ 2 * e0_ ^ 2 / (21.2 ^ 2 * tm_))
        """

        # we'll use numpy arrays since we can do the math without loops
        z = np.asarray(z)
        a = np.asarray(a)
        t = np.asarray(t)
        X_t_o = np.asarray(X_t_o)
        position = np.asarray(position)

        rel_t = t / X_t_o
        if (self.verbose == 3):
            print("Relative thickness X/Xo: {}".format(rel_t))

        C_prime = 0.157 * z*(z+4) * t / (a * np.power(E_o_MeV,2))
        C_prime = 1e-300 + np.sum(C_prime)

        if (self.verbose == 3):
            print("C_prime {}".format(C_prime))
        # calculate weighted position of the sum of "targets"
        
        sum_thick = np.sum(rel_t) + 1e-300
        pos = np.sum(rel_t * position) / sum_thick
        # determine theta as the larger of the two theta variables
        theta_1 = np.arctan(crit_distance_cm / (crit_distance_cm - pos))
        theta_2 = np.arctan(crit_radius_cm / crit_distance_cm)

        
        theta = theta_1 if theta_1 < theta_2 else theta_2

        if (self.verbose == 3):
            print("Theta {}, theta1 {}, theta2 {}".format(theta, theta_1, theta_2))
        
        # ok, now we can actually do the computation of the scatter production
        # we are going right from the technote on this, which differs from Excel
        multiplier = C_prime / (np.power(theta, 2) + C_prime)

        production = (E_o_MeV * current_uA* (multiplier + (1 - multiplier) * np.exp(-np.power(theta*E_o_MeV,2) / (449.4*sum_thick))))
        production = f_scatter * production

        return production
        

    def PhotoProduction(self, a, t, E_o_MeV, X_t_o, X_fe_o):
        """
        Calculated the photoneutron production - validated on Cu radiator to spreadsheet 3/5/2021 AJS
        a - atomic number of the element used to determin factor
        t - target thickness 
        E_o_MeV - photon energy (MeV)
        X_t_o - target radiation length (g/cm2)
        X_fe_o - iron radiation length (g/cm2)

        Formula:
        rel_t - target thicknesss (X/X_o)
        n = rel_t^2 * X_t_o /2 * integral(dk/k) from 100 to E_o
        d = 0.572 * X_fe_o * Eo * integral (dk / k^2) from 100 to E_o

        result = n / d 
        """
        photo_production = 0

        rel_t = t / X_t_o
        if (self.verbose == 3):
            print("Relative thickness X/Xo: {}".format(rel_t))

        # this is still a bit strange, but we sum the relative thicknesses then multiply 
        # by a factor
        # note factor for A <= 2 should be 0.5, we use 1 as default
        fac = np.ones(len(a))
        fac[a < 2] = 0.5
        numerator = (np.sum(rel_t) * 
                        np.sum(rel_t * X_t_o / 2 * 
                            np.log(E_o_MeV / self.integral_low_energy) * fac
                            )
                        )

        denominator =  - 0.572 * X_fe_o * E_o_MeV * (1/E_o_MeV - 1 / self.integral_low_energy)

        if (self.verbose == 3):
            print("Nump: {:2.4e}".format(numerator))        
            print("Denp: {:2.4e}".format(denominator))   
        
        photo_production = numerator / denominator 

        return photo_production


    def ElectroProduction(self, a, t, E_o_MeV, X_t_o, X_fe_o, scatter=1):
            """
            Calculated the electroneutron production - validated on Cu radiator to spreadsheet 3/5/2021 AJS
            a - atomic number of the target, used to determine fac
            t - tgt_thick_mg_cm2 / 1000
            X_t_o - target radiation length (g/cm2)
            X_fe_o - iron radiation length (g/cm2)
            E_o_MeV - beam energy (MeV)
            scatter - scatter factor

            Formula:

            result = numerator / denominator - numpy array
            """

            electro_production = 0
            fac = np.ones(len(a))

            # we don't include hydrogen in the electroprodution based on the excel calculation
            # there was a step that has:
            # if att(i) == 1 then goto 11
            # this skips the loselec for hydrogen, we'll just set fac =0 for this which makes it easy
            # interestingly the loselec equation in Excel (f3()) contained a statement 
            # if (a < 2) then fac = 0.5, but this couldn't happen due to program flow described already
            fac[a < 2] = 0           

            # Note: we modified the equation by noticing that rel_t * X_t_o = t
            numerator = self.tni * (t) * fac
            denominator = - 0.572 * X_fe_o * E_o_MeV * (1/E_o_MeV - 1 / self.integral_low_energy)
            if (self.verbose > 0):
                
                print("num_e: {}".format(np.sum(numerator)))
                print(numerator)
                print("den_p: {}".format(denominator))

            electro_production = numerator / denominator 

            return electro_production

    def get_element_prop(self, elem="H", fle="target_za_data.csv"):
        if (self.verbose > 0):
            print("Reading element data from {}".format(fle))
        tgt_za = pd.read_csv(os.path.join(self.data_dir, fle))

        if (self.verbose == 3):
            print(tgt_za)

        z = tgt_za.loc[tgt_za["material"] == elem, "Z_data"].values[0]
        a = tgt_za.loc[tgt_za["material"] == elem, "A_data"].values[0]

        if (self.verbose > 0):
            print("Returning {} values z: {}, a: {}".format(elem, z, a))

        return z, a

    def get_material_prop(self, elem="H", fle="mat_radlength_density.csv"):

        if (self.verbose > 0):
            print("Reading element data from {}".format(fle))

        mat_prop = pd.read_csv(os.path.join(self.data_dir, fle))

        if (self.verbose == 3):
            print(mat_prop)

        elm_Xo_g_cm2 = mat_prop.loc[mat_prop["material"] == elem, "Xo (g/cm2)"].values[0]
        elm_rho_g_cm3 = mat_prop.loc[mat_prop["material"] == elem, "density (g/cm3)"].values[0]

        if (self.verbose > 0):
            print("Returning {} values Xo: {}, rho: {}".format(elem, elm_Xo_g_cm2, elm_rho_g_cm3))

        return elm_Xo_g_cm2, elm_rho_g_cm3

    def get_hall_data(self, fle="hall_data.csv"):

        hall_df = pd.read_csv(os.path.join(self.data_dir, fle))
        self.hall_df = hall_df

        if (self.verbose == 3):
            print(hall_df)

    def get_data_params(self, fle="data_parameters.json"):

        f = open(os.path.join(self.data_dir, fle),'r')
        j_data = json.load(f)
        f.close()

        self.j_data = j_data

        if (self.verbose == 3):
            print(j_data)

    def run_lookup_params(self, energy):
        j_data = self.j_data

        self.tni = self.linear_interp(energy, j_data["ee0"], j_data["tni"])
        self.d = self.linear_interp(energy, j_data["e01"], j_data["d"])
        self.lam = self.linear_interp(energy, j_data["e02"], j_data["lam"])

        if (self.verbose > 0):

            print("interpolated tni: {:2.2e}".format(self.tni))
            print("interpolated d: {:2.2e}".format(self.d))
            print("interpolated lam: {:2.2e}".format(self.lam))

    def linear_interp(self, key, key_array, val_array):
        """
        Simple function for linear interpolation - 
        
        key: intermediate value at which the value should be interpolated to
        key_array: array of keys defining the independent variable
        val_array: array of values defining the dependent variable
        
        Note: could use built-in functions but we'll keep it explicit for now
        """

        # we'll ensure numpy arrays first
        key_array = np.asarray(key_array)
        val_array = np.asarray(val_array)

        if (key > key_array.max()):

            print ("BAD Key in interpolation --> Setting to maximum value for output")
            print( "\tNote this is not an error for energies > 6GeV, code was modified to accept max value")
            print("\tkey: {}, key_max: {}".format(key, key_array.max()))
            out_val = val_array.max()

        else:

            # find all entries where key is greater than array entries and get the max
            low_idx = np.argwhere(key > key_array).max()
            
            # find all entries where key is less than array entries and get min
            high_idx = np.argwhere(key < key_array).min()

            if (self.verbose == 3):
                print("Interpolated indicies {} - {}".format(low_idx, high_idx))

            slope = (val_array[high_idx] - val_array[low_idx]) / (key_array[high_idx] - key_array[low_idx])
            base = val_array[low_idx]

            # calculate the value based on the linear equation
            out_val = slope * (key - key_array[low_idx]) + base

        return out_val



    def RunForHall(self, hall="hall_a_rr",
            setup_config_name = "1",
            position_m = [0, 0],
            tgt_thick_mg_cm2 = [774, 774],
            current_uA = 50,
            energy_GeV = 2.2,
            element = ["Pb", "Pb"],
            critical_radius_cm = 3.125,
            critical_distance_m = 1.25,
            scatter_factor = 0.5,
            rbm_run = True,
            start_dist = 180,
            end_dist = 1000,
            hall_loc_x = 0,
            hall_loc_y = 0,
            verbose = 0
            ):



        # setup data common across all targets that could be installed
        self.get_hall_data()
        self.get_data_params()

        # perform interpoloations up front
        self.run_lookup_params(energy_GeV * 1000)


        # get information for the shielding thickness around the hall and limits of integration (b)
        hall_df = self.hall_df
        tr = hall_df.loc[hall_df["name"] == "tr", hall].values[0]
        td = hall_df.loc[hall_df["name"] == "td", hall].values[0]
        ri = hall_df.loc[hall_df["name"] == "ri", hall].values[0]
        b = hall_df.loc[hall_df["name"] == "b", hall].values[0]

        # calculate the shielding coefficients with the energy sent as MeV
        self.SheildingCoefficients(energy_GeV * 1000)

        # the original workbook implmented trapezoidal integration (TRAPZD function), based on 
        # Numerical Recipes in Fortran 77, Willam Press et al., Cambridge Press, 1986
        
        # the reference indicates an improvement by utilizing QTRAPZD.  We utilize
        # the python equivalent of Fortran QUADPACK which integrates the function
        # RoofThickness over the range 0 to b with argument tr, td, ri, which correspond
        # to hall dimensions
        q_integral, abserr = integrate.quad(self.RoofThickness, 0, b, args=(tr, td, ri))

        if (self.verbose > 0):
            print("Integral calculation yeilded: {}, abs error = {}".format(q_integral, abserr))

        q_integral = q_integral * 2 * np.pi / self.d

        # individual target calculations
        z = []
        a = []
        tgt_Xo = []
        frac_scatter = []

        for elm in element:
            z_, a_ = self.get_element_prop(elem=elm)
            z.append(z_)
            a.append(a_)

            # this is used in photo and electro production scatter factor
            if (a_ <= 2):
                frac_scatter.append(0.5)
            else:
                frac_scatter.append(1)

        # if we don't supply the radiation length in material prop we calculate it
            try:
                tgt_Xo_, tgt_rho_ = self.get_material_prop(elem=elm)
                tgt_Xo.append(tgt_Xo_)
            except:
                tgt_Xo_ = self.RadiationLength(z_,a_)
                tgt_Xo.append(tgt_Xo_)

        if (self.verbose ==3):
            print("tgt_Xo: {}".format(tgt_Xo))

        # ok, now we'll convert all of our lists to numpy arrays to make calculations without loops
        z = np.asarray(z)
        a = np.asarray(a)
        tgt_Xo = np.asarray(tgt_Xo)
        frac_scatter = np.asarray(frac_scatter)
        tgt_thick_mg_cm2 = np.asarray(tgt_thick_mg_cm2)
        position_m = np.asarray(position_m)

        # only need one iron result and we know we supply this
        fe_Xo, fe_rho = self.get_material_prop(elem="Fe")



        # note photo and electroproduction take numpy arrays, photoproduction returns single value
        # also, note frac_scatter is a per element term, while the scatter_factor term applies 
        # to the scattering production
        photo_n = self.PhotoProduction(a, tgt_thick_mg_cm2 / 1000, energy_GeV * 1000, tgt_Xo, fe_Xo)

        #electroproduction returns a numpy array which needs to be accumulated
        elect_n = self.ElectroProduction(a, tgt_thick_mg_cm2 / 1000, energy_GeV * 1000, tgt_Xo, fe_Xo)
        elect_n = np.sum(elect_n)

        # Note scatter production takes the full list of targets at once to compute result
        scat = self.ScatterProduction(z, 
                                        a, 
                                        tgt_thick_mg_cm2 / 1000, 
                                        tgt_Xo,
                                        position_m, 
                                        energy_GeV * 1000, 
                                        current_uA,
                                        critical_radius_cm,
                                        critical_distance_m * 100,
                                        f_scatter=scatter_factor
                                    )
        
        # the current and energy must be multiplied to obtain the total los, scatter already has this included
        los_dir = (photo_n + elect_n) * energy_GeV * 1000 * current_uA
        total_loss = los_dir + scat

        if (self.verbose > 0):
            print("Photo loss {}".format((photo_n)))
            print("Elec loss {}".format(np.sum(elect_n)))
            print("Scattering loss: {}".format(scat))
            print("Total loss: {}".format(total_loss))

        # last thing we need is the dose --> calculation follows exactly from PyELEC5b with update to perform at
        # multiple distances (~continuous) if rbm_run is set false


        # first determine rbm doses at given locations from the 
        location = []
        rbm_doses = []
        r = []


        for name in hall_df[hall_df["type"]=="location"]["name"].unique():
            #name = "RBM " + str(i)
            location.append(name)
            rr = hall_df.loc[hall_df['name'] == name, hall].values[0]

            r.append(rr)
            rbm_dose = q_integral * total_loss * 2e-15  * np.exp(-rr / self.lam) / np.power(rr+40, 2)
            rbm_doses.append(rbm_dose) 
            if (self.verbose > 0):
                print("Dose at {}: {} ".format(name, rbm_dose))
        r = np.asarray(r)
        rbm_doses = np.asarray(rbm_doses)



        # now we'll calculated a dense mapping for the dose rate in uR/hr to allow display if desired.
        dose_map = []
        x = np.linspace(start_dist, end_dist)
        y = np.linspace(start_dist, end_dist)
        ux, vy = np.meshgrid(x,y)
        dist = np.sqrt((ux - hall_loc_x)**2 + (vy-hall_loc_y)**2)
        dose_map = q_integral * total_loss * 2e-15  * np.exp(-dist / self.lam) / np.power(dist+40, 2) * 1e8


        # we say doses but really in dose rate! uR/hr after applying 1e-8 multiplier
        rbm_doses = rbm_doses * 1e8
        dose_df = pd.DataFrame({"location": location, "r":r, "dose_rate-"+str(setup_config_name):rbm_doses})


        return dose_df, dose_map
        

def main_test():

    # initialize the class
    pylec = PyElec()

    pylec.verbose = 3    

    x_default = 100
    y_default = 200


    hall_a_locx, hall_a_locy = (x_default, pylec.hall_a_w/2+y_default)

    # create some area around the hall to compute the doses for rbm_run = false
    start_dist = -100
    end_dist = 500



    a_dose_df, a_dose_map = pylec.RunForHall(hall="hall_a_rr", 
                                position_m = [0, -0.4, -0.4, 0, 0, 0, 0.5],#, 0, 0],
                                tgt_thick_mg_cm2 = [90, 46.9, 1.4, 26.2, 39.2, 44.8, 93.9],#, 220, 220],
                                current_uA = 60,
                                energy_GeV = 4.4,
                                element = ["He[3]", "Be", "Al", "N", "Si", "O", "Be"],#, "He", "He"],
                                critical_radius_cm = 4.45,
                                critical_distance_m = 1.33,
                                scatter_factor = 0.5,
                                start_dist=start_dist, 
                                end_dist=end_dist, 
                                hall_loc_y=56.4/2+200, 
                                hall_loc_x=100,
                                verbose=3,
                                rbm_run=True)    

    a_dose_df.to_csv("Hall_A-redo.csv")    
        



def main(config):

    pylec = PyElec()
    pylec.verbose = config["verbose"]

    x_default = 100
    y_default = 100


    hall_a_locx, hall_a_locy = (x_default, pylec.hall_a_w/2+y_default)

    # create some area around the hall to compute the doses for rbm_run = false
    start_dist = -2000
    end_dist = 2000

    res_dose_df = None
    total_map = None
    for i in range(len(config["setups"])):

        dose_df, dose_map = pylec.RunForHall(hall=config["setups"][i]["hall"]+ "_rr", 
                                            setup_config_name =config["setups"][i]["setup_name"],
                                            position_m =config["setups"][i]["position_m"],
                                            tgt_thick_mg_cm2 =config["setups"][i]["tgt_thick_mg_cm2"],
                                            current_uA =config["setups"][i]["current_uA"],
                                            energy_GeV =config["setups"][i]["energy_GeV"],
                                            element =config["setups"][i]['target_elem'],
                                            critical_radius_cm =config["setups"][i]["critical_radius_cm"],
                                            critical_distance_m =config["setups"][i]["critical_distance_m"],
                                            scatter_factor =config["setups"][i]["scatter_factor"],
                                            start_dist=start_dist, 
                                            end_dist=end_dist, 
                                            hall_loc_y=x_default, 
                                            hall_loc_x=y_default,
                                            verbose=3,
                                            rbm_run=True)     

        if (res_dose_df is None):
            res_dose_df = dose_df
            total_map = dose_map
        else:
            res_dose_df["dose_rate-" +config["setups"][i]["setup_name"]] = dose_df["dose_rate-" +config["setups"][i]["setup_name"]]
            total_map += dose_map
    
    dose_limit_ur = config["dose_limit_ur"]
    if (config["plot_map"]):
        
        fig, ax = plt.subplots()
        im = ax.imshow(np.log(np.asarray(total_map)), cmap='bone')
        # plt.contour((np.asarray(total_map)), 
        #         levels=dose_vals, cmap=cm.jet, origin='lower')
        lim_diff = (ax.get_xlim()[1] - ax.get_xlim()[0]) 
        circle_loc = x_default/(end_dist - start_dist) * (lim_diff) +lim_diff/2

        color = iter(cm.jet(np.linspace(0, 1, len(res_dose_df))))
        for r in res_dose_df.values:
            c = next(color)
            rad = r[1] * (lim_diff)/(end_dist - start_dist)  
            circle2 = plt.Circle((circle_loc, circle_loc), rad, color=c, fill=False, label=r[0])
            ax.add_patch(circle2)
        ax.text(circle_loc, circle_loc, "Hall A", color="tab:red", 
                bbox=dict(fill=True, edgecolor='red', linewidth=2, facecolor='white', boxstyle='round'))
        ax.legend(loc='upper left')
        print(ax.get_xlim())
        ax.set_xticks([])
        ax.set_yticks([])

        fig.colorbar(im)
  
        
        ax.set_title("Dose Rate Map (uR) Total - {}".format(config["plot_base_name"]))
        plt.savefig(config["plot_base_name"] + "-rbm_map_plot.png")

    if (config["plot_rbm"]):
        fontsize= 14
        plot_df = res_dose_df.drop("r", axis=1)

        fig, ax = plt.subplots()

        plot_df.plot.bar(x="location", ax=ax)
        ax.set_ylabel("Dose Rate $(\mu R/hr)$", fontsize=fontsize)
        plt.savefig(config["plot_base_name"] + "-rbm_bar_plot.png")
    
    print(res_dose_df)
    res_dose_df.to_csv(config["out_file"])


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Run PyElec with configuration arguments')

    parser.add_argument('--config', dest='config',
                        default=os.path.join(os.getcwd(), "Run_Configuration", "setup_single.json"),
                        help='sum the integers (default: find the max)')

    args = parser.parse_args()
    config_path = args.config

    try:
        f = open(config_path, 'r')
        
        run_config = json.load(f)
        f.close()
    except:
        print("Couldn't load json file, check your path or the file using a JSON linter")
    print(run_config.keys())

    main(run_config)




    
    