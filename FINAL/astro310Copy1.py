# Python routines for Astronomy 310 / Stellar Astrophysics
#
# (c) 2017 Rich Townsend, Jacqueline Goldstein
#
# Version 1: original
# Version 2: removed usetex, font specifications

# Imports

import numpy as np
import re
import zipfile as zp

#### EZ-Web routines ####

# Read an EZ-Web history file

def read_history (filename) :

    """
    Read history data from an EZ-Web output zip file
    
    Parameters
    ----------
    
    filename : string giving name of zip file
    
    Returns
    -------

    data: dictionary containing the history data. The following keys/values 
          appear in the dictionary:

        i         step index (gives index of corresponding profile)
        t         age (years)
        M         mass (Msun)
        log_L     log10(luminosity / Lsun)
        log_R     lop10(radius / Rsun)
        logT_s    log10(surface temperature / K)
        log_T_c   log10(central temprature / K)
        log_rho_c log10(central density / kg/m^3)
        log_P_c   log10(central pressure / N/m^2)
        Psi_c     central electron degeneracy parameter
        X_c       central hydrogen mass fraction
        Y_c       central helium mass fraction
        X_Cc      central carbon mass fraction
        X_Nc      central nitrogen mass fraction
        X_Oc      central oxgyen mass fraction
        tau_dyn   dynamical timescale (s)
        tau_KH    Kelvin-Helmholtz timescale (years)
        tau_nuc   nuclear timescale (years)
        L_PP      luminosity from PP chain (Lsun)
        L_CNO     luminosity from CNO cycle (Lsun)
        L_3a      luminosity from triple-alpha reactions (Lsun)
        L_Z       luminosity from metal burning (Lsun)
        L_nu      luminosity from neutrino losses (Lsun)
        M_He      mass of helium core (Msun)
        M_C       mass of carbon core (Msun)
        M_O       mass of oxygen core (Msun)
        R_He      radius of helium core (Rsun)
        R_C       radius of carbon core (Rsun)
        R_O       radius of oxygen core (Rsun)

    """
    
    col_names = ['i', 't', 'M', 'log_L', 'log_R', 'log_T_s', 
                 'log_T_c', 'log_rho_c', 'log_P_c', 'Psi_c', 
                 'X_c', 'Y_c', 'X_Cc', 'X_Nc', 'X_Oc', 
                 'tau_dyn', 'tau_KH', 'tau_nuc', 'L_PP',
                 'L_CNO', 'L_3a', 'L_Z', 'L_nu', 'M_He',
                 'M_C', 'M_O', 'R_He', 'R_C', 'R_O']

    return read_ez_zip(filename, 'summary.txt', col_names)

# Read an EZ-Web profile

def read_profile (filename, index):

    """
    Read profile data from an EZ-Web output zip file
    
    Parameters
    ----------
    
    filename : string giving name of zip file
    index    : integer step index
    
    Returns
    -------

    data: dictionary containing the profile data. The following keys/values 
          appear in the dictionary:

        m         mass coordinate (Msun)
        r         radius coordinate (Rsun)
        F         interior luminosity (Lsun)
        P         pressure (N/m^2)
        rho       density (kg/m^3)
        T         temperature (K)
        u         specific internal energy (J/kg)
        s         specific entropy (J/K/kg)
        c_P       specific heat at constant pressure (J/K/kg)
        gamma_ad  adiabatic exponent
        nabla_ad  adiabatic temperature gradient
        mu        mean molecular weight
        n_e       electron number density (/m^3)
        P_e       electron partial pressure (N/m^2)
        P_rad     radiation partial pressure (N/m^2)
        nabla_rad radiative temperature gradient
        nabla     material temperature gradient
        v_c       convective velocity (m/s)
        kappa     Rosseland mean opacity (m^2/kg)
        q_nuc     specific nuclear energy release rate (J/s/kg)
        q_PP      specific nuclear energy release rate from PP chain (J/s/kg)
        q_CNO     specific nuclear energy release rate from CNO cycle (J/s/kg)
        q_3a      specific nuclear energy release rate from triple-alpha reactions (J/s/kg)
        q_nunuc   specific nuclear neutrino energy loss rate (J/s/kg)
        q_nu      specific non-nuclear neutrino energy loss rate (J/s/kg)
        q_grav    specific energy release rate from gravitational contraction (J/s/kg)
        X         hydrogen mass fraction (all ionization stages)
        X_mol     molecular hydrogen mass fraction
        X_+       singly-ionized hydrogen mass fraction
        Y         helium mass fraction (all ionization stages)
        Y_+       singly-ionized helium mass fraction
        Y_++      doubly-ionized helium mass fraction
        X_C       carbon mass fraction
        X_N       nitrogen mass fraction
        X_O       oxygen mass fraction
        Psi       electron degeneracy parameter
    """
    
    col_names = ['m', 'r', 'F', 'P', 'rho', 'T', 'u', 's', 
                 'c_P', 'gamma_ad', 'nabla_ad', 'mu', 'n_e', 
                 'P_e', 'P_rad', 'nabla_rad', 'nabla', 'v_c',
                 'kappa', 'q_nuc', 'q_PP', 
                 'q_CNO', 'q_3a', 'q_nunuc',
                 'q_nu', 'q_grav', 'X', 'X_mol',
                 'X_+', 'Y', 'Y_+', 'Y_++', 'X_C', 'X_N', 
                 'X_O', 'Psi']

    return read_ez_zip(filename, 'structure_{:05d}.txt'.format(index), col_names)

# Read data from an EZ-Web zipfile

def read_ez_zip(zip_filename, data_filename, col_names):

    # Open the zipfile

    z = zp.ZipFile(zip_filename, 'r')

    # Open the data file and read data

    with z.open(data_filename, 'r') as f:
        table = np.loadtxt(f)

    data = dict()

    for i in range(0,len(col_names)):
        data[col_names[i]] = table[:,i]

    return data


#### Plotting Routines ####

import matplotlib as mpl

def plot_defaults(width=6, height=6, fontsize=12, legend_fontsize=12):

    """
    Set up plot parameters
    
    Parameters
    ----------
    
    width           : width of figure (inches)
    height          : height of figure (inches)
    fontsize        : font size for axis labels (points)
    legend_fontsize : font size for legend (points)
    
    Returns
    -------

    Nothing

    """

    params = {'backend': 'pdf',
              'figure.figsize': [width, height],
              'font.size': fontsize,
              'axes.titlesize': 'medium',
              'axes.labelsize': 'medium',
              'legend.fontsize': legend_fontsize,
              'legend.frameon' : False,
              'figure.dpi': 600,
              'lines.markersize': 4,
              'lines.linewidth': 1,
              'lines.antialiased': False,
              'path.simplify': False }
              
    mpl.rcParams.update(params)
