#----------------------------------------------------------------------------------------
#
# This file is part of CosmicFish.
#
# Copyright (C) 2015-2016 by the CosmicFish authors
#
# The CosmicFish code is free software;
# You can use it, redistribute it, and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# The full text of the license can be found in the file LICENSE at
# the top level of the CosmicFish distribution.
#
#----------------------------------------------------------------------------------------

"""

Simple Python code to perform analysis of Information Gain for several experiments

Developed by Marco Raveri (mraveri@sissa.it) and
Matteo Martinelli (m.martinelli@thphys.uni-heidelberg.de) for the CosmicFish code.

"""

# ***************************************************************************************

# import first dependencies:
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines    as mlines
import numpy as np
import argparse
import math
import sys
import os
import copy
import itertools as it
import ConfigParser

# get the path of the application and the CosmicFish library:
here          = os.path.dirname(os.path.abspath(__file__))
cosmicfish_dir = os.environ.get('COSMICFISH_DIR')
# check:
if not cosmicfish_dir:
    print 'Please export the environment variable COSMICFISH_DIR for the script to work correctly.'
    exit(1)
# add to the python path the library:
cosmicfish_pylib_path = cosmicfish_dir+'/python'
sys.path.insert(0, os.path.normpath(cosmicfish_pylib_path))

# import the CosmicFish pylib
import cosmicfish_pylib.utilities            as fu
import cosmicfish_pylib.colors               as fc
import cosmicfish_pylib.fisher_matrix        as fm
import cosmicfish_pylib.fisher_derived       as fd
import cosmicfish_pylib.fisher_operations    as fo
import cosmicfish_pylib.fisher_plot_settings as fps
import cosmicfish_pylib.fisher_plot_analysis as fpa
import cosmicfish_pylib.fisher_plot          as fp

# ***************************************************************************************

""" Hard coded options """

bit_to_GB     = 1.25e-10
main_fontsize = 10
x_size        = 18.0
y_size        = 29.7
dot_size      = 4
dostat        = False
filling       = True
alpha         = 0.2
wide_flat_prior_coefficient = 0.0

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

data_dir        = here+'/../raw_results'
additional_data = here+'/../additional_fishers'

# ***************************************************************************************

fu.CosmicFish_write_header(' LCDM Information Gain plotter')
# parameters considered:
param_names   = [ 'omegabh2', 'omegach2', 'h', 'logA', 'ns', 'tau' ]
# define model names:
model_names   = [ 'mnu', 'r', 'w0wa', 'Linder', 'EFT', 'Horava', 'EFT2' ]
# define model param names:
model_param_names = [ ['omeganuh2'],
                      [ 'r' ],
                      [ 'w0_ppf', 'wa_ppf' ],
                      [ 'Linder_gamma'],
                      [ 'EFTOmega0'],
                      [ 'Horava_lambda', 'Horava_eta' ],
                      [ 'EFTOmega0', 'EFTGamma10', 'EFTGamma20', 'EFTGamma30' ],
                      ]
# write the header:
fu.CosmicFish_write_header(' Dark Energy and Gravity Information Gain plotter')
# create two dictionaries with the parameter names:
model_param_names = dict( zip(model_names, model_param_names) )
param_names_total = dict( zip(model_names, [ param_names + model_param_names[names] for names in model_names ] ) )
# import the fisher matrices:
fishers = fpa.CosmicFish_FisherAnalysis( fisher_path=[data_dir,additional_data], search_fisher_guess=True, with_derived=False )
# model fishers marginalized over other parameters:
model_fishers = dict( zip(model_names, [ fishers.marginalise( params=model_param_names[names], names=[ fish for fish in fishers.get_fisher_name_list() if '_'+names+'_' in fish ] ) for names in model_names] ) )
# total fisher marginalizing over other parameters:
model_fishers = dict( zip(model_names, [ fishers.marginalise( params=model_param_names[names], names=[ fish for fish in fishers.get_fisher_name_list() if '_'+names+'_' in fish ] ) for names in model_names] ) )
full_fishers  = dict( zip(model_names, [ fishers.marginalise( params=param_names_total[names], names=[ fish for fish in fishers.get_fisher_name_list() if '_'+names+'_' in fish ] ) for names in model_names] ) )
# remove model fishers:
for name in model_names:
    fishers.delete_fisher_matrix( names=[ fish for fish in fishers.get_fisher_name_list() if '_'+name+'_' in fish ]  )
# marginalize the left over fisher matrices:
fishers      = fishers.marginalise( params=param_names )
fisher_prior = [ fish for fish in fishers.get_fisher_list() if '_Prior_' in fish.name ]
if len(fisher_prior) == 0:
    print 'Warning no base prior Fisher matrix found'
    exit(1)
elif len(fisher_prior) > 1:
    print 'Warning too many prior Fisher matrices found'
    exit(1)
fisher_prior = fisher_prior[0]
fisher_wide_flat_prior = fm.fisher_matrix( fisher_matrix=wide_flat_prior_coefficient*fisher_prior.get_fisher_matrix(),
                                           param_names=fisher_prior.get_param_names(),
                                           param_names_latex=fisher_prior.get_param_names_latex(),
                                           fiducial=fisher_prior.get_param_fiducial() )
# get the prior Fisher matrices:
full_fisher_prior = []
for name in model_names:
    fisher_temp = [ fish for fish in full_fishers[name].get_fisher_list() if '_Prior_' in fish.name ]
    if len(fisher_temp) == 0:
        print 'Warning no prior Fisher matrix found for model: ', name
        exit(1)
    elif len(fisher_temp) > 1:
        print 'Warning too many prior Fisher matrices found for model: ', name
        exit(1)
    full_fisher_prior.append(fisher_temp[0])
# get the prior dictionary:
model_fisher_prior = dict( zip(model_names, full_fisher_prior ) )
full_fisher_prior  = dict( zip(model_names, [ fish + fisher_prior for fish in full_fisher_prior ] ) )
# create the flat wide prior:
model_wide_flat_prior = dict( zip(model_names, [ fm.fisher_matrix( fisher_matrix=wide_flat_prior_coefficient*model_fisher_prior[name].get_fisher_matrix(),
                                                  param_names=model_fisher_prior[name].get_param_names(),
                                                  param_names_latex=model_fisher_prior[name].get_param_names_latex(),
                                                  fiducial=model_fisher_prior[name].get_param_fiducial() )
                                                  for name in model_names ]))
full_wide_flat_prior  = dict( zip(model_names, [ fm.fisher_matrix( fisher_matrix=wide_flat_prior_coefficient*full_fisher_prior[name].get_fisher_matrix(),
                                                  param_names=full_fisher_prior[name].get_param_names(),
                                                  param_names_latex=full_fisher_prior[name].get_param_names_latex(),
                                                  fiducial=full_fisher_prior[name].get_param_fiducial() )
                                                  for name in model_names ]))

# do CMB:

labels_CMB = ['COBE','WMAP','Planck 2015','Simons Array','CORE']

years_CMB = [ 1992, 2005, 2014.5, 2018, 2030 ]

names = [ '1_COBE_fisher_matrix_cls_marginal',
          '2_WMAP_9yr_fisher_matrix_cls_marginal',
          '3_Planck_2015_fisher_matrix_cls_marginal',
          '5_Simons_fisher_matrix_cls_marginal',
          '4_CORE_fisher_matrix_cls_marginal' ]

names_models = {}

names_models['mnu'] = [ '1_COBE_mnu_fisher_matrix_cls_marginal',
              '2_WMAP_9yr_mnu_fisher_matrix_cls_marginal',
              '3_Planck_2015_mnu_fisher_matrix_cls_marginal',
              '5_Simons_mnu_fisher_matrix_cls_marginal',
              '4_CORE_mnu_fisher_matrix_cls_marginal' ]

names_models['r'] = [ '1_COBE_r_fisher_matrix_cls_marginal',
            '2_WMAP_9yr_r_fisher_matrix_cls_marginal',
            '3_Planck_2015_r_fisher_matrix_cls_marginal',
            '5_Simons_r_fisher_matrix_cls_marginal',
            '4_CORE_r_fisher_matrix_cls_marginal' ]

names_models['w0wa'] = [ '1_COBE_w0wa_fisher_matrix_cls_marginal',
              '2_WMAP_9yr_w0wa_fisher_matrix_cls_marginal',
              '3_Planck_2015_w0wa_fisher_matrix_cls_marginal',
              '5_Simons_w0wa_fisher_matrix_cls_marginal',
              '4_CORE_w0wa_fisher_matrix_cls_marginal' ]

names_models['Linder'] = [ '1_COBE_Linder_fisher_matrix_cls_marginal',
            '2_WMAP_9yr_Linder_fisher_matrix_cls_marginal',
            '3_Planck_2015_Linder_fisher_matrix_cls_marginal',
            '5_Simons_Linder_fisher_matrix_cls_marginal',
            '4_CORE_Linder_fisher_matrix_cls_marginal' ]

names_models['EFT'] = [ '1_COBE_EFT_fisher_matrix_cls_marginal',
            '2_WMAP_9yr_EFT_fisher_matrix_cls_marginal',
            '3_Planck_2015_EFT_fisher_matrix_cls_marginal',
            '5_Simons_EFT_fisher_matrix_cls_marginal',
            '4_CORE_EFT_fisher_matrix_cls_marginal' ]

names_models['Horava'] = [ '1_COBE_Horava_fisher_matrix_cls_marginal',
              '2_WMAP_9yr_Horava_fisher_matrix_cls_marginal',
              '3_Planck_2015_Horava_fisher_matrix_cls_marginal',
              '5_Simons_Horava_fisher_matrix_cls_marginal',
              '4_CORE_Horava_fisher_matrix_cls_marginal' ]

names_models['EFT2'] = [ '1_COBE_EFT2_fisher_matrix_cls_marginal',
            '2_WMAP_9yr_EFT2_fisher_matrix_cls_marginal',
            '3_Planck_2015_EFT2_fisher_matrix_cls_marginal',
            '5_Simons_EFT2_fisher_matrix_cls_marginal',
            '4_CORE_EFT2_fisher_matrix_cls_marginal' ]


data_from_prior_CMB     = np.array( [ fo.information_gain( fishers.get_fisher_matrix(name)[0] , fisher_prior, fisher_wide_flat_prior, stat=dostat ) for name in names ] )

data_from_prior_CMB_alt = {}
for m_name in model_names:
    data_from_prior_CMB_alt[m_name] = np.array( [ fo.information_gain( full_fishers[m_name].get_fisher_matrix(name)[0] , full_fisher_prior[m_name], full_wide_flat_prior[m_name], stat=dostat ) for name in names_models[m_name] ] )
model_data_from_prior_CMB_alt = {}
for m_name in model_names:
    model_data_from_prior_CMB_alt[m_name] = np.array( [ fo.information_gain( model_fishers[m_name].get_fisher_matrix(name)[0] , model_fisher_prior[m_name], model_wide_flat_prior[m_name], stat=dostat ) for name in names_models[m_name] ] )

# do SN:

labels_SN = ['Low z SN','Hi z SN','HST SN','DETF SIV SN ground','DETF SIV SN space']

years_SN = [ 1996, 1998, 2003, 2022, 2035 ]

data_from_prior_SN = np.array( [
         fo.information_gain( fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]   , fisher_prior, fisher_wide_flat_prior, stat=dostat ),
         fo.information_gain( fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0] + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0] , fisher_prior, fisher_wide_flat_prior, stat=dostat ),
         fo.information_gain( fishers.get_fisher_matrix('9_SN_HST_fisher_matrix_SN_marginal')[0]    , fisher_prior, fisher_wide_flat_prior, stat=dostat ),
         fo.information_gain( fishers.get_fisher_matrix('11_SN_DETFIV_ground_fisher_matrix_SN_marginal')[0]  , fisher_prior, fisher_wide_flat_prior, stat=dostat ),
         fo.information_gain( fishers.get_fisher_matrix('10_SN_DETFIV_space_fisher_matrix_SN_marginal')[0]   , fisher_prior, fisher_wide_flat_prior, stat=dostat ),
         ] )

data_from_prior_SN_alt = {}
for m_name in model_names:
    data_from_prior_SN_alt[m_name] = np.array( [
             fo.information_gain( full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]   , full_fisher_prior[m_name], full_wide_flat_prior[m_name], stat=dostat ),
             fo.information_gain( full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0] + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0] , full_fisher_prior[m_name], full_wide_flat_prior[m_name], stat=dostat ),
             fo.information_gain( full_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]    , full_fisher_prior[m_name], full_wide_flat_prior[m_name], stat=dostat ),
             fo.information_gain( full_fishers[m_name].get_fisher_matrix('11_SN_DETFIV_ground_'+m_name+'_fisher_matrix_SN_marginal')[0]  , full_fisher_prior[m_name], full_wide_flat_prior[m_name], stat=dostat ),
             fo.information_gain( full_fishers[m_name].get_fisher_matrix('10_SN_DETFIV_space_'+m_name+'_fisher_matrix_SN_marginal')[0]   , full_fisher_prior[m_name], full_wide_flat_prior[m_name], stat=dostat ),
             ] )
model_data_from_prior_SN_alt = {}
for m_name in model_names:
    model_data_from_prior_SN_alt[m_name] = np.array( [
             fo.information_gain( model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]   , model_fisher_prior[m_name], model_wide_flat_prior[m_name], stat=dostat ),
             fo.information_gain( model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0] + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0] , model_fisher_prior[m_name], model_wide_flat_prior[m_name], stat=dostat ),
             fo.information_gain( model_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]    , model_fisher_prior[m_name], model_wide_flat_prior[m_name], stat=dostat ),
             fo.information_gain( model_fishers[m_name].get_fisher_matrix('11_SN_DETFIV_ground_'+m_name+'_fisher_matrix_SN_marginal')[0]  , model_fisher_prior[m_name], model_wide_flat_prior[m_name], stat=dostat ),
             fo.information_gain( model_fishers[m_name].get_fisher_matrix('10_SN_DETFIV_space_'+m_name+'_fisher_matrix_SN_marginal')[0]   , model_fisher_prior[m_name], model_wide_flat_prior[m_name], stat=dostat ),
             ] )

# 3 Hubble:

labels_Hubble = ['Hubble $H_0$', 'HST $H_0$', '$H_0$ 2016']

years_Hubble  = [ 1927 , 2001, 2016 ]

data_from_prior_Hubble = np.array( [
         fo.information_gain( fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0]      , fisher_prior, fisher_wide_flat_prior, stat=dostat ),
         fo.information_gain( fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]  , fisher_prior, fisher_wide_flat_prior, stat=dostat ),
         fo.information_gain( fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0] , fisher_prior, fisher_wide_flat_prior, stat=dostat ),
         ] )

# 4 LSS
labels_LSS = ['DETF SIV LSS space','DETF SIV LSS ground']

years_LSS = [ 2020, 2028 ]

data_from_prior_LSS = np.array( [
        fo.information_gain( fishers.get_fisher_matrix('14_DETFIV_GC_WL_space_opt_Simons_fisher_matrix_cls_marginal')[0]  , fisher_prior, fisher_wide_flat_prior, stat=dostat ),
        fo.information_gain( fishers.get_fisher_matrix('13_DETFIV_GC_WL_ground_opt_Simons_fisher_matrix_cls_marginal')[0] , fisher_prior, fisher_wide_flat_prior, stat=dostat ),
        ])

data_from_prior_LSS_alt = {}
for m_name in model_names:
    data_from_prior_LSS_alt[m_name] = np.array( [
            fo.information_gain( full_fishers[m_name].get_fisher_matrix('14_DETFIV_GC_WL_space_opt_Simons_'+m_name+'_fisher_matrix_cls_marginal')[0]  , full_fisher_prior[m_name], full_wide_flat_prior[m_name], stat=dostat ),
            fo.information_gain( full_fishers[m_name].get_fisher_matrix('13_DETFIV_GC_WL_ground_opt_Simons_'+m_name+'_fisher_matrix_cls_marginal')[0] , full_fisher_prior[m_name], full_wide_flat_prior[m_name], stat=dostat ),
             ] )
model_data_from_prior_LSS_alt = {}
for m_name in model_names:
    model_data_from_prior_LSS_alt[m_name] = np.array( [
             fo.information_gain( model_fishers[m_name].get_fisher_matrix('14_DETFIV_GC_WL_space_opt_Simons_'+m_name+'_fisher_matrix_cls_marginal')[0]  , model_fisher_prior[m_name], model_wide_flat_prior[m_name], stat=dostat ),
             fo.information_gain( model_fishers[m_name].get_fisher_matrix('13_DETFIV_GC_WL_ground_opt_Simons_'+m_name+'_fisher_matrix_cls_marginal')[0] , model_fisher_prior[m_name], model_wide_flat_prior[m_name], stat=dostat ),
             ] )

# 5 RD

labels_RD = ['E-ELT RD']

years_RD = [ 2032 ]

data_from_prior_RD = np.array( [
        fo.information_gain( fishers.get_fisher_matrix('12_RD_ELT_fisher_matrix_RD_marginal')[0] , fisher_prior, fisher_wide_flat_prior, stat=dostat ),
])

data_from_prior_RD_alt = {}
for m_name in model_names:
    data_from_prior_RD_alt[m_name] = np.array( [
            fo.information_gain( full_fishers[m_name].get_fisher_matrix('12_RD_ELT_'+m_name+'_fisher_matrix_RD_marginal')[0] , full_fisher_prior[m_name], full_wide_flat_prior[m_name], stat=dostat ),
            ] )
model_data_from_prior_RD_alt = {}
for m_name in model_names:
    model_data_from_prior_RD_alt[m_name] = np.array( [
             fo.information_gain( model_fishers[m_name].get_fisher_matrix('12_RD_ELT_'+m_name+'_fisher_matrix_RD_marginal')[0] , model_fisher_prior[m_name], model_wide_flat_prior[m_name], stat=dostat ),
             ] )

# add all together:
label = labels_CMB + labels_SN + labels_Hubble + labels_LSS + labels_RD
years = years_CMB + years_SN + years_Hubble + years_LSS + years_RD

label = [ y for (x,y) in sorted( zip(years, label) ) ]

data_from_prior     = np.concatenate( (data_from_prior_CMB, data_from_prior_SN, data_from_prior_Hubble, data_from_prior_LSS, data_from_prior_RD) )
data_from_prior     = [ y for (x,y) in sorted( zip(years, data_from_prior) ) ]
data_from_prior_alt       = {}
model_data_from_prior_alt = {}
for m_name in model_names:
    temp1  =  np.concatenate( (data_from_prior_CMB_alt[m_name], data_from_prior_SN_alt[m_name], data_from_prior_Hubble, data_from_prior_LSS_alt[m_name], data_from_prior_RD_alt[m_name]) )
    temp2  =  np.concatenate( (model_data_from_prior_CMB_alt[m_name], model_data_from_prior_SN_alt[m_name], data_from_prior_Hubble, model_data_from_prior_LSS_alt[m_name], model_data_from_prior_RD_alt[m_name]) )
    data_from_prior_alt[m_name]       = [ y for (x,y) in sorted( zip(years, temp1) ) ]
    model_data_from_prior_alt[m_name] = [ y for (x,y) in sorted( zip(years, temp2) ) ]

years = sorted( years )

# do cumulative Information Gain:

fisher_array = [
                # Start with Hubble H0:
                fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0] ,
                # Add COBE:
                fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('1_COBE_fisher_matrix_cls_marginal')[0] ,
                # Add LowZ SN:
                fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('1_COBE_fisher_matrix_cls_marginal')[0] ,
                # Add HiZ SN:
                fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('1_COBE_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0] ,
                # Add HST H0:
                fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('1_COBE_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0] ,
                # Add HST SN:
                fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('1_COBE_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('9_SN_HST_fisher_matrix_SN_marginal')[0] ,
                # Add WMAP:
                fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('2_WMAP_9yr_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('9_SN_HST_fisher_matrix_SN_marginal')[0] ,
                # Add Planck:
                fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('3_Planck_2015_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('9_SN_HST_fisher_matrix_SN_marginal')[0] ,
                # Add Riess Hubble:
                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('3_Planck_2015_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('9_SN_HST_fisher_matrix_SN_marginal')[0] ,
                # Add Simons array:
                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('5_Simons_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('9_SN_HST_fisher_matrix_SN_marginal')[0] ,
                # Add DETF LSS:
                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('9_SN_HST_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('14_DETFIV_GC_WL_space_opt_Simons_fisher_matrix_cls_marginal')[0] ,
                # Add DETF SN:
                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('9_SN_HST_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('11_SN_DETFIV_ground_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('14_DETFIV_GC_WL_space_opt_Simons_fisher_matrix_cls_marginal')[0] ,
                # Add Space DETF:
                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('9_SN_HST_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('11_SN_DETFIV_ground_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('13_DETFIV_GC_WL_ground_opt_Simons_fisher_matrix_cls_marginal')[0] ,
                # Add CORE:
                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('9_SN_HST_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('11_SN_DETFIV_ground_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('15_DETFIV_GC_WL_ground_opt_CORE_fisher_matrix_cls_marginal')[0] ,
                # Add RD:
                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('9_SN_HST_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('11_SN_DETFIV_ground_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('12_RD_ELT_fisher_matrix_RD_marginal')[0]
                + fishers.get_fisher_matrix('15_DETFIV_GC_WL_ground_opt_CORE_fisher_matrix_cls_marginal')[0] ,
                # Add Space DETF SN:
                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                + fishers.get_fisher_matrix('6_SN_lowZ_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('7_SN_SDSS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('8_SN_SNLS_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('9_SN_HST_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('11_SN_DETFIV_ground_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('10_SN_DETFIV_space_fisher_matrix_SN_marginal')[0]
                + fishers.get_fisher_matrix('12_RD_ELT_fisher_matrix_RD_marginal')[0]
                + fishers.get_fisher_matrix('15_DETFIV_GC_WL_ground_opt_CORE_fisher_matrix_cls_marginal')[0] ,
]

fisher_array_alt = {}
for m_name in model_names:
    fisher_array_alt[m_name] =[
                                # Start with Hubble H0:
                                fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0] ,
                                # Add COBE:
                                fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('1_COBE_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add LowZ SN:
                                fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('1_COBE_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add HiZ SN:
                                fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('1_COBE_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add HST H0:
                                fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('1_COBE_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add HST SN:
                                fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('1_COBE_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add WMAP:
                                fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('2_WMAP_9yr_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add Planck:
                                fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('3_Planck_2015_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add Riess Hubble:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('3_Planck_2015_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add Simons array:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('5_Simons_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add DETF LSS:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('14_DETFIV_GC_WL_space_opt_Simons_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add DETF SN:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('11_SN_DETFIV_ground_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('14_DETFIV_GC_WL_space_opt_Simons_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add Space DETF:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('11_SN_DETFIV_ground_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('13_DETFIV_GC_WL_ground_opt_Simons_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add CORE:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('11_SN_DETFIV_ground_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('15_DETFIV_GC_WL_ground_opt_CORE_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add RD:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('11_SN_DETFIV_ground_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('12_RD_ELT_'+m_name+'_fisher_matrix_RD_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('15_DETFIV_GC_WL_ground_opt_CORE_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add Space DETF SN:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('11_SN_DETFIV_ground_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('10_SN_DETFIV_space_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('12_RD_ELT_'+m_name+'_fisher_matrix_RD_marginal')[0]
                                + full_fishers[m_name].get_fisher_matrix('15_DETFIV_GC_WL_ground_opt_CORE_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                ]

model_fisher_array_alt = {}
for m_name in model_names:
    model_fisher_array_alt[m_name] =[
                                # Start with Hubble H0:
                                fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0] ,
                                # Add COBE:
                                fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('1_COBE_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add LowZ SN:
                                fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('1_COBE_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add HiZ SN:
                                fishers.get_fisher_matrix('1_Hubble_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('1_COBE_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add HST H0:
                                fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('1_COBE_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add HST SN:
                                fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('1_COBE_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add WMAP:
                                fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('2_WMAP_9yr_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add Planck:
                                fishers.get_fisher_matrix('2_Hubble_HST_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('3_Planck_2015_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add Riess Hubble:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('3_Planck_2015_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add Simons array:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('5_Simons_'+m_name+'_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0] ,
                                # Add DETF LSS:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('14_DETFIV_GC_WL_space_opt_Simons_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add DETF SN:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('11_SN_DETFIV_ground_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('14_DETFIV_GC_WL_space_opt_Simons_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add Space DETF:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('11_SN_DETFIV_ground_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('13_DETFIV_GC_WL_ground_opt_Simons_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add CORE:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('11_SN_DETFIV_ground_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('15_DETFIV_GC_WL_ground_opt_CORE_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add RD:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('11_SN_DETFIV_ground_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('12_RD_ELT_'+m_name+'_fisher_matrix_RD_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('15_DETFIV_GC_WL_ground_opt_CORE_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                # Add Space DETF SN:
                                fishers.get_fisher_matrix('3_Hubble_2016_fisher_matrix_cls_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('6_SN_lowZ_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('7_SN_SDSS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('8_SN_SNLS_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('9_SN_HST_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('11_SN_DETFIV_ground_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('10_SN_DETFIV_space_'+m_name+'_fisher_matrix_SN_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('12_RD_ELT_'+m_name+'_fisher_matrix_RD_marginal')[0]
                                + model_fishers[m_name].get_fisher_matrix('15_DETFIV_GC_WL_ground_opt_CORE_'+m_name+'_fisher_matrix_cls_marginal')[0] ,
                                ]

# get cumulative data:

years_cumulative  = years #[ 1927, 1992, 1996, 1998, 2001, 2003, 2005, 2015, 2016, 2018, 2020, 2022, 2028, 2030, 2031, 2035 ]
cumulative_info_gain = [ fo.information_gain( fish, fisher_prior, fisher_wide_flat_prior, stat=dostat ) for fish in fisher_array ]

cumulative_info_gain_alt = {}
for m_name in model_names:
    cumulative_info_gain_alt[m_name] = [ fo.information_gain( fish, full_fisher_prior[m_name], full_wide_flat_prior[m_name], stat=dostat ) for fish in fisher_array_alt[m_name] ]
model_cumulative_info_gain_alt = {}
for m_name in model_names:
    model_cumulative_info_gain_alt[m_name] = [ fo.information_gain( fish, model_fisher_prior[m_name], model_wide_flat_prior[m_name], stat=dostat ) for fish in model_fisher_array_alt[m_name] ]

# get differential data:
years_differential  = years

differential_info_gain = [ fo.information_gain( fisher_array[0], fisher_prior, fisher_prior, stat=dostat ) ]
for ind in xrange(1,len(fisher_array)) :
    differential_info_gain.append( fo.information_gain( fisher_array[ind], fisher_array[ind-1], fisher_prior, stat=dostat ) )

differential_info_gain_alt = {}
for m_name in model_names:
    differential_temp = [ fo.information_gain( fisher_array_alt[m_name][0], full_fisher_prior[m_name], full_fisher_prior[m_name], stat=dostat ) ]
    for ind in xrange(1,len(fisher_array_alt[m_name])):
        differential_temp.append( fo.information_gain( fisher_array_alt[m_name][ind], fisher_array_alt[m_name][ind-1], full_fisher_prior[m_name], stat=dostat ) )
    differential_info_gain_alt[m_name] = differential_temp

model_differential_info_gain_alt = {}
for m_name in model_names:
    differential_temp = [ fo.information_gain( model_fisher_array_alt[m_name][0], full_fisher_prior[m_name], full_fisher_prior[m_name], stat=dostat ) ]
    for ind in xrange(1,len(fisher_array_alt[m_name])):
        differential_temp.append( fo.information_gain( model_fisher_array_alt[m_name][ind], model_fisher_array_alt[m_name][ind-1], full_fisher_prior[m_name], stat=dostat ) )
    model_differential_info_gain_alt[m_name] = differential_temp
