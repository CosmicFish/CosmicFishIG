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

# load all common stuff:
from load_fisher import *

# prepare:
fig = matplotlib.pyplot.gcf()
fig.set_size_inches( x_size/2.54, y_size/2.54 )

gs = gridspec.GridSpec(3, 5)

ax_3_p1 = plt.subplot(gs[2, 0])
ax_3_p2 = plt.subplot(gs[2, 1:], sharey=ax_3_p1 )

ax_2_p1 = plt.subplot(gs[1, 0], sharex=ax_3_p1)
ax_2_p2 = plt.subplot(gs[1, 1:], sharey=ax_2_p1, sharex=ax_3_p2 )

ax_1_p1 = plt.subplot(gs[0, 0], sharex=ax_3_p1)
ax_1_p2 = plt.subplot(gs[0, 1:], sharey=ax_1_p1, sharex=ax_3_p2 )


# vertical lines corresponding to data:
for plot in [ax_1_p1,ax_1_p2,ax_2_p1,ax_2_p2,ax_3_p1,ax_3_p2]:
    for year in years:
        dashes = [2,2,2,2]  # 10 points on, 5 off, 100 on, 5 off
        line = plot.axvline( year, color='k', linestyle='--', linewidth=0.5)
        line.set_dashes(dashes)
    dashes_2 = [2,2,2,2]  # 10 points on, 5 off, 100 on, 5 off
    line_2 = plot.axvline( 2016, color='red', linestyle='--', linewidth=0.5)
    line_2.set_dashes(dashes_2)

# plot single experiment prior information gain:
# 2- mnu:
ax_1_p1.plot( years, model_data_from_prior_alt['mnu'], '-o', markersize=dot_size, color=fc.nice_colors(1) )
ax_1_p2.plot( years, model_data_from_prior_alt['mnu'], '-o', markersize=dot_size, color=fc.nice_colors(1) )
# 3- r:
ax_1_p1.plot( years, model_data_from_prior_alt['r'], '-o', markersize=dot_size, color=fc.nice_colors(2) )
ax_1_p2.plot( years, model_data_from_prior_alt['r'], '-o', markersize=dot_size, color=fc.nice_colors(2) )

# plot cumulative experiment prior information gain:
# 2- mnu:
ax_2_p1.plot( years_cumulative, model_cumulative_info_gain_alt['mnu'], '-o', markersize=dot_size, color=fc.nice_colors(1) )
ax_2_p2.plot( years_cumulative, model_cumulative_info_gain_alt['mnu'], '-o', markersize=dot_size, color=fc.nice_colors(1) )
# 3- r:
ax_2_p1.plot( years_cumulative, model_cumulative_info_gain_alt['r'], '-o', markersize=dot_size, color=fc.nice_colors(2) )
ax_2_p2.plot( years_cumulative, model_cumulative_info_gain_alt['r'], '-o', markersize=dot_size, color=fc.nice_colors(2) )

# plot differential experiment prior information gain:
# 2- mnu:
ax_3_p1.plot( years_differential, model_differential_info_gain_alt['mnu'], '-o', markersize=dot_size, color=fc.nice_colors(1) )
ax_3_p2.plot( years_differential, model_differential_info_gain_alt['mnu'], '-o', markersize=dot_size, color=fc.nice_colors(1) )
# 3- r:
ax_3_p1.plot( years_differential, model_differential_info_gain_alt['r'], '-o', markersize=dot_size, color=fc.nice_colors(2) )
ax_3_p2.plot( years_differential, model_differential_info_gain_alt['r'], '-o', markersize=dot_size, color=fc.nice_colors(2) )

ax_1_p1.set_ylabel('bits', fontsize=main_fontsize)
ax_2_p1.set_ylabel('bits', fontsize=main_fontsize)
ax_3_p1.set_ylabel('bits', fontsize=main_fontsize)

# set log plot:
ax_1_p1.set_yscale('log')
ax_1_p2.set_yscale('log')
ax_2_p1.set_yscale('log')
ax_2_p2.set_yscale('log')
ax_3_p1.set_yscale('log')
ax_3_p2.set_yscale('log')

# limits:
ax_1_p1.set_xlim( 1920, 1940 )
ax_1_p2.set_xlim( 1990, 2040 )

ax_2_p1.set_xlim( 1920, 1940 )
ax_2_p2.set_xlim( 1990, 2040 )

ax_3_p1.set_xlim( 1920, 1940 )
ax_3_p2.set_xlim( 1990, 2040 )

# ticks with experiments:
ax_1_p1_temp = ax_1_p1.twiny()
ax_1_p1_temp.set_xlim( ax_1_p1.get_xlim() )
ax_1_p1_temp.spines['right'].set_visible(False)
ax_1_p1_temp.set_xticks( [1920, 1940 ] )

ax_1_p2_temp = ax_1_p2.twiny()
ax_1_p2_temp.set_xlim( ax_1_p2.get_xlim() )
ax_1_p2_temp.spines['left'].set_visible(False)

fig.canvas.draw()
ax_1_p1_temp.set_xticklabels( [item.get_text() for item in ax_1_p1_temp.get_xticklabels()], fontsize=0.9*main_fontsize )
ax_1_p2_temp.set_xticklabels( [item.get_text() for item in ax_1_p2_temp.get_xticklabels()], fontsize=0.9*main_fontsize )

ax_1_p1_temp.xaxis.get_majorticklabels()[0].set_horizontalalignment('left')
ax_1_p1_temp.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right')

ax_1_p2_temp.xaxis.get_majorticklabels()[0].set_horizontalalignment('left')
ax_1_p2_temp.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right')

# hide what needs to be hidden:

# the vertical axes:
ax_1_p1.spines['right'].set_visible(False)
ax_1_p2.spines['left'].set_visible(False)
ax_2_p1.spines['right'].set_visible(False)
ax_2_p2.spines['left'].set_visible(False)
ax_3_p1.spines['right'].set_visible(False)
ax_3_p2.spines['left'].set_visible(False)

# move ticks to the proper position
ax_1_p1.yaxis.tick_left()
ax_1_p2.yaxis.tick_right()
ax_2_p1.yaxis.tick_left()
ax_2_p2.yaxis.tick_right()
ax_3_p1.yaxis.tick_left()
ax_3_p2.yaxis.tick_right()

ax_1_p1.tick_params(labelbottom='off')
ax_1_p2.tick_params(labelright='off',labelbottom='off')
ax_2_p1.tick_params(labelbottom='off')
ax_2_p2.tick_params(labelright='off',labelbottom='off')
ax_3_p1.tick_params()
ax_3_p2.tick_params(labelright='off')

# set the ticks and ticks labels:

ax_3_p1.set_xticks( [years[0]] )
ax_3_p2.set_xticks( years[1:] )

ax_3_p1.set_xticklabels( [label[0]], fontsize=0.9*main_fontsize, rotation=45, horizontalalignment='right' )
ax_3_p2.set_xticklabels( label[1:], fontsize=0.9*main_fontsize, rotation=45, horizontalalignment='right' )

# draw the canvas to adjust the ticks:
fig.canvas.draw()
ax_1_p1.set_yticklabels( [item.get_text() for item in ax_1_p1.get_yticklabels()], fontsize=0.9*main_fontsize )
ax_2_p1.set_yticklabels( [item.get_text() for item in ax_2_p1.get_yticklabels()], fontsize=0.9*main_fontsize )
ax_3_p1.set_yticklabels( [item.get_text() for item in ax_3_p1.get_yticklabels()], fontsize=0.9*main_fontsize )

ax_1_p1.get_yticklabels()[0].set_verticalalignment('bottom')
ax_1_p1.get_yticklabels()[-1].set_verticalalignment('top')
ax_2_p1.get_yticklabels()[0].set_verticalalignment('bottom')
ax_2_p1.get_yticklabels()[-1].set_verticalalignment('top')
ax_3_p1.get_yticklabels()[0].set_verticalalignment('bottom')
ax_3_p1.get_yticklabels()[-1].set_verticalalignment('top')

# update structure and sizes:
gs.update( bottom=0.12, top=0.9, left=0.1, right=0.98, wspace=0.05, hspace=0.10)

ax_1_p2.set_title('Prior Model Information Gain', fontsize=main_fontsize, loc='right',y=1.15 )
ax_2_p2.set_title('Cumulative Prior Model Information Gain', fontsize=main_fontsize, loc='right')
ax_3_p2.set_title('Differential Model Information Gain', fontsize=main_fontsize, loc='right')

# legend:
names  = [ '$\\Lambda$CDM', '$\\Lambda$CDM +mnu', '$\\Lambda$CDM +r' ]
colors = [ fc.nice_colors(0), fc.nice_colors(1), fc.nice_colors(2) ]

leg_handlers = []
for name, color in zip( names, colors ):
    leg_handlers.append( mlines.Line2D([], [], color = color ) )

fig.legend( handles  = leg_handlers,
            labels   = names,
            fontsize = main_fontsize,
            frameon  = True,
            fancybox = True,
            ncol     = 3,
            borderaxespad = 0.0,
            loc = 'lower center',
            )

# export:
plt.savefig(here+'/../results/4_model_information_gain_LCDM.pdf', transparent=True)
plt.clf()

# do the single plotting:
for index in xrange(3):
    # prepare:
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches( x_size/2.54, y_size/2.54/2.5 )

    gs = gridspec.GridSpec(1, 5)

    ax_1_p1 = plt.subplot(gs[0, 0])
    ax_1_p2 = plt.subplot(gs[0, 1:], sharey=ax_1_p1 )

    # vertical lines corresponding to data:
    for plot in [ax_1_p1,ax_1_p2]:
        for year in years:
            dashes = [2,2,2,2]  # 10 points on, 5 off, 100 on, 5 off
            line = plot.axvline( year, color='k', linestyle='--', linewidth=0.5)
            line.set_dashes(dashes)
        dashes_2 = [2,2,2,2]  # 10 points on, 5 off, 100 on, 5 off
        line_2 = plot.axvline( 2016, color='red', linestyle='--', linewidth=0.5)
        line_2.set_dashes(dashes_2)

    # plot single experiment prior information gain:
    if index==0:
        # 2- mnu:
        ax_1_p1.plot( years, model_data_from_prior_alt['mnu'], '-o', markersize=dot_size, color=fc.nice_colors(1) )
        ax_1_p2.plot( years, model_data_from_prior_alt['mnu'], '-o', markersize=dot_size, color=fc.nice_colors(1) )
        # 3- r:
        ax_1_p1.plot( years, model_data_from_prior_alt['r'], '-o', markersize=dot_size, color=fc.nice_colors(2) )
        ax_1_p2.plot( years, model_data_from_prior_alt['r'], '-o', markersize=dot_size, color=fc.nice_colors(2) )

    # plot cumulative experiment prior information gain:
    if index==1:
        # 2- mnu:
        ax_1_p1.plot( years_cumulative, model_cumulative_info_gain_alt['mnu'], '-o', markersize=dot_size, color=fc.nice_colors(1) )
        ax_1_p2.plot( years_cumulative, model_cumulative_info_gain_alt['mnu'], '-o', markersize=dot_size, color=fc.nice_colors(1) )
        # 3- r:
        ax_1_p1.plot( years_cumulative, model_cumulative_info_gain_alt['r'], '-o', markersize=dot_size, color=fc.nice_colors(2) )
        ax_1_p2.plot( years_cumulative, model_cumulative_info_gain_alt['r'], '-o', markersize=dot_size, color=fc.nice_colors(2) )

    # plot differential experiment prior information gain:
    if index==2:
        # 2- mnu:
        ax_1_p1.plot( years_differential, model_differential_info_gain_alt['mnu'], '-o', markersize=dot_size, color=fc.nice_colors(1) )
        ax_1_p2.plot( years_differential, model_differential_info_gain_alt['mnu'], '-o', markersize=dot_size, color=fc.nice_colors(1) )
        # 3- r:
        ax_1_p1.plot( years_differential, model_differential_info_gain_alt['r'], '-o', markersize=dot_size, color=fc.nice_colors(2) )
        ax_1_p2.plot( years_differential, model_differential_info_gain_alt['r'], '-o', markersize=dot_size, color=fc.nice_colors(2) )
        # filling:
        if filling:
            ax_1_p1.fill_between( years, model_differential_info_gain_alt['mnu'], 1.e-16, color=fc.nice_colors(1), alpha = alpha )
            ax_1_p2.fill_between( years, model_differential_info_gain_alt['mnu'], 1.e-16, color=fc.nice_colors(1), alpha = alpha )
            ax_1_p1.fill_between( years, model_differential_info_gain_alt['r'], 1.e-16, color=fc.nice_colors(2), alpha = alpha )
            ax_1_p2.fill_between( years, model_differential_info_gain_alt['r'], 1.e-16, color=fc.nice_colors(2), alpha = alpha )

    ax_1_p1.set_ylabel('bits', fontsize=main_fontsize)

    # set log plot:
    ax_1_p1.set_yscale('log')
    ax_1_p2.set_yscale('log')

    # limits:
    ax_1_p1.set_xlim( 1920, 1940 )
    ax_1_p2.set_xlim( 1990, 2040 )

    # ticks with experiments:
    ax_1_p1_temp = ax_1_p1.twiny()
    ax_1_p1_temp.set_xlim( ax_1_p1.get_xlim() )
    ax_1_p1_temp.spines['right'].set_visible(False)
    ax_1_p1_temp.set_xticks( [1920, 1940 ] )

    ax_1_p2_temp = ax_1_p2.twiny()
    ax_1_p2_temp.set_xlim( ax_1_p2.get_xlim() )
    ax_1_p2_temp.spines['left'].set_visible(False)

    fig.canvas.draw()
    ax_1_p1_temp.set_xticklabels( [item.get_text() for item in ax_1_p1_temp.get_xticklabels()], fontsize=0.9*main_fontsize )
    ax_1_p2_temp.set_xticklabels( [item.get_text() for item in ax_1_p2_temp.get_xticklabels()], fontsize=0.9*main_fontsize )

    ax_1_p1_temp.xaxis.get_majorticklabels()[0].set_horizontalalignment('left')
    ax_1_p1_temp.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right')

    ax_1_p2_temp.xaxis.get_majorticklabels()[0].set_horizontalalignment('left')
    ax_1_p2_temp.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right')

    # hide what needs to be hidden:

    # the vertical axes:
    ax_1_p1.spines['right'].set_visible(False)
    ax_1_p2.spines['left'].set_visible(False)

    # move ticks to the proper position
    ax_1_p1.yaxis.tick_left()
    ax_1_p2.yaxis.tick_right()

    ax_1_p1.tick_params()
    ax_1_p2.tick_params(labelright='off')

    # set the ticks and ticks labels:

    ax_1_p1.set_xticks( [years[0]] )
    ax_1_p2.set_xticks( years[1:] )

    ax_1_p1.set_xticklabels( [label[0]], fontsize=0.9*main_fontsize, rotation=45, horizontalalignment='right' )
    ax_1_p2.set_xticklabels( label[1:], fontsize=0.9*main_fontsize, rotation=45, horizontalalignment='right' )

    # draw the canvas to adjust the ticks:
    fig.canvas.draw()
    ax_1_p1.set_yticklabels( [item.get_text() for item in ax_1_p1.get_yticklabels()], fontsize=0.9*main_fontsize )

    ax_1_p1.get_yticklabels()[0].set_verticalalignment('bottom')
    ax_1_p1.get_yticklabels()[-1].set_verticalalignment('top')

    # update structure and sizes:
    gs.update( bottom=0.3, top=0.85, left=0.1, right=0.98, wspace=0.05, hspace=0.10)

    if index==0:
        ax_1_p2.set_title('Prior Model Information Gain', fontsize=main_fontsize, loc='right',y=1.15 )
    elif index==1:
        ax_1_p2.set_title('Cumulative Prior Model Information Gain', fontsize=main_fontsize, loc='right',y=1.15 )
    elif index==2:
        ax_1_p2.set_title('Differential Model Information Gain', fontsize=main_fontsize, loc='right',y=1.15 )

    # legend:
    names  = [ '$\\Lambda$CDM +mnu', '$\\Lambda$CDM +r' ]
    colors = [ fc.nice_colors(1), fc.nice_colors(2) ]

    leg_handlers = []
    for name, color in zip( names, colors ):
        leg_handlers.append( mlines.Line2D([], [], color = color ) )

    fig.legend( handles  = leg_handlers,
                labels   = names,
                fontsize = main_fontsize,
                frameon  = True,
                fancybox = True,
                ncol     = 2,
                borderaxespad = 0.0,
                loc = 'lower center',
                )

    # export:
    if index==0:
        plt.savefig(here+'/../results/single/4_model_information_gain_LCDM.pdf', transparent=True)
    elif index==1:
        plt.savefig(here+'/../results/single/4_model_cumulative_information_gain_LCDM.pdf', transparent=True)
    elif index==2:
        plt.savefig(here+'/../results/single/4_model_Differential_information_gain_LCDM.pdf', transparent=True)

    plt.clf()

exit(0)
