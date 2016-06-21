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

color_1 = fc.nice_colors(0)
color_2 = fc.nice_colors(3)
color_3 = fc.nice_colors(2)
color_4 = fc.nice_colors(1)
color_5 = fc.nice_colors(4)

y_size  = 9.0

# prepare the plot
fig = matplotlib.pyplot.gcf()
fig.set_size_inches( x_size/2.54, y_size/2.54 )

gs = gridspec.GridSpec(1, 5)

ax_1_p1 = plt.subplot(gs[0, 0:])

# vertical lines corresponding to data:
for plot in [ax_1_p1]:
    for year in years:
        dashes = [2,2,2,2]  # 10 points on, 5 off, 100 on, 5 off
        line = plot.axvline( year, color='k', linestyle='--', linewidth=0.5)
        line.set_dashes(dashes)
    dashes_2 = [2,2,2,2]  # 10 points on, 5 off, 100 on, 5 off
    line_2 = plot.axvline( 2016, color='red', linestyle='--', linewidth=1.0)
    line_2.set_dashes(dashes_2)

# 1- std:
ax_1_p1.plot( years[1:], model_cumulative_info_gain_alt['mnu'][1:], '-o', markersize=dot_size, color=color_1 )
ax_1_p1.plot( years[1:], model_cumulative_info_gain_alt['r'][1:], '-o', markersize=dot_size, color=color_2 )
ax_1_p1.plot( years[1:], model_cumulative_info_gain_alt['w0wa'][1:], '-o', markersize=dot_size, color=color_3 )
ax_1_p1.plot( years[1:], model_cumulative_info_gain_alt['Horava'][1:], '-o', markersize=dot_size, color=color_4 )
ax_1_p1.plot( years[1:], model_cumulative_info_gain_alt['Linder'][1:], '-o', markersize=dot_size, color=color_5 )

# filling:
if filling:
    ax_1_p1.fill_between( years[1:], model_cumulative_info_gain_alt['Horava'][1:], model_cumulative_info_gain_alt['w0wa'][1:], color=color_4, alpha = alpha )
    ax_1_p1.fill_between( years[1:], model_cumulative_info_gain_alt['r'][1:], model_cumulative_info_gain_alt['w0wa'][1:], color=color_2, alpha = alpha )
    ax_1_p1.fill_between( years[1:], model_cumulative_info_gain_alt['mnu'][1:], model_cumulative_info_gain_alt['w0wa'][1:], color=color_1, alpha = alpha )
    ax_1_p1.fill_between( years[1:], model_cumulative_info_gain_alt['Linder'][1:], model_cumulative_info_gain_alt['w0wa'][1:], color=color_5, alpha = alpha )

ax_1_p1.set_ylabel('bits', fontsize=main_fontsize)

# set log plot:
ax_1_p1.set_yscale('log')

# limits:
ax_1_p1.set_xlim( 1990, 2040 )

# ticks with experiments:
ax_1_p1_temp = ax_1_p1.twiny()
ax_1_p1_temp.set_xlim( ax_1_p1.get_xlim() )
ax_1_p1_temp.spines['right'].set_visible(False)

fig.canvas.draw()
ax_1_p1_temp.set_xticklabels( [item.get_text() for item in ax_1_p1_temp.get_xticklabels()], fontsize=0.9*main_fontsize )
ax_1_p1_temp.xaxis.get_majorticklabels()[0].set_horizontalalignment('left')
ax_1_p1_temp.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right')

# set the ticks and ticks labels:
ax_1_p1.set_xticks( years[1:] )
ax_1_p1.set_xticklabels( label[1:], fontsize=0.8*main_fontsize, rotation=45, horizontalalignment='right' )

# draw the canvas to adjust the ticks:
fig.canvas.draw()
ax_1_p1.set_yticklabels( [item.get_text() for item in ax_1_p1.get_yticklabels()], fontsize=0.9*main_fontsize )
ax_1_p1.get_yticklabels()[0].set_verticalalignment('bottom')
ax_1_p1.get_yticklabels()[-1].set_verticalalignment('top')

# update structure and sizes:
gs.update( bottom=0.28, top=0.88, left=0.08, right=0.99, wspace=0.05, hspace=0.10)

ax_1_p1.set_title('Model Information Gain:', fontsize=main_fontsize, loc='left',y=1.10 )

# legend:
names  = [ 'Massive Neutrinos', 'Tensor', 'CPL', 'Ho\\v rava', 'Growth $\\gamma$' ]
colors = [ color_1, color_2, color_3, color_4, color_5 ]

leg_handlers = []
for name, color in zip( names, colors ):
    leg_handlers.append( mlines.Line2D([], [], color = color ) )

ax_1_p1.legend( handles  = leg_handlers,
            labels   = names,
            fontsize = main_fontsize,
            frameon  = False,
            fancybox = False,
            ncol     = len(names),
            borderaxespad = 0.0,
            columnspacing = 0.4,
            handlelength = 1.5,
            loc = 'lower right',
            bbox_to_anchor=(1.0, 1.09)
            )

# export:
plt.savefig(here+'/../results/figure_2.pdf', transparent=True)
plt.clf()

exit(0)
