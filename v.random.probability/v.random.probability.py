#!/usr/bin/env python
# -*- coding: utf-8 -*-

############################################################################
#
# MODULE:     v.random.probability
# AUTHOR(S):  Anna Petrasova
# PURPOSE:    Generates random vector points based on probability raster
# COPYRIGHT: (C) 2013 by the GRASS Development Team
#
#                This program is free software under the GNU General Public
#                License (>=v2). Read the file COPYING that comes with GRASS
#                for details.
#
#############################################################################

#%module
#% description: Generates random vector points based on probability raster
#% keywords: vector
#% keywords: random
#%end
#%option
#% type: string
#% gisprompt: old,cell,raster
#% key: probability
#% label: Name of the input probability raster map
#% description: Probability of generating points
#% required : yes
#%end
#%option
#% type: string
#% gisprompt: new,vector,vector
#% key: output
#% description: Name for output vector map
#% required : yes
#%end
#%option
#% key: count
#% type: integer
#% description: Approximate number of particles
#% required : yes
#%end

import os
import atexit
from math import sqrt

from grass.script import core as gcore
from grass.script import raster as grast


TMP_RAST = []
TMP_VECT = []


def main():
    options, flags = gcore.parser()
    probability = options['probability']
    output = options['output']
    count = int(options['count'])

    gcore.use_temp_region()

    # probability map
    probab_01 = 'probability_01_' + str(os.getpid())
    TMP_RAST.append(probab_01)
    info = grast.raster_info(probability)
    gcore.write_command('r.recode', flags='d', input=probability, output=probab_01,
                        title="Recoded probability map to 0 to 1",
                        rules='-', stdin='{minim}:{maxim}:0:1'.format(minim=info['min'], maxim=info['max']))
    mean = gcore.parse_key_val(gcore.read_command('r.univar', map=probab_01, flags='g'),
                               val_type=float)['mean']
    resolution = count / (mean * (info['north'] - info['south'] + info['east'] - info['west']))
    resolution = sqrt((mean * (info['north'] - info['south']) * (info['east'] - info['west'])) / count)
    gcore.run_command('g.region', res=resolution)

    random_name = 'random_' + str(os.getpid())
    point_map = 'points_' + str(os.getpid())
    point_grid = 'points_' + str(os.getpid())
    TMP_RAST.append(random_name)
    TMP_RAST.append(point_map)
    TMP_VECT.append(point_grid)

    gcore.run_command('r.surf.random', output=random_name, min=0, max=1)
    grast.mapcalc(exp='{point_map} = if({rand} <= {prob}, 1, null())'.format(rand=random_name,
                                                                             prob=probab_01,
                                                                             point_map=point_map))
    gcore.run_command('r.to.vect', flags='t', input=point_map, output=point_grid, type='point')
    gcore.run_command('v.perturb', input=point_grid, output=output,
                      parameter=resolution / 2., seed=os.getpid())


def cleanup():
    if len(TMP_RAST + TMP_VECT):
        gcore.info(_("Cleaning %d temporary maps...") % len(TMP_RAST + TMP_VECT))

    gcore.run_command('g.remove', rast=','.join(TMP_RAST), quiet=True)
    gcore.run_command('g.remove', vect=','.join(TMP_VECT), quiet=True)
    gcore.del_temp_region()


if __name__ == '__main__':
    atexit.register(cleanup)
    main()
