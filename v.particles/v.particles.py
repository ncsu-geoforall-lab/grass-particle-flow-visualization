#!/usr/bin/env python
# -*- coding: utf-8 -*-

############################################################################
#
# MODULE:     v.particles
# AUTHOR(S):  Anna Petrasova
# PURPOSE:    Creates vector maps for particle flow visualization
# COPYRIGHT: (C) 2013 by the GRASS Development Team
#
#                This program is free software under the GNU General Public
#                License (>=v2). Read the file COPYING that comes with GRASS
#                for details.
#
#############################################################################

#%module
#% description: Creates vector maps for particle flow visualization
#% keywords: visualization
#% keywords: vector
#%end
#%option
#% type: string
#% gisprompt: old,cell,raster
#% key: aspect
#% description: Name of the input aspect raster map
#% required : yes
#%end
#%option
#% type: string
#% gisprompt: old,cell,raster
#% key: speed
#% description: Name of the input speed raster map
#% required : yes
#%end
#%option
#% type: string
#% gisprompt: new,vector,vector
#% key: particle_base
#% label: Basename of the particle output maps
#% description: Separator _ will be used
#% required : yes
#%end
#%option
#% type: string
#% gisprompt: string
#% key: particles
#% description: Name of the output vector map
#% required: yes
#% guisection: comet-like visualization
#%end
#%option
#% key: min_size
#% type: double
#% description: Minimum size of particle
#% required: no
#% answer: 1
#% guisection: comet-like visualization
#%end
#%option
#% key: max_size
#% type: double
#% description: Maximum size of particle
#% required: no
#% answer: 8
#% guisection: comet-like visualization
#%end
#%option
#% key: comet_length
#% type: integer
#% description: Number of steps in the comet
#% required: no
#% answer: 10
#% guisection: comet-like visualization
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
#% key: total_time
#% type: integer
#% description: Time of the simulation in seconds
#% required : yes
#%end
#%option
#% key: step
#% type: integer
#% description: Time step of the simulation to output vector
#% required : yes
#%end
#%option
#% key: age
#% type: integer
#% description: Maximum particle age in seconds
#% required : yes
#%end
#%option
#% key: count
#% type: integer
#% description: Number of particles
#% required : yes
#%end

import os
import atexit
import tempfile

from grass.script import core as gcore
from grass.script import raster as grast
from grass.script import vector as gvect
import grass.temporal as tgis


TMP_VECT = []
TMP_RAST = []


def main():
    options, flags = gcore.parser()
    aspect = options['aspect']
    speed = options['speed']
    probability = options['probability']
    if options['particle_base']:
        particle_base = options['particle_base'] + '_'
    else:
        particle_base = None
    if options['particles']:
        particles = options['particles']
        min_size = float(options['min_size'])
        max_size = float(options['max_size'])
        comet_length = int(options['comet_length'])
    else:
        particles = min_size = max_size = comet_length = None
    try:
        total_time = int(options['total_time'])
        step = int(options['step'])
        age = int(options['age'])
        count = int(options['count'])
    except ValueError:
        gcore.fatal(_("Parameter should be integer"))

    gcore.use_temp_region()

    # create aspect in x and y direction
    aspect_x = 'aspect_x_' + str(os.getpid())
    aspect_y = 'aspect_y_' + str(os.getpid())
    xshift_tmp = 'xshift_tmp_' + str(os.getpid())
    yshift_tmp = 'yshift_tmp_' + str(os.getpid())
    TMP_RAST.append(aspect_x)
    TMP_RAST.append(aspect_y)
    grast.mapcalc(exp="{aspect_x} = cos({aspect})".format(aspect_x=aspect_x, aspect=aspect))
    grast.mapcalc(exp="{aspect_y} = sin({aspect})".format(aspect_y=aspect_y, aspect=aspect))
    grast.mapcalc(exp="{xshift} = {aspect_x}*{speed}*{t}".format(xshift=xshift_tmp, t=step, speed=speed,
                                                                 aspect_x=aspect_x), overwrite=True)
    grast.mapcalc(exp="{yshift} = {aspect_y}*{speed}*{t}".format(yshift=yshift_tmp, t=step, speed=speed,
                                                                 aspect_y=aspect_y), overwrite=True)

    # initialize
    vector_tmp1 = 'vector_tmp1_' + str(os.getpid())
    vector_tmp2 = 'vector_tmp2_' + str(os.getpid())
    vector_tmp3 = 'vector_tmp3_' + str(os.getpid())
    vector_region = 'vector_region_' + str(os.getpid())
    TMP_VECT.extend([vector_tmp1, vector_tmp2, vector_tmp3, vector_region])
    random_tmp = 'random_tmp_' + str(os.getpid())
    TMP_RAST.extend([xshift_tmp, yshift_tmp, random_tmp])
    gcore.run_command('v.in.region', output=vector_region, type='area')

    loop = 0
    vector_1 = particle_base + "{0:03d}".format(loop)
    generate_points(name=vector_1, probability_map=probability, count=count)

    grast.mapcalc(exp="{random} = int(rand(1, {maxt}))".format(random=random_tmp, maxt=age + 1))
    gcore.run_command('v.what.rast', map=vector_1, raster=random_tmp, column='t')
    write_vect_history('v.particles', options, flags, vector_1)
    vector_names = [vector_1, ]
    for time in range(0, total_time + step, step):
        vector_1 = particle_base + "{0:03d}".format(loop)
        vector_2 = particle_base + "{0:03d}".format(loop + 1)
        vector_names.append(vector_2)

        gcore.run_command('v.what.rast', map=vector_1, raster=xshift_tmp, column='xshift')
        gcore.run_command('v.what.rast', map=vector_1, raster=yshift_tmp, column='yshift')
        gcore.run_command('v.transform', layer=1, input=vector_1, output=vector_2,
                          columns='xshift:xshift,yshift:yshift', quiet=True)
        # increase age
        gcore.info("Increasing age...")
        sql = 'UPDATE {table} SET t=t+1;'.format(table=vector_2)
        gcore.run_command('db.execute', sql=sql)

        # remove old points
        gcore.info("Removing old points...")
        gcore.run_command('v.select', overwrite=True, ainput=vector_2, atype='point',
                          binput=vector_region, btype='area', operator='within', output=vector_tmp1)
        gcore.run_command('v.extract', input=vector_tmp1, layer=1, type='point',
                          where="t <= " + str(age) + " AND xshift IS NOT NULL", output=vector_tmp2, overwrite=True)

        # generate new points
        gcore.info("Generating new points...")
        count_to_generate = count - gvect.vector_info(vector_tmp2)['points']
        if count_to_generate > 0:
            generate_points(name=vector_tmp3, probability_map=probability,
                            count=count_to_generate, overwrite=True)

            gcore.info("Patchig new and old points...")
            gcore.run_command('v.patch', flags='e', input=[vector_tmp2, vector_tmp3],
                              output=vector_2, overwrite=True)
            sql = 'UPDATE {table} SET t={t} WHERE t IS NULL;'.format(table=vector_2, t=0)
            gcore.run_command('db.execute', sql=sql)
            
        write_vect_history('v.particles', options, flags, vector_2)

        loop += 1
    # Make sure the temporal database exists
    tgis.init()

    tgis.open_new_space_time_dataset(particle_base[:-1], type='stvds',
                                     temporaltype='relative',
                                     title="title", descr='desc',
                                     semantic='mean', dbif=None,
                                     overwrite=gcore.overwrite())
    # TODO: we must start from 1 because there is a bug in register_maps_in_space_time_dataset
    tgis.register_maps_in_space_time_dataset(
        type='vect', name=particle_base[:-1], maps=','.join(vector_names),
        start=str(1), end=None, unit='seconds', increment=step,
        interval=False, dbif=None)
        
    # create one vector map with multiple layers
    fd, path = tempfile.mkstemp(text=True)
    tmpfile = open(path, 'w')
    k = 0
    for vector in vector_names:
        k += 1
        layers = [x for x in range(k - comet_length + 1, k + 1) if x > 0]
        categories = list(range(len(layers), 0, -1))
        text = ''
        for layer, cat in zip(layers, categories):
            text += '{l} {c}\n'.format(l=layer, c=cat)
        coords = gcore.read_command('v.to.db', flags='p', quiet=True, map=vector,
                                    type='point', option='coor', separator=" ").strip()
        for coord in coords.split('\n'):
            coord = coord.split()
            tmpfile.write('P 1 {n_cat}\n{x} {y}\n'.format(n_cat=len(categories), x=coord[1], y=coord[2]))
            tmpfile.write(text)
    tmpfile.close()

    gcore.run_command('v.in.ascii', flags='n', overwrite=True, input=path, output=particles,
                      format='standard', separator=" ")
    os.close(fd)
    os.remove(path)
    k = 0
    sql = []
    sizes = get_sizes(max_size, min_size, comet_length)
    temporal_maps = []
    for vector in vector_names:
        k += 1
        table = 't' + str(k)
        gcore.run_command('v.db.addtable', map=particles, table=table, layer=k,
                          column="width double precision")
        temporal_maps.append(particles + ':' + str(k))
        for i in range(comet_length):
            sql.append("UPDATE {table} SET width={w:.1f} WHERE cat={c}".format(table=table,
                                                                               w=sizes[i][1], c=sizes[i][0]))
    gcore.write_command('db.execute', input='-', stdin=';\n'.join(sql))

    tgis.open_new_space_time_dataset(particles, type='stvds',
                                     temporaltype='relative',
                                     title="title", descr='desc',
                                     semantic='mean', dbif=None,
                                     overwrite=True)
    # TODO: we must start from 1 because there is a bug in register_maps_in_space_time_dataset
    tgis.register_maps_in_space_time_dataset(
        type='vect', name=particles, maps=','.join(temporal_maps),
        start=str(1), end=None, unit='seconds', increment=step,
        interval=False, dbif=None)

    write_vect_history('v.particles', options, flags, particles)


def generate_points(name, probability_map, count, overwrite=False):
    gcore.run_command('v.random.probability', output=name, probability=probability_map,
                      count=count, overwrite=overwrite)
    gcore.run_command('v.db.addtable', map=name,
                      columns="xshift double precision,"
                              "yshift double precision,"
                              "t integer")


def get_sizes(max_size, min_size, count):
    step = (max_size - min_size) / float(count - 1)
    cats = range(1, count + 1)
    sizes = []
    for i in range(count):
        sizes.append(min_size + i * step)
    return zip(cats, sizes)


def write_vect_history(module, options, flags, map_name):
    def check_string(s):
        try:
            float(s)
            return s
        except ValueError:
            return '"{}"'.format(s)

    command = module
    for key, value in options.iteritems():
        command += ' {key}={value}'.format(key=key, value=check_string(value))
    for key, value in flags.iteritems():
        if value:
            command += ' -{flag}'.format(flag=key)
    gcore.run_command('v.support', map=map_name, cmdhist=command)


def cleanup():
    if len(TMP_VECT + TMP_RAST):
        gcore.info(_("Cleaning %d temporary vector maps...") % len(TMP_VECT + TMP_RAST))

    gcore.run_command('g.remove', vect=','.join(TMP_VECT), quiet=True)
    gcore.run_command('g.remove', rast=','.join(TMP_RAST), quiet=True)


if __name__ == '__main__':
    atexit.register(cleanup)
    main()
