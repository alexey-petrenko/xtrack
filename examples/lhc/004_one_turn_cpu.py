from pathlib import Path
import numpy as np
from scipy.special import factorial

import xtrack as xt
import xobjects as xo
import sixtracktools
import pysixtrack

short_test = True # Short line (5 elements)

api_conf = {'prepointer': ' /*gpuglmem*/ '}

context = xo.ContextCpu()
context = xo.ContextCupy()
#context = xo.ContextPyopencl('0.0')

six = sixtracktools.SixInput(".")
pyst_line = pysixtrack.Line.from_sixinput(six)

if short_test:
    new_elements = []
    new_names = []
    found_types = []
    for ee, nn in zip(pyst_line.elements, pyst_line.element_names):
        if ee.__class__ not in found_types:
            new_elements.append(ee)
            new_names.append(nn)
            found_types.append(ee.__class__)
    pyst_line.elements = new_elements
    pyst_line.element_names = new_names

    pyst_line.elements[0] = pysixtrack.elements.Drift(length=77.)

sixdump = sixtracktools.SixDump101("res/dump3.dat")

# TODO: The two particles look identical, to be checked
part0_pyst = pysixtrack.Particles(**sixdump[0::2][0].get_minimal_beam())
part1_pyst = pysixtrack.Particles(**sixdump[1::2][0].get_minimal_beam())
pysixtrack_particles = [part0_pyst, part1_pyst]

particles = xt.Particles(pysixtrack_particles=[part0_pyst, part1_pyst],
                         _context=context)

print('Creating line...')
xtline = xt.Line(_context=context, sequence=pyst_line)

print('Build capi')
sources = []
kernels = {}
cdefs = []

sources.append(xt._pkg_root.joinpath('headers/constants.h'))

# Particles
source_particles, kernels_particles, cdefs_particles = (
                            xt.Particles.XoStruct._gen_c_api(conf=api_conf))
sources.append(source_particles)
kernels.update(kernels_particles)
cdefs += cdefs_particles.split('\n')

# Local particles
sources.append(xt.particles.gen_local_particle_api())

# Elements
element_classes = xtline._ElementRefClass._rtypes
for cc in element_classes:
    ss, kk, dd = cc._gen_c_api(conf=api_conf)
    sources.append(ss)
    kernels.update(kk)
    cdefs += dd.split('\n')

    sources.append(cc.track_function_source)

cdefs_norep=[]
for cc in cdefs:
    if cc not in cdefs_norep:
        cdefs_norep.append(cc)

src_lines = []
src_lines.append(r'''
    /*gpukern*/
    void track_line(
        /*gpuglmem*/ int8_t* buffer,
        /*gpuglmem*/ int64_t* ele_offsets,
        /*gpuglmem*/ int64_t* ele_typeids,
                     ParticlesData particles,
                     int ele_start,
                     int num_ele_track){


    LocalParticle lpart;

    int64_t part_id = 0;                    //only_for_context cpu_serial cpu_openmp
    int64_t part_id = blockDim.x * blockIdx.x + threadIdx.x; //only_for_context cuda
    int64_t part_id = get_global_id(0);                    //only_for_context opencl

    int64_t n_part = ParticlesData_get_num_particles(particles);
    if (part_id<n_part){
    Particles_to_LocalParticle(particles, &lpart, part_id);


    for (int64_t ee=ele_start; ee<ele_start+num_ele_track; ee++){
        /*gpuglmem*/ int8_t* el = buffer + ele_offsets[ee];
        int64_t ee_type = ele_typeids[ee];

        switch(ee_type){
''')

for ii, cc in enumerate(element_classes):
    ccnn = cc.__name__.replace('Data', '')
    src_lines.append(f'''
            case {ii}:
                {ccnn}_track_local_particle(({ccnn}Data) el, &lpart);
                break;''')

src_lines.append('''
        } //switch
    } //for
    }//if
}//kernel
''')

source_track = '\n'.join(src_lines)
sources.append(source_track)

kernel_descriptions = {
    "track_line": xo.Kernel(
        args=[
            xo.Arg(xo.Int8, pointer=True, name='buffer'),
            xo.Arg(xo.Int64, pointer=True, name='ele_offsets'),
            xo.Arg(xo.Int64, pointer=True, name='ele_typeids'),
            xo.Arg(xt.particles.ParticlesData, name='particles'),
            xo.Arg(xo.Int32, name='ele_start'),
            xo.Arg(xo.Int32, name='num_ele_track'),
        ],
    )
}

# Internal API can be exposed only on CPU
if not isinstance(context, xo.ContextCpu):
    kernels = {}
kernels.update(kernel_descriptions)

# Compile!
context.add_kernels(sources, kernels, extra_cdef='\n\n'.join(cdefs_norep),
                    save_source_as='source.c',
                    specialize=True)

print('Start check')
ele_offsets = np.array([ee._offset for ee in xtline.elements], dtype=np.int64)
ele_typeids = np.array(
        [element_classes.index(ee._xobject.__class__) for ee in xtline.elements],
        dtype=np.int64)
ele_offsets_dev = context.nparray_to_context_array(ele_offsets)
ele_typeids_dev = context.nparray_to_context_array(ele_typeids)

ip_check = 1
pyst_part = pysixtrack_particles[ip_check].copy()
vars_to_check = ['x', 'px', 'y', 'py', 'zeta', 'delta', 's']
problem_found = False
for ii, (eepyst, nn) in enumerate(zip(pyst_line.elements, pyst_line.element_names)):
    print(f'\nelement {nn}')
    vars_before = {vv :getattr(pyst_part, vv) for vv in vars_to_check}
    particles.set_one_particle_from_pysixtrack(ip_check, pyst_part)

    context.kernels.track_line.description.n_threads = particles.num_particles
    context.kernels.track_line(buffer=xtline._buffer.buffer,
                               ele_offsets=ele_offsets_dev,
                               ele_typeids=ele_typeids_dev,
                               particles=particles._xobject,
                               ele_start=ii,
                               num_ele_track=1)

    eepyst.track(pyst_part)
    for vv in vars_to_check:
        pyst_change = getattr(pyst_part, vv) - vars_before[vv]
        xt_change = getattr(particles, vv)[ip_check] -vars_before[vv]
        passed = np.isclose(xt_change, pyst_change, rtol=1e-10, atol=1e-14)
        if not passed:
            problem_found = True
            print(f'Not passend on var {vv}!\n'
                  f'    pyst:   {pyst_change: .7e}\n'
                  f'    xtrack: {xt_change: .7e}\n')
            break

    if not passed:
        break
    else:
        print("Check passed!")


if not problem_found:
    print('All passed on context:')
    print(context)
