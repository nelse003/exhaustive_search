from __future__ import division
from __future__ import print_function

import os
from cStringIO import StringIO

import cctbx.eltbx
import iotbx.pdb
import matplotlib.lines as mlines
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mmtbx.utils
import numpy as np
from iotbx import reflection_file_utils
from matplotlib.legend import Legend
from mmtbx.utils import data_and_flags_master_params

from exhaustive import master_phil

"""
Exploring Bhat 1989: R value for atoms

NOT TESTED after refactor
"""


def dump(obj):
    for attr in dir(obj):
        print("obj.%s = %r" % (attr, getattr(obj, attr)))


def calc_bhat_r(occ, atom_type):
    r_top_sum = np.zeros(700)
    r_bot_sum = np.zeros(700)

    for stol_sq in fmodel.f_obs_work().sin_theta_over_lambda_sq().data():
        f0 = cctbx.eltbx.xray_scattering.best_approximation(atom_type).at_stol_sq(stol_sq)

        b = np.linspace(0, 70, 700)

        r_top = np.abs(occ * f0 - f0 * np.exp(-b * stol_sq))
        r_bot = occ * f0
        r_top_sum += r_top
        r_bot_sum += r_bot

    r = r_top_sum / r_bot_sum

    return b, r

if __name__ == '__main__':

    params = master_phil.extract()
    params.input.xtal_name = "FALZA-x0085"
    params.input.in_path = os.path.join(os.path.realpath(
        "./exhaustive/test/resources"), params.input.xtal_name)
    params.validate.input.base_mtz = os.path.join(params.input.in_path,
                                                  "FALZA-x0085.free.mtz")
    params.input.mtz = os.path.join(params.input.in_path,
                                    "FALZA-x0085.free.mtz")
    params.input.pdb = os.path.join(params.input.in_path, "refine.pdb")

    args = [params.input.pdb, params.input.mtz]
    inputs = mmtbx.utils.process_command_line_args(args=args)

    rfs = reflection_file_utils.reflection_file_server(
        crystal_symmetry=inputs.crystal_symmetry,
        force_symmetry=True,
        reflection_files=inputs.reflection_files,
        err=StringIO())

    pdb_inp = iotbx.pdb.input(file_name=inputs.pdb_file_names[0])
    ph = pdb_inp.construct_hierarchy()

    xrs = ph.extract_xray_structure(
        crystal_symmetry=inputs.crystal_symmetry)
    xrs.show_summary()

    data_flags_params = data_and_flags_master_params().extract()
    data_flags_params.labels = params.exhaustive.options.column_type

    determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
        reflection_file_server=rfs,
        parameters=data_flags_params,
        keep_going=True,
        log=StringIO())

    f_obs = determined_data_and_flags.f_obs
    r_free_flags = determined_data_and_flags.r_free_flags

    sites_frac = xrs.sites_frac()

    f_obs.show_summary()

    mask_params = mmtbx.masks.mask_master_params.extract()
    mask_params.ignore_hydrogens = False
    mask_params.ignore_zero_occupancy_atoms = False

    fmodel = mmtbx.f_model.manager(
        f_obs=f_obs,
        r_free_flags=r_free_flags,
        mask_params=mask_params,
        xray_structure=xrs)

    d_spacings = fmodel.f_obs_work().d_spacings().data()

    s_sqrd = fmodel.f_obs_work().sin_theta_over_lambda_sq().data()

    print(len(d_spacings))
    print(len(s_sqrd))

    for i in xrange(0, 100):
        print(d_spacings[i], s_sqrd[i])

    print(fmodel.f_obs_work().d_max_min())

    print(xrs.scattering_type_registry())

    dump(xrs.scattering_type_registry())

    print("----------------------------------------------")

    xrs.scattering_type_registry().show_summary()

    print("----------------------------------------------")

    print(xrs.scattering_dictionary_as_string())

    print("----------------------------------------------")

    dump(xrs)

    print("----------------------------------------------")

    dump(xrs.structure_factors)

    print("----------------------------------------------")

    print(xrs.scattering_type_registry().last_table())

    print("----------------------------------------------")

    print(cctbx.eltbx.xray_scattering.best_approximation("O"))
    dump(cctbx.eltbx.xray_scattering.best_approximation("O"))

    print("----------------------------------------------")

    print(cctbx.eltbx.xray_scattering.best_approximation("O").show())

    jet = plt.get_cmap('jet')
    colors = iter(jet(np.linspace(0.1, 1, 10)))

    fig, ax = plt.subplots(figsize=(8, 12))
    ax_top = plt.axes([0.4, 0.4, 0.35, 0.35])
    ax_top.set_xlim(0, 20)
    ax_top.set_ylim(0, 1)
    ax_top.set_xlabel("B")
    ax_top.set_ylabel("R")
    ax.set_xlabel("B")
    ax.set_ylabel("R")
    ax.set_title("Bhat 1989: R value for atoms")

    for occ in np.linspace(0.1, 1, 10):
        colour = next(colors)

        B_ox, R_ox = calc_bhat_R(occ=occ, atom_type="O")
        ax.plot(B_ox, R_ox, c=colour, label='{}'.format(occ))
        ax_top.plot(B_ox, R_ox, c=colour)

        B_cl, R_cl = calc_bhat_R(occ=occ, atom_type="CL")
        ax.plot(B_cl, R_cl, c=colour, linestyle='dashed')
        ax_top.plot(B_cl, R_cl, c=colour, linestyle='dashed')

        B_h, R_h = calc_bhat_R(occ=occ, atom_type="H")
        ax.plot(B_h, R_h, c=colour, linestyle=':')
        ax_top.plot(B_h, R_h, c=colour, linestyle=':')

    dashed = mlines.Line2D([], [], color='k', linestyle='dashed', label="CL")
    solid = mlines.Line2D([], [], color='k', label="O")
    dotted = mlines.Line2D([], [], color='k', linestyle=':', label="H")

    ax.legend(loc='best', title="Occupancy", frameon=False)

    leg = Legend(parent=ax,
                 handles=[dashed, solid, dotted],
                 labels=["CL", "O", "H"],
                 loc='upper center',
                 frameon=False)

    ax.add_artist(leg)
    plt.savefig("bhat_o_cl_h.png", dpi=300)

    # print(fmodel.f_obs_work().structure_factors_from_scatterers(xray_structure=xrs))

    # print(fmodel.f_obs_work().structure_factors_from_scatterers(
    #     xray_structure=xrs).f_calc().size())
