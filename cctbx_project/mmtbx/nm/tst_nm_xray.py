from __future__ import division
from mmtbx.bulk_solvent import bulk_solvent_and_scaling
from mmtbx import masks
from iotbx import reflection_file_reader
from iotbx import pdb
from mmtbx.nm import tools
from mmtbx_nm_ext import *
import mmtbx.f_model
import mmtbx.model
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
from libtbx.utils import user_plus_sys_time
from pprint import pprint
from cStringIO import StringIO
import time
import libtbx.load_env
import iotbx.pdb
import random
import sys, os
import struct

def run(args):
    #    if (len(args) == 0 ):
#        raise RuntimeError("Please specify one or more pdb file names.")
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
#   pdb_file = libtbx.env.find_in_repositories(
#          relative_path="phenix_regression/pdb/1A32.pdb",
#          test=os.path.isfile)
    pdb_file = "./4KI8.pdb"
#   print pdb_file
    processed_pdb_file = monomer_library.pdb_interpretation.process(
           mon_lib_srv               = mon_lib_srv,
           ener_lib                  = ener_lib,
           file_name                 = pdb_file,
           raw_records               = None,
           force_symmetry            = True)
    xray_structure = processed_pdb_file.xray_structure()
    pdb_inp = iotbx.pdb.input(file_name=pdb_file)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
#    xray_structure.scattering_type_registry(table = "wk1995")
#    reflection_file = reflection_file_reader.any_reflection_file(file_name = "./4ksc-sf.cif")
#    miller_arrays = reflection_file.as_miller_arrays()
#    for miller_array in miller_arrays:
#        if (miller_array.is_xray_amplitude_array()):
#            f_obs = miller_array
    dummy = xray_structure.structure_factors(d_min = 3.0).f_calc()
    f_obs = abs(dummy.structure_factors_from_scatterers(
                xray_structure = xray_structure).f_calc())
    flags = f_obs.generate_r_free_flags(fraction=0.01)
#    f_calc = f_obs.structure_factors_from_scatterers(
#            xray_structure = xray_structure).f_calc()
    fmodel = mmtbx.f_model.manager(xray_structure = xray_structure,
                                   f_obs          = f_obs,
                                   r_free_flags   = flags)
#    bulk_solvent_manager = bulk_solvent_and_scaling.bulk_solvent_and_scales(
#            fmodel_kbu = fmodel.fmodel_kbu)
#    f_bulk = bulk_solvent_manager.f_bulk()
#    f_bulk.set_info("Bulk solvent structure factors")
#    f_bulk.show_summary()


#    xray_structure.convert_to_isotropic()
#    u_iso_start = xray_structure.extract_u_iso_or_u_equiv()
#    xray_structure.set_b_iso(value = 100.0)
#    xray_structure.convert_to_anisotropic()
    xray_structure.show_u_statistics(
            text = "After set u_iso_start")
    atoms = pdb_hierarchy.atoms()
    atomsxyz = atoms.extract_xyz()
    sites_cart = xray_structure.sites_cart()
    atomic_weights = xray_structure.atomic_weights()
    u_cart = xray_structure.scatterers().extract_u_cart(xray_structure.unit_cell())
#    nm_init_manager = nm_init(filename = "4ki8_evec.dat",
#            n_modes = 50,
#            atoms = atoms,
#            zero_mode_input_flag = False,
#            zero_mode_flag = True)
#    eigenvector = nm_init_manager.return_modes(1)
    selections = []
    selection_strings = ["chain A", "chain B", "chain C", "chain D", "chain E", "chain F", "chain G"]
#    selection_strings = ["chain A"]
    for string in selection_strings:
        selections.append(processed_pdb_file.all_chain_proxies.selection(
                                                                string = string))
    geometry = processed_pdb_file.geometry_restraints_manager(
                                            show_energies = False,
                                            plain_pairs_radius = 5.0)
    restraints_manager = mmtbx.restraints.manager(geometry = geometry)
    modes = tools.generate_evec(selections = selections,
                                xray_structure = xray_structure,
                                pdb_hierarchy = pdb_hierarchy,
                                filename = "4ki8_evec.dat",
                                n_modes = 20)
    uanisos = flex.sym_mat3_double(xray_structure.sites_cart().size(), [0,0,0,0,0,0])
    nmval = tools.read_nmval_file(evalin = "4ki8_eval.dat",
                                  n_modes = 20,
                                  zero_mode_input_flag = False,
                                  zero_mode_flag = True)
    xs_initial = []
    adp_nmas = []
    for selection in selections:
        x = tools.init_nm_para(nmval = nmval,
                               n_modes = 20)
        modes1d = tools.selected_modes_to_1D(modes = modes, n_modes = 20, selection = selection)
        weights_selected = xray_structure.atomic_weights().select(selection)
        print "The number of selected atoms is %d." % len(weights_selected)
        u_cart_selected = u_cart.select(selection)
        adp_nma = init_nm_adp(modes = modes1d,
                              weights = weights_selected,
                              n_modes = 20,
                              zero_mode_flag = True)
        adp_nmas.append(adp_nma)
        uaniso_from_s_manager = uaniso_from_s(x = x,                            
                                              adp_nma = adp_nma, 
                                              n_modes = 20,
                                              zero_mode_flag = True)
        u = uaniso_from_s_manager.u_cart()
        x_scaled = scale_x(x = x,
                           uanisos = u_cart_selected,
                           adp_all = u,
                           n_modes = 20,
                           zero_mode_flag = True)
        xs_initial.append(x_scaled)
#    t1 = time.time()
    sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
    sfg_params.cos_sin_table = False
    tools.update_xray_structure_with_nm(xray_structure = xray_structure,
                                        adp_nmas = adp_nmas,
                                        selections = selections,
                                        xs = xs_initial,
                                        n_modes = 20,
                                        zero_mode_flag = True)
    del adp_nmas
    class refinement_flags: pass
    refinement_flags.adp_tls = selections
    model = mmtbx.model.manager(
            refinement_flags = refinement_flags,
            restraints_manager = restraints_manager,
            xray_structure = xray_structure,
            pdb_hierarchy = pdb_hierarchy)
    fmodel.update_xray_structure(xray_structure = xray_structure,
                                 update_f_calc = True)
    fmodel_cp = fmodel.deep_copy()
    nm_refinement_manager = tools.nm_refinement(
            fmodel = fmodel_cp,
            model  = model,
            selections = selections,
            selections_1d = None,
            n_modes = 20,
            number_of_macro_cycles = 30,
            max_number_of_iterations = 10000,
            evecin = "4ki8_evec.dat",
            evalin = "4ki8_eval.dat",
            start_xs_value = xs_initial,
            run_finite_differences_test = False)
    tools.show_time()

if (__name__ == "__main__"):
    run(sys.argv[1:])
