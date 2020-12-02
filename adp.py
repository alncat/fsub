from __future__ import division
from mmtbx.refinement import print_statistics
from mmtbx.refinement import adp_refinement
import scitbx.lbfgs
import mmtbx.refinement.minimization

def run(
      model,
      params,
      fmodels,
      save_optimize_adp_weight,
      log,
      macro_cycle,
      prefix,
      target_weights):
  target_weights_params = params.target_weights
  h_params = params.hydrogens
  print_statistics.make_header(prefix, out = log)
  if(save_optimize_adp_weight):
    r_work = fmodels.fmodel_xray().r_work()
    r_free = fmodels.fmodel_xray().r_free()
    if ((r_free < r_work or (r_free-r_work)<0.01) and
        (not target_weights_params.force_optimize_weights)) :
      target_weights_params.optimize_adp_weight = False
    else:
      target_weights_params.optimize_adp_weight=save_optimize_adp_weight
  fmodels.fmodel_xray().xray_structure.scatterers().flags_set_grads(
    state = False)
  # Exclude H from ADP refinement
  if ((params.hydrogens.refine == "riding") or
      (params.hydrogens.force_riding_adp)) :
    excl_h = model.xray_structure.hd_selection().count(True) > 0
    if(excl_h and model.refinement_flags.adp_individual_iso is not None):
      model.refinement_flags.adp_individual_iso.set_selected(
        model.xray_structure.hd_selection(), False)
    if(excl_h and model.refinement_flags.adp_individual_aniso is not None):
      model.refinement_flags.adp_individual_aniso.set_selected(
        model.xray_structure.hd_selection(), False)
  #
  save_xray_structure = \
    fmodels.fmodel_xray().xray_structure.deep_copy_scatterers() # XXX ?
  ###> Make ADP of H/D sites equal
  model.reset_adp_of_hd_sites_to_be_equal()
  fmodels.update_xray_structure(xray_structure = model.xray_structure,
                                update_f_calc  = True)
  tidy_us(params=params, fmodels=fmodels, model=model, log=log)
  adp_refinement_manager = adp_refinement.manager(
    fmodels                = fmodels,
    model                  = model,
    group_adp_selections   = model.refinement_flags.adp_group,
    group_adp_selections_h = model.refinement_flags.group_h,
    group_adp_params       = params.group_b_iso,
    tls_selections         = model.refinement_flags.adp_tls,
    nm_selections          = model.refinement_flags.adp_nm,
    all_params             = params,
    tls_params             = params.tls,
    nm_params              = params.nm,
    individual_adp_params  = params.adp,
    adp_restraints_params  = params.adp_restraints,
    refine_adp_individual  = model.refinement_flags.individual_adp,
    refine_adp_group       = model.refinement_flags.group_adp,
    refine_tls             = model.refinement_flags.tls,
    refine_nm              = model.refinement_flags.nm,
    tan_b_iso_max          = 0,
    restraints_manager     = model.restraints_manager,
    macro_cycle            = macro_cycle,
    target_weights         = target_weights,
    log                    = log,
    h_params               = h_params)
  assert fmodels.fmodel_xray().xray_structure.min_u_cart_eigenvalue()>=0
  print >> log
  model.show_adp_statistics(padded = True)
  print >> log
  fmodels.update_xray_structure(update_f_calc = True)
  ###> Make ADP of H/D sites equal
  model.reset_adp_of_hd_sites_to_be_equal()
  fmodels.update_xray_structure(xray_structure = model.xray_structure,
                                update_f_calc  = True)
  assert fmodels.fmodel_xray().xray_structure.min_u_cart_eigenvalue()>=0

def refine_adp_for_ias_only(fmodels, model, params, prefix, log):
  print_statistics.make_header(prefix, out = log)
  xrs_start = model.xray_structure.deep_copy_scatterers()
  fmodels.show_targets(log = log, text = "start")
  fmodels.fmodel_xray().xray_structure.scatterers().flags_set_grads(
    state = False)
  i_selection = model.ias_selection.iselection()
  scatterers = fmodels.fmodel_xray().xray_structure.scatterers()
  for i_sel in i_selection:
    if(scatterers[i_sel].flags.use_u_iso()):
      scatterers[i_sel].flags.set_grad_u_iso(True)
    if(scatterers[i_sel].flags.use_u_aniso()):
      scatterers[i_sel].flags.set_grad_u_aniso(True)
  #
  lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
    max_iterations = params.main.max_number_of_iterations)
  minimized = mmtbx.refinement.minimization.lbfgs(
    restraints_manager       = None,
    fmodels                  = fmodels,
    model                    = model,
    refine_adp               = True,
    lbfgs_termination_params = lbfgs_termination_params)
  fmodels.show_targets(log = log, text = "final")

def refine_adp_for_h_only(fmodels, model, params, prefix, log):
  print_statistics.make_header(prefix, out = log)
  xrs_start = model.xray_structure.deep_copy_scatterers()
  fmodels.show_targets(log = log, text = "start")
  fmodels.fmodel_xray().xray_structure.scatterers().flags_set_grads(
    state = False)
  i_selection = model.xray_structure.hd_selection().iselection()
  fmodels.fmodel_xray().xray_structure.scatterers(
    ).flags_set_grad_u_iso(iselection = i_selection)
  lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
    max_iterations = params.main.max_number_of_iterations)
  is_neutron_scat_table = False
  if(params.main.scattering_table == "neutron"):
    is_neutron_scat_table = True
  minimized = mmtbx.refinement.minimization.lbfgs(
    restraints_manager       = None,
    is_neutron_scat_table    = is_neutron_scat_table,
    fmodels                  = fmodels,
    model                    = model,
    refine_adp               = True,
    lbfgs_termination_params = lbfgs_termination_params)
  fmodels.show_targets(log = log, text = "final")
  # sanity check


def tidy_us(params, fmodels, model, log):
  if(("individual_adp" in params.refine.strategy or
     "tls" in params.refine.strategy or "nm" in params.refine.strategy) and
     fmodels.fmodel_neutron() is None): # XXX does not work in Joint XN refinement
    r_work1 = fmodels.fmodel_xray().r_work()
    fmodels.fmodel_xray().xray_structure.tidy_us()
    fmodels.update_xray_structure(update_f_calc=True)
    r_work2 = fmodels.fmodel_xray().r_work()
    delta = (r_work2-r_work1)*100
    if(delta > 1.):
      print >> log, \
        "Tidying ADP resulted in r_work increase by %6.4f percent."%delta
      print >> log, "  r_work = %6.4f r_free = %6.4f"%(r_work2,
        fmodels.fmodel_xray().r_free())
    model.xray_structure=fmodels.fmodel_xray().xray_structure
