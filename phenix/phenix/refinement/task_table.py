from __future__ import division
import inspect

class run(object):
  """
  Create dictionary of True/False switches for running known refinement strategies
  """
  def __init__(self,
               params,
               macro_cycle,
               model,
               fmodel,
               neutron_refinement,
               neutron_refinement2):
    self.params              = params
    self.macro_cycle         = macro_cycle
    self.fmodel              = fmodel
    self.neutron_refinement  = neutron_refinement  #XXX REMOVE later
    self.neutron_refinement2 = neutron_refinement2 #XXX REMOVE later
    self.model               = model
    self.hd_present = self.model.xray_structure.hd_selection().count(True)>0
    self.solvent_present = self.model.solvent_selection().count(True)>0
    self.table = {}
    self.counts = {}
    self._initialize_counts()
    self.__call__()

  def _initialize_counts(self):
    for method_name in dir(self):
      if(method_name.startswith("run_")):
        self.counts.setdefault(method_name[4:], 0)

  def __call__(self):
    for method_name in dir(self):
      if(method_name.startswith("run_")):
        self.__getattribute__(method_name)()

  def _get_name(self): return inspect.stack()[1][3][4:]

  def status(self, func):
    self.__call__()
    if(not type(func)==type("")): func = func.__name__
    result = self.table[func]
    if(result): self.counts[func] += 1
    return result

  def update(
        self,
        params,
        macro_cycle,
        model,
        fmodel,
        neutron_refinement,
        neutron_refinement2):
    self.params              = params
    self.macro_cycle         = macro_cycle
    self.fmodel              = fmodel
    self.neutron_refinement  = neutron_refinement  #XXX REMOVE later
    self.neutron_refinement2 = neutron_refinement2 #XXX REMOVE later
    self.model               = model
    self.hd_present = self.model.xray_structure.hd_selection().count(True)>0
    self.solvent_present = self.model.solvent_selection().count(True)>0
    self.__call__()

  def _need_to_run(self, mode, nmc, mc):
    flag = False
    if(mode == "every_macro_cycle"): flag = True
    elif(str(mode).lower()=="none"): flag = False
    elif(mode == "every_macro_cycle_after_first"):
      if(mc>1 or nmc==1): flag = True
    elif(mode == "second_and_before_last"):
      if(mc==2 or mc==nmc-1 or nmc-1==0): flag=True
    elif(mode == "once"):
      if(int((nmc+1)/2) == mc): flag = True
    elif(mode == "first"):
      if(mc == 1): flag = True
    elif(mode == "first_half"):
      if(nmc/2 >= mc): flag = True
    elif(mode == "second_half"):
      if(nmc/2 <= mc): flag = True
    elif(mode == "filter_only"): flag = True # XXX REMOVE PVA: when replacing ordered_solvent.py->water.py
    else: raise RuntimeError("Not a valid mode: %s"%mode)
    return flag

  def run_bss(self):
    self.table[self._get_name()] = self.params.main.bulk_solvent_and_scale

  def run_sotr(self):
    fl = self.macro_cycle == 1 and \
         ("individual_sites" in self.params.refine.strategy or
          "individual_sites_real_space" in self.params.refine.strategy) and \
         (self.params.main.ncs or self.params.main.reference_model_restraints)
    self.table[self._get_name()] = fl

  def run_ncsrot(self):
    fl = self.macro_cycle > 0 and \
         ("individual_sites" in self.params.refine.strategy or
          "individual_sites_real_space" in self.params.refine.strategy) and \
         self.params.main.ncs
    self.table[self._get_name()] = fl

  def run_sol(self):
    d_min = self.fmodel.f_obs().d_min()
    fl = self.params.main.ordered_solvent and \
         self.params.ordered_solvent.low_resolution > d_min and \
         self._need_to_run(
           mode = self.params.ordered_solvent.mode,
           nmc  = self.params.main.number_of_macro_cycles,
           mc   = self.macro_cycle)
    self.table[self._get_name()] = fl
    return fl

  def run_ion(self):
    fl = self.params.main.place_ions and \
         self.run_sol()
    self.table[self._get_name()] = fl

  def run_fwh(self):
    self.table[self._get_name()] = self.params.main.find_and_add_hydrogens

  def run_ias(self):
    fl = self.params.main.ias and self.macro_cycle in [1,2,5]
    self.table[self._get_name()] = fl

  def run_rbr(self):
    fl = "rigid_body" in self.params.refine.strategy and \
         ((self.params.rigid_body.mode=="first_macro_cycle_only" and
          self.macro_cycle==1) or
          self.params.rigid_body.mode=="every_macro_cycle")
    self.table[self._get_name()] = fl

  def run_addcbetar(self):
    fl = (self.run_realsrl() or self.fmodel.f_obs().d_min() > 1.5) and \
         self.counts["addcbetar"] == 0 and \
         self.model.restraints_manager is not None
    self.table[self._get_name()] = fl

  def run_realsrl(self):
    mc = self.macro_cycle
    nmc = self.params.main.number_of_macro_cycles
    fl = (not self.neutron_refinement and not self.neutron_refinement2 and \
         "individual_sites_real_space" in self.params.refine.strategy and \
         ((mc > 1 and nmc>1) or (mc==1 and nmc==1)) and not \
         self.hd_present) \
         or \
         (len(self.params.refine.strategy)==1 and
          "individual_sites_real_space" in self.params.refine.strategy)
    self.table[self._get_name()] = fl
    return fl

  def run_realsrg(self):
    fl = (not self.neutron_refinement and not \
         self.neutron_refinement2 and \
         self.fmodel.f_obs().d_min() > 1.5 and \
         self.fmodel.r_work() > 0.25 and \
         self.model.refinement_flags.sites_individual is not None and \
         self.model.refinement_flags.sites_individual.count(True)>0 and \
         "individual_sites_real_space" in self.params.refine.strategy) \
         or \
         (len(self.params.refine.strategy)==1 and
          "individual_sites_real_space" in self.params.refine.strategy)
    self.table[self._get_name()] = fl

  def run_weight(self):
    fl = self.params.main.use_geometry_restraints and \
         [self.model.refinement_flags.individual_sites,
          self.model.refinement_flags.torsion_angles,
          self.model.refinement_flags.individual_adp,
          self.params.main.simulated_annealing].count(True) > 0 or \
         self.params.main.target_weights_only
    self.table[self._get_name()] = fl

  def run_den(self):
    self.table[self._get_name()] = "den" in self.params.refine.strategy

  def run_tardy(self):
    fl = self.params.main.simulated_annealing_torsion and \
         self._need_to_run(
           mode = self.params.tardy.mode,
           nmc  = self.params.main.number_of_macro_cycles,
           mc   = self.macro_cycle)
    self.table[self._get_name()] = fl

  def run_sacart(self):
    fl = self.params.main.simulated_annealing and \
         "individual_sites" in self.params.refine.strategy and \
          self._need_to_run(
            mode = self.params.simulated_annealing.mode,
            nmc  = self.params.main.number_of_macro_cycles,
            mc   = self.macro_cycle)
    self.table[self._get_name()] = fl

  def run_nqh(self):
    fl = self.macro_cycle > 0 and \
         "individual_sites" in self.params.refine.strategy and \
         self.params.main.nqh_flips and self.params.main.use_molprobity
    self.table[self._get_name()] = fl

  def run_fitrh(self):
    fl = "individual_sites" in self.params.refine.strategy and \
         self.hd_present and not \
         self.params.main.ias and \
         self.params.hydrogens.real_space_optimize_x_h_orientation and \
         ((self.neutron_refinement or self.neutron_refinement2) or
          (self.fmodel.f_obs().d_min()<1.2))
    self.table[self._get_name()] = fl

  def run_morph(self):
    pass

  def run_xyzrec(self):
    fl = (self.model.refinement_flags.individual_sites and
         self.params.main.max_number_of_iterations > 0 and
         self.model.refinement_flags.sites_individual.count(True) > 0 and
         'den' not in self.params.refine.strategy) or \
         (self.model.refinement_flags.individual_sites and
          self.model.refinement_flags.sites_individual.count(True) > 0 and
          'den' in self.params.refine.strategy and
          self.params.den.final_refinement_cycle)
    self.table[self._get_name()] = fl

  def run_regHxyz(self):
    fl = self.hd_present and \
         self.params.hydrogens.refine == "riding" and \
         "individual_sites" in self.params.refine.strategy and \
         self.model is not None
    self.table[self._get_name()] = fl

  def run_adp(self):
    fl = [self.model.refinement_flags.group_adp,
          self.model.refinement_flags.individual_adp,
          self.model.refinement_flags.tls,
          self.model.refinement_flags.nm].count(True) > 0 and \
         'den' not in self.params.refine.strategy
    self.table[self._get_name()] = fl

  def run_regHadp(self):
    fl = self.model is not None and \
         (self.params.hydrogens.refine == "riding" or
          self.params.hydrogens.force_riding_adp) and \
         "individual_adp" in self.params.refine.strategy and \
         self.hd_present
    self.table[self._get_name()] = fl

  def run_adph(self):
    fl = self.hd_present and not \
         self.neutron_refinement and not \
         self.neutron_refinement2 and \
         self.params.hydrogens.contribute_to_f_calc and \
         len(self.params.refine.strategy) > 0 and \
         self.model.refinement_flags.individual_adp and \
         self.params.main.scattering_table != "neutron" and \
         self.params.hydrogens.refine == "individual" and \
         not self.params.hydrogens.force_riding_adp
    self.table[self._get_name()] = fl

  def run_occh(self):
    fl = self.hd_present and not \
         self.neutron_refinement and not \
         self.neutron_refinement2 and \
         self.params.hydrogens.contribute_to_f_calc and \
         len(self.params.refine.strategy) > 0 and \
         self.model.refinement_flags.occupancies and \
         self.params.main.scattering_table != "neutron" and \
         self.params.hydrogens.refine == "individual"
    self.table[self._get_name()] = fl

  def run_adpias(self):
    fl = self.model.ias_selection is not None and \
         self.model.ias_selection.count(True) > 0 and \
         len(self.params.refine.strategy) > 0 and \
         self.model.refinement_flags.individual_adp
    self.table[self._get_name()] = fl

  def run_occias(self):
    fl = self.model.ias_selection is not None and \
         self.model.ias_selection.count(True) > 0 and \
         len(self.params.refine.strategy) > 0 and \
         self.model.refinement_flags.occupancies
    self.table[self._get_name()] = fl

  def run_occ(self):
    fl = self.model.refinement_flags.s_occupancies is not None and \
         self.model.refinement_flags.occupancies and \
         self.params.main.max_number_of_iterations > 0
    self.table[self._get_name()] = fl

  def run_fp_fdp(self):
    self.table[self._get_name()] = self.model.refinement_flags.group_anomalous

  def run_updatecdl(self):
    fl = self.params.pdb_interpretation.cdl and self.macro_cycle!=1
    self.table[self._get_name()] = fl

  def run_settarget(self):
    self.table[self._get_name()] = True
