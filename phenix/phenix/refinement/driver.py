from __future__ import division
from phenix import refinement
import phenix.refinement.strategies
from mmtbx.refinement import rigid_body
from mmtbx.torsion_restraints.reference_model import reference_model
from mmtbx.torsion_restraints.torsion_ncs import torsion_ncs
import mmtbx.torsion_restraints.utils as torsion_utils
import mmtbx.refinement.print_statistics
import mmtbx.secondary_structure
from mmtbx import ncs
import mmtbx.maps
from copy import deepcopy
import mmtbx.restraints
import mmtbx.ncs.restraints
import iotbx.phil
from cctbx import adptbx
import iotbx.xplor.map
import iotbx.mtz
from iotbx import pdb
from cctbx.array_family import flex
from libtbx.utils import plural_s
from libtbx.utils import Sorry, date_and_time, multi_out, null_out
from libtbx.str_utils import show_string
import copy
import sys, os
from cctbx import xray
from cctbx import geometry_restraints
import mmtbx.model
from mmtbx.tls import tools
from mmtbx.nm import tools as nmtools
from mmtbx.refinement import print_statistics
from phenix.command_line import simple_ncs_from_pdb
from mmtbx import utils
from mmtbx import pdbtools
from mmtbx import model_statistics
from phenix.refinement import pdb_header
import libtbx.callbacks # import dependency
from libtbx.str_utils import format_value, make_header
from libtbx import easy_pickle
from libtbx import easy_mp
from libtbx import Auto
import copy
import iotbx.cif
import cStringIO
from libtbx.test_utils import approx_equal

def overwrite_info():
  return """\
  Options for resolving this problem:
    - Rename or remove the file(s) from previous runs.
    - Add --overwrite to the command line.
    - To permanently disable this safety check, edit your
      .cshrc or .bashrc file and add:
        tcsh: setenv PHENIX_OVERWRITE_ALL true
        bash: export PHENIX_OVERWRITE_ALL=true"""

# XXX make sure there are no parameter conflicts before running
def final_parameter_check (params) :
  if (params.main.use_statistical_model_for_missing_atoms) :
    raise Sorry("Modeling of missing atoms is currently unavailable.")
  if (params.main.ias) :
    if ("rigid_body" in params.refine.strategy) :
      raise Sorry("rigid_body strategy not allowed when main.ias=True.")
  if (params.main.harmonic_restraints) :
    if (params.main.reference_model_restraints) :
      raise Sorry("Reference model restraints may not be combined with "+
        "harmonic restraints to the starting model.")
  if (params.main.place_ions) :
    from mmtbx.ions import check_supported
    check_supported(params.ion_placement.elements)
  if (params.twinning.twin_law not in [None, Auto]) :
    from cctbx import sgtbx
    try :
      rt_mx = sgtbx.rt_mx(symbol=params.twinning.twin_law, r_den=12, t_den=144)
    except ValueError :
      raise Sorry(("Invalid twin_law value '%s'.  If you want to run twinned "+
        "refinement, please double-check the input, or run Xtriage "+
        "(phenix.xtriage) to get a list of possible twin laws for this "+
        "lattice.") % params.twinning.twin_law)
    map_params = params.electron_density_maps.map_coefficients + \
                 params.electron_density_maps.map
    for mp in map_params :
      if (mp.ncs_average) :
        raise Sorry("NCS averaging of output maps not available when "+
          "performing twinned refinement.")

def warn_if_bad_strategy (params, d_min, raise_sorry=False) :
  def _warn (msg) :
    if (raise_sorry) :
      raise Sorry(msg)
    else :
      libtbx.warn(msg)
  if (params.pdb_interpretation.peptide_link.ramachandran_restraints):
    if (d_min < 3.0) :
      _warn("Use of Ramachandran restraints at resolutions better "+
        "than 3.0A is strongly discouraged, as they severely bias "+
        "validation results.  Even at low resolution, you need to check "+
        "the fit of the model to electron density after refining with "+
        "these restraints.")
    else :
      _warn("You have enabled Ramachandran restraints; even at low "+
        "resolution, you need to check the fit of the model to electron "+
        "density after refining with these restraints.")


class refine:

  def __init__(
            self,
            inputs,
            overwrite=False,
            log=None,
            master_params=None):
    self.overwrite = overwrite
    self.inputs = inputs
    self.master_params = master_params
    if (master_params is None) :
      self.master_params = refinement.master_params()
    if(log is None): log = inputs.processed_pdb_file.log
    if(log is None):
       log = sys.stdout
       inputs.processed_pdb_file.log = log
       if(inputs.processed_pdb_file_reference is not None):
          inputs.processed_pdb_file_reference.log = log
    self.log = log
    if(inputs.pdb_inp is not None and inputs.pdb_inp.model_ids().size()>1):
      raise Sorry("Refinement of multiple models not implemented.")
    nproc = easy_mp.enable_multiprocessing_if_possible(
      nproc=inputs.params.refinement.main.nproc,
      log=self.log)
    inputs.params.refinement.main.nproc = nproc
    self.processed_pdb_file = inputs.processed_pdb_file
    self.processed_pdb_file_reference = inputs.processed_pdb_file_reference
    self.pdb_inp = inputs.pdb_inp
    self.params = inputs.params.refinement
    final_parameter_check(self.params)
    self.neutron_r_free_flags = None
    self.f_obs_neutron = None
    self.reference_model_manager = None
    self.den_manager = None
    self.ncs_object = None
    self.cif_block = None
    self.external_energy_manager = None # XXX Rosetta
    warn_if_bad_strategy(
      params=self.params,
      d_min=inputs.f_obs.d_min())
    #
    xsfppf = mmtbx.utils.xray_structures_from_processed_pdb_file(
      processed_pdb_file = self.processed_pdb_file,
      scattering_table   = self.params.main.scattering_table,
      d_min              = inputs.f_obs.d_min(),
      log                = log)
    xray_structure = xsfppf.xray_structures[0].deep_copy_scatterers()
    self.neutron_scattering_dict = xsfppf.neutron_scattering_dict
    self.xray_scattering_dict = xsfppf.xray_scattering_dict
    #
    all_chain_proxies  = self.processed_pdb_file.all_chain_proxies
    self.apply_cif_links = all_chain_proxies.apply_cif_links
    #
    anomalous_scatterer_groups = self.process_anomalous_scatterer_groups()
    potential_non_zero_f_double_prime = False
    for group in anomalous_scatterer_groups:
      group.copy_to_scatterers_in_place(scatterers=xray_structure.scatterers())
      if (group.f_double_prime != 0 or group.refine_f_double_prime):
        potential_non_zero_f_double_prime = True
    #
    sctr_keys = \
           xray_structure.scattering_type_registry().type_count_dict().keys()
    has_hd = "H" in sctr_keys or "D" in sctr_keys
    self.sec_str = mmtbx.secondary_structure.process_structure(
      params=self.params.secondary_structure,
      processed_pdb_file=self.processed_pdb_file,
      tmp_dir=os.getcwd(),
      log=log,
      assume_hydrogens_all_missing=(not has_hd))
    restraints_manager = None
    reference_sites_cart = None
    selection_moving = None
    amber_mdgx_structs = use_amber = None
    if(self.params.main.use_geometry_restraints):
      h_bond_proxies = None
      custom_nb_excl = None
      if self.params.main.secondary_structure_restraints :
        mmtbx.refinement.print_statistics.make_header(
          "Secondary structure", out=log)
        self.sec_str.initialize(log=log)
        build_proxies = self.sec_str.create_hbond_proxies(
          log=log,
          hbond_params=self.params.hydrogen_bonding)
        h_bond_proxies = build_proxies.proxies
        if (self.params.hydrogen_bonding.exclude_nonbonded) :
          custom_nb_excl = build_proxies.exclude_nb_list
      # XXX temporarily disabled 101011
      if (self.params.main.hydrogen_bonds) :
        if (h_bond_proxies is not None) :
          raise Sorry("Secondary structure and generic hydrogen-bond "+
            "restraints may not be used in combination.")
        mmtbx.refinement.print_statistics.make_header("Hydrogen bonds",out=log)
        from mmtbx.geometry_restraints import hbond
        h_bond_proxies = hbond.find_implicit_hydrogen_bonds(
          pdb_hierarchy=self.processed_pdb_file.all_chain_proxies.pdb_hierarchy,
          xray_structure=xray_structure,
          params=self.params.hydrogen_bonding,
          log=log).proxies
      # TODO Amber hooks
      if hasattr(self.params, "amber") :
        use_amber = self.params.amber.use_amber
        if (use_amber) :
          amber_params = self.params.amber
          import amber_adaptbx
          make_header("Initializing AMBER", out=log)
          print >> log, "  topology: %s" % amber_params.topology_file_name
          amber_mdgx_structs = amber_adaptbx.get_amber_structs(
            prmtop=amber_params.topology_file_name,
            ambcrd=amber_params.coordinate_file_name)
          if (self.params.hydrogens.refine.lower() == "auto") :
            self.params.hydrogens.refine = "individual"
            self.params.hydrogens.force_riding_adp = True
      rama_selection = None
      if (self.params.pdb_interpretation.peptide_link.ramachandran_restraints):
        from mmtbx.geometry_restraints import ramachandran as rama_restraints
        rama_selection = rama_restraints.process_refinement_settings(
          params=self.params.pdb_interpretation.peptide_link,
          pdb_hierarchy=self.processed_pdb_file.all_chain_proxies.pdb_hierarchy,
          secondary_structure_manager=self.sec_str,
          log=log)
      if ("den" in self.params.refine.strategy) :
        from mmtbx.den import den_restraints
        if self.params.den.reference_file is not None:
          pdb_io_ref = pdb.input(self.params.den.reference_file)
          pdb_hierarchy_ref = pdb_io_ref.construct_hierarchy()
          pdb_hierarchy_ref.atoms().reset_i_seq()
        else: #restrain model to starting coordinates
          pdb_hierarchy_ref = None
        den_manager = den_restraints(
          pdb_hierarchy=\
            self.processed_pdb_file.all_chain_proxies.pdb_hierarchy,
          pdb_hierarchy_ref=pdb_hierarchy_ref,
          params=self.params.den,
          log=self.log)
        self.den_manager = den_manager
      mmtbx.refinement.print_statistics.make_header(
        "Summary of geometry restraints",
        out = log)
      assert ((self.external_energy_manager is None) or
              (hasattr(self.external_energy_manager, "__call__")))
      geometry = self.processed_pdb_file.geometry_restraints_manager(
        show_energies      = False,
        plain_pairs_radius = self.params.adp_restraints.iso.plain_pairs_radius,
        params_edits       = self.params.geometry_restraints.edits,
        params_remove      = self.params.geometry_restraints.remove,
        hydrogen_bond_proxies = h_bond_proxies,
        hydrogen_bond_params = self.params.hydrogen_bonding,
        custom_nonbonded_exclusions = custom_nb_excl,
        external_energy_function=self.external_energy_manager,
        assume_hydrogens_all_missing = not has_hd,
        ramachandran_atom_selection=rama_selection,
        den_manager=self.den_manager)

      restraints_manager = mmtbx.restraints.manager(
        geometry      = geometry,
        normalization = self.params.main.use_normalized_geometry_target,
        use_amber     = use_amber,
        amber_mdgx_structs = amber_mdgx_structs)
####> Fixing bad ADP in input model
    if("individual_adp" in self.params.refine.strategy or
       "tls" in self.params.refine.strategy or
       "nm" in self.params.refine.strategy or
       "group_adp" in self.params.refine.strategy):
      print_statistics.make_header(
                             "Fixing bad ADP in input model (if any)", out = log)
      xray_structure.tidy_us()
####> Extract f_obs and r-free flags (xray data)
    self.f_obs = inputs.f_obs
    self.r_free_flags = inputs.r_free_flags
    self.test_flag_value = inputs.test_flag_value
    self.r_free_flags_md5_hexdigest = inputs.r_free_flags_md5_hexdigest
####> Extract experimental phases (xray data)
    extract_experimental_phases_obj = extract_experimental_phases(
                              experimental_phases = inputs.experimental_phases,
                              f_obs               = self.f_obs,
                              log                 = log)
    self.experimental_phases = \
                          extract_experimental_phases_obj.experimental_phases()
####> Extract f_obs and r-free flags (neutron data)
    if(inputs.f_obs_neutron is not None):
       self.f_obs_neutron = inputs.f_obs_neutron
       self.neutron_r_free_flags = inputs.r_free_flags_neutron
       self.neutron_test_flag_value = inputs.test_flag_value_neutron
       self.neutron_r_free_flags = inputs.r_free_flags_neutron
       self.neutron_test_flag_value = inputs.test_flag_value_neutron
    self.neutron_refinement = (self.neutron_r_free_flags is not None and
      self.f_obs_neutron is not None) or \
      self.params.main.scattering_table == "neutron"
    if(self.neutron_refinement):
      self.neutron_scattering_dict = xray_structure.deep_copy_scatterers().\
        switch_to_neutron_scattering_dictionary()
    #
    self.set_hydrogens_refine()
####> Check point; START
    #f_calc_start = xray_structure.structure_factors(
    #                                            anomalous_flag = None,
    #                                            d_min          = 2.0,
    #                                            algorithm      = "fft",
    #                                            cos_sin_table  = False,
    #                                            quality_factor = None,
    #                                            u_base         = None,
    #                                            b_base         = None,
    #                                            wing_cutoff    = None).f_calc()
####> Create fake f_obs (replace original f_obs in-place)
    modify_f_obs_result = modify_f_obs(
      params       = self.params,
      f_obs        = self.f_obs,
      r_free_flags = self.r_free_flags,
      log = log)
    self.f_obs = modify_f_obs_result.f_obs
    self.r_free_flags = modify_f_obs_result.r_free_flags
    #
    if(self.params.main.fake_f_obs):
       print_statistics.make_header("Generating fake Fobs", out = log)
       if(self.params.input.xray_data.outliers_rejection):
         print >> log
         print >> log, \
           "Automatic adjustment: input.xray_data.outliers_rejection = False"
         print >> log
         self.params.input.xray_data.outliers_rejection = False
       if(self.params.input.neutron_data.outliers_rejection):
         print >> log
         print >> log, \
           "Automatic adjustment: input.neutron_data.outliers_rejection = False"
         print >> log
         self.params.input.neutron_data.outliers_rejection = False
       fake_f_obs = abs(utils.fmodel_from_xray_structure(
         xray_structure = xray_structure,
         f_obs          = self.f_obs,
         params         = self.params.fake_f_obs).f_model)
       fake_f_obs.set_observation_type_xray_amplitude()
       self.f_obs = self.f_obs.customized_copy(data = fake_f_obs.data())
####> Modify initial model
    print_statistics.make_header("Modifying start model if requested",out= log)
    sites_cart_start = xray_structure.sites_cart()
    modify_model_obj = pdbtools.modify(
                            xray_structure    = xray_structure,
                            params            = self.params.modify_start_model,
                            all_chain_proxies = all_chain_proxies,
                            log               = log)
    xray_structure = modify_model_obj.xray_structure
    sites_cart_final = xray_structure.sites_cart()
    if(restraints_manager is not None and
       restraints_manager.geometry is not None and
       flex.max((sites_cart_start-sites_cart_final).dot())>0):
      restraints_manager.geometry.replace_site_symmetry(new_site_symmetry_table =
        xray_structure.site_symmetry_table())
    ####> omit selection
    if (self.params.modify_start_model.omit_selection is not None) :
      print_statistics.make_header("Processing selection for omit map", out=log)
      all_chain_proxies = self.processed_pdb_file.all_chain_proxies
      omit_selection = all_chain_proxies.selection(
        self.params.modify_start_model.omit_selection)
      if (omit_selection.count(True) == 0) :
        raise Sorry("The atom selection to be omitted does not match any "+
          "atoms in the input model.")
      print >> log, "  %d atoms will have their occupancies set to zero." % \
        omit_selection.count(True)
      xray_structure.set_occupancies(
        value=0,
        selection=omit_selection)
    ####> Fixing bad ADP in input model
    if(self.params.refine.strategy in ["individual_adp","tls","nm","group_adp"]):
      print_statistics.make_header(
                             "Fixing bad ADP in input model (if any)", out = log)
      xray_structure.tidy_us()
    self.extract_tls_selections_from_pdb_file_header()
    if ("tls" in self.params.refine.strategy) :
      if ((self.params.tls.find_automatically == True) or
          ((self.params.tls.find_automatically is Auto) and
           (len(self.params.refine.adp.tls) == 0))) :
        print_statistics.make_header("Identifying TLS groups", out=log)
        from mmtbx.command_line import find_tls_groups
        tls_params = find_tls_groups.master_phil.fetch().extract()
        tls_params.nproc = self.params.main.nproc
        tls_selections = find_tls_groups.find_tls(
          params=tls_params,
          pdb_inp=all_chain_proxies.pdb_inp,
          pdb_hierarchy=all_chain_proxies.pdb_hierarchy,
          xray_structure=xray_structure,
          return_as_list=True,
          out=cStringIO.StringIO())
        print >> log
        for sele_str in tls_selections :
          print >> log, "    %s" % sele_str
        self.params.refine.adp.tls = tls_selections
      if (self.params.main.ordered_solvent) :
        tools.check_tls_selections_for_waters(
          tls_selections=self.params.refine.adp.tls,
          pdb_hierarchy=all_chain_proxies.pdb_hierarchy)
    if ("rigid_body" in self.params.refine.strategy) :
      if (len(self.params.refine.sites.rigid_body) == 0) :
        print_statistics.make_header(
          "Separating polymer chains into rigid bodies", out=log)
        rigid_body_selections = rigid_body.rigid_groups_from_pdb_chains(
          pdb_hierarchy=all_chain_proxies.pdb_hierarchy)
        if (len(rigid_body_selections) == 0) :
          print >> log, "  No rigid groups found."
        else :
          for group_sele in rigid_body_selections :
            print >> log, "  %s" % group_sele
        self.params.refine.sites.rigid_body = rigid_body_selections
    hd_sel = xray_structure.hd_selection()
    # XXX make sure we don't have any anisotropic hydrogen atoms if running
    # conventional X-ray refinement
    if ((hd_sel.count(True) > 0) and
        (inputs.params.refinement.main.scattering_table != "neutron") and
        (inputs.params.refinement.input.neutron_data.file_name is None) and
        ("individual_adp" in inputs.params.refinement.refine.strategy) and
        (inputs.params.refinement.main.number_of_macro_cycles > 0)) :
      xray_structure.convert_to_isotropic(selection=hd_sel.iselection())
####> Extract refinement strategy and set up selections
    self.refinement_strategy_and_selections_obj = \
      extract_refinement_strategy_and_selections(
         params            = inputs.params.refinement,
         xray_structure    = xray_structure,
         all_chain_proxies = all_chain_proxies,
         neutron_refinement= self.neutron_refinement,
         log               = log)
    refinement_flags = \
                 self.refinement_strategy_and_selections_obj.refinement_flags()
    print >> log
    xray_structure.show_scatterer_flags_summary(out = log)
####> Check refinement_flags for sanity:
    refinement_flags.check_all()
####> Automatic TLS+ind ADP refinement parameters adjustment
    if(refinement_flags.tls and refinement_flags.individual_adp):
       self.adjust_params_for_TLSplusADP_refinement()
####> Set up model object
    if(refinement_flags.tls and self.params.main.number_of_macro_cycles > 0):
      tls_groups = tools.tls_groups(selection_strings =
        self.params.refine.adp.tls)
      if(refinement_flags is not None and [refinement_flags,
         refinement_flags.adp_tls].count(None)==0):
        tlsos = tools.generate_tlsos(
          selections     = refinement_flags.adp_tls,
          xray_structure = xray_structure,
          value          = 0.0)
        tls_groups.tlsos = tlsos
    else: tls_groups = None
    if(refinement_flags.nm and self.params.main.number_of_macro_cycles > 0):
      nm_groups = nmtools.nm_groups(selection_strings =
        self.params.refine.adp.nm)
      if(refinement_flags is not None and [refinement_flags,
         refinement_flags.adp_nm].count(None)==0):
        if(self.params.nm.start_xs_value is not None):
          xs = nmtools.nm_params_reader(filename = self.params.nm.start_xs_value)
        else:
          xs = None
        nm_groups.xs = xs
# to be modified
    else: nm_groups = None

    self.model = mmtbx.model.manager(
      processed_pdb_files_srv = inputs.processed_pdb_files_srv,
      refinement_flags = refinement_flags,
      restraints_manager = restraints_manager,
      xray_structure = xray_structure,
      pdb_hierarchy = self.processed_pdb_file.all_chain_proxies.pdb_hierarchy,
      tls_groups = tls_groups,
      nm_groups = nm_groups,
      reference_sites_cart = reference_sites_cart,
      selection_moving = selection_moving,
      anomalous_scatterer_groups = anomalous_scatterer_groups,
#      all_chain_proxies = inputs.processed_pdb_file.all_chain_proxies,
      log = log)
####> Get model_ref
    if(inputs.processed_pdb_file_reference is not None):
       geometry = inputs.processed_pdb_file_reference.geometry_restraints_manager(
        show_energies      = False,
        plain_pairs_radius = self.params.adp_restraints.iso.plain_pairs_radius)
       restraints_manager = mmtbx.restraints.manager(
               geometry      = geometry,
               normalization = self.params.main.use_normalized_geometry_target)
       ippfr = inputs.processed_pdb_file_reference
       self.model_ref = mmtbx.model.manager(
         processed_pdb_files_srv = inputs.processed_pdb_files_srv,
         refinement_flags = refinement_flags,
         restraints_manager = restraints_manager,
         xray_structure = ippfr.xray_structure(),
         pdb_hierarchy = ippfr.all_chain_proxies.pdb_hierarchy,
         log = log)

    else:
      self.model_ref = None
####> Add reference dihedral restraints if requested
    if (self.params.main.reference_model_restraints):
      if (self.params.reference_model.use_starting_model_as_reference and \
          len(self.params.reference_model.file) > 0):
        raise Sorry("Cannot not restrain working model to self and a "+
                    "reference model simultaneously")
      if(not self.params.reference_model.use_distance_based_target):
        self.reference_file_list = []
        if self.params.reference_model.use_starting_model_as_reference:
          self.reference_file_list.append(self.params.input.pdb.file_name[0])
          print >> self.log, \
            "*** Restraining model to starting coordinates ***"
        else:
          for file_name in self.params.reference_model.file:
            self.reference_file_list.append(file_name)
        print >> self.log, "*** Adding Reference Model Restraints ***"
        #test for inserted TER cards in working model
        ter_indices = self.processed_pdb_file.all_chain_proxies.\
          pdb_inp.ter_indices()
        if ter_indices is not None:
          torsion_utils.check_for_internal_chain_ter_records(
            pdb_hierarchy=
              self.processed_pdb_file.all_chain_proxies.pdb_hierarchy,
            ter_indices=ter_indices,
            file_name=self.params.input.pdb.file_name[0])
        rm = reference_model(
          pdb_hierarchy=
            self.processed_pdb_file.all_chain_proxies.pdb_hierarchy,
          reference_file_list=self.reference_file_list,
          has_hd=has_hd,
          params=self.params.reference_model,
          selection=refinement_flags.sites_individual,
          log=self.log)
        rm.add_reference_dihedral_proxies(geometry=geometry)
        rm.show_reference_summary(log=self.log)
        self.reference_model_manager = rm
        if geometry.generic_restraints_manager.\
             c_beta_dihedral_proxies is None:
          torsion_utils.add_c_beta_restraints(
            geometry=geometry,
            pdb_hierarchy=
              self.processed_pdb_file.all_chain_proxies.pdb_hierarchy,
            log=self.log)
####> Harmonic restraints on starting coordinates
    if (self.params.main.harmonic_restraints) :
      assert (not self.params.main.reference_model_restraints)
      sites_start = xray_structure.sites_cart()
      restraints_selection = flex.bool(sites_start.size(), True)
      if (self.params.harmonic_restraints.selection is not None) :
        acp = self.processed_pdb_file.all_chain_proxies
        restraints_selection = acp.selection(
          self.params.harmonic_restraints.selection)
        if (restraints_selection.count(True) == 0) :
          raise Sorry(("No atoms selected for harmonic restraints (input "+
            "selection string: %s)") %
            self.params.harmonic_restraints.selection)
      print >> self.log, "*** Restraining %d atoms to initial coordinates ***"\
        % restraints_selection.count(True)
      isel = restraints_selection.iselection()
      geometry.generic_restraints_manager.reference_manager.\
        add_coordinate_restraints(
          sites_cart=sites_start.select(isel),
          selection=isel,
          sigma=self.params.harmonic_restraints.sigma)
####> Stop refinement if the resolution of input data is too high
    f_obs_d_min = self.f_obs.d_min()
    max_d_min = self.params.main.max_d_min
    if(f_obs_d_min < max_d_min):
       raise Sorry("""\
Proper refinement at a resolution this high is currently not supported.
  Resolution of the data:    %8.4f
  refinement.main.max_d_min: %8.4f
  To proceed: 1) cut the data at a lower resolution, e.g.:
                   main.high_resolution=%.1f
           or 2) change the threshold for this message, e.g.:
                   main.max_d_min=%.1f"""
          % (f_obs_d_min, max_d_min, max_d_min+0.1, f_obs_d_min-0.1))
    if (refinement_flags.group_anomalous):
      self.check_anomalous_consistency()
    self.set_ignore_hydrogens_in_mask_calculation_based_on_data_type()
    self.check_hydrogen_refinement()
    #
    self.fmodel = None
    #
####>
    self.assert_occupancies_not_all_zero()
####> misc setup
    developer_params = getattr(self.params, "developer", None)
    if (self.params.main.place_ions) :
      # more permissive settings for water picking when finding ions
      self.params.main.ordered_solvent = True
      self.params.ordered_solvent.b_iso_min = 0.0
      self.params.ordered_solvent.h_bond_min_mac = 1.0
      self.params.ordered_solvent.h_bond_min_sol = 1.0
      self.params.ordered_solvent.h_bond_max = 6.0
      self.params.ordered_solvent.mode = "every_macro_cycle_after_first"
####> Write out modified model to file
    if(self.params.modify_start_model.output.file_name is not None):
       print_statistics.make_header("Write out modified model", out = log)
       print >> log, self.params.modify_start_model.output.file_name
       ofn = open(self.params.modify_start_model.output.file_name,"w")
       self.model.write_pdb_file(out = ofn)
       ofn.close()
####> Check point: END
    #if(not modified_flag):
    #   f_calc_end = xray_structure.structure_factors(
    #                                            anomalous_flag = None,
    #                                            d_min          = 2.0,
    #                                            algorithm      = "fft",
    #                                            cos_sin_table  = False,
    #                                            quality_factor = None,
    #                                            u_base         = None,
    #                                            b_base         = None,
    #                                            wing_cutoff    = None).f_calc()
    #   bug = abs(flex.max(flex.abs(f_calc_start.data()-f_calc_end.data())))
    #   if(bug > 1.e-6):
    #      print >> log, "Bug value: ", bug
    #      raise Sorry("\nBug encountered. Report it to PAfonine@lbl.gov\n")
    #else:
    #   xrs = self.processed_pdb_file.xray_structure(show_summary = False)
    #   f_calc_end = xrs.structure_factors(anomalous_flag = None,
    #                                      d_min          = 2.0,
    #                                      algorithm      = "fft",
    #                                      cos_sin_table  = False,
    #                                      quality_factor = None,
    #                                      u_base         = None,
    #                                      b_base         = None,
    #                                      wing_cutoff    = None).f_calc()
    #   bug = abs(flex.max(flex.abs(f_calc_start.data()-f_calc_end.data())))
    #   if(bug < 0.1):
    #      print >> log, "Bug value: ", bug
    #      raise Sorry("\nBug encountered. Report it to PAfonine@lbl.gov\n")
####> Write initial parameters into .eff file
#    print_statistics.make_header("Write initial parameters into .eff file", out = log)
#    self.write_eff_file()
####> Analyze / find NCS
    # XXX bug: not sensitive to input params, like weights
    # XXX
    print_statistics.make_header(
                            "Process input NCS or/and find new NCS", out = log)
    if(not self.params.main.ncs):
       print >> log, "Using existing and finding new NCS is disabled."
       print >> log, "Use refinement.main.ncs=true to activate it."
       print >> log, "Look at refinement.ncs for more NCS related parameters."
    elif self.params.ncs.type == "torsion":
      if (len(self.params.ncs.restraint_group) > 0) :
        libtbx.warn("Global NCS restraint groups are defined (parameter scope "+
          "refinement.ncs.restraint_group), but torsion NCS is selected.")
      #test for inserted TER cards
      ter_indices = self.processed_pdb_file.all_chain_proxies.\
        pdb_inp.ter_indices()
      if ter_indices is not None:
        torsion_utils.check_for_internal_chain_ter_records(
          pdb_hierarchy=
            self.processed_pdb_file.all_chain_proxies.pdb_hierarchy,
          ter_indices=ter_indices,
          file_name=self.params.input.pdb.file_name[0])
      if not self.params.ncs.restrain_b_factors:
        print >> self.log, "Not restraining NCS-related b-factors:"
        print >> self.log, "refinement.ncs.b_factor_weight = 0.0\n"
        self.params.ncs.b_factor_weight = 0.0
      grm_ncs = None
      if(self.params.main.use_geometry_restraints):
        geometry.generic_restraints_manager.ncs_manager = torsion_ncs(
          pdb_hierarchy= \
          self.processed_pdb_file.all_chain_proxies.pdb_hierarchy,
          fmodel=self.fmodel,
          params=self.params.ncs.torsion,
          b_factor_weight=self.params.ncs.b_factor_weight,
          coordinate_sigma=self.params.ncs.coordinate_sigma,
          selection=refinement_flags.sites_individual,
          log=log)
        if geometry.generic_restraints_manager.\
           ncs_manager.ncs_dihedral_proxies is None:
          geometry.generic_restraints_manager.\
            ncs_manager = None
        grm_ncs = geometry.generic_restraints_manager.ncs_manager
      if grm_ncs is not None:
        grm_ncs.add_ncs_dihedral_proxies(geometry=geometry)
        if grm_ncs.remove_conflicting_torsion_restraints:
          grm_ncs.sync_dihedral_restraints(geometry=geometry)
        if grm_ncs.found_ncs != None:
          print >> log, "\nNew torsion NCS groups: "
          found_ncs_as_phil = iotbx.phil.parse(grm_ncs.found_ncs)
          found_ncs_as_phil.show(out = self.log)
          import phenix.refinement
          self.params.ncs.torsion.restraint_group = \
            self.master_params.\
            fetch(found_ncs_as_phil).extract().refinement.\
            ncs.torsion.restraint_group
        if([self.params.ncs.coordinate_sigma,
            self.params.ncs.b_factor_weight].count(None) != 2) and \
            self.params.ncs.restrain_b_factors:
          for rg in self.params.ncs.torsion.restraint_group:
            if(self.params.ncs.coordinate_sigma is not None):
              rg.coordinate_sigma = self.params.ncs.coordinate_sigma
            if(self.params.ncs.b_factor_weight is not None):
              rg.b_factor_weight = self.params.ncs.b_factor_weight
        #if self.params.ncs.restrain_b_factors:
        grm_ncs.process_ncs_restraint_groups(
          model=self.model,
          processed_pdb_file=self.processed_pdb_file)
        if self.model.restraints_manager.torsion_ncs_groups is None:
          print >> log, \
            "Torsion NCS pairs could not be determined...skipping."
        if self.reference_model_manager != None and \
           self.params.reference_model.auto_shutoff_for_ncs:
          self.reference_model_manager.remove_restraints_with_ncs_matches(
            ncs_dihedral_proxies=grm_ncs.ncs_dihedral_proxies,
            ncs_match_hash=grm_ncs.ncs_match_hash)
        if geometry.generic_restraints_manager.c_beta_dihedral_proxies is None:
          torsion_utils.add_c_beta_restraints(
            geometry=geometry,
            pdb_hierarchy=
              self.processed_pdb_file.all_chain_proxies.pdb_hierarchy,
            log=self.log)
    else:
       n_ncs_groups = 0
       for i_seq, group in enumerate(self.params.ncs.restraint_group):
           if(group.reference is not None):
              n_ncs_groups += 1
       if(n_ncs_groups == 0):
          print >> log, "NCS groups found in parameter file:", n_ncs_groups
       else:
          print >> log, "NCS groups found in parameter file:", n_ncs_groups
          for i_seq, group in enumerate(self.params.ncs.restraint_group):
              if(group.reference is not None):
                 print >> log, "NCS group %d:"%i_seq
                 print >> log, "  reference=", group.reference
                 for ncs_selection in group.selection:
                     print >> log, "  selection=", ncs_selection
       print >> log
       if(self.params.ncs.find_automatically):
          print >> log, "Automatic NCS search:"
          try :
            self.ncs_object = self.find_cartesian_ncs_groups(log=log)
          except Exception, e:
            print >> log, "*"*30
            print >> log, "\nAutomatic NCS search failed:"
            print >> log, str(e)
            print >> log, "*"*30
            print >> log, "\nphenix.refine continues anyway ..."
          else :
            found_ncs = self.ncs_object.format_all_for_phenix_refine(
                                                 log=None,quiet=True,out=None)
            print >> log, "\nUpdated/new NCS groups: "
            found_ncs_as_phil = iotbx.phil.parse(found_ncs)
            found_ncs_as_phil.show(out = self.log)
            import phenix.refinement
            self.params.ncs.restraint_group = self.master_params.\
              fetch(found_ncs_as_phil).extract().refinement.ncs.restraint_group
            print >> log
       if not self.params.ncs.restrain_b_factors:
         print >> self.log, "Not restraining NCS-related b-factors:"
         print >> self.log, "refinement.ncs.b_factor_weight = 0.0"
         self.params.ncs.b_factor_weight = 0.0
       if([self.params.ncs.coordinate_sigma,
                             self.params.ncs.b_factor_weight].count(None) != 2):
            for rg in self.params.ncs.restraint_group:
                if(self.params.ncs.coordinate_sigma is not None):
                   rg.coordinate_sigma = self.params.ncs.coordinate_sigma
                if(self.params.ncs.b_factor_weight is not None):
                   rg.b_factor_weight = self.params.ncs.b_factor_weight
       self.process_ncs_restraint_groups()
    #
####>
    if(len(self.params.refine.strategy)==1 and
      self.params.refine.strategy[0]=="rigid_body" and
      self.params.rigid_body.mode == "first_macro_cycle_only" and
      self.params.main.number_of_macro_cycles > 1):
      print_statistics.make_header("Automatic adjustment", out = self.log)
      print >> log, "Reason: refine.strategy=rigid_body and"
      print >> log, "        rigid_body.mode=first_macro_cycle_only"
      print >> log, "Old value main.number_of_macro_cycles=",\
                                        self.params.main.number_of_macro_cycles
      self.params.main.number_of_macro_cycles = 1
      print >> log, "New value main.number_of_macro_cycles=", \
                                        self.params.main.number_of_macro_cycles
####> Write initial parameters into .eff file
    print_statistics.make_header("Write initial parameters into .eff file", out = log)
    if (self.params.output.write_eff_file): self.write_eff_file()
####>
#    if (not("nm" in self.params.refine.strategy)):
    self.assert_no_anisotropic_atoms_at_low_resolution()
####> MUST BE THE LAST LINE IN THIS CONSTRUCTOR:
    check_input_params(params         = self.params,
                       f_obs          = self.f_obs,
                       xray_structure = self.model.xray_structure)

  def extract_tls_selections_from_pdb_file_header(self):
    if(len(self.params.refine.adp.tls)==0):
      tls_selections = []
      acp = self.processed_pdb_file.all_chain_proxies
      pdb_inp_tls = acp.pdb_inp.extract_tls_params(acp.pdb_hierarchy)
      if(pdb_inp_tls.tls_present):
        print_statistics.make_header(
          "TLS group selections from PDB file header", out=self.log)
        print >> self.log, "TLS group selections:"
        atom_counts = []
        for t in pdb_inp_tls.tls_params:
          n_atoms = acp.pdb_hierarchy.atom_selection_cache().selection(
            string = t.selection_string).count(True)
          print >> self.log, "  selection string:"
          print >> self.log, "    %s"%t.selection_string
          print >> self.log, "    selects %d atoms"%n_atoms
          tls_selections.append(t.selection_string)
          atom_counts.append(n_atoms)
        if(pdb_inp_tls.tls_present):
          if(pdb_inp_tls.error_string is not None):
            print >> self.log, "  %s"%pdb_inp_tls.error_string
            tls_selections = []
        if(0 in atom_counts):
          msg="""
  One of TLS selections is an empty selection: skipping TLS infromation found in
  PDB file header.
"""
          print >> self.log, msg
          return
        self.params.refine.adp.tls = tls_selections
        return
        # XXX temporarily disabled for 1.7.3 release
        if(len(tls_selections)>0 and
           not "tls" in self.params.refine.strategy and
           len(self.params.refine.strategy)>0):
          msg = """
Automatic adjustment:
  TLS refinement is enabled because valid TLS group selections found in input
  PDB file header. To disable TLS refinement: remove TLS related section from
  PDB file header.
"""
          print >> self.log, msg
          self.params.refine.strategy.append("tls")

  def assert_no_anisotropic_atoms_at_low_resolution(self):
    if(self.f_obs.d_min() > self.params.main.switch_to_isotropic_high_res_limit and
       "individual_adp" in self.params.refine.strategy):
      convert_to_u_iso = self.model.xray_structure.use_u_aniso()
      if(convert_to_u_iso.count(True) > 0):
        rf = self.model.refinement_flags
        if(self.params.refine.adp.individual.anisotropic is not None):
          adp_individual_aniso = utils.get_atom_selections(
             iselection          = False,
             all_chain_proxies   = self.processed_pdb_file.all_chain_proxies,
             selection_strings   = self.params.refine.adp.individual.anisotropic,
             xray_structure      = self.model.xray_structure,
             one_selection_array = True,
             parameter_name="refinement.refine.adp.individual.anisotropic")
          convert_to_u_iso.set_selected(adp_individual_aniso, False)
        if("tls" in self.params.refine.strategy):
          for tls_sel in rf.adp_tls:
            tls_sel = flex.bool(convert_to_u_iso.size(), tls_sel)
            convert_to_u_iso.set_selected(tls_sel, False)
        if("nm" in self.params.refine.strategy):
          for tls_sel in rf.adp_nm:
            tls_sel = flex.bool(convert_to_u_iso.size(), tls_sel)
            convert_to_u_iso.set_selected(tls_sel, False)
        self.model.xray_structure.convert_to_isotropic(selection =
          convert_to_u_iso.iselection())
        if(rf.adp_individual_iso is not None):
          rf.adp_individual_iso.set_selected(convert_to_u_iso, True)
        if(rf.adp_individual_aniso is not None):
          rf.adp_individual_aniso.set_selected(convert_to_u_iso, False)

  def assert_occupancies_not_all_zero(self):
    occ = self.model.xray_structure.scatterers().extract_occupancies()
    if(self.model.xray_structure.scatterers().extract_occupancies().all_eq(0)):
      raise Sorry("All occupancies in input PDB file are zero.")

  def adjust_params_for_TLSplusADP_refinement(self):
    print_statistics.make_header("Automatic adjustment due to TLS refinement.",
                                                                out = self.log)
    print >> self.log, "Automatic adjustment due to TLS refinement. "
    print >> self.log, "  adp_restraints.iso.use_u_local_only   = True"
    print >> self.log, "  adp_restraints.iso.sphere_radius      = 1.55"
    print >> self.log, "  adp_restraints.iso.distance_power     = 0.0 "
    print >> self.log, "  adp_restraints.iso.average_power      = 0.0 "
    print >> self.log, "  adp_restraints.iso.plain_pairs_radius = 5   "
    self.params.adp_restraints.iso.use_u_local_only   = True
    self.params.adp_restraints.iso.sphere_radius      = 1.55
    self.params.adp_restraints.iso.distance_power     = 0.0
    self.params.adp_restraints.iso.average_power      = 0.0
    self.params.adp_restraints.iso.plain_pairs_radius = 5
    d_min = self.f_obs.d_min()
    if(self.params.target_weights.wxu_scale == 1.0):
      if  (d_min < 1.6):                  wxu_scale =     3.00
      elif(d_min >= 1.6 and d_min < 2.0): wxu_scale =     2.72
      elif(d_min >= 2.0 and d_min < 2.4): wxu_scale =     2.16
      elif(d_min >= 2.4 and d_min < 2.8): wxu_scale =     1.81
      elif(d_min >= 2.8 and d_min < 3.2): wxu_scale =     1.66
      else:                               wxu_scale = 0.50
      print >> self.log, "  target_weights.wxu_scale              = %4.2f"%wxu_scale
      self.params.target_weights.wxu_scale = wxu_scale

  def check_anomalous_consistency(self):
    if (not self.f_obs.anomalous_flag()):
      n_groups = 0
      n_scatterers = 0
      for group in self.model.anomalous_scatterer_groups:
        if (group.refine_f_double_prime):
          n_groups += 1
          n_scatterers += group.iselection.size()
      if (n_scatterers != 0):
        raise Sorry(
          "Invalid combination: F-obs not anomalous, but refinement of"
          + " %d f_double_prime value%s requested" % plural_s(n_groups)
          + " (%d scatterer%s):\n" % plural_s(n_scatterers)
          + "  Please provide anomalous reflection data or\n"
          + "  fix all f_double_prime values (see\n"
          + "    refinement.refine.anomalous_scatterer.group\n"
          + "  input blocks).")
      if ("group_anomalous" in self.params.refine.strategy):
        raise Sorry(
          "Invalid combination: F-obs not anomalous, but"
          + " refinement.refine.strategy=group_anomalous:\n"
          + "  Please provide anomalous reflection data or\n"
          + '  remove the "group_anomalous" option from the strategy list.')
    #
    if ("group_anomalous" in self.params.refine.strategy):
      asgs = self.model.anomalous_scatterer_groups
      if (asgs is None): asgs = []
      n_both_fixed = 0
      for group in asgs:
        if (    not group.refine_f_prime
            and not group.refine_f_double_prime):
          n_both_fixed += 1
      if (n_both_fixed == len(asgs)):
        if (len(asgs) == 0):
          raise Sorry(
            "Invalid combination: refinement.refine.strategy=group_anomalous"
            + " but no refinement.refine.anomalous_scatterers.group:\n"
            + "  Please supply one or more"
            + " refinement.refine.anomalous_scatterers.group scopes or\n"
            + '  remove the "group_anomalous" option from the strategy list.')
        raise Sorry(
          "Invalid combination: refinement.refine.strategy=group_anomalous"
          + " but no refinable f_prime and f_double_prime parameters:\n"
          + "  Please enable refinement of one or more f_prime or"
          + " f_double_prime in the\n"
          + "  refinement.refine.anomalous_scatterers.group scopes or\n"
          + '  remove the "group_anomalous" option from the strategy list.')

  def set_ignore_hydrogens_in_mask_calculation_based_on_data_type(self):
    log = self.log
    params = self.params
    if(params.main.scattering_table == "neutron" and
       params.mask.ignore_hydrogens):
      params.mask.ignore_hydrogens=False
      print >> log, "*"*79
      print >> log, "Automatic adjustment due to refinement against neutron data:"
      print >> log, "  params.mask.ignore_hydrogens=False"
      print >> log, "*"*79
      print >> log

  def set_hydrogens_refine(self):
    log = self.log
    params = self.params
    # XXX not sure if it's the right way of using Auto
    if(params.hydrogens.refine is None or
       params.hydrogens.refine.lower() == "auto"):
      if(params.main.scattering_table == "neutron"):
        if(self.f_obs.d_min()>1.5):
          params.hydrogens.refine="riding"
        else:
          params.hydrogens.refine="individual"
      else:
        if(self.f_obs.d_min()<0.9):
          params.hydrogens.refine="individual"
        else:
          params.hydrogens.refine="riding"
      print >> log, "*"*79
      print >> log, "Automatic adjustment:"
      print >> log, "  hydrogens.refine=%s"%params.hydrogens.refine
      print >> log, "*"*79
      print >> log

  # XXX this seems like a somewhat arbitrary limitation...
  def check_hydrogen_refinement (self) :
    if ((self.params.hydrogens.refine == "individual") and
        (self.params.main.scattering_table != "neutron") and
        (self.params.input.neutron_data.file_name is None)) :
      if (self.params.refine.adp.individual.anisotropic is not None):
        adp_individual_aniso = utils.get_atom_selections(
           iselection          = False,
           all_chain_proxies   = self.processed_pdb_file.all_chain_proxies,
           selection_strings   = self.params.refine.adp.individual.anisotropic,
           xray_structure      = self.model.xray_structure,
           one_selection_array = True,
           parameter_name="refinement.refine.adp.individual.anisotropic")
        hd_selection = self.model.xray_structure.hd_selection()
        if ((adp_individual_aniso & hd_selection).count(True) > 0) :
          raise Sorry("Your refinement strategy includes anisotropic ADPs "+
            "for hydrogen atoms refined individually; hydrogens should be "+
            "left as isotropic even for individual refinement.")

  def process_ncs_restraint_groups(self):
    log = self.log
    mmtbx.refinement.print_statistics.make_sub_header(
                                            "Building NCS restraints", out=log)
    ncs_groups = ncs.restraints.groups()
    sites_cart = None
    for param_group in self.params.ncs.restraint_group:
      print >> log, "NCS restraint group %d:" % (
        len(ncs_groups.members)+1)
      group = ncs.restraints.group.from_atom_selections(
        processed_pdb              = self.processed_pdb_file,
        reference_selection_string = param_group.reference,
        selection_strings          = param_group.selection,
        coordinate_sigma           = param_group.coordinate_sigma,
        b_factor_weight            = param_group.b_factor_weight,
        special_position_warnings_only
          = self.params.ncs.special_position_warnings_only,
        log = log)
      if (sites_cart is None):
        sites_cart = self.model.xray_structure.sites_cart()
      ncs_operators = group.operators(sites_cart=sites_cart)
      ncs_operators.show(sites_cart=sites_cart, out=log, prefix="  ")
      ncs_groups.members.append(group)
      print >> log
    if (len(ncs_groups.members) == 0):
      print >> log, "No NCS restraint groups specified."
      print >> log
    else:
      self.model.restraints_manager.ncs_groups = ncs_groups

  def active_anomalous_scatterer_group_scopes(params):
    for group in params.refine.anomalous_scatterers.group:
      if (   group.selection is not None
          or group.f_prime != 0
          or group.f_double_prime != 0
          or len(group.refine) != 2):
        yield group
  active_anomalous_scatterer_group_scopes = staticmethod(
    active_anomalous_scatterer_group_scopes)

  def process_anomalous_scatterer_groups(self):
    log = self.log
    mmtbx.refinement.print_statistics.make_header(
      "Anomalous scatterer groups", out=log)
    result = []
    n_anomalous_total = 0
    if (self.params.group_anomalous.find_automatically) :
      from mmtbx.refinement import anomalous_scatterer_groups
      result = anomalous_scatterer_groups.find_anomalous_scatterer_groups(
        pdb_atoms=self.processed_pdb_file.all_chain_proxies.pdb_atoms,
        xray_structure=self.processed_pdb_file.xray_structure(),
        group_same_element=False,
        out=self.log)
      for group in result :
        n_anomalous_total += group.iselection.size()
      if (n_anomalous_total == 0) :
        if ("group_anomalous" in self.params.refine.strategy) :
          self.params.refine.strategy.remove("group_anomalous")
    else :
      groups = list(self.active_anomalous_scatterer_group_scopes(
        params=self.params))
      if (len(groups) != 0):
        chain_proxies = self.processed_pdb_file.all_chain_proxies
        sel_cache = chain_proxies.pdb_hierarchy.atom_selection_cache()
        for group in groups:
          if (group.f_prime is None): group.f_prime = 0
          if (group.f_double_prime is None): group.f_double_prime = 0
          aag = xray.anomalous_scatterer_group(
            iselection=chain_proxies.phil_atom_selection(
              cache=sel_cache,
              scope_extract=group,
              attr="selection",
              raise_if_empty_selection=True).iselection(),
            f_prime=group.f_prime,
            f_double_prime=group.f_double_prime,
            refine=group.refine,
            selection_string=group.selection)
          aag.show_summary(out=log)
          result.append(aag)
          n_anomalous_total += aag.iselection.size()
          print >> log
    if (len(result) == 0):
      print >> log, "All atoms refined with f_prime=0 and f_double_prime=0."
    else:
      print >> log, "Total number of atoms in anomalous groups:", \
        n_anomalous_total
    return result

  def copy_anomalous_scatterer_groups_to_defaults(self, defaults):
    scatterers = self.fmodel.xray_structure.scatterers()
    scopes = list(self.active_anomalous_scatterer_group_scopes(
      params=defaults.refinement))
    if (self.model.anomalous_scatterer_groups is None):
      assert (len(scopes) == 0)
      return
    elif (self.params.group_anomalous.find_automatically) :
      return # XXX should we be saving them?
    assert len(self.model.anomalous_scatterer_groups) == len(scopes)
    for group,params in zip(self.model.anomalous_scatterer_groups, scopes):
      assert group.selection_string == params.selection
      group.extract_from_scatterers_in_place(scatterers=scatterers)
      params.f_prime = group.f_prime
      params.f_double_prime = group.f_double_prime

  def open_output_file_helper(self,
        output_type,
        ext,
        message=None,
        increment_serial=False):
    assert self.params.output.prefix is not None
    assert self.params.output.serial is not None
    serial = self.params.output.serial
    if (increment_serial): serial += 1
    file_name = os.path.abspath(
      ("%s_"+self.params.output.serial_format+ext) % (
        self.params.output.prefix, serial))
    if (os.path.exists(file_name) and not self.overwrite):
      raise Sorry("%s file exists already: %s\n%s" % (
        output_type, file_name, overwrite_info()))
    if (message is not None):
      print >> self.log, "%s:\n  %s" % (message, file_name)
    return file_name

  def open_output_file(self,
        output_type,
        ext,
        message=None,
        increment_serial=False):
    file_name = self.open_output_file_helper(
      output_type = output_type,
      ext = ext,
      message = message,
      increment_serial = increment_serial)
    return open(file_name, "w")

  def start_log_file(self):
    self.log.replace_stringio(
      old_label="log_buffer",
      new_label="log",
      new_file_object=self.open_output_file(output_type="Log", ext=".log"))
    if ("stdout" not in self.log.labels):
      err_log = multi_out()
      err_log.register(label="stdout", file_object=sys.stdout)
      err_log.register(label="log", file_object=self.log)
      sys.stderr = err_log

  def write_params_file(f, title_line, params, master_params=None):
    print >> f, "#", title_line
    print >> f, "#", date_and_time()
    print >> f
    print >> f, "# Command to extract only non-defaults:"
    print >> f, "#   phenix.refine --diff-params %s" % show_string(
      os.path.basename(f.name))
    print >> f
    if master_params is None:
      master_params = refinement.master_params()
    master_params.format(python_object=params).show(out=f)
  write_params_file = staticmethod(write_params_file)

  def write_eff_file(self):
    f = self.open_output_file(
      output_type="Parameter",
      ext=".eff",
      message="Writing effective parameters to file")
    self.write_params_file(
      f=f,
      title_line="Effective refinement parameters",
      params=self.inputs.params,
      master_params=self.master_params)
    f.close()
    print >> self.log

  def write_geo_file(self):
    if(self.model.restraints_manager is None): return
    f = self.open_output_file(
      output_type="Geometry restraints",
      ext=".geo",
      message="Writing geometry restraints to file")
    print >> f, "# Geometry restraints before refinement"
    print >> f, "#", date_and_time()
    print >> f
    xray_structure = self.model.xray_structure
    sites_cart = xray_structure.sites_cart()
    site_labels = xray_structure.scatterers().extract_labels()
    self.model.restraints_manager.geometry.show_sorted(
      sites_cart=sites_cart,
      site_labels=site_labels,
      f=f)
    n_excessive = 0
    if (self.model.restraints_manager.ncs_groups is not None):
      n_excessive = self.model.restraints_manager.ncs_groups \
        .show_sites_distances_to_average(
          sites_cart=sites_cart,
          site_labels=site_labels,
          excessive_distance_limit=self.params.ncs.excessive_distance_limit,
          out=f)
      print >> f
      self.model.restraints_manager.ncs_groups \
        .show_adp_iso_differences_to_average(
          u_isos=xray_structure.extract_u_iso_or_u_equiv(),
          site_labels=site_labels,
          out=f)
      print >> f
      self.processed_pdb_file.show_atoms_without_ncs_restraints(
        ncs_restraints_groups=self.model.restraints_manager.ncs_groups,
        out=f)
      print >> f
    f.close()
    if (n_excessive != 0):
      raise Sorry("Excessive distances to NCS averages:\n"
        + "  Please inspect the file\n"
        + "    %s\n" % show_string(f.name)
        + "  for a full listing of the distances to the NCS averages.\n"
        + '  Look for the word "EXCESSIVE".\n'
        + "  The current limit is:\n"
        + "    refinement.ncs.excessive_distance_limit=%.6g\n"
            % self.params.ncs.excessive_distance_limit
        + "  The number of distances exceeding this limit is: %d\n"
            % n_excessive
        + "  Please correct your model or redefine the limit, e.g. with:\n"
        + "    refinement.ncs.excessive_distance_limit=%.2g\n"
            % abs(2*self.params.ncs.excessive_distance_limit)
        + "  To disable this message completely define:\n"
        + "    refinement.ncs.excessive_distance_limit=None")
    print >> self.log

  def write_final_geo_file(self):
    if(self.model.restraints_manager is None): return
    f = self.open_output_file(
      output_type="Geometry restraints",
      ext="_final.geo",
      message="Writing final geometry restraints to file")
    print >> f, "# Geometry restraints after refinement"
    print >> f, "#", date_and_time()
    print >> f
    xray_structure = self.model.xray_structure
    sites_cart = xray_structure.sites_cart()
    site_labels = xray_structure.scatterers().extract_labels()
    self.model.restraints_manager.geometry.show_sorted(
      sites_cart=sites_cart,
      site_labels=site_labels,
      f=f)
    n_excessive = 0
    if (self.model.restraints_manager.ncs_groups is not None):
      n_excessive = self.model.restraints_manager.ncs_groups \
        .show_sites_distances_to_average(
          sites_cart=sites_cart,
          site_labels=site_labels,
          excessive_distance_limit=self.params.ncs.excessive_distance_limit,
          out=f)
      print >> f
      self.model.restraints_manager.ncs_groups \
        .show_adp_iso_differences_to_average(
          u_isos=xray_structure.extract_u_iso_or_u_equiv(),
          site_labels=site_labels,
          out=f)
      print >> f
      self.processed_pdb_file.show_atoms_without_ncs_restraints(
        ncs_restraints_groups=self.model.restraints_manager.ncs_groups,
        out=f)
      print >> f
    f.close()
    if (n_excessive != 0):
      raise Sorry("Excessive distances to NCS averages:\n"
        + "  Please inspect the file\n"
        + "    %s\n" % show_string(f.name)
        + "  for a full listing of the distances to the NCS averages.\n"
        + '  Look for the word "EXCESSIVE".\n'
        + "  The current limit is:\n"
        + "    refinement.ncs.excessive_distance_limit=%.6g\n"
            % self.params.ncs.excessive_distance_limit
        + "  The number of distances exceeding this limit is: %d\n"
            % n_excessive
        + "  Please correct your model or redefine the limit, e.g. with:\n"
        + "    refinement.ncs.excessive_distance_limit=%.2g\n"
            % abs(2*self.params.ncs.excessive_distance_limit)
        + "  To disable this message completely define:\n"
        + "    refinement.ncs.excessive_distance_limit=None")
    print >> self.log

  def open_defaults_file(self):
    return self.open_output_file(
      output_type="Parameter default",
      ext=".def",
      increment_serial=True)

  def write_defaults_file(self, f, output_mtz_file_name, output_mtz_labels):
    print >> self.log, \
      "Writing default parameters for subsequent refinement:\n  %s" % f.name
    defaults = self.master_params.clone(
      python_object=self.inputs.params,
      converter_registry=iotbx.phil.default_converter_registry)
    defaults.refinement.output.serial += 1
    ri = defaults.refinement.input
    ri.pdb.file_name = [os.path.abspath(
      ("%s_"+self.params.output.serial_format+".pdb") % (
        self.params.output.prefix, self.params.output.serial))]
    if (output_mtz_file_name is not None):
      ri.xray_data.file_name = output_mtz_file_name
      ri.xray_data.labels = output_mtz_labels.f_obs_xray
      ri.xray_data.r_free_flags.file_name = output_mtz_file_name
      ri.xray_data.r_free_flags.label = output_mtz_labels.r_free_flags_xray
      if(self.f_obs_neutron is not None):
        ri.neutron_data.file_name = output_mtz_file_name
        ri.neutron_data.r_free_flags.file_name = output_mtz_file_name
        ri.neutron_data.labels = output_mtz_labels.f_obs_neutron
        ri.neutron_data.r_free_flags.label = output_mtz_labels.r_free_flags_neutron
      if (ri.experimental_phases.file_name is not None):
        ri.experimental_phases.file_name = output_mtz_file_name
      ri.experimental_phases.labels = None
    defaults.refinement.main.random_seed = int(
      (defaults.refinement.main.random_seed+92365) % (2**31-1))
    defaults.refinement.input.xray_data.r_free_flags.generate = False
    defaults.refinement.input.neutron_data.r_free_flags.generate = False
    self.copy_anomalous_scatterer_groups_to_defaults(defaults=defaults)
    self.write_params_file(
      f=f,
      title_line="Default parameters for subsequent refinement",
      params=defaults,
      master_params=self.master_params)
    f.close()
    print >> self.log

  def open_final_f_model(self, format):
    return self.open_output_file(
      output_type=format.upper(), ext="_f_model.%s" % format)

  def export_final_f_model(self, f, format):
    print >> self.log, \
      "Writing refined fmodel to %s file:\n  %s" % (format.upper(), f.name)
    self.fmodel.export(out=f, format=format)
    print >> self.log

  def open_refined_pdb(self):
    return self.open_output_file(output_type="PDB", ext=".pdb")

  def write_fmodel_pickle(self):
    file_name = self.open_output_file_helper(
        output_type = "fmodel pickle file",
        ext = "_fmodel.pickle",
        message="Writing fmodel pickle file",
        increment_serial=False)
    easy_pickle.dump(file_name, self.fmodel) # XXX joint XN refinement

  def write_stats_pickle (self) :
    file_name = self.open_output_file_helper(
      output_type="monitor statistics pickle",
      ext="_stats.pkl",
      message="Writing statistics pickle file",
      increment_serial=False)
    self.monitor.dump_statistics(file_name)

  def open_refined_model_cif(self):
    return self.open_output_file(output_type="CIF_model", ext=".cif")

  def add_phenix_refine_citation_to_cif_block(self):
    import iotbx.cif.model
    from phenix import phenix_info
    phenix_version, tag = phenix_info.version_and_release_tag()
    if phenix_version is not None:
      if tag is not None:
        phenix_version = "_".join((phenix_version, tag))
    else:
      phenix_version = "?"
    if (libtbx.env.is_development_environment()):
      phenix_version += "+SVN"
    self.cif_block['_computing.structure_refinement'] \
      = 'PHENIX (phenix.refine: %s)' %phenix_version
    software_loop = iotbx.cif.model.loop(header=(
      '_software.pdbx_ordinal',
      '_software.name',
      '_software.version',
      #'_software.date',
      '_software.type',
      '_software.contact_author',
      '_software.contact_author_email',
      '_software.location',
      '_software.classification',
      '_software.citation_id',
      '_software.language',
    ))
    software_loop.add_row((
      '1', 'phenix.refine', phenix_version, 'program',
      'Paul D. Adams', 'pdadams@lbl.gov', 'https://www.phenix-online.org/',
      'refinement', 'phenix.refine', 'Python/C++',))
    software_loop.add_row((
      '1', 'Phenix', phenix_version, 'program',
      'Paul D. Adams', 'pdadams@lbl.gov', 'https://www.phenix-online.org/',
      'refinement', 'phenix', 'Python/C++',))
    self.cif_block.add_loop(software_loop)
    from phenix.utilities import citations
    citations.citations_as_cif_block(
      ('phenix.refine', 'phenix'), cif_block=self.cif_block)

  def write_model_cif_file(self, f):
    print >> self.log, "Writing refined structure to CIF file:\n  %s" % f.name
    cif = iotbx.cif.model.cif()
    cif_block = iotbx.cif.model.block()
    cif_block.add_loop(iotbx.cif.atom_type_cif_loop(
      self.processed_pdb_file.xray_structure(), format="mmcif"))
    if (not self.inputs.params.refinement.main.wavelength in [None, Auto]) :
      cif_block['_diffrn_source.pdbx_wavelength'] = \
        "%g" % self.inputs.params.refinement.main.wavelength
    exptl_method = []
    if self.xray_scattering_dict is not None: exptl_method.append(
      "X-RAY DIFFRACTION")
    if self.neutron_scattering_dict is not None: exptl_method.append(
      "NEUTRON DIFFRACTION")
    assert len(exptl_method) > 0
    if len(exptl_method) == 1:
      exptl_method = exptl_method[0]
    cif_block["_exptl.method"] = exptl_method
    cif_block.update(self.model.pdb_hierarchy().as_cif_block(
      crystal_symmetry=self.fmodel.f_obs().crystal_symmetry()))
    if self.cif_block is not None:
      self.add_phenix_refine_citation_to_cif_block()
      self.cif_block.update(cif_block)
      cif_block = self.cif_block
    cif[self.params.output.prefix] = cif_block
    f = open(f.name, 'w')
    print >> f, cif

  def write_refined_pdb(self, f):
    joint_xn_flag = [self.neutron_scattering_dict,
      self.xray_scattering_dict].count(None)==0
    print >> self.log, "Writing refined structure to PDB file:\n  %s" % f.name
    inp = self.inputs.params.refinement.input
    pr = "REMARK "
    print >>f,"REMARK", date_and_time()
    print >>f,pr+"PHENIX refinement"
    print >>f,pr
    print >>f,pr+"****************** INPUT FILES AND LABELS ******************************"
    if(joint_xn_flag):
      print >>f,pr+"X-RAY DATA."
    print >>f,pr+"Reflections:"
    print >>f,pr+"  file name      : %-s"%inp.xray_data.file_name
    print >>f,pr+"  labels         : %-s"%inp.xray_data.labels
    print >>f,pr+"R-free flags:"
    print >>f,pr+"  file name      : %-s"%inp.xray_data.r_free_flags.file_name
    print >>f,pr+"  label          : %-s"%inp.xray_data.r_free_flags.label
    print >>f,pr+"  test_flag_value: %-d"%inp.xray_data.r_free_flags.test_flag_value
    phin = inp.experimental_phases
    if((phin.file_name,phin.labels) != (None, None)):
      print >>f,pr+"Experimental phase information:"
      print >>f,pr+"  file name      : %-s"%phin.file_name
      print >>f,pr+"  label          : %-s"%phin.labels
    if(joint_xn_flag):
      print >>f,pr+"NEUTRON DATA."
      inp_n = self.inputs.params.refinement.input
      print >>f,pr+"Reflections:"
      print >>f,pr+"  file name      : %-s"%inp_n.neutron_data.file_name
      print >>f,pr+"  labels         : %-s"%inp_n.neutron_data.labels
      print >>f,pr+"R-free flags:"
      print >>f,pr+"  file name      : %-s"%inp_n.neutron_data.r_free_flags.file_name
      print >>f,pr+"  label          : %-s"%inp_n.neutron_data.r_free_flags.label
      print >>f,pr+"  test_flag_value: %-d"%inp_n.neutron_data.r_free_flags.test_flag_value
    print >>f,pr+"Model file name(s): "
    for pdb_file_name in inp.pdb.file_name:
      print >>f,pr+"  "+pdb_file_name
    print >>f,pr
    print >>f,pr+"******************** REFINEMENT SUMMARY: QUICK FACTS *******************"
    rw_start = self.monitor.r_works[0]
    rf_start = self.monitor.r_frees[0]
    if(len(self.monitor.steps) > 1):
       if(self.monitor.steps[1].count("bss")):
          rw_start = self.monitor.r_works[1]
          rf_start = self.monitor.r_frees[1]
    rw_final = self.monitor.r_works[len(self.monitor.r_works)-1]
    rf_final = self.monitor.r_frees[len(self.monitor.r_frees)-1]
    rw, rf = self.fmodel.r_work(), self.fmodel.r_free()
    assert approx_equal(rw, rw_final), [rw, rw_final]
    assert approx_equal(rf, rf_final), [rf, rf_final]
    if(joint_xn_flag):
      print >>f,pr+"X-ray data:"
    print >>f,pr+\
      "Start: r_work = %6.4f r_free = %6.4f bonds = %s angles = %s" % \
      (rw_start,rf_start,format_value("%5.3f",self.monitor.bond_start),
       format_value("%5.3f",self.monitor.angle_start))
    print >>f,pr+\
      "Final: r_work = %6.4f r_free = %6.4f bonds = %s angles = %s" % \
      (rw_final,rf_final,format_value("%5.3f",self.monitor.bond_final),
       format_value("%5.3f",self.monitor.angle_final))
    #
    #if(flex.max(abs(self.fmodels.fmodel_xray().f_part1()).data())!=0):
    #  raise RuntimeError("Not immplemented")
      #fmxdc = self.fmodels.fmodel_xray().deep_copy()
      #fmxdc.update_core(f_part)
      #fmxdc.update_solvent_and_scale()
      #print >>f,pr+\
      #"Final: r_work = %6.4f r_free = %6.4f (no F_part contribution)" % \
      #(fmxdc.r_work(),fmxdc.r_free())
      #del fmxdc
    #
    if(joint_xn_flag):
      print >>f,pr+"Neutron data:"
      rw_start = self.monitor_neutron.r_works[0]
      rf_start = self.monitor_neutron.r_frees[0]
      if(len(self.monitor_neutron.steps) > 1):
         if(self.monitor_neutron.steps[1].count("bss")):
            rw_start = self.monitor_neutron.r_works[1]
            rf_start = self.monitor_neutron.r_frees[1]
      rw_final = \
        self.monitor_neutron.r_works[len(self.monitor_neutron.r_works)-1]
      rf_final = \
        self.monitor_neutron.r_frees[len(self.monitor_neutron.r_frees)-1]
      print >>f,pr+"Start: r_work = %6.4f r_free = %6.4f"%(rw_start,rf_start)
      print >>f,pr+"Final: r_work = %6.4f r_free = %6.4f"%(rw_final,rf_final)
    print >>f,pr+"*"*72
    print >>f,pr
    if("rigid_body" in self.params.refine.strategy):
      print >>f,pr+"Rigid body refinement target: %-s" % \
        self.params.rigid_body.target
    self.monitor.show(out=f, remark=pr)
    model_statistics.model_content(model = self.model).show(out = f,
      prefix = pr, pdb_deposition = False)
    print >> f, pr+"-"*71
    if(joint_xn_flag):
      self.monitor_neutron.show(out=f, remark=pr)
    print >>f,pr+"r_free_flags.md5.hexdigest",self.r_free_flags_md5_hexdigest
    print >>f,pr
    self.model.xray_structure.show_scatterer_flags_summary(out = self.log)
    print >>f,pr+\
      "IF THIS FILE IS FOR PDB DEPOSITION: REMOVE ALL FROM THIS LINE UP."
    m_info =model_statistics.info(model = self.model, fmodel_x = self.fmodel,
      fmodel_n = self.fmodel_neutron, refinement_params = self.params,
      ignore_hd = not self.neutron_refinement,
      use_molprobity=self.inputs.params.refinement.main.use_molprobity)
    pdb_header.write_remark_3(info = m_info, out = f)
    if self.params.output.write_model_cif_file:
      self.cif_block = m_info.as_cif_block()
    if(not self.params.main.ias): # XXX selection is broken for ias option: bug
      m = self.model.deep_copy()
    else: m = self.model
    fl = not joint_xn_flag and not self.neutron_refinement and \
       self.params.hydrogens.optimize_scattering_contribution and not \
       self.params.main.ias
    if(fl):
      m.reset_occupancies_for_hydrogens()
      if(self.params.main.truncate_at_pdb_format_precision):
        m.xray_structure.truncate_at_pdb_format_precision()
    ###DEBUG
    #from libtbx import easy_pickle
    #easy_pickle.dump("fmodel_pr_%s.pickle"%self.params.output.prefix, self.fmodel)
    ####
    m.write_pdb_file(out = f)
    f.close()
    print >> self.log

  def open_map_files(self):
    return map_manager(refine_object=self)

  def box_around_selection(self, iselection, buffer):
    structure = self.processed_pdb_file.xray_structure() # XXX self.model.xray_structure
    sites_cart = structure.sites_cart()
    if (iselection is not None):
      sites_cart = sites_cart.select(iselection)
    return structure.unit_cell().box_frac_around_sites(
      sites_cart=sites_cart, buffer=buffer)

  def find_cartesian_ncs_groups (self, log, hierarchy=None) :
    scattering_types = \
      self.model.xray_structure.scatterers().extract_scattering_types()
    number_of_d = (scattering_types == "D").count(True)
    number_of_h = (scattering_types == "H").count(True)
    if number_of_d > 0:
        exclude_d=True
    else:
        exclude_d=False
    if number_of_h > 0:
        exclude_h=True
    else:
        exclude_h=False

    ncs_search = simple_ncs_from_pdb.simple_ncs_from_pdb(
      params               = self.params.ncs,
      suggested_ncs_groups = self.params.ncs.restraint_group,
      pdb_inp              = self.pdb_inp,
      source_info          = str(self.params.input.pdb.file_name),
      all_chain_proxies    = self.processed_pdb_file.all_chain_proxies,
      hierarchy            = hierarchy,
      suppress_print       = True,
      exclude_h            = exclude_h,
      exclude_d            = exclude_d,
      log                  = log)
    return ncs_search.ncs_object

  def run(self, call_back_handler=None):
    log = self.log
    mmtbx.refinement.print_statistics.make_header(
      "Non-default parameters", out = log)
    print >> log, "A complete record of all parameters was written to the" \
      " .eff file above."
    print >> log, "Below are only the non-defaults."
    print >> log
    print >> log, "#phil __ON__"
    mp = self.master_params
    mp.fetch_diff(source=mp.format(
      python_object=self.inputs.params)).show(out=log)
    print >> log, "#phil __OFF__"
    print >> log
    params = self.params
    flags = self.r_free_flags.data()
    if(self.neutron_r_free_flags is not None):
      flags_neutron_data = self.neutron_r_free_flags.data()
    else:
      flags_neutron_data = None
    r_free_flags = self.f_obs.array(data = flags)
    del flags
    self.f_obs = self.f_obs.sort(reverse=True, by_value="packed_indices")
    r_free_flags = r_free_flags.sort(reverse=True, by_value="packed_indices")
    if(self.experimental_phases):
      self.experimental_phases = self.experimental_phases.sort(reverse=True, by_value="packed_indices")
    ###
    result = refinement.strategies.run(
      model                   = self.model,
      verbose                 = 1,
      params                  = params,
      f_obs                   = self.f_obs,
      flags                   = r_free_flags,
      model_ref               = self.model_ref,
      f_obs_neutron           = self.f_obs_neutron,
      flags_neutron           = flags_neutron_data,
      neutron_scattering_dict = self.neutron_scattering_dict,
      xray_scattering_dict    = self.xray_scattering_dict,
      abcd                    = self.experimental_phases,
      pdb_inp                 = self.pdb_inp,
      apply_cif_links         = self.apply_cif_links,
      log                     = log,
      call_back_handler       = call_back_handler,
      reference_model_manager = self.reference_model_manager).result
    self.model, self.fmodels, self.fmodel, self.fmodel_neutron, self.monitor, \
         self.monitor_neutron, self.tlsos = result[0], result[1],\
           result[1].fmodel_xray(), result[1].fmodel_neutron(), \
           result[2].monitor_xray, result[2].monitor_neutron, result[3]
    mmtbx.refinement.print_statistics.make_header(
                        "overall refinement statistics: step by step", out=log)
    self.monitor.show()
    if(self.monitor_neutron is not None):
       self.monitor_neutron.show()
    if (self.model.restraints_manager is not None and
        self.model.restraints_manager.ncs_groups is not None):
      mmtbx.refinement.print_statistics.make_header(
        "NCS restraints at end of refinement", out=log)
      self.model.restraints_manager.ncs_groups.show_operators(
                        sites_cart = self.model.xray_structure.sites_cart(),
                        out        = log)

def combine_data_and_map_coeffs(
      original_data_mtz_dataset,
      map_coeff_dataset,
      fmodels,
      file_name,
      wavelength=None):
  class labels_decorator:
    def __init__(self, amplitudes_label, phases_label):
      self._amplitudes = amplitudes_label
      self._phases = phases_label
    def amplitudes(self):
      return self._amplitudes
    def phases(self, root_label, anomalous_sign=None):
      assert anomalous_sign is None or not anomalous_sign
      return self._phases
  if (wavelength in [None, Auto]) :
    wavelength = 1
  xray_suffix="_xray"
  neutron_suffix="_neutron"
  if(fmodels.fmodel_neutron() is None):
    xray_suffix=""
    neutron_suffix=""
  original_data_mtz_dataset.set_name("Original-experimental-data-mapped-to-asu")
  new_dataset = original_data_mtz_dataset.mtz_crystal().add_dataset(
    name = "Experimental-data-used-in-refinement",
    wavelength=wavelength)
  new_dataset.add_miller_array(
    miller_array = fmodels.fmodel_xray().arrays.f_obs,
    column_root_label="F-obs-filtered"+xray_suffix)
  another_dataset = new_dataset.mtz_crystal().add_dataset(
    name = "Model-structure-factors-(bulk-solvent-and-all-scales-included)",
    wavelength=wavelength)
  another_dataset.add_miller_array(
    miller_array = fmodels.fmodel_xray().f_model_scaled_with_k1_composite_work_free(),
    column_root_label="F-model"+xray_suffix)
  if(fmodels.fmodel_neutron() is not None):
    another_dataset.add_miller_array(
      miller_array = fmodels.fmodel_neutron().arrays.f_obs,
      column_root_label="F-obs-filtered"+neutron_suffix)
    another_dataset.add_miller_array(
      miller_array = fmodels.fmodel_neutron().f_model_scaled_with_k1_composite_work_free(),
      column_root_label="F-model"+neutron_suffix)
  yet_another_dataset = another_dataset.mtz_crystal().add_dataset(
    name = "Fourier-map-coefficients", wavelength=1)
  if(map_coeff_dataset is not None):
    for ma in map_coeff_dataset.mtz_object().as_miller_arrays():
      labels=ma.info().labels
      ld = labels_decorator(amplitudes_label=labels[0], phases_label=labels[1])
      yet_another_dataset.add_miller_array(
        miller_array      = ma,
        column_root_label = labels[0],
        label_decorator   = ld)
  yet_another_dataset.mtz_object().write(file_name = file_name)

class map_manager(object):

  def __init__(self, refine_object):
    self.write_xplor = False
    self.write_mtz = False
    self.write_phs = False
    self.write_cif = False
    self.refine_object = refine_object
    self.check_map_parameters()
    self._callback = None
    self.map_params = refine_object.params.electron_density_maps
    self.setup_ncs_averaging()
    if(len(self.map_params.map_coefficients)>0 and
       refine_object.params.output.write_map_coefficients):
      for i_seq, mtz_mc in enumerate(self.map_params.map_coefficients):
        if(mtz_mc.map_type is not None and "mtz" in mtz_mc.format):
          self.write_mtz = True
          break
      for i_seq, mtz_mc in enumerate(self.map_params.map_coefficients):
        if(mtz_mc.map_type is not None and "phs" in mtz_mc.format):
          self.write_phs = True
          break
    if(len(self.map_params.map)>0 and refine_object.params.output.write_maps):
      for mtz_mc in self.map_params.map:
        if(mtz_mc.map_type is not None):
          self.write_xplor = True
          break
    if refine_object.params.output.write_reflection_cif_file and self.write_mtz:
      # XXX currently we can only write the CIF if we are also writing the mtz
      self.write_cif = True

  def setup_ncs_averaging (self) :
    ncs_average = False
    for params in self.map_params.map_coefficients+self.map_params.map :
      if (params.ncs_average) :
        ncs_average = True
        break
    if (ncs_average) :
      from mmtbx import map_tools
      ncs_object = self.refine_object.find_cartesian_ncs_groups(
        log=null_out(),
        hierarchy=self.refine_object.model.pdb_hierarchy(True))
      self._callback = map_tools.ncs_averager(
        ncs_object=ncs_object,
        params=self.refine_object.params.ncs.map_averaging,
        log=self.refine_object.log)
      print >> self.refine_object.log, \
        "NCS averaging will be applied to one or more maps."

  def check_map_parameters(self):
    mcp = self.refine_object.params.electron_density_maps.map_coefficients
    if(len(mcp)>0):
      all_labels = []
      for mp in mcp:
        all_labels.append(mp.mtz_label_amplitudes)
      new_params = []
      counter = 0
      for mp in mcp:
        if(all_labels.count(mp.mtz_label_amplitudes)>1):
          counter += 1
          mla = mp.mtz_label_amplitudes
          mlp = mp.mtz_label_phases
          cntr = "_%s"%str(counter)
          if(mla is None):
            mla = mp.map_type.replace("+","").replace("-","")
            mlp = "P"+mla
            cntr=""
          mp.mtz_label_amplitudes = mla + cntr
          mp.mtz_label_phases = mlp + cntr
        new_params.append(mp)
      self.refine_object.params.electron_density_maps.map_coefficients = new_params

  def write_mtz_file(self):
    if(self.refine_object.fmodels.fmodel_neutron() is not None):
      cmos = []
      for i_suf, suffix in enumerate(["_xray", "_neutron"]):
        mcp_dc = copy.deepcopy(self.map_params.map_coefficients)
        for i_seq, mcp in enumerate(mcp_dc):
          mcp_dc[i_seq].mtz_label_amplitudes = mcp.mtz_label_amplitudes+suffix
          mcp_dc[i_seq].mtz_label_phases     = mcp.mtz_label_phases+suffix
        if(suffix == "_xray"):
          fmodel = self.refine_object.fmodels.fmodel_xray()
        if(suffix == "_neutron"):
          fmodel = self.refine_object.fmodels.fmodel_neutron()
        if(i_suf == 0):
          cmo = mmtbx.maps.compute_map_coefficients(fmodel = fmodel,
            params = mcp_dc,
            log=self.refine_object.log)
        else:
          cmo = mmtbx.maps.compute_map_coefficients(fmodel = fmodel,
            params = mcp_dc, mtz_dataset = cmos[0].mtz_dataset,
            log=self.refine_object.log)
        cmos.append(mmtbx.maps.compute_map_coefficients(
          fmodel = fmodel, params = mcp_dc))
      del cmos
    else:
      cmo = mmtbx.maps.compute_map_coefficients(
        fmodel = self.refine_object.fmodels.fmodel_xray(),
        params = self.map_params.map_coefficients,
        post_processing_callback=self._callback,
        log=self.refine_object.log,
        pdb_hierarchy=self.refine_object.model.pdb_hierarchy())
    return cmo.mtz_dataset


  def write_files(self):
    if(not (self.write_xplor or self.write_mtz or self.write_phs)): return
    if(self.write_mtz):
      mtz_file_object = self.refine_object.open_output_file(
        output_type="MTZ file: data and map coefficients", ext=".mtz")
      if (not self.refine_object.params.main.wavelength in [None, Auto]) :
        self.refine_object.inputs.mtz_dataset_orig.set_wavelength(
          self.refine_object.params.main.wavelength)
      combine_data_and_map_coeffs(
        original_data_mtz_dataset=self.refine_object.inputs.mtz_dataset_orig,
        map_coeff_dataset=self.write_mtz_file(),
        fmodels=self.refine_object.fmodels,
        file_name=mtz_file_object.name,
        wavelength=self.refine_object.params.main.wavelength)
    if(self.write_phs):
      for i_seq, mtz_mc in enumerate(self.map_params.map_coefficients):
        if(mtz_mc.map_type is not None):
          phs_file_object = self.refine_object.open_output_file(
            output_type="PHS map coefficients", ext="_"+mtz_mc.map_type+".phs")
          coeffs = mmtbx.maps.map_coefficients_from_fmodel(
            fmodel = self.refine_object.fmodels.fmodel_xray(),
            params = mtz_mc)
          if (coeffs is None) : continue
          print >> self.refine_object.log, "Writing %s coefficients to PHS file:\n  %s" % (
            mtz_mc.map_type, phs_file_object.name)
          coeffs.as_phases_phs(out=phs_file_object)
          phs_file_object.close()
    if(self.write_xplor):
      map_type_list = []
      for mp in self.map_params.map:
        if(mp.map_type is not None):
          coeffs = mmtbx.maps.map_coefficients_from_fmodel(
            fmodel=self.refine_object.fmodels.fmodel_xray(),
            params=mp,
            post_processing_callback=self._callback)
          if (coeffs is None) : continue
          addl_ins=""
          if mp.format == "xplor" :
            format_ext = ".map"
          else :
            format_ext = ".ccp4"
          if(mp.fill_missing_f_obs):
            addl_ins += "_filled"
          if(self.refine_object.fmodels.fmodel_neutron() is not None):
            addl_ins += "_xray"
          basic_map_type = mp.map_type + addl_ins
          if (basic_map_type in map_type_list) :
            addl_ins += "-%d" % (map_type_list.count(basic_map_type) + 1)
          map_type_list.append(basic_map_type)
          set_map_file_name = (mp.file_name is None)
          if (set_map_file_name) :
            mp.file_name = self.refine_object.open_output_file(
              output_type = "XPLOR map",
              ext="_"+mp.map_type+addl_ins+format_ext).name
          mmtbx.maps.write_xplor_map_file(
            params = mp,
            coeffs = coeffs,
            xray_structure = self.refine_object.fmodels.fmodel_xray().xray_structure)
          if(self.refine_object.fmodels.fmodel_neutron() is not None):
            coeffs = mmtbx.maps.map_coefficients_from_fmodel(
              fmodel = self.refine_object.fmodels.fmodel_neutron(), params = mp)
            addl_ins=""
            if(mp.fill_missing_f_obs):
              addl_ins = "_filled"
            mp.file_name = self.refine_object.open_output_file(
              output_type = "XPLOR map",
              ext="_"+mp.map_type+addl_ins+"_neutron"+format_ext).name
            mmtbx.maps.write_xplor_map_file(
              params = mp,
              coeffs = coeffs,
              xray_structure = self.refine_object.fmodels.fmodel_neutron().xray_structure)
          if (set_map_file_name) :
            mp.file_name = None
    if self.write_cif:
      import iotbx.cif.model

      from iotbx.command_line import mtz_as_cif

      mtz_dataset_orig = self.refine_object.inputs.mtz_dataset_orig
      mtz_object_orig = mtz_dataset_orig.mtz_object()
      cif_blocks = mtz_as_cif.mtz_as_cif_blocks(
        mtz_object_orig, log=self.refine_object.log,
        test_flag_value=self.refine_object.params.input.xray_data.r_free_flags.test_flag_value).cif_blocks

      data_types = set(["xray"])
      if cif_blocks['neutron'] is not None:
        data_types.add("neutron")

      for data_type in data_types:
        cif_blocks[data_type].add_miller_array(
          array=getattr(self.refine_object.fmodels, "fmodel_%s" %data_type)()\
            .f_model_scaled_with_k1_composite_work_free(),
          column_names=("_refln.pdbx_F_calc_with_solvent",
                        "_refln.pdbx_phase_calc_with_solvent"))

      cif_model = iotbx.cif.model.cif()
      cif_model[self.refine_object.params.output.prefix] = cif_blocks["xray"].cif_block
      if self.refine_object.neutron_refinement:
        cif_model[self.refine_object.params.output.prefix+"_neutron"] = cif_blocks["neutron"].cif_block
      cif_file_object = self.refine_object.open_output_file(
        output_type="CIF file: data and map coefficients", ext=".reflections.cif")
      print >> self.refine_object.log, "Writing data and map coefficients to CIF file:\n  %s" % (
        cif_file_object.name)
      print >> cif_file_object, cif_model
      cif_file_object.close()


class extract_refinement_strategy_and_selections(object):
  def __init__(self, params, xray_structure, all_chain_proxies,
                     neutron_refinement, log):
    # ! Order of events below is important !
    print_statistics.make_header("Extract refinement strategy and selections",
                                                                     out = log)
    self.params = params.refine
    self.params_nm = params.nm
    self.neutron_refinement = neutron_refinement
    self.params_all = params
    self.xray_structure = xray_structure
    self.strategy = self.params.strategy
    self.all_chain_proxies = all_chain_proxies
    self.log = log
    self.options = ["individual_sites",
                    "individual_sites_real_space",
                    "rigid_body",
                    "individual_adp",
                    "group_adp",
                    "tls",
                    "nm",
                    "occupancies",
                    "group_anomalous",
                    "rosetta",
                    "den"]
    self.check_strategy()

  def check_strategy(self):
    for item in self.strategy:
      assert item in self.options
    if (not self.params_all.main.use_geometry_restraints) :
      if (self.params_all.target_weights.optimize_xyz_weight) :
        raise Sorry("Optimization of X-ray/stereochemistry weight is not "+
          "applicable when geometry restraints are turned off.")

  def get_anisotropic_flags(self):
    iso_sel = None
    aniso_sel = None
    if(self.params.adp.individual.isotropic is not None):
       iso_sel = utils.get_atom_selections(
         iselection        = False,
         all_chain_proxies = self.all_chain_proxies,
         selection_strings = self.params.adp.individual.isotropic,
         xray_structure    = self.xray_structure,
         parameter_name    = "refinement.refine.adp.individual.isotropic")[0]
       self.xray_structure.convert_to_isotropic(selection=iso_sel.iselection())
    if(self.params.adp.individual.anisotropic is not None):
       aniso_sel = utils.get_atom_selections(
          iselection        = False,
          all_chain_proxies = self.all_chain_proxies,
          selection_strings = self.params.adp.individual.anisotropic,
          xray_structure    = self.xray_structure,
          parameter_name = "refinement.refine.adp.individual.anisotropic")[0]
       self.xray_structure.convert_to_anisotropic(selection = aniso_sel)
    if(self.adp_nm is not None):
      uc = self.xray_structure.unit_cell()
      nm_sel_st = flex.size_t()
      for nm_sel in self.adp_nm:
         nm_sel_st.extend(nm_sel)
      modes = nmtools.generate_evec(selections = self.adp_nm,
                                    selections_1d = nm_sel_st,
                                    xray_structure = self.xray_structure,
                                    pdb_hierarchy = self.all_chain_proxies.pdb_hierarchy,
                                    filename = self.params_nm.evecin,
                                    n_modes = self.params_nm.n_modes)
      for nm_sel in self.adp_nm:
#        self.xray_structure.convert_to_anisotropic(selection = nm_sel)
        for item in nm_sel:
          sc = self.xray_structure.scatterers()[item]
          if(sc.u_iso == -1.0):
            assert sc.u_star != (-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)
            assert sc.flags.use_u_aniso()
            u_cart = adptbx.u_star_as_u_cart(uc, sc.u_star)
            # Subtract maximal possible amount from u_cart to leave it
            # positive definite and having anisotropy > 0.25
            # Ideally, do subtraction from diagonalized u_cart
            u_min = 0
            evalues = adptbx.eigensystem(u_cart).values()
            if(max(evalues) and min(evalues)/max(evalues) > 0.25):
              u_min = (min(evalues)-0.25*max(evalues))*0.75
              if(u_min<0): u_min=0
            u_cart = [u_cart[0]-u_min,
                      u_cart[1]-u_min,
                      u_cart[2]-u_min,u_cart[3],u_cart[4],u_cart[5]]
            sc.u_star = adptbx.u_cart_as_u_star(uc, u_cart)
            sc.u_iso = u_min
            sc.flags.set_use_u_iso(True)
          if(sc.u_star == (-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)):
            assert sc.u_iso != -1.0
            assert sc.flags.use_u_iso()
            sc.u_star = (0.0,0.0,0.0,0.0,0.0,0.0)
            sc.flags.set_use_u_aniso(True)

# convert all the atoms in nm refinement selections into anisotropic
    if(iso_sel is not None and aniso_sel is not None ):
       if(not (iso_sel&aniso_sel).all_eq(False)):
          raise Sorry("\nDuplicated selections in adp.individual.anisotropic "
                                              "and adp.individual.isotropic\n")

    if(self.adp_tls is not None or self.adp_group is not None):
       uc = self.xray_structure.unit_cell()
       is_adp_tls, is_adp_group = False, False
       if(self.adp_tls is not None):
          selections = self.adp_tls
          is_adp_tls = True
       else:
          selections = self.adp_group
          is_adp_group = True
       for sel_ in selections:
           for item in sel_:
               sc = self.xray_structure.scatterers()[item]
               if(sc.u_iso == -1.0):
                  assert sc.u_star != (-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)
                  assert sc.flags.use_u_aniso()
                  u_cart = adptbx.u_star_as_u_cart(uc, sc.u_star)
                  # Subtract maximal possible amount from u_cart to leave it
                  # positive definite and having anisotropy > 0.25
                  # Ideally, do subtraction from diagonalized u_cart
                  u_min = 0
                  evalues = adptbx.eigensystem(u_cart).values()
                  if(max(evalues) and min(evalues)/max(evalues) > 0.25):
                    u_min = (min(evalues)-0.25*max(evalues))*0.75
                    if(u_min<0): u_min=0
                  u_cart = [u_cart[0]-u_min,
                            u_cart[1]-u_min,
                            u_cart[2]-u_min,u_cart[3],u_cart[4],u_cart[5]]
                  sc.u_star = adptbx.u_cart_as_u_star(uc, u_cart)
                  sc.u_iso = u_min
                  sc.flags.set_use_u_iso(True)
               if(sc.u_star == (-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)):
                  assert sc.u_iso != -1.0
                  assert sc.flags.use_u_iso()
                  if(is_adp_tls):
                     sc.u_star = (0.0,0.0,0.0,0.0,0.0,0.0)
                     sc.flags.set_use_u_aniso(True)

  def refinement_flags(self):
    sites_individual           = None
    sites_torsion_angles       = None
    sites_rigid_body           = None
    adp_individual_iso         = None
    adp_individual_aniso       = None
    self.adp_group             = None
    adp_group_h                = None
    self.adp_tls               = None
    self.adp_nm                = None
    s_occupancies              = None
    hd_selection = self.xray_structure.hd_selection()
    if(self.xray_structure.hd_selection().count(True) > 0):
       adp_group_h = utils.get_atom_selections(
                                all_chain_proxies     = self.all_chain_proxies,
                                one_group_per_residue = True,
                                hydrogens_only        = True,
                                xray_structure        = self.xray_structure)
    if("individual_sites" in self.strategy):
      sites_individual = utils.get_atom_selections(
        iselection          = False,
        all_chain_proxies   = self.all_chain_proxies,
        selection_strings   = self.params.sites.individual,
        xray_structure      = self.xray_structure,
        one_selection_array = True,
        parameter_name      = "refinement.refine.sites.individual")
    if(self.params_all.main.simulated_annealing_torsion):
      sites_torsion_angles = utils.get_atom_selections(
        iselection          = False,
        all_chain_proxies   = self.all_chain_proxies,
        selection_strings   = self.params.sites.torsion_angles,
        xray_structure      = self.xray_structure,
        one_selection_array = True,
        parameter_name      = "refinement.refine.sites.torsion_angles")
    if("rigid_body" in self.strategy):
      sites_rigid_body = utils.get_atom_selections(
        all_chain_proxies = self.all_chain_proxies,
        selection_strings = self.params.sites.rigid_body,
        xray_structure    = self.xray_structure,
        parameter_name = "refinement.refine.sites.rigid_body")
      self.rigid_body_selections_and_atoms_at_special_positions(
                                    rigid_body_selections = sites_rigid_body)
    if("individual_adp" in self.strategy):
       if(self.params.adp.individual.isotropic is not None):
         adp_individual_iso = utils.get_atom_selections(
           iselection          = False,
           all_chain_proxies   = self.all_chain_proxies,
           selection_strings   = self.params.adp.individual.isotropic,
           allow_empty_selection = None,
           xray_structure      = self.xray_structure,
           one_selection_array = True,
           parameter_name      = "refinement.refine.adp.individual.isotropic")
         if (adp_individual_iso is None):
           self.params.adp.individual.isotropic = None
       if(self.params.adp.individual.anisotropic is not None):
         adp_individual_aniso = utils.get_atom_selections(
           iselection          = False,
           all_chain_proxies   = self.all_chain_proxies,
           selection_strings   = self.params.adp.individual.anisotropic,
           allow_empty_selection = None,
           xray_structure      = self.xray_structure,
           one_selection_array = True,
           parameter_name      = "refinement.refine.adp.individual.anisotropic")
         if (adp_individual_aniso is None):
           self.params.adp.individual.anisotropic = None
       if(self.params.adp.individual.isotropic is None and
          self.params.adp.individual.anisotropic is None):
         adp_individual_iso = self.xray_structure.use_u_iso()
         adp_individual_aniso = self.xray_structure.use_u_aniso()
    if("group_adp" in self.strategy or "tls" in self.strategy or "nm" in self.strategy):
      if(self.params.adp.group_adp_refinement_mode == "group_selection"):
        if(len(self.params.adp.group) == 0):
          raise Sorry("No selections for group B-factor refinement is given.")
        self.adp_group = utils.get_atom_selections(
          all_chain_proxies     = self.all_chain_proxies,
          selection_strings     = self.params.adp.group,
          one_group_per_residue = False,
          xray_structure        = self.xray_structure,
          parameter_name        = "refinement.refine.adp.group")
      else:
        if("group_adp" in self.strategy and len(self.params.adp.group) > 0):
          msg1 = "Selections for group B-factor refinement are given, "
          msg2 = "but no refinement requested."
          raise Sorry(msg1+msg2)
      if(self.params.adp.group_adp_refinement_mode ==
         "one_adp_group_per_residue"):
        self.adp_group = utils.get_atom_selections(
          all_chain_proxies     = self.all_chain_proxies,
          selection_strings     = None,
          one_group_per_residue = True,
          xray_structure        = self.xray_structure)
      elif(self.params.adp.group_adp_refinement_mode ==
         "two_adp_groups_per_residue"):
        ogpr = utils.get_atom_selections(
          all_chain_proxies     = self.all_chain_proxies,
          selection_strings     = self.params.adp.group,
          one_group_per_residue = True,
          xray_structure        = self.xray_structure,
          parameter_name        = "refinement.refine.adp.group")
        sc = utils.get_atom_selections(
          iselection            = True,
          one_selection_array   = True,
          all_chain_proxies     = self.all_chain_proxies,
          selection_strings     = ["sidechain"],
          one_group_per_residue = False,
          xray_structure        = self.xray_structure)
        bb = utils.get_atom_selections(
          iselection            = True,
          one_selection_array   = True,
          all_chain_proxies     = self.all_chain_proxies,
          selection_strings     = ["backbone"],
          one_group_per_residue = False,
          xray_structure        = self.xray_structure)
        self.adp_group = []
        for g in ogpr:
          s1 = g.intersection(sc)
          s2 = g.intersection(bb)
          s3 = s1.intersection(s2)
          assert s3.size() == 0
          if(s1.size() > 0):
            self.adp_group.append( s1 )
          if(s2.size() > 0):
            self.adp_group.append( s2 )
    if("tls" in self.strategy or "nm" in self.strategy):
      self.adp_tls = utils.get_atom_selections(
        all_chain_proxies = self.all_chain_proxies,
        selection_strings = self.params.adp.tls,
        xray_structure    = self.xray_structure,
        parameter_name    = "refinement.refine.adp.tls")
      for sa, ss in zip(self.adp_tls, self.params.adp.tls):
        if(sa.size()==0):
          raise Sorry("Empty TLS selection: %s"%str(ss))
        if(sa.size()<self.params_all.tls.min_tls_group_size):
          small_tls_group_message = """
TLS group selected with '%s' contains %d atoms (less than
refinement.tls.min_tls_group_size=%d. Please either review TLS groups selection
or decrease the refinement.tls.min_tls_group_size parameter.
"""
          raise Sorry(small_tls_group_message%(ss,sa.size(),
            self.params_all.tls.min_tls_group_size))
      #
      # XXX Force atoms in TLS group to NOT refine as individual ANISOtropic.
      # Temporary fix start
      if(adp_individual_aniso is not None):
        if(adp_individual_iso is None):
          adp_individual_iso = flex.bool(adp_individual_aniso.size(), False)
        for tls_i in self.adp_tls:
          if(adp_individual_aniso.select(tls_i).count(True) > 0):
            adp_individual_aniso.set_selected(tls_i, False)
            adp_individual_iso.set_selected(tls_i, True)
      # Temporary fix end
    if("nm" in self.strategy):
      self.adp_nm = utils.get_atom_selections(
        all_chain_proxies = self.all_chain_proxies,
        selection_strings = self.params.adp.nm,
        xray_structure    = self.xray_structure,
        parameter_name    = "refinement.refine.adp.nm")
      for sa, ss in zip(self.adp_nm, self.params.adp.nm):
        if(sa.size()==0):
          raise Sorry("Empty NM selection: %s"%str(ss))
        if(sa.size()<self.params_all.nm.min_nm_group_size):
          small_nm_group_message = """
NM group selected with '%s' contains %d atoms (less than
refinement.nm.min_nm_group_size=%d. Please either review NM groups selection
or decrease the refinement.nm.min_nm_group_size parameter.
"""
          raise Sorry(small_nm_group_message%(ss,sa.size(),
            self.params_all.nm.min_tls_group_size))
      #
      # XXX Force atoms in NM group to NOT refine as individual ANISOtropic.
      # Temporary fix start
      if(adp_individual_aniso is not None):
        if(adp_individual_iso is None):
          adp_individual_iso = flex.bool(adp_individual_aniso.size(), False)
        for nm_i in self.adp_nm:
          if(adp_individual_aniso.select(nm_i).count(True) > 0):
            adp_individual_aniso.set_selected(nm_i, False)
            adp_individual_iso.set_selected(nm_i, False)
      # Temporary fix end

    if("occupancies" in self.strategy):
      s_occupancies = utils.occupancy_selections(
        all_chain_proxies   = self.all_chain_proxies,
        xray_structure      = self.xray_structure,
        add_water           = self.params_all.ordered_solvent.refine_occupancies,
        other_constrained_groups = self.params.occupancies.constrained_group,
        other_individual_selection_strings = self.params.occupancies.individual,
        remove_selection    = self.params.occupancies.remove_selection,
        as_flex_arrays      = True)
    if("tls" in self.strategy and adp_group_h is not None): #XXX EXPERIMENTAL
    #if("tls" in self.strategy and adp_group_h is not None and not self.neutron_refinement): #XXX EXPERIMENTAL
      self.adp_tls = self.exclude_hydrogens(hd_selection = hd_selection,
                                            other        = self.adp_tls)
    if("group_adp" in self.strategy and adp_group_h is not None):
      self.adp_group = self.exclude_hydrogens(hd_selection = hd_selection,
                                              other        = self.adp_group)
    if("nm" in self.strategy and adp_group_h is not None): #XXX EXPERIMENTAL
    #if("tls" in self.strategy and adp_group_h is not None and not self.neutron_refinement): #XXX EXPERIMENTAL
      self.adp_nm = self.exclude_hydrogens(hd_selection = hd_selection,
                                           other        = self.adp_nm)
    self.get_anisotropic_flags()
    from mmtbx.refinement import refinement_flags as reffl
    refinement_flags_obj = reffl.manager(
      sites_individual       = sites_individual,
      sites_torsion_angles   = sites_torsion_angles,
      sites_rigid_body       = sites_rigid_body,
      adp_individual_iso     = adp_individual_iso,
      adp_individual_aniso   = adp_individual_aniso,
      adp_group              = self.adp_group,
      group_h                = adp_group_h,
      adp_tls                = self.adp_tls,
      adp_nm                 = self.adp_nm,
      s_occupancies          = s_occupancies)
    for key in refinement_flags_obj.__dict__.keys():
      for item in self.strategy:
        if(key == item):
          refinement_flags_obj.__dict__[key]=True
    if(self.params_all.main.simulated_annealing_torsion):
      refinement_flags_obj.torsion_angles=True
    refinement_flags_obj.show(log = self.log)
    refinement_flags_obj.check_all()
    return refinement_flags_obj

  def exclude_hydrogens(self, hd_selection, other):
    size = hd_selection.size()
    result = []
    for sel in other:
      sel_b = flex.bool(size, sel)
      sel_b.set_selected(hd_selection, False)
      result.append(sel_b.iselection())
    return result

  def rigid_body_selections_and_atoms_at_special_positions(self,
                                                        rigid_body_selections):
    special_position_indices = self.xray_structure.special_position_indices()
    atoms_at_sp_and_in_rigid_body_selections = []
    if(special_position_indices.size() > 0):
       for rb_sel in rigid_body_selections:
           for i_seq in special_position_indices:
               if i_seq in rb_sel:
                  atoms = self.all_chain_proxies.pdb_atoms
                  line = "    %s" % atoms[i_seq].quote()
                  atoms_at_sp_and_in_rigid_body_selections.append(line)
    if(len(atoms_at_sp_and_in_rigid_body_selections) > 0):
       print >> self.log, \
            "\nSome atoms selected for rigid body refinement are in "+\
            "special positions (listed\n below). Please remove these "+\
            "atoms from the rigid body selection and run\n refinement again.\n"
       for item in atoms_at_sp_and_in_rigid_body_selections:
           print >> self.log, item
       print >> self.log
       raise Sorry("\nAtoms at special positions are within rigid groups.\n")

class extract_experimental_phases(object):
  def __init__(self, experimental_phases, f_obs, log):
    if(experimental_phases is not None):
       mmtbx.refinement.print_statistics.make_header(
                                              "Experimental phases", out = log)
       experimental_phases.show_comprehensive_summary(f = log)
       print >> log
       if(not f_obs.anomalous_flag()):
          if(experimental_phases.anomalous_flag()):
             print >> log, "Averaging Bijvoet mates of experimental phases."
             experimental_phases = experimental_phases.average_bijvoet_mates()
             print >> log
       elif(not experimental_phases.anomalous_flag()):
          print >> log, "Generating Bijvoet mates of experimental phases."
          experimental_phases = experimental_phases.generate_bijvoet_mates()
          print >> log
       print >> log, "Average figures of merit by resolution:"
       figures_of_merit = abs(experimental_phases.phase_integrals())
       figures_of_merit.setup_binner(n_bins=10)
       legend_len = figures_of_merit.mean(use_binning=True).show(
                     data_fmt="%6.3f", show_unused=False, f = log, prefix="  ")
       print >> log, " ", ("%%%ds"%legend_len)%"overall", \
                                                "%6.3f"%figures_of_merit.mean()
       print >> log, "  Note: Figures of merit are determined by integration of"
       print >> log, "        Hendrickson-Lattman coefficients."
       print >> log
       experimental_phases = experimental_phases.map_to_asu().matching_set(
                                      other = f_obs, data_substitute=(0,0,0,0))
       print >> log, \
          "Number and fraction of available experimental phases by resolution:"
       experimental_phases.setup_binner(n_bins=10)
       experimental_phases.count_and_fraction_in_bins(
         data_value_to_count = (0,0,0,0),
         count_not_equal= True).show(show_unused = False, f = log, prefix="  ")
       print >> log
    self.experimental_phases_ = experimental_phases

  def experimental_phases(self):
    return self.experimental_phases_

def check_input_params(params, f_obs, xray_structure):
  # IAS refinement
  if(params.main.ias):
    if(params.main.ordered_solvent or
       params.target_weights.wc != 0 or
       f_obs.d_min()>=1.0):
      msg = """
Cannot do IAS refinement. To proceed make sure that:
  refinement.main.ordered_solvent=False
  refinement.target_weights.wc=0
  and data resolution is higher than 1.0A (at least).
"""
      raise Sorry(msg)

class modify_f_obs(object):

  def __init__(self, params, f_obs, r_free_flags, log):
    p = params.modify_f_obs
    self.f_obs = f_obs
    self.r_free_flags = r_free_flags
    self.log = log
    self.f_obs_r = None
    if(p.remove is None): return
    elif(p.remove == "random"):
      sel = flex.random_double(size=f_obs.data().size())<=(1.-p.remove_fraction)
      #
      self.f_obs_r = self.f_obs.select(selection = ~sel)
      #
      self.f_obs = self.f_obs.select(selection = sel)
      self.r_free_flags = self.r_free_flags.select(selection = sel)
    elif(p.remove == "strong"):
      self.selection_helper(reverse         = False,
                            remove_fraction = p.remove_fraction,
                            by_value        = "data")
    elif(p.remove == "weak"):
      self.selection_helper(reverse         = True,
                            remove_fraction = p.remove_fraction,
                            by_value        = "data")
    elif(p.remove == "strong_and_weak"):
      self.selection_helper(reverse         = False,
                            remove_fraction = p.remove_fraction/2.,
                            by_value        = "data")
      self.selection_helper(reverse         = True,
                            remove_fraction = p.remove_fraction/2.,
                            by_value        = "data")
    elif(p.remove == "low"):
      self.selection_helper(reverse         = False,
                            remove_fraction = p.remove_fraction,
                            by_value        = "resolution")
    elif(p.remove == "other"):
      R=0.001
      while 1:
        sel = flex.bool(self.f_obs.indices().size(), True)
        for i_seq, mi in enumerate(self.f_obs.indices()):
          if(mi[0]**2+mi[1]**2 < R*mi[2]**2):
            sel[i_seq] = False
        if(sel.count(False)/sel.size() >= p.remove_fraction): break
        R += 0.001
      print "ELLIPTICAL REMOVAL:", sel.count(False)/sel.size(), p.remove_fraction

      #
      self.f_obs_r = self.f_obs.select(selection = ~sel)
      #
      self.f_obs = self.f_obs.select(selection = sel)
      self.r_free_flags = self.r_free_flags.select(selection = sel)
    else:
      raise RuntimeError
    #
    d = self.f_obs.d_spacings().data()
    fo = self.f_obs.data()
    d_r = self.f_obs_r.d_spacings().data()
    fo_r = self.f_obs_r.data()
    #for i in xrange(fo.size()):
    #  if(i < d_r.size()):
    #    print >> self.log, "%8.4f %15.3f %8.4f %15.3f" % (d[i], fo[i], d_r[i], fo_r[i])
    #  else:
    #    print >> self.log, "%8.4f %15.3f" % (d[i], fo[i])
    #
    mtz_history_buffer = flex.std_string()
    file_name = params.output.prefix+"_data_manipulated.mtz"
    mtz_dataset = self.f_obs.as_mtz_dataset(
      column_root_label="F-obs",
      wavelength=params.main.wavelength)
    mtz_dataset.add_miller_array(
      miller_array = self.r_free_flags, column_root_label="R-free-flags")
    ha = mtz_history_buffer.append
    ha(date_and_time())
    ha("file name: %s" % os.path.basename(file_name))
    ha("directory: %s" % os.path.dirname(file_name))
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.add_history(lines=mtz_history_buffer)
    mtz_object.write(file_name = file_name)

  def selection_helper(self, reverse, remove_fraction, by_value):
    sel = self.f_obs.sort_permutation(by_value = by_value, reverse = reverse)
    self.f_obs = self.f_obs.select(selection = sel)
    self.r_free_flags = self.r_free_flags.select(selection = sel)
    n_remove = int(self.f_obs.data().size()*remove_fraction)
    sel = flex.bool(n_remove, False).concatenate(
      flex.bool(self.f_obs.data().size()-n_remove, True))
    #
    if(self.f_obs_r is None):
      self.f_obs_r = self.f_obs.select(selection = ~sel)
    else:
       self.f_obs_r = self.f_obs_r.concatenate(self.f_obs.select(selection = ~sel))
    #
    self.f_obs = self.f_obs.select(selection = sel)
    self.r_free_flags = self.r_free_flags.select(selection = sel)
