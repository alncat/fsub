from __future__ import division
from mmtbx.tls import tools
import mmtbx.f_model
import mmtbx.model
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import libtbx.load_env
import iotbx.pdb
import random
import sys, os

def run(args):
#    if (len(args) == 0 ):
#        raise RuntimeError("Please specify one or more pdb file names.")
   mon_lib_srv = monomer_library.server.server()
   ener_lib = monomer_library.server.ener_lib()
#   pdb_file = libtbx.env.find_in_repositories(
#          relative_path="phenix_regression/pdb/1A32.pdb",
#          test=os.path.isfile)
   pdb_file = "./4KSC.pdb"
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
   xray_structure.scattering_type_registry(table = "wk1995")
   xray_structure.convert_to_isotropic()
   u_iso_start = xray_structure.extract_u_iso_or_u_equiv()
   xray_structure.convert_to_anisotropic()
   atoms = pdb_hierarchy.atoms()
   atomsxyz = atoms.extract_xyz()
   sites_cart = xray_structure.sites_cart()
   for i in range(len(sites_cart)):
      print atoms[i].parent().resname, atoms[i].name
      for j in range(3):
           if sites_cart[i][j] - atomsxyz[i][j] <= 1.e-9:
               pass
           else:
               print i, j
               raise RuntimeError("don't match!")
#   for file_name in args:
#       pdb_obj = iotbx.pdb.hierarchy.input(file_name=file_name)
#        pdb_obj.hierarchy.overall_counts().show()
#        for model in pdb_obj.hierarchy.models():
#            for chain in model.chains():
#                for rg in chain.residue_groups():
#                    print 'resid: "%s"' % rg.resid()
#                    for ag in rg.atom_groups():
#                        print ' altloc: "%s", resname: "%s"' % (ag.altloc, ag.resname)
#                        for atom in ag.atoms():
#                            print '    ', atom.name
        

if (__name__ == "__main__"):
    run(sys.argv[1:])
