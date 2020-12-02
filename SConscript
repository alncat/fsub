import libtbx.load_env
import os
Import("env_etc")

env_etc.mmtbx_dist = libtbx.env.dist_path("mmtbx")
env_etc.mmtbx_include = os.path.dirname(env_etc.mmtbx_dist)

env_etc.mmtbx_common_includes = [
  env_etc.libtbx_include,
  env_etc.mmtbx_include,
  env_etc.cctbx_include,
  env_etc.scitbx_include,
  env_etc.boost_include,
]

if (not env_etc.no_boost_python):
  Import("env_scitbx_boost_python_ext")
  env_bpl = env_scitbx_boost_python_ext.Clone()
  env_etc.include_registry.append(
    env=env_bpl,
    paths=[env_etc.mmtbx_include] + env_etc.scitbx_common_includes)
  env_bpl.SharedLibrary(
    target="#lib/mmtbx_alignment_ext",
    source=["alignment_ext.cpp"])

SConscript("masks/SConscript")
SConscript("ncs/SConscript")
SConscript("max_lik/SConscript")
SConscript("dynamics/SConscript")
SConscript("bulk_solvent/SConscript")
SConscript("tls/SConscript")
SConscript("nm/SConscript")
SConscript("scaling/SConscript")
SConscript("f_model/SConscript")
SConscript("utils/SConscript")
SConscript("rsr/SConscript")
SConscript("geometry_restraints/SConscript")
SConscript("secondary_structure/SConscript")
SConscript("den/SConscript")
SConscript("cablam/SConscript")
SConscript("geometry/SConscript")
