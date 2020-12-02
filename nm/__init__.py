from __future__ import division
import cctbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("mmtbx_nm_ext")
from mmtbx_nm_ext import *
