from __future__ import division
from iotbx import pdb
from cctbx.array_family import flex
import math, random
import sys, time
from scitbx import lbfgs
from mmtbx_nm_ext import *
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
from cctbx import adptbx
from cctbx import xray
from libtbx.utils import user_plus_sys_time, Sorry
from libtbx.str_utils import line_breaker
from libtbx import group_args
from stdlib import math
import mmtbx.utils
import iotbx
import os
import subprocess

time_generate_evec = 0.0
time_convert_modes = 0.0
time_nm_from_uanisos = 0.0
time_update_xray_structure_with_nm = 0.0
time_u_cart_from_nm = 0.0
time_nm_total = 0.0
time_make_nm_compatible_with_u_positive_definite = 0.0
time_split_u = 0.0
def show_time(out = None):
    if(out is None): out = sys.stdout
    total = time_generate_evec +\
            time_convert_modes +\
            time_nm_from_uanisos +\
            time_update_xray_structure_with_nm +\
            time_make_nm_compatible_with_u_positive_definite +\
            time_nm_total +\
            time_split_u
    if(total > 0.01):
        print >> out, "NM refinement:"
        print >> out, " time_generate_evec                               = %-7.2f"% time_generate_evec
        print >> out, " time_convert_modes                               = %-7.2f"% time_convert_modes
        print >> out, " time_nm_from_uanisos                             = %-7.2f"% time_nm_from_uanisos
        print >> out, " time_update_xray_structure_with_nm               = %-7.2f"% time_update_xray_structure_with_nm
        print >> out, " time_u_cart_from_nm                              = %-7.2f"% time_u_cart_from_nm
        print >> out, " time_make_nm_compatible_with_u_positive_definite = %-7.2f"% time_make_nm_compatible_with_u_positive_definite
        print >> out, " time_nm_total                                    = %-7.2f"% time_nm_total
        print >> out, " time_split_u                                     = $-7.2f"% time_split_u
    return total

def split_u(u_eq, xray_structure, tls_selections, offset):
    global time_split_u
    timer = user_plus_sys_time()
    uc = xray_structure.unit_cell()
    for tls_selection in tls_selections:
        if(offset):
            offset_ = adptbx.b_as_u(5.0)
        else: offset_ = 0.0
        for i_seq in tls_selection:
            sc = xray_structure.scatterers()[i_seq]
            u_nm = sc.u_star
            u_nm_iso = adptbx.u_star_as_u_iso(uc, u_nm)
            u_iso_new = u_eq[i_seq] - u_nm_iso
            if u_iso_new < 0.0:
                u_iso_new = 0.0
            assert u_iso_new >= 0.0
            sc.u_iso = u_iso_new
            assert sc.flags.use_u_aniso()
            assert sc.flags.use_u_iso()
    u_iso = xray_structure.scatterers().extract_u_iso().select(
          xray_structure.use_u_iso())
    assert (u_iso < 0.0).count(True) == 0
    time_split_u += timer.elapsed()

class nm_groups(object):
    def __init__(self, xs = None, selection_strings = None):
        self.xs, self.selection_strings = xs, selection_strings

class show_nm(object):
    def __init__(self, xs, n_modes, text="", zero_mode_flag=True, out=None):
        if(out is None): out = sys.stdout
        counter = 0
        print >> out, text
        for item in xs:
            counter += 1
            print >> out, "NM group number %d: " % counter
            nt = 0
            if zero_mode_flag == True:
                for i in range(6):
                    for j in range(i+1):
                        print >> out, '{0:6.2f}'.format(item[nt]),
                        nt += 1
                    print >> out, ""
                for i in range(6, n_modes):
                    for j in range(i+1):
                        if j < 6:
                            print >> out, " "*6,
                        else:
                            print >> out, '{0:6.2f}'.format(item[nt]),
                            nt += 1
                    print >> out, ""
            else:
#		print len(item)
                for i in range(n_modes):
                    for j in range(i+1):
                        print >> out, '{0:6.2f}'.format(item[nt]),
                        nt += 1
                    print >> out, ""

def read_nmval_file(evalin = "./eigenvalues.dat",
                    n_modes = 20,
                    zero_mode_input_flag = True,
                    zero_mode_flag = False):
    assert os.path.isfile(evalin) is True, "Cannot find: %s" % evalin
    nmval = []
    if zero_mode_flag and not zero_mode_input_flag:
        nmode_start = 6
    else:
        nmode_start = 0
    print "reading eigenvalues..."
    for i in range(nmode_start):
        nmval.append(0.0)
    i = 0
    with open(evalin, 'r') as file_eval:
        for line in file_eval:
            i += 1
            if i <= n_modes:
                nmval.append(float(line.split()[1]))
            else:
                break
    if zero_mode_flag:
        for i in range(6):
            if abs(nmval[i]) > 1.E-4:
                print "check modes: zero_modes may not be provided"
                nmval[i] = 0.0
    return nmval


def generate_evec(selections,
		  selections_1d,
                  xray_structure,
                  pdb_hierarchy,
                  filename,
                  n_modes,
                  zero_mode_input_flag = False,
                  zero_mode_flag = True):
    global time_generate_evec
    t1 = time.time()
    atoms = pdb_hierarchy.atoms()
    print "reading eigenvectors from", filename
    nm_init_manager = nm_init(filename = filename,
                            n_modes = n_modes,
                            atoms = atoms,
                            zero_mode_input_flag = zero_mode_input_flag,
                            zero_mode_flag = zero_mode_flag)
    modes = []
    for i in range(n_modes):
        modes.append(nm_init_manager.return_modes(i))
    if zero_mode_flag == True and zero_mode_input_flag == False:
        count = 0
        for selection in selections:
#	    selection is size_t array
# 	    print len(selections)
#	    for selection_ in selection:
#		print selection_
	    bool_selection = flex.bool(xray_structure.sites_cart().size(), False)
	    bool_selection.set_selected(selection, True)
            sites_cart_selected = xray_structure.sites_cart().select(selection)
            atomic_weights_selected = xray_structure.atomic_weights().select(selection)
            nm_init_manager.gen_zero_modes(sites_cart_selected, atomic_weights_selected, bool_selection)
            padd = 6*count
            for i in range(6):
                selected_zero_modes = nm_init_manager.return_zero_modes(i+padd)
                modes[i].set_selected(selection, selected_zero_modes)
            count += 1
	new_sel_st = flex.size_t()
        for i in xrange(count):
            new_bool_selection = nm_init_manager.return_new_selection(i)
	    new_selection = []
            for selection_ in selections[i]:
	        if new_bool_selection[selection_] != False:
		    new_selection.append(selection_)
	    new_selection = flex.size_t(new_selection)
	    new_sel_st.extend(new_selection)
            selections[i] = new_selection
        selections_1d = new_sel_st
#    scatterers = xray_structure.scatterers()
#    for selection_ in selections:
#        for nm_selection__ in selection_:
#            scatterers[nm_selection__].flags.set_use_u_aniso(True)

# selection will be modified in this step, we will only keep those atoms with normal modes being assigned.
    nm_init_manager.normalize(n_modes)
    t2 = time.time()
    time_generate_evec += (t2 - t1)
    return modes

def selected_modes_to_1D(modes,
                        n_modes,
                        selection):
    global time_convert_modes
    t1 = time.time()
    assert len(modes) == n_modes
    for i in range(n_modes):
        assert len(modes[i]) == len(modes[0])
    len_selected_modes = len(modes[0].select(selection))
#    modes1d = flex.vec3_double(len_selected_modes*n_modes, [0,0,0])
    modes1d = []
    for i in range(n_modes):
        modes_i_selected = modes[i].select(selection)
        normalize_mode(modes_i_selected)
#        for j in range(len_selected_modes):
#            modes1d.append(modes_i_selected[j])
        modes1d.extend(modes_i_selected)
    modes1d = flex.vec3_double(modes1d)
    t2 = time.time()
    time_convert_modes += (t2 - t1)
    return modes1d

def init_nm_para(nmval, n_modes, zero_mode_flag = True):
    freq = []
    for i in range(n_modes):
        freq.append(0.0)
    if zero_mode_flag:
        nstart = 6
    else:
        nstart = 0
    for i in range(nstart, n_modes):
        if nmval[i] < 0.0:
            print "warning: non-positive modes"
            nmval[i] = abs(nmval[i])
        freq[i] = max(math.sqrt(nmval[i]), 1.e-2)

    ave_freq = 0.0
    for i in range(nstart, n_modes):
        ave_freq += freq[i]
    if nstart != n_modes:
        ave_freq = ave_freq/(n_modes - nstart)
    s = flex.double(n_modes*n_modes, 0.0)
    if zero_mode_flag:
        for i in range(6, n_modes):
            s[i+i*n_modes] = 1.0/freq[i]
        if n_modes < 6:
            mag = 1.0
        else:
            mag = .3/ave_freq
        for i in range(6):
            for j in range(i+1):
                s[j+i*n_modes] = mag*random.uniform(-1, 1)
        for i in range(7, n_modes):
            for j in range(6, i):
                s[j+i*n_modes] = mag*random.uniform(-1, 1)
    else:
        for i in range(n_modes):
            s[i+i*n_modes] = 1.0/freq[i]
        mag = .3/ave_freq
        for i in range(1, n_modes):
            for j in range(i):
                s[j+i*n_modes] = mag*random.uniform(-1, 1)
    if zero_mode_flag:
        n_nmpars = 21 + (n_modes - 5)*(n_modes - 6)/2
    else:
        n_nmpars = (n_modes + 1)*n_modes/2
    x = flex.double(int(n_nmpars), 0.0)
    s2x(s = s, x = x, n_modes = n_modes, zero_mode_flag = zero_mode_flag)
    return x

def s2x(s, x, n_modes, zero_mode_flag = True):
    nt = 0
    if zero_mode_flag:
        for i in range(6):
            for j in range(i+1):
                x[nt] = s[j + i*n_modes]
                nt += 1
        for i in range(6, n_modes):
            for j in range(6, i+1):
                x[nt] = s[j + i*n_modes]
                nt += 1
    else:
        for i in range(n_modes):
            for j in range(i+1):
                x[nt] = s[j + i*n_modes]
                nt += 1

class nm_from_uaniso_minimizer(object):
    def __init__(self,
                 uaniso,
                 x_initial,
                 adp_nma,
#                 weights,
                 n_modes,
                 zero_mode_flag,
                 max_iterations):
        adopt_init_args(self, locals())
#        assert self.uaniso.size() == self.weights.size()
        self.x     = self.x_initial
        self.x_min = self.x_initial
#        self.n = self.x.size()
        t1 = time.time()
        self.minimizer = lbfgs.run(
                                  target_evaluator = self,
                                  termination_params = lbfgs.termination_parameters(
                                      max_iterations = max_iterations,
                                      max_calls = int(max_iterations*1.5)),
                                  exception_handling_params =
                                  lbfgs.exception_handling_parameters(
                                      ignore_line_search_failed_step_at_lower_bound = True,
                                      ignore_line_search_failed_step_at_upper_bound = True,
                                      ignore_line_search_failed_maxfev              = True)
                                  )
        self.compute_functional_and_gradients()
        self.x_min = self.x
        del self.x
        t2 = time.time()
#        print t2 - t1

    def compute_functional_and_gradients(self):
        manager = nm_from_uaniso_target_and_grads(self.x,
#                                                 self.weights,
                                                  self.adp_nma,
                                                  self.uaniso,
                                                  self.n_modes,
                                                  self.zero_mode_flag)
        self.f = manager.target()
        self.g = manager.grad_nm()
        return self.f, self.g

def nm_from_uanisos(xray_structure,
                    selections,
                    modes,
                    xs_initial,
                    n_modes,
                    adp_nmas               = None,
                    number_of_macro_cycles = 1,
                    max_iterations         = 100,
                    zero_mode_flag         = True,
                    verbose                = -1,
                    out                    = None):
    global time_nm_from_uanisos
    t1 = time.time()
    if(out is None): out = sys.stdout
    if(verbose > 0):
        show_nm(xs = xs_initial,
                n_modes = n_modes,
                text = "NM from ADP: start NM values",
                zero_mode_flag = zero_mode_flag,
                out = out)
    u_cart = xray_structure.scatterers().extract_u_cart(xray_structure.unit_cell())
    xs_min = []
    group_counter = 0
    if( modes is not None and adp_nmas is None ):
        for x_initial, selection in zip(xs_initial, selections):
            group_counter += 1
            stop_flag = 0
            target_stop = -1.0
            weights_selected = xray_structure.atomic_weights().select(selection)
            u_cart_selected = u_cart.select(selection)
            modes1d_selected = selected_modes_to_1D(modes = modes, n_modes = n_modes,
                                                    selection = selection)
            assert len(modes1d_selected)/n_modes == len(weights_selected)
            adp_nma_selected = init_nm_adp(modes = modes1d_selected,
                                  weights = weights_selected,
                                  n_modes = n_modes,
                                  zero_mode_flag = zero_mode_flag)
            for i in range(1, number_of_macro_cycles+1):
                target_start = target_stop
                minimized = nm_from_uaniso_minimizer(uaniso = u_cart_selected,
                                                     x_initial = x_initial,
                                                     adp_nma = adp_nma_selected,
    #                                                 weights = weights_selected,
                                                     n_modes = n_modes,
                                                     zero_mode_flag = zero_mode_flag,
                                                     max_iterations = max_iterations)
                x_initial = minimized.x_min
            if(verbose > 0):
                print >> out, "NM group %d: minimized target = " %(group_counter), minimized.f
            x_min_ = minimized.x_min
            xs_min.append(x_min_)
    elif adp_nmas is not None and modes is None:
        for x_initial, selection, adp_nma_selected in zip(xs_initial, selections, adp_nmas):
            u_cart_selected = u_cart.select(selection)
            group_counter += 1
            stop_flag = 0
            target_stop = -1.0
            for i in range(1, number_of_macro_cycles+1):
                target_start = target_stop
                minimized = nm_from_uaniso_minimizer(uaniso = u_cart_selected,
                                                     x_initial = x_initial,
                                                     adp_nma = adp_nma_selected,
    #                                                 weights = weights_selected,
                                                     n_modes = n_modes,
                                                     zero_mode_flag = zero_mode_flag,
                                                     max_iterations = max_iterations)
                x_initial = minimized.x_min
            if(verbose > 0):
                print >> out, "NM group %d: minimized target = " %(group_counter), minimized.f
                print >> out, "Macrocycle %d finished!" %i
            x_min_ = minimized.x_min
            xs_min.append(x_min_)
    if(verbose > 0):
        show_nm(xs = xs_min,
                n_modes = n_modes,
                text = "NM from ADP: final NM values",
                zero_mode_flag = zero_mode_flag,
                out = out)
    t2 = time.time()
    time_nm_from_uanisos += (t2 - t1)
    return xs_min

class nm_xray_grads(object):

    def __init__(self, gradients, adp_nmas, selections, xs, n_modes, zero_mode_flag):
        self.grad_nm = []
#        d_target_d_uaniso = target_result.gradients_wrt_atomic_parameters(
#            u_aniso=True)
        for adp_nma, sel, x in zip(adp_nmas, selections, xs):
            d_target_d_nm_manager = d_target_d_nm(
                adp_nma = adp_nma,
                d_target_d_uaniso = gradients.select(sel),
                x       = x,
                n_modes = n_modes,
                zero_mode_flag = zero_mode_flag)
            self.grad_nm.append(list(d_target_d_nm_manager.grad_nm()))

class nm_xray_target_minimizer(object):
    def __init__(self,
                 fmodel,
                 model,
                 xs_initial,
                 adp_nmas,
                 fsub_adp,
                 selections,
                 selections_1d,
                 max_iterations,
                 n_modes,
                 weight_nmre = 1.0,
                 run_finite_differences_test = False,
                 correct_adp = True,
                 zero_mode_flag = True):
        adopt_init_args(self, locals())
        fmodel.xray_structure.scatterers().flags_set_grads(state=False)
        xray.set_scatterer_grad_flags(scatterers = fmodel.xray_structure.scatterers(),
                                      u_aniso = True)
        if(self.run_finite_differences_test): self.correct_adp = False
        self.fmodel_copy = self.fmodel.deep_copy()
        self.target_functor = self.fmodel_copy.target_functor()
        self.run_finite_differences_test_counter = 0
        self.counter = 0
        self.n_groups = len(self.xs_initial)
        self.dim_x    = len(self.xs_initial[0])
        self.xs_min    = self.xs_initial
        self.x = self.pack(self.xs_min)
#        for adp_nma_selected, selection in zip(adp_nmas, selections):
#            weights_selected = self.fmodel.xray_structure.atomic_weights.select(selection)
#            modes1d_selected = selected_modes_to_1D(modes = self.modes, n_modes = self.n_modes,
#                                                selection = self.selection)
#            assert len(modes1d_selected)/n_modes == len(weights_selected)
#            adp_nma_selected = init_nm_adp(modes = modes1d_selected,
#                              weights = weights_selected,
#                              n_modes = n_modes,
#                              zero_mode_flag = zero_mode_flag)
#            self.adp_nmas.append(adp_nma_selected)
        self.minimizer = lbfgs.run(
                target_evaluator = self,
                core_params = lbfgs.core_parameters(maxfev = 10),
                termination_params = lbfgs.termination_parameters(
                    min_iterations = max_iterations,
                    max_calls = int(max_iterations*1.5)),
                exception_handling_params = lbfgs.exception_handling_parameters(
                    ignore_line_search_failed_step_at_lower_bound = True,
                    ignore_line_search_failed_step_at_upper_bound = True,
                    ignore_line_search_failed_maxfev              = True))
        self.compute_functional_and_gradients()
        del self.x
        self.xs_result = self.xs_min

    def pack(self, xs):
        v = []
        for x in xs:
            v += list(x)
        return flex.double(tuple(v))

    def unpack(self):
        i = 0
        for j in xrange(self.n_groups):
            self.xs_min[j] = flex.double(tuple(self.x)[i:i+self.dim_x])
            i += self.dim_x

    def compute_functional_and_gradients(self):
        self.counter += 1
        self.unpack()
        update_xray_structure_with_nm(
                xray_structure = self.fmodel_copy.xray_structure,
                adp_nmas       = self.adp_nmas,
                selections     = self.selections,
                xs             = self.xs_min,
                n_modes        = self.n_modes,
                selections_1d  = self.selections_1d,
                correct_adp    = self.correct_adp,
                zero_mode_flag = self.zero_mode_flag)
        self.fmodel_copy.update_xray_structure(update_f_calc=True)
        t_r = self.target_functor(compute_gradients=True)
        self.f = t_r.target_work()
        if (self.weight_nmre != 0.0  and self.fsub_adp is not None):
#            energies_adp_aniso = self.model.restraints_manager.energies_adp_aniso(
#                    xray_structure = self.fmodel_copy.xray_structure,
#                    compute_gradients = True,
#                    selection      = None)
#            resi_target = energies_adp_aniso.target
#	    print resi_target
#            u_aniso_gradients = energies_adp_aniso.gradients_aniso_cart
#           u_aniso_gradients_sel = u_aniso_gradients.select(selections_1d)
#            d_xray_d_uaniso = t_r.gradients_wrt_atomic_parameters(
#                          u_aniso=True)
#           d_xray_d_uaniso_sel = d_target_d_uaniso.select(selections_1d)
            rest_target = 0.0
            rest_gradients = []
            for x, adp_nma, selection in zip(self.xs_min, self.adp_nmas, self.selections):
                fsub_adp_selected = self.fsub_adp.select(selection)
                restraints_manager = nm_from_uaniso_target_and_grads(x = x,
                                        adp_nma = adp_nma,
                                        uanisos = fsub_adp_selected,
                                        n_modes = self.n_modes,
                                        zero_mode_flag = self.zero_mode_flag)
                rest_target += restraints_manager.target()
                rest_gradients.append(restraints_manager.grad_nm())
#	    if(self.counter < 5):
#                rest_resi = grad_resi(gradients = u_aniso_gradients)
#                xray_resi = grad_resi(gradients = d_xray_d_uaniso)
#                rxtogeom = math.sqrt(rest_resi/xray_resi)
#	        if rxtogeom == 0.0:
#		    self.scalrest = self.weight_nmre
#                else:
#		    self.scalrest = self.weight_nmre/rxtogeom
#           self.f += resi_target*self.scalrest
#           d_target_d_uaniso = add_grads(gradients1 = d_xray_d_uaniso,
#                                      gradients2 = u_aniso_gradients,
#                                      weight = self.scalrest)
        d_xray_d_uaniso = t_r.gradients_wrt_atomic_parameters(
                            u_aniso=True)
        grad_manager = nm_xray_grads(
                gradients     = d_xray_d_uaniso,
                adp_nmas      = self.adp_nmas,
                selections    = self.selections,
                xs            = self.xs_min,
                n_modes       = self.n_modes,
                zero_mode_flag= self.zero_mode_flag)
        if (self.weight_nmre != 0.0  and self.fsub_adp is not None):
            d_xray_d_nm = grad_manager.grad_nm
            rest_resi = 0.
            xray_resi = 0.
            for xray_grad, rest_grad in zip(d_xray_d_nm, rest_gradients):
                rest_resi += grad_resi(gradients = rest_grad)
                xray_resi += grad_resi(gradients = xray_grad)
            rxtogeom = math.sqrt(rest_resi/xray_resi)
            if rxtogeom == 0.:
    	        self.scalrest = self.weight_nmre
            else:
    	        self.scalrest = self.weight_nmre/rxtogeom
            xray_grad = self.pack(d_xray_d_nm)
            rest_grad = self.pack(rest_gradients)
            self.g = add_grads(gradients1 = xray_grad,
                               gradients2 = rest_grad,
                               weight     = self.scalrest)
            self.f += rest_target*self.scalrest
        else:
            self.g = self.pack(grad_manager.grad_nm)
#        self.g = self.pack(grad_manager.grad_nm)
#        sumg = flex.mean_sq(u_aniso_gradients)
#        sumx = flex.mean_sq(self.g)
#        rxtogeom = math.sqrt(sumg/sumx)
#        scalrest = self.weight_nmre/rxtogeom

        if(self.run_finite_differences_test and
           self.run_finite_differences_test_counter < 2):
            tolerance = 1.e-3
            self.run_finite_differences_test_counter += 1
            grad_x = finite_differences_grads_of_xray_target_wrt_nm(
                    target_functor = self.target_functor,
                    xs = self.xs_min,
                    adp_nmas = self.adp_nmas,
                    selections = self.selections,
                    n_modes = self.n_modes,
                    zero_mode_flag = self.zero_mode_flag,
                    delta = 0.00001)
            GX = self.pack(grad_x)
            for m1, m2 in zip(self.g, GX):
                assert approx_equal(m1, m2, tolerance)
        return self.f, self.g

def update_xray_structure_with_nm(xray_structure,
                                 adp_nmas,
                                 selections,
                                 xs,
                                 n_modes,
                                 selections_1d = None,
                                 correct_adp = True,
                                 zero_mode_flag = True):
    global time_update_xray_structure_with_nm
    timer = user_plus_sys_time()
    total_size = xray_structure.sites_cart().size()
    u_cart_from_nm_ = u_cart_from_nm(adp_nmas = adp_nmas,
                                     selections = selections,
                                     xs = xs,
                                     n_modes = n_modes,
				     total_size = total_size,
                                     zero_mode_flag = zero_mode_flag)
    xray_structure.set_u_cart(u_cart=u_cart_from_nm_, selection = selections_1d)
    if(correct_adp): xray_structure.tidy_us(u_min = 1.e-6)
    time_update_xray_structure_with_nm += timer.elapsed()


def u_cart_from_nm(adp_nmas, selections, xs, n_modes, total_size, zero_mode_flag = True):
    global time_u_cart_from_nm
    t1 = time.time()
    uanisos = flex.sym_mat3_double(total_size, [0,0,0,0,0,0])
# where the problem is!!
#fixed!
    for adp_nma, selection, x in zip(adp_nmas, selections, xs):
        uaniso_form_s_manager = uaniso_from_s(x = x,
                          adp_nma = adp_nma,
                          n_modes = n_modes,
                          zero_mode_flag = zero_mode_flag)
        u = uaniso_form_s_manager.u_cart()
        uanisos.set_selected(selection, u)
    t2 = time.time()
    time_u_cart_from_nm += (t2 - t1)
    return uanisos

#def nm_from_u_cart(xray_structure,
#                   xs_initial,
#                   selections,
#                   n_modes,
#                   zero_mode_flag = True,
#                   number_of_macro_cycles = 100,
#                   max_iterations         = 100):
#    global time_nm_from_u_cart
#    timer = user_plus_sys_time()
#    uc = xray_structure.unit_cell()
#    xray_structure.tidy_us(u_min = 1.e-6)
#    ueq = xray_structure.extract_u_iso_or_u_equiv()
#    assert (ueq < 0.0).count(True) == 0
#    u_cart = xray_structure.scatteres().extract_u_cart(uc)
#    for selection in selections:
#        u_cart_selected = u_cart.select(selection)
#        assert adptbx.is_positive_definite(u_cart_selected, 1.e-6).count(False)==0
#    xray_structure.tidy_us(u_min = 1.e-6)


def finite_differences_grads_of_xray_target_wrt_nm(target_functor,
                                                   xs,
                                                   adp_nmas,
                                                   selections,
                                                   n_modes,
                                                   zero_mode_flag = True,
                                                   delta=0.00001):
    fmodel = target_functor.manager
    derivative_xs = []
    for j in xrange(len(xs)):
        dx = []
        for i in xrange(len(xs[j])):
            target_values = []
            for d_sign in (-1, 1):
                x_ = []
                for item in xs:
                    x_.append(list(item))
                d = d_sign*delta
                x_[j][i] += d
                update_xray_structure_with_nm(
                        xray_structure = fmodel.xray_structure,
                        adp_nmas       = adp_nmas,
                        selections     = selections,
                        xs             = x_,
                        n_modes        = n_modes,
                        correct_adp    = False,
                        zero_mode_flag = zero_mode_flag)
                fmodel.update_xray_structure(update_f_calc=True)
                t_w = target_functor(compute_gradients=False).target_work()

                target_values.append(t_w)
            derivative = (target_values[1] - target_values[0]) / (2 * delta)
            dx.append(derivative)
        derivative_xs.append(dx)

    return derivative_xs

def make_nm_compatible_with_u_positive_definite(
        xray_structure,
        adp_nmas,
        selections,
        max_iterations,
        number_of_u_nonpositive_definite,
        eps,
        number_of_macro_cycles_for_nm_from_uanisos,
        n_modes,
        xs,
	selections1d = None,
        zero_mode_flag = True,
        out = None):
    global time_make_nm_compatible_with_u_positive_definite
    t1 = time.time()
    if(out is None): out = sys.stdout
    for i in range(1, max_iterations+1):
        update_xray_structure_with_nm(
                xray_structure = xray_structure,
                adp_nmas       = adp_nmas,
                selections     = selections,
                xs             = xs,
                n_modes        = n_modes,
		selections_1d  = selections_1d,
                zero_mode_flag = zero_mode_flag)
        ipad_1 = xray_structure.is_positive_definite_u()
        if(i == 1 or i == max_iterations):
            xray_structure.show_u_statistics(out = out)
        xray_structure.tidy_us(u_min = 1.e-6)
        xs = nm_from_uanisos(
                xray_structure = xray_structure,
                selections     = selections,
                modes          = None,
                xs_initial     = xs,
                n_modes        = n_modes,
                adp_nmas       = adp_nmas,
                number_of_macro_cycles = number_of_macro_cycles_for_nm_from_uanisos,
                max_iterations = 100,
                zero_mode_flag = zero_mode_flag)
        if(i == max_iterations): xray_structure.show_u_statistics(out = out)
        if(ipad_1.count(False) == number_of_u_nonpositive_definite):
            break
    assert xray_structure.is_positive_definite_u().count(False) == 0
    t2 = time.time()
    time_make_nm_compatible_with_u_positive_definite += (t2 - t1)
    return xs

class nm_refinement(object):
    def __init__(self,
                 fmodel,
                 model,
                 selections,
                 selections_1d,
                 n_modes,
                 number_of_macro_cycles,
                 max_number_of_iterations,
                 weight_nmre = 1.0,
                 evecin = "evec.dat",
                 evalin = "eval.dat",
                 fsub_evecin = "proj.dat",
		 nm_params_filename = "nm_params.dat",
                 zero_mode_flag = True,
                 zero_mode_input_flag = False,
  		 update_evec = False,
	         enmcalc_path = "enmcalc_fix.csh",
                 start_xs_value = None,
                 zero_mode_corr = True,
                 fsub_ref = False,
                 fsub_n_modes = 20,
                 run_finite_differences_test = False,
                 eps = 1.e-6,
                 out = None,
                 macro_cycle = None,
                 verbose = True):
        adopt_init_args(self, locals())
        global time_nm_total
        timer = user_plus_sys_time()
        if(out is None): out = sys.stdout
        prefix = "NM refinement:"
        fmodel.info().show_targets(text = prefix + " start model", out = out)
        fmodel.xray_structure.show_u_statistics(text = prefix+" start model",
                                                out = out)
        xrs = fmodel.xray_structure
        pdb_hierarchy = model.pdb_hierarchy()
#pdb_hierarchy misssing a parameter, I am not sure whether it is all right.
#        xrs.tidy_us(u_min = 1.e-6)
	if(update_evec):
	    f_pdb = open("tmp.pdb", "w")
	    m = model.deep_copy()
	    m.write_pdb_file(out = f_pdb)
	    f_pdb.close()
	    try:
		retcode = subprocess.call(enmcalc_path)
 		if retcode < 0:
		    print >>out, "Enmcalc was terminated by signal", -retcode
		else:
		    print >>out, "Enmcalc completed", retcode
		    evalin = "eval.dat"
		    evecin = "evec.dat"
                    fsub_evecin = "proj.dat"
		    print >>out, "evalin was changed to", evalin
		    print >>out, "evecin was changed to", evecin
                    print >>out, "fsub_evecin was change to", fsub_evecin
	    except OSError as e:
		print >>out, "Executing enmcalc failed:", e
        if(fsub_ref):
            self.total_evec = evecin
            evecin = fsub_evecin
        modes = generate_evec(selections = selections,
			      selections_1d = selections_1d,
                              xray_structure = xrs,
                              pdb_hierarchy = pdb_hierarchy,
                              filename = evecin,
                              n_modes = n_modes,
                              zero_mode_input_flag = zero_mode_input_flag,
                              zero_mode_flag = zero_mode_flag)
#selections will be modified during this step
#        for nm_sel in selections:
#            fmodel.xray_structure.convert_to_anisotropic(selection = nm_sel)
        if(zero_mode_corr):
            curr_zero_mode_flag = zero_mode_flag
            zero_mode_flag = False
        adp_nmas = []
        for selection in selections:
            modes1d = selected_modes_to_1D(modes = modes, n_modes = n_modes, selection = selection)
            weights_selected = xrs.atomic_weights().select(selection)
            print "The number of selected atoms is %d." % len(weights_selected)
            adp_nma = init_nm_adp(modes = modes1d,
                                  weights = weights_selected,
                                  n_modes = n_modes,
                                  zero_mode_flag = zero_mode_flag)
            adp_nmas.append(adp_nma)
	u_cart = xrs.scatterers().extract_u_cart(xrs.unit_cell())
        if(fsub_ref):
            fsub_ref_manager = self.fsub_total_refinement()
            fsub_adp = xrs.scatterers().extract_u_cart(xrs.unit_cell())
        else:
            fsub_adp = None
        if(start_xs_value is not None):
            xs = nm_params_reader(filename = start_xs_value)
        elif(macro_cycle == 1):
            xs_initial = []
            nmval = read_nmval_file(evalin = evalin,
                                    n_modes = n_modes,
                                    zero_mode_input_flag = zero_mode_input_flag,
                                    zero_mode_flag = zero_mode_flag)
            for adp_nma, selection in zip(adp_nmas, selections):
                x = init_nm_para(nmval = nmval,
                                 n_modes = n_modes,
				 zero_mode_flag = zero_mode_flag)
                u_cart_selected = u_cart.select(selection)
                uaniso_from_s_manager = uaniso_from_s(x = x,
                                                      adp_nma = adp_nma,
                                                      n_modes = n_modes,
                                                      zero_mode_flag = zero_mode_flag)
                u = uaniso_from_s_manager.u_cart()
                x_scaled = scale_x(x = x,
                                   uanisos = u_cart_selected,
                                   adp_all = u,
                                   n_modes = n_modes,
                                   zero_mode_flag = zero_mode_flag)
                xs_initial.append(x_scaled)
            xs = nm_from_uanisos(xray_structure = xrs,
                                 selections = selections,
                                 modes      = None,
                                 xs_initial = xs_initial,
                                 n_modes    = n_modes,
                                 adp_nmas   = adp_nmas,
                                 number_of_macro_cycles = 2,
                                 max_iterations = self.max_number_of_iterations,
				 zero_mode_flag = zero_mode_flag,
                                 verbose    = -1)
        else:
            if(os.path.isfile("fsub_params.dat") is True and fsub_ref is False):
                xs_initial = nm_params_reader(filename = "fsub_params.dat")
            else:
                xs_initial = model.nm_groups.xs
            xs = nm_from_uanisos(xray_structure = xrs,
                                 selections = selections,
                                 modes      = None,
                                 xs_initial = xs_initial,
                                 n_modes    = n_modes,
                                 adp_nmas   = adp_nmas,
                                 number_of_macro_cycles = 3,
                                 max_iterations = self.max_number_of_iterations,
				 zero_mode_flag = zero_mode_flag,
                                 verbose    = -1)
        if(verbose):
            show_nm(xs = xs,
                    n_modes = n_modes,
                    text = prefix + " start parameters",
                    zero_mode_flag = zero_mode_flag,
                    out = out)
        for macro_cycle in range(1, number_of_macro_cycles+1):
            print >> out
            prefix = "NM refinement: after macrocycle "+str(macro_cycle)
            minimized = nm_xray_target_minimizer(
                    fmodel                      = fmodel,
		    model			= model,
                    xs_initial                  = xs,
                    adp_nmas                    = adp_nmas,
                    fsub_adp                    = fsub_adp,
                    selections                  = selections,
                    selections_1d               = selections_1d,
                    max_iterations              = max_number_of_iterations,
                    n_modes                     = n_modes,
                    weight_nmre                 = weight_nmre,
                    run_finite_differences_test = run_finite_differences_test,
                    zero_mode_flag              = zero_mode_flag)
            xrs = minimized.fmodel_copy.xray_structure
            xrs.show_u_statistics(text = prefix, out = out)
            if(verbose):
                show_nm(xs = minimized.xs_result,
                        n_modes = n_modes,
                        text = prefix,
                        zero_mode_flag = zero_mode_flag,
                        out = out)
            fmodel.update_xray_structure(xray_structure = xrs,
                                         update_f_calc  = True)
            fmodel.info().show_targets(text = prefix, out = out)
            if(xrs.is_positive_definite_u().count(False) > 0):
                xrs.tidy_us(u_min = 1.e-6)
                xrs.show_u_statistics(
                        text = prefix+": after making positive definite",
                        out  = out)
                fmodel.update_xray_structure(xray_structure = xrs,
                                             update_f_calc  = True)
                fmodel.info().show_targets(text = prefix+": after making positive definite",
                                           out  = out)
                xs = make_nm_compatible_with_u_positive_definite(
                        xray_structure = xray_structure,
                        adp_nmas       = adp_nmas,
                        selections     = selections,
                        max_iterations = 100,
                        number_of_u_nonpositive_definite = 0,
                        eps            = eps,
                        number_of_macro_cycles_for_nm_from_uanisos = 2,
                        n_modes        = n_modes,
                        xs             = minimize.xs_result,
			selections_1d  = selections_1d,
                        zero_mode_flag = zero_mode_flag,
                        out            = out)
            else: xs = minimized.xs_result
        if(verbose):
            show_nm(xs = xs,
                    n_modes = n_modes,
                    text = "NM refinement: final values",
                    zero_mode_flag = zero_mode_flag,
                    out = out)
        self.xs = xs
        filename = nm_params_filename
        nm_params_writer(xs = minimized.xs_result,
                         n_modes = n_modes,
                         zero_mode_flag = zero_mode_flag,
                         filename = filename)
        if(os.path.isfile("fsub_params.dat") is False or fsub_ref is True):
            model.nm_groups.xs = self.xs
        self.fmodel = fmodel
        time_nm_total += timer.elapsed()
    def fsub_total_refinement(self):
        selections = []
        selections.append(self.selections_1d)
        fsub_refinement = self.__class__(
                 fmodel = self.fmodel,
                 model = self.model,
                 selections = selections,
                 selections_1d = self.selections_1d,
                 n_modes = self.fsub_n_modes,
                 number_of_macro_cycles = self.number_of_macro_cycles,
                 max_number_of_iterations = self.max_number_of_iterations,
                 weight_nmre = 0.0,
                 evecin = self.total_evec,
                 evalin = self.evalin,
                 fsub_evecin = self.fsub_evecin,
		 nm_params_filename = "fsub_params.dat",
                 zero_mode_flag = True,
                 zero_mode_input_flag = False,
  		 update_evec = False,
	         enmcalc_path = "mgr.csh",
                 start_xs_value = None,
                 zero_mode_corr = self.zero_mode_corr,
                 fsub_ref = False,
                 fsub_n_modes = self.fsub_n_modes,
                 run_finite_differences_test = False,
                 eps = 1.e-6,
                 out = self.out,
                 macro_cycle = self.macro_cycle,
                 verbose = self.verbose)


def nm_params_writer(xs,
                    n_modes,
                    zero_mode_flag,
                    filename = "nm_params.dat"):
    with open(filename, 'w') as f:
        f.write("n_modes " + str(n_modes) + "\n")
        if(zero_mode_flag):
            f.write("zero_mode_flag True\n")
        else:
            f.write("zero_mode_flag False\n")
        f.write("n_groups "+str(len(xs))+"\n")
        for x in xs:
            for x_ in x:
                f.write(str(x_)+"\n")

def nm_params_reader(filename = "nm_params.dat"):
    with open(filename, 'r') as f:
        line = f.readline()
        n_modes = int(line.split()[1])
        line = f.readline()
        if(line.split()[1] == "True"):
            n_params = 21 + (n_modes - 6)*(n_modes - 5)/2
        elif(line.split()[1] == "False"):
            n_params = (n_modes + 1)*n_modes/2
        line = f.readline()
        n_groups = int(line.split()[1])
        xs = []
        for i in range(n_groups):
            x = []
            for j in range(int(n_params)):
                line = f.readline()
                x.append(float(line))
            x = flex.double(tuple(x))
            xs.append(x)
    return xs
