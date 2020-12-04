#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/nm/nm.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python.hpp>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(mmtbx::nm::common)

namespace mmtbx { namespace nm {
    namespace bp = boost::python;

namespace {
    void init_module()
    {
        using namespace boost::python;
        using boost::python::arg;

        class_<nm_init>("nm_init",
                    init<const char*,
                        std::size_t,
                        af::shared<iotbx::pdb::hierarchy::atom> const&,
                        bool,
                        bool>(
                            (arg("filename"),
                             arg("n_modes"),
                             arg("atoms"),
                             arg("zero_mode_input_flag"),
                             arg("zero_mode_flag"))))
        .def("normalize", &nm_init::normalize)
        .def("return_modes", &nm_init::return_modes)
        .def("print_eigenvector", &nm_init::print_eigenvector)
        .def("gen_zero_modes", &nm_init::gen_zero_modes)
        .def("return_zero_modes", &nm_init::return_zero_modes)
        .def("return_new_selection", &nm_init::return_new_selection)
        ;
        def("normalize_mode", (void (*)(af::shared<vec3<double> >)) normalize_mode, (arg("mode")))
        ;
        def("init_nm_adp", (af::shared<sym_mat3<double> >(*)(af::shared<vec3<double> >,
                                                            af::shared<double>,
                                                            std::size_t,
                                                            bool)) init_nm_adp,
                            (arg("modes"), arg("weights"), arg("n_modes"), arg("zero_mode_flag")))
        ;
        def("unpack_x", (af::versa<double, af::c_grid<2> >(*)(af::shared<double>, std::size_t, bool))
                                                                    unpack_x,
                        (arg("x"), arg("n_modes"), arg("zero_mode_flag")))
        ;
        def("scale_x", (af::shared<double>(*)(af::shared<double>, af::shared<sym_mat3<double> >,
                                              af::shared<sym_mat3<double> >,
                                              std::size_t,
                                              bool)) scale_x,
                       (arg("x"), arg("uanisos"), arg("adp_all"), arg("n_modes"), arg("zero_mode_flag")))
        ;
        def("grad_resi", (double(*)(af::shared<double>)) grad_resi,
                       (arg("gradients")))
        ;
        def("add_grads", (af::shared<double>(*)(af::shared<double>, af::shared<double>, double)) add_grads,
                       (arg("gradients1"), arg("gradients2"), arg("weight")))
        ;
        class_<uaniso_from_s>("uaniso_from_s",
                              init<af::shared<double> const&,
                              af::shared<sym_mat3<double> > const&,
                              std::size_t,
                              bool>(
                                  (arg("x"),
                                   arg("adp_nma"),
                                   arg("n_modes"),
                                   arg("zero_mode_flag"))))
        .def("u_cart", &uaniso_from_s::u_cart)
        ;
        class_<d_target_d_nm>("d_target_d_nm",
                              init<af::shared<sym_mat3<double> > const&,
                              af::shared<sym_mat3<double> > const&,
                              af::shared<double> const&,
                              std::size_t,
                              bool>(
                                  (arg("adp_nma"),
                                   arg("d_target_d_uaniso"),
                                   arg("x"),
                                   arg("n_modes"),
                                   arg("zero_mode_flag"))))
        .def("grad_nm", &d_target_d_nm::grad_nm)
        ;
        class_<nm_from_uaniso_target_and_grads>("nm_from_uaniso_target_and_grads",
                                                init<af::shared<double> const&,
                                                af::shared<sym_mat3<double> > const&,
                                                af::shared<sym_mat3<double> > const&,
                                                std::size_t,
                                                bool>(
                                                    (arg("x"),
                                                     arg("adp_nma"),
                                                     arg("uanisos"),
                                                     arg("n_modes"),
                                                     arg("zero_mode_flag"))))
         .def("target", &nm_from_uaniso_target_and_grads::target)
         .def("grad_nm", &nm_from_uaniso_target_and_grads::grad_nm)
         ;

    }
}
}}

BOOST_PYTHON_MODULE(mmtbx_nm_ext)
{
    mmtbx::nm::init_module();
}
