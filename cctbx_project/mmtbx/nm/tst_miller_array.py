from iotbx import reflection_file_reader

reflection_file = reflection_file_reader.any_reflection_file(file_name = "./4ksc-sf.cif")
miller_arrays = reflection_file.as_miller_arrays()
for miller_array in miller_arrays:
    if (miller_array.is_xray_amplitude_array()):
        f_obs = miller_array
f_obs.show_comprehensive_summary()
