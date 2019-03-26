__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.3.0"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

from numpy import array

SET_INIT = {'Rhub': 0.006, 'R': 0.45, 'B': 3, 'turbine_name': 'norvezan_karlsen_ntnu_njegove_krivulje', 'r': array([0.0075, 0.0225, 0.049 , 0.055 , 0.0675, 0.0825, 0.0975, 0.1125,
       0.1275, 0.1425, 0.1575, 0.1725, 0.1875, 0.2025, 0.2175, 0.2325,
       0.2475, 0.2775, 0.2925, 0.3075, 0.3225, 0.3375, 0.3675, 0.3825,
       0.3975, 0.4125, 0.4275, 0.4425]), 'c': array([0.0135    , 0.0135    , 0.0135    , 0.0495    , 0.08143281,
       0.08011114, 0.07701157, 0.07312555, 0.06900783, 0.064952  ,
       0.06110231, 0.05752027, 0.0542229 , 0.05120426, 0.04844738,
       0.04593074, 0.04363179, 0.03960108, 0.03783062, 0.03620078,
       0.03469692, 0.03201685, 0.03081913, 0.02970404, 0.02866371,
       0.0276912 , 0.02678034, 0.02592566]), 'theta': array([ 1.200000e+02,  1.200000e+02,  1.200000e+02,  3.800000e+01,
        3.705499e+01,  3.254410e+01,  2.867674e+01,  2.526190e+01,
        2.243033e+01,  1.998782e+01,  1.803445e+01,  1.634885e+01,
        1.466325e+01,  1.306715e+01,  1.182907e+01,  1.075265e+01,
        8.882724e+00,  7.987743e+00,  7.252693e+00,  6.564977e+00,
        5.918680e+00,  4.718537e+00,  4.131617e+00,  3.543894e+00,
        2.943336e+00,  2.218509e+00,  1.097028e-01, -7.167448e-01]), 'foils': ['s826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826'], 'tip_loss': False, 'hub_loss': False, 'new_tip_loss': False, 'new_hub_loss': False, 'cascade_correction': False, 'rotational_augmentation_correction': False, 'rotational_augmentation_correction_method': 0, 'max_iterations': 100.0, 'convergence_limit': 0.001, 'rho': 1.225, 'method': 7, 'fix_reynolds':False, 'reynolds':50000, 'linspace_interp': False, 'num_interp': 5.0, 'v_min': 5.0, 'v_max': 15.0, 'v_num': 10.0, 'rpm_min': 500.0, 'rpm_max': 3000.0, 'rpm_num': 10.0, 'relaxation_factor': 0.3, 'print_all': False, 'print_out': False, 'target_rpm': 1500.0, 'target_speed': 10.0, 'delta_start': 5.0, 'decrease_factor': 0.1, 'min_delta': 0.001, 'min_add_angle': -45.0, 'max_add_angle': 45.0, 'angle_step': 0.5, 'curves': {'s826': {'x': [1.0, 0.99664, 0.987056, 0.972314, 0.953332, 0.930531, 0.903794, 0.872736, 0.837646, 0.799007, 0.75733, 0.71316, 0.667064, 0.619624, 0.571426, 0.523051, 0.475063, 0.428007, 0.382245, 0.337323, 0.293568, 0.251455, 0.211424, 0.173909, 0.139301, 0.107962, 0.08022, 0.056312, 0.036491, 0.020843, 0.009548, 0.002548, 0.000169, 1e-06, 0.0, 0.000198, 0.000946, 0.002165, 0.003667, 0.013639, 0.029183, 0.049938, 0.075768, 0.106324, 0.141288, 0.179884, 0.220074, 0.261526, 0.304974, 0.350096, 0.397636, 0.447637, 0.500087, 0.554624, 0.610364, 0.666289, 0.721305, 0.774276, 0.82406, 0.869551, 0.909558, 0.9427, 0.968241, 0.986101, 0.996568, 1.0], 'y': [0.0, 0.000985, 0.004265, 0.009834, 0.017037, 0.024949, 0.032768, 0.040403, 0.0481, 0.055842, 0.063544, 0.071084, 0.078312, 0.085056, 0.091118, 0.096277, 0.10028, 0.102811, 0.103258, 0.101581, 0.0983, 0.093644, 0.087785, 0.080881, 0.073078, 0.064532, 0.055386, 0.04581, 0.035989, 0.026134, 0.016553, 0.007584, 0.001654, 0.000102, 6e-06, -0.001417, -0.002691, -0.003992, -0.005223, -0.010342, -0.015169, -0.019578, -0.023565, -0.027219, -0.03072, -0.034729, -0.038653, -0.040769, -0.040625, -0.038079, -0.032821, -0.025646, -0.017114, -0.008267, 4.9e-05, 0.007219, 0.012793, 0.016474, 0.018113, 0.017689, 0.015218, 0.011278, 0.006966, 0.003275, 0.00085, 0.0], 'max_thickness': 0.5}}, 'dr': array([0.    , 0.015 , 0.0265, 0.006 , 0.0125, 0.015 , 0.015 , 0.015 ,
       0.015 , 0.015 , 0.015 , 0.015 , 0.015 , 0.015 , 0.015 , 0.015 ,
       0.015 , 0.03  , 0.015 , 0.015 , 0.015 , 0.015 , 0.03  , 0.015 ,
       0.015 , 0.015 , 0.015 , 0.015 ]), 'r_in': array([0.0075, 0.0225, 0.049 , 0.055 , 0.0675, 0.0825, 0.0975, 0.1125,
       0.1275, 0.1425, 0.1575, 0.1725, 0.1875, 0.2025, 0.2175, 0.2325,
       0.2475, 0.2775, 0.2925, 0.3075, 0.3225, 0.3375, 0.3675, 0.3825,
       0.3975, 0.4125, 0.4275, 0.4425]), 'c_in': array([0.0135    , 0.0135    , 0.0135    , 0.0495    , 0.08143281,
       0.08011114, 0.07701157, 0.07312555, 0.06900783, 0.064952  ,
       0.06110231, 0.05752027, 0.0542229 , 0.05120426, 0.04844738,
       0.04593074, 0.04363179, 0.03960108, 0.03783062, 0.03620078,
       0.03469692, 0.03201685, 0.03081913, 0.02970404, 0.02866371,
       0.0276912 , 0.02678034, 0.02592566]), 'theta_in': array([ 1.200000e+02,  1.200000e+02,  1.200000e+02,  3.800000e+01,
        3.705499e+01,  3.254410e+01,  2.867674e+01,  2.526190e+01,
        2.243033e+01,  1.998782e+01,  1.803445e+01,  1.634885e+01,
        1.466325e+01,  1.306715e+01,  1.182907e+01,  1.075265e+01,
        8.882724e+00,  7.987743e+00,  7.252693e+00,  6.564977e+00,
        5.918680e+00,  4.718537e+00,  4.131617e+00,  3.543894e+00,
        2.943336e+00,  2.218509e+00,  1.097028e-01, -7.167448e-01]), 'foils_in': ['s826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826', 's826']}
        