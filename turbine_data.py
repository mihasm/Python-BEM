__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.2.6"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

from numpy import array

SET_INIT = {
    "Rhub": 0.1,
    "R": 0.776,
    "B": 5,
    "r": array(
        [
            0.15,
            0.18,
            0.21,
            0.24,
            0.27,
            0.3,
            0.33,
            0.36,
            0.39,
            0.42,
            0.45,
            0.48,
            0.51,
            0.54,
            0.57,
            0.6,
            0.63,
            0.66,
            0.69,
            0.72,
            0.75,
            0.76,
        ]
    ),
    "c": array(
        [
            0.12429,
            0.14448,
            0.15224,
            0.14968,
            0.14213,
            0.13129,
            0.12198,
            0.1136,
            0.10666,
            0.10015,
            0.09487,
            0.08959,
            0.08441,
            0.08029,
            0.07728,
            0.07371,
            0.07097,
            0.06797,
            0.06637,
            0.06337,
            0.06072,
            0.05984,
        ]
    ),
    "theta": array(
        [
            29.24,
            33.41,
            34.9,
            34.6,
            32.39,
            20.96,
            21.71,
            13.33,
            11.04,
            8.18,
            7.05,
            5.88,
            3.41,
            1.09,
            0.57,
            -0.2,
            -0.66,
            -2.05,
            -2.29,
            -2.6,
            -2.92,
            -2.37,
        ]
    ),
    "AoA_cL": [
        -0.271260997,
        0.542521994,
        1.085043988,
        1.627565982,
        2.170087977,
        2.712609971,
        3.255131965,
        4.340175953,
        4.882697947,
        6.239002933,
        7.052785924,
        8.409090909,
        10.03665689,
        10.85043988,
        11.66422287,
        12.47800587,
        13.56304985,
        14.10557185,
        14.37683284,
        14.64809384,
        14.64809384,
        14.64809384,
        14.91935484,
        16.27565982,
        18.4457478,
        20.61583578,
        22.24340176,
        23.87096774,
        25.49853372,
        27.3973607,
        30.65249267,
        33.09384164,
        34.99266862,
        37.97653959,
        40.68914956,
        44.48680352,
        47.47067449,
        49.91202346,
        53.16715543,
        54.52346041,
        57.23607038,
        59.40615836,
        62.39002933,
        64.2888563,
        66.45894428,
        69.17155425,
        72.69794721,
        74.59677419,
        78.12316716,
        81.64956012,
        85.44721408,
        88.97360704,
        -86.68621701,
        -82.88856305,
        -80.17595308,
        -75.2932551,
        -70.9530792,
        -67.6979472,
        -64.9853372,
        -61.730205299999994,
        -59.288856300000006,
        -55.2199413,
        -50.60850439999999,
        -46.2683284,
        -44.640762499999994,
        -42.4706745,
        -40.0293255,
        -38.13049849999999,
        -35.96041059999999,
        -33.247800600000005,
        -30.26392960000001,
        -27.5513196,
        -25.3812317,
        -22.66862169999999,
        -20.227272699999986,
        -17.514662799999996,
        -14.802052800000013,
        -12.0894428,
        -10.461876799999999,
        -9.648093799999998,
        -9.376832799999988,
        -8.56304990000001,
        -7.478005899999999,
        -6.664222899999999,
        -5.579178899999988,
        -4.494134900000006,
        -3.4090908999999954,
        -2.8665689000000043,
        -2.3240469000000132,
        -1.7815248999999937,
        -1.2390029000000027,
        -0.6964809000000116,
        -0.42521990000000187,
    ],
    "AoA_cD": [
        0.00059512,
        1.050981948,
        2.62507439,
        6.299940488,
        8.925014878,
        12.33802817,
        15.22614561,
        16.27772267,
        17.59710375,
        21.0208292,
        23.65721087,
        27.60761754,
        30.24042849,
        31.82642333,
        34.4568538,
        36.03987304,
        40.25570323,
        43.67942868,
        47.89763936,
        52.89724261,
        58.433049,
        63.43086689,
        68.42808966,
        73.68359452,
        78.94326523,
        84.19400913,
        89.97083912,
        -84.77841698,
        -79.7925015,
        -74.0228129,
        -69.0470145,
        -63.8051974,
        -58.302122600000004,
        -53.322753399999996,
        -47.5608014,
        -42.3267209,
        -37.09442569999999,
        -31.86629640000001,
        -26.6399524,
        -21.935528699999992,
        -16.183693699999992,
        -5.182305100000008,
    ],
    "cL": [
        0.009150327,
        0.04725273,
        0.107149611,
        0.183386362,
        0.243283244,
        0.308626748,
        0.412096614,
        0.499210639,
        0.613573751,
        0.722466282,
        0.798695047,
        0.880354462,
        0.967452514,
        1.000108294,
        1.021870827,
        1.049079984,
        1.059941285,
        1.032692197,
        0.989111226,
        0.896510647,
        0.733111954,
        0.83115117,
        0.667744491,
        0.613238329,
        0.602281193,
        0.607663926,
        0.640295747,
        0.672927568,
        0.721899258,
        0.792649455,
        0.901486082,
        0.966773682,
        1.00484414,
        1.031989407,
        1.048249414,
        1.048137606,
        1.042603134,
        1.020744766,
        0.982522569,
        0.960696146,
        0.911596675,
        0.873406424,
        0.807959098,
        0.764330209,
        0.704353465,
        0.644360749,
        0.546217712,
        0.480802331,
        0.404445786,
        0.306302749,
        0.208151726,
        0.110008689,
        0.00639507,
        -0.097202576,
        -0.162641916,
        -0.282611376,
        -0.402564864,
        -0.495253292,
        -0.566139255,
        -0.642487813,
        -0.69702592,
        -0.784291683,
        -0.860680173,
        -0.92616743,
        -0.958895086,
        -0.975298846,
        -0.975370722,
        -0.959086756,
        -0.915577661,
        -0.850298047,
        -0.774133172,
        -0.697960312,
        -0.665344463,
        -0.63819121,
        -0.621923217,
        -0.649236195,
        -0.681995796,
        -0.752881759,
        -0.807395907,
        -0.856439474,
        -0.894573822,
        -0.861918042,
        -0.818377002,
        -0.76393473,
        -0.704053821,
        -0.627833043,
        -0.513485903,
        -0.404569413,
        -0.295652924,
        -0.192183058,
        -0.088713191,
        0.009310052,
        0.14002102,
    ],
    "cD": [
        -0.00500206,
        -6.71417e-05,
        -0.005136344,
        -0.000335709,
        0.004518639,
        0.014321334,
        0.024150886,
        0.039063067,
        0.098859506,
        0.198457571,
        0.298095921,
        0.412633025,
        0.482339585,
        0.577043016,
        0.62679505,
        0.696555323,
        0.836022157,
        0.935620222,
        1.095041582,
        1.20453634,
        1.408788244,
        1.503317107,
        1.592857339,
        1.647463721,
        1.736990524,
        1.751687852,
        1.776335587,
        1.791032915,
        1.785789144,
        1.750573299,
        1.66052279,
        1.600390643,
        1.530267804,
        1.470149085,
        1.370081028,
        1.245096669,
        1.105146415,
        0.930275739,
        0.740439168,
        0.575572612,
        0.390697816,
        0.210543085,
    ],
    "tip_loss": False,
    "hub_loss": False,
    "new_tip_loss": False,
    "new_hub_loss": False,
    "cascade_correction": False,
    "rotational_augmentation_correction": False,
    "rotational_augmentation_correction_method": 1,
    "max_iterations": 100.0,
    "convergence_limit": 0.001,
    "rho": 1.225,
    "method": 10,
    "linspace_interp": False,
    "num_interp": 25.0,
    "v_min": 3.0,
    "v_max": 20.0,
    "v_num": 10.0,
    "rpm_min": 100.0,
    "rpm_max": 3000.0,
    "rpm_num": 10.0,
    "relaxation_factor": 0.3,
    "print_all": False,
    "print_out": False,
    "target_rpm": 1500.0,
    "target_speed": 10.0,
    "delta_start": 5.0,
    "decrease_factor": 0.1,
    "min_delta": 0.001,
    "min_add_angle": -45.0,
    "max_add_angle": 45.0,
    "angle_step": 0.5,
    "dr": array(
        [
            0.,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.01,
        ]
    ),
    "r_in": array(
        [
            0.15,
            0.18,
            0.21,
            0.24,
            0.27,
            0.3,
            0.33,
            0.36,
            0.39,
            0.42,
            0.45,
            0.48,
            0.51,
            0.54,
            0.57,
            0.6,
            0.63,
            0.66,
            0.69,
            0.72,
            0.75,
            0.76,
        ]
    ),
    "c_in": array(
        [
            0.12429,
            0.14448,
            0.15224,
            0.14968,
            0.14213,
            0.13129,
            0.12198,
            0.1136,
            0.10666,
            0.10015,
            0.09487,
            0.08959,
            0.08441,
            0.08029,
            0.07728,
            0.07371,
            0.07097,
            0.06797,
            0.06637,
            0.06337,
            0.06072,
            0.05984,
        ]
    ),
    "theta_in": array(
        [
            29.24,
            33.41,
            34.9,
            34.6,
            32.39,
            20.96,
            21.71,
            13.33,
            11.04,
            8.18,
            7.05,
            5.88,
            3.41,
            1.09,
            0.57,
            -0.2,
            -0.66,
            -2.05,
            -2.29,
            -2.6,
            -2.92,
            -2.37,
        ]
    ),
}
