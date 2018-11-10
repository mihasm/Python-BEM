__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.3.0"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

from numpy import array

SET_INIT = {
    "Rhub": 0.006,
    "R": 0.45,
    "B": 3,
    "turbine_name": "norvezan_karlsen_ntnu_njegove_krivulje",
    "r": array(
        [
            0.0075,
            0.0225,
            0.049,
            0.055,
            0.0675,
            0.0825,
            0.0975,
            0.1125,
            0.1275,
            0.1425,
            0.1575,
            0.1725,
            0.1875,
            0.2025,
            0.2175,
            0.2325,
            0.2475,
            0.2775,
            0.2925,
            0.3075,
            0.3225,
            0.3375,
            0.3675,
            0.3825,
            0.3975,
            0.4125,
            0.4275,
            0.4425,
        ]
    ),
    "c": array(
        [
            0.0135,
            0.0135,
            0.0135,
            0.0495,
            0.08143281,
            0.08011114,
            0.07701157,
            0.07312555,
            0.06900783,
            0.064952,
            0.06110231,
            0.05752027,
            0.0542229,
            0.05120426,
            0.04844738,
            0.04593074,
            0.04363179,
            0.03960108,
            0.03783062,
            0.03620078,
            0.03469692,
            0.03201685,
            0.03081913,
            0.02970404,
            0.02866371,
            0.0276912,
            0.02678034,
            0.02592566,
        ]
    ),
    "theta": array(
        [
            1.200000e+02,
            1.200000e+02,
            1.200000e+02,
            3.800000e+01,
            3.705499e+01,
            3.254410e+01,
            2.867674e+01,
            2.526190e+01,
            2.243033e+01,
            1.998782e+01,
            1.803445e+01,
            1.634885e+01,
            1.466325e+01,
            1.306715e+01,
            1.182907e+01,
            1.075265e+01,
            8.882724e+00,
            7.987743e+00,
            7.252693e+00,
            6.564977e+00,
            5.918680e+00,
            4.718537e+00,
            4.131617e+00,
            3.543894e+00,
            2.943336e+00,
            2.218509e+00,
            1.097028e-01,
            -7.167448e-01,
        ]
    ),
    "foils": [
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
    ],
    "tip_loss": False,
    "hub_loss": False,
    "new_tip_loss": False,
    "new_hub_loss": False,
    "cascade_correction": False,
    "rotational_augmentation_correction": False,
    "rotational_augmentation_correction_method": 0,
    "max_iterations": 100.0,
    "convergence_limit": 0.001,
    "rho": 1.225,
    "method": 7,
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
    "curves": {
        "s826": {
            "max_thickness":0.144027,
            "AoA_cL": [
                -4.7462686567,
                -3.3987347215,
                -1.2626497144,
                -0.7286284626,
                0.4626497144,
                0.791278177,
                1.5717707758,
                1.8779436153,
                2.2906455377,
                4.2624163135,
                5.6385480007,
                7.3227688717,
                7.9663329443,
                9.0480683005,
                9.5820895522,
                12.3343529267,
                13.5223880597,
                14.4704379338,
                16.565444383,
                18.5372151588,
                20.1803574719,
                21.9018242123,
                25.1097844113,
                28.8068546158,
            ],
            "AoA_cD": [
                -9.9402985075,
                -9.5564400221,
                -9.2072722806,
                -8.0365333825,
                -7.2355015048,
                -5.9209876543,
                -5.4691235182,
                -4.633859509,
                -4.1956882255,
                -3.1481850009,
                -2.7168601437,
                -2.2239174498,
                -1.3475748828,
                -0.9094035993,
                0.2407960199,
                0.8980529452,
                1.8428597752,
                2.6826880686,
                3.1162950679,
                3.9173269455,
                4.9032123334,
                6.2998832996,
                7.141993735,
                7.8403292181,
                8.6208218168,
                9.8531785517,
                10.5515140348,
                11.1060745654,
                11.6606350961,
                12.3589705792,
                12.8792989784,
                14.2075056815,
                14.9058411645,
                15.8506479946,
                16.7905247425,
                17.6991830969,
                17.9905656901,
                18.2126650697,
                19.8352681039,
                20.4309071924,
                21.6427246484,
                22.690227873,
                23.3269455193,
                24.2443666441,
                24.9700878324,
                26.038130336,
                26.599537293,
                27.3526441865,
                27.8455868804,
                28.3796081322,
                28.7903937105,
                29.6941219827,
                30.0597014925,
            ],
            "cL": [
                0.0044091711,
                0.1559866745,
                0.39812456,
                0.4546744471,
                0.5911954478,
                0.627601048,
                0.6913108484,
                0.7252487643,
                0.7656389488,
                0.9461500497,
                1.0644682503,
                1.1645836509,
                1.2151469845,
                1.2619483093,
                1.2829018515,
                1.3625899072,
                1.3800705467,
                1.375432752,
                1.3254149701,
                1.2401455934,
                1.1731189351,
                0.9092182522,
                0.9265718786,
                0.9520761208,
            ],
            "cD": [
                0.0940140845,
                0.0849443575,
                0.0773465484,
                0.0555077378,
                0.0427882107,
                0.0297052686,
                0.0282516084,
                0.02501478,
                0.0234933055,
                0.022209181,
                0.0215840723,
                0.0208911494,
                0.0214290268,
                0.0216518866,
                0.0225725961,
                0.0227664754,
                0.0233672405,
                0.0250890666,
                0.0246852721,
                0.0264345331,
                0.0272291775,
                0.0312728221,
                0.0327482177,
                0.0344053208,
                0.0369735698,
                0.0449687011,
                0.0471491914,
                0.0518735872,
                0.0572569988,
                0.0649565293,
                0.0689540949,
                0.0856711876,
                0.0980272996,
                0.1169248826,
                0.1428484901,
                0.1710737263,
                0.1830664232,
                0.1917883846,
                0.2557494349,
                0.269195792,
                0.28772996,
                0.3080812033,
                0.3182568249,
                0.3352161943,
                0.3512064568,
                0.3727690836,
                0.3812487683,
                0.3974813076,
                0.4003886281,
                0.4156520605,
                0.425100852,
                0.4410911146,
                0.448943662,
            ],
        }
    },
    "dr": array(
        [
            0.,
            0.015,
            0.0265,
            0.006,
            0.0125,
            0.015,
            0.015,
            0.015,
            0.015,
            0.015,
            0.015,
            0.015,
            0.015,
            0.015,
            0.015,
            0.015,
            0.015,
            0.03,
            0.015,
            0.015,
            0.015,
            0.015,
            0.03,
            0.015,
            0.015,
            0.015,
            0.015,
            0.015,
        ]
    ),
    "r_in": array(
        [
            0.0075,
            0.0225,
            0.049,
            0.055,
            0.0675,
            0.0825,
            0.0975,
            0.1125,
            0.1275,
            0.1425,
            0.1575,
            0.1725,
            0.1875,
            0.2025,
            0.2175,
            0.2325,
            0.2475,
            0.2775,
            0.2925,
            0.3075,
            0.3225,
            0.3375,
            0.3675,
            0.3825,
            0.3975,
            0.4125,
            0.4275,
            0.4425,
        ]
    ),
    "c_in": array(
        [
            0.0135,
            0.0135,
            0.0135,
            0.0495,
            0.08143281,
            0.08011114,
            0.07701157,
            0.07312555,
            0.06900783,
            0.064952,
            0.06110231,
            0.05752027,
            0.0542229,
            0.05120426,
            0.04844738,
            0.04593074,
            0.04363179,
            0.03960108,
            0.03783062,
            0.03620078,
            0.03469692,
            0.03201685,
            0.03081913,
            0.02970404,
            0.02866371,
            0.0276912,
            0.02678034,
            0.02592566,
        ]
    ),
    "theta_in": array(
        [
            1.200000e+02,
            1.200000e+02,
            1.200000e+02,
            3.800000e+01,
            3.705499e+01,
            3.254410e+01,
            2.867674e+01,
            2.526190e+01,
            2.243033e+01,
            1.998782e+01,
            1.803445e+01,
            1.634885e+01,
            1.466325e+01,
            1.306715e+01,
            1.182907e+01,
            1.075265e+01,
            8.882724e+00,
            7.987743e+00,
            7.252693e+00,
            6.564977e+00,
            5.918680e+00,
            4.718537e+00,
            4.131617e+00,
            3.543894e+00,
            2.943336e+00,
            2.218509e+00,
            1.097028e-01,
            -7.167448e-01,
        ]
    ),
    "foils_in": [
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
        "s826",
    ],
}