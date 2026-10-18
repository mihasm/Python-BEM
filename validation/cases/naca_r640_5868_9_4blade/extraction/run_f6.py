from figdig import *
B = (390, 1640, 2200, 3005)
# seeds: dark-line peaks on row y=2910 whose up-left trace ends on the CP axis (selection printed in the log)
xs = [738, 871, 1009, 1145, 1300, 1519, 1776]
seeds = {a: (x, 2910, -112, 12) for a, x in zip([15, 20, 25, 30, 35, 40, 45], xs)}
digitize('f6_cp_2b', 8, B, (0.0, 0.1), (0.0, 0.01), seeds, (390, 1600, 2195, 3000), 2998.5, (398, 2193),
         vthr=238, norm=True, tkw=dict(min_score=0.04, max_gap=16), dkw=dict(atten=0.6))
