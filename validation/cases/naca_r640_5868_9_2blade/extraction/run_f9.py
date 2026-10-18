from figdig import *
B = (285, 310, 2270, 1610)
# seeds: dark-line peaks on row y=1330 whose up-left trace runs to the CT axis (see probe output in log)
xs = [718, 873, 1035, 1208, 1395, 1611, 1916]
seeds = {a: (x, 1330, -115, 12) for a, x in zip([15, 20, 25, 30, 35, 40, 45], xs)}
digitize('f9_ct_3b', 10, B, (0.0, 0.1), (0.0, 0.01), seeds, (295, 320, 2258, 1606), 1604, (293, 2258),
         tkw=dict(min_score=0.04, max_gap=16), dkw=dict(atten=0.6))
