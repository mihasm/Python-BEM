from figdig import *
B = (246, 1995, 2170, 3030)
runs = auto_seeds(9, B, 2921, 260, 2160)
assert len(runs) == 7, runs
seeds = {a: (x, 2921, 58, 12) for a, (x, w) in zip([15, 20, 25, 30, 35, 40, 45], runs)}
guides = {35: [(254, 2219), (330, 2210), (420, 2207), (500, 2210), (560, 2232), (620, 2262), (690, 2291), (790, 2318), (830, 2325)]}
digitize('f8_ct_2b', 9, B, (0.0, 0.1), (0.0, 0.01), seeds, (240, 1990, 2190, 3040), 3023.5, (253, 2161), guides=guides, tkw=dict(min_score=0.04, max_gap=16), dkw=dict(atten=0.6))
