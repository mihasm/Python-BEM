from figdig import *
B = (295, 330, 2290, 1060)
xs = [849, 1005, 1180, 1347, 1540, 1781, 2077]
seeds = {a: (x, 912, -92, 12) for a, x in zip([15, 20, 25, 30, 35, 40, 45], xs)}
digitize('f5_eta_2b', 8, B, (0.0, 0.1, 302.0), (0.0, 0.1, 1046.0), seeds, (300, 335, 2285, 1050), 1046, (302, 2279),
         tkw=dict(min_score=0.04, max_gap=16), dkw=dict(atten=0.6), extend=False)
