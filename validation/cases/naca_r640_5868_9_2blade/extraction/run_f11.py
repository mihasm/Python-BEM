from figdig import *
B = (290, 300, 2210, 2275)
xs = [651, 776, 909, 1049, 1200, 1406, 1637]
seeds = {a: (x, 2160, -110, 12) for a, x in zip([15, 20, 25, 30, 35, 40, 45], xs)}
# 45 deg: tracker stops at the '.10' C_T label (x~1296); upper part follows hand waypoints read on f11_top.png
guides = {45: [(1300, 1070), (1210, 860), (1150, 760), (1100, 675), (1050, 640), (1000, 635), (925, 612), (850, 580),
               (790, 522), (745, 475), (700, 432), (650, 385), (600, 345), (550, 325), (500, 318), (450, 315)]}
digitize('f11_cp_3b', 11, B, (0.0, 0.1), (0.0, 0.01), seeds, (295, 305, 2205, 2270), 2262, (301, 2195),
         tkw=dict(min_score=0.04, max_gap=16), dkw=dict(atten=0.6), guides=guides, tol=4)
