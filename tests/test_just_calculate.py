import sys, os
sys.path.append("..")
from calculation_runner import calculate_power_3d
import json

results = calculate_power_3d(json.loads(open(os.path.join("..","karlsen.bem")).read()))