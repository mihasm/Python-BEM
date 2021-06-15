import scipy.optimize as optimize
bounds = [(-5,5),(-5,5)]
def func(x):
	return x[0]**2+x[1]**2
result = optimize.fsolve(func,(0.5,0.5))
print(result)