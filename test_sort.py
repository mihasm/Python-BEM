from scipy import optimize

def func(x):
    return x**2

a = optimize.fixed_point(func,50)
print(a)