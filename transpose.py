a = [["a",1,2],
	 ["b",3,4]]

o = [["a","b"],
	 [1,3],
	 [2,4]]

def transpose(a):
	o=[]
	for i in range(len(a[0])):
		o.append([])
		for r in a:
			o[i].append(r[i])
	return o

