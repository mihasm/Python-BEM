import requests
from bs4 import BeautifulSoup as bs

def get_polars(link="http://airfoiltools.com/airfoil/details?airfoil=s826-nr"):
	results = {}

	r = requests.get(link)
	txt = r.text
	soup = bs(txt,features='html.parser')
	#print(txt)
	table_of_polars = soup.find_all("table",{"class":"polar"})[0]
	rows = table_of_polars.find_all("tr")[1:-1]
	links = []
	for r in rows:
		data = r.find_all("td")
		checkbox,name,reynolds,ncrit,maxclcd,description,source,details = [d for d in data]
		details_link = "http://airfoiltools.com"+details.find_all("a")[0].get('href')
		links.append(details_link)
		#print(details_link)
		#print(name,reynolds,details)

	csv_links = []

	for l in links:
		r = requests.get(l)
		txt = r.text
		soup = bs(txt,features='html.parser')
		details_table = soup.find_all("table",{"class":"details"})[0]
		td1 = details_table.find_all("td",{"class":"cell1"})
		#print(len(details_table))
		links = details_table.find_all("a")
		#print(links)
		for l in links:
			#print(l.get("href"))
			if "csv?polar" in l.get("href"):
				csv_links.append("http://airfoiltools.com"+l.get("href"))
				break
	
	alpha_dict = {}
	for l in csv_links:
		r = requests.get(l)
		text = r.text
		lines = text.splitlines()
		start = lines.index("Alpha,Cl,Cd,Cdp,Cm,Top_Xtr,Bot_Xtr")
		reynolds_number = [float(i.split(",")[1]) for i in lines if "Reynolds number," in i][0]
		ncrit = [float(i.split(",")[1]) for i in lines if "Ncrit," in i][0]
		mach = [float(i.split(",")[1]) for i in lines if "Mach," in i][0]
		#print(reynolds_number)
		for line_num in range(start+1,len(lines)):
			alpha,cl,cd = lines[line_num].split(",")[:3]
			alpha,cl,cd = float(alpha),float(cl),float(cd)
			alpha_dict[alpha] = {"cl":cl,"cd":cd}
		results[reynolds_number] = {ncrit:alpha_dict}
	return results

get_polars()