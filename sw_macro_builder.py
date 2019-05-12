import os
import sys
from turbine_data import SET_INIT

def create_macro_text(list_of_files):
	template_start = """
Dim swApp As Object

Dim part As Object
Dim boolstatus As Boolean
Dim longstatus As Long, longwarnings As Long

Sub main()

Set swApp = Application.SldWorks

template = swApp.GetUserPreferenceStringValue(swDefaultTemplatepart)
Set part = swApp.NewDocument(template, 0, 0, 0)

Dim myModelView As Object
"""

	template_end = """
End Sub
"""
	str_between = ""
	for f in list_of_files:
		str_between += 'boolstatus = part.InsertCurveFile("%s")\n' % f
	return template_start+"\n"+str_between+"\n"+template_end

import win32com.client
import pythoncom
import pywintypes
import os
import win32job

print(dir(pywintypes))

cwd = os.getcwd()

#swYearLastDigit = 5
#sw = win32com.client.DispatchEx("SldWorks.Application")  # e.g. 20 is SW2012,  23 is SW2015

#sw.Visible = 1
#part = sw.Newpart


foil_x = SET_INIT["airfoils"]["s826"]["x"]
foil_y = SET_INIT["airfoils"]["s826"]["y"]
z = [0]*len(foil_x)
list_of_coordinates = zip(foil_x,foil_y,z)



class SolidworksApp:
	def __init__(self,sw_com_name="SldWorks.Application",Visible=True):
		self.sw = win32com.client.DispatchEx(sw_com_name)
		self.sw.Visible=Visible

	def create_part(self):
		self.part = self.sw.Newpart

	def create_curve(self,list_of_coordinates):
		self.part.InsertCurveFileBegin()
		for x,y,z in list_of_coordinates:
			self.part.InsertCurveFilePoint(x,y,z)
		self.part.InsertCurveFileEnd

	def create_point(self,x,y,z):
		self.part.CreatePoint(1,0,0)
		self.part.CreatePoint

	def select_by_id(self,name,_type,x=0,y=0,z=0):
		arg1 = win32com.client.VARIANT(pythoncom.VT_DISPATCH, None)
		self.part.Extension.SelectByID2(name, _type, 0, 0, 0, False, 0, arg1, 0)

	def insert_sketch(self):
		self.select_by_id("Front Plane","PLANE")
		self.part.SketchManager.InsertSketch(True)
		


app = SolidworksApp()
app.create_part()
app.create_curve(list_of_coordinates)
app.insert_sketch()
app.create_point(1,0,0)
