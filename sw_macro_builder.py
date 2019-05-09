import os
import sys

def create_macro_text(list_of_files):
	template_start = """
Dim swApp As Object

Dim Part As Object
Dim boolstatus As Boolean
Dim longstatus As Long, longwarnings As Long

Sub main()

Set swApp = Application.SldWorks

template = swApp.GetUserPreferenceStringValue(swDefaultTemplatePart)
Set Part = swApp.NewDocument(template, 0, 0, 0)

Dim myModelView As Object
"""

	template_end = """
End Sub
"""
	str_between = ""
	for f in list_of_files:
		str_between += 'boolstatus = Part.InsertCurveFile("%s")\n' % f
	return template_start+"\n"+str_between+"\n"+template_end

