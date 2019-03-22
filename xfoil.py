#!/usr/bin/env python3
import os
import time
from subprocess import Popen, PIPE, DEVNULL, STDOUT
import re
from threading import Timer
import sys

xfoil_path = "xfoil.exe"
if sys.platform.startswith("darwin"):
    xfoil_path = "xfoil"

def get_coefficients_from_output(output_str):
    lines = output_str.splitlines()
    if "VISCAL:  Convergence failed" in output_str:
        return False
    if "TRCHEK2: N2 convergence failed." in output_str:
        return False
    i=-1
    cur_line = 0
    _,__,a,CL,Cm,CD,CDf,CDp = [None]*8
    for l in lines:
        a = re.match("(.+)a =(.+)CL =(.+)",l)
        _,a,CL = a.groups() if a != None else [_,a,CL]
        b = re.match("(.+)Cm =(.+)CD =(.+)=>.+CDf =(.+)CDp =(.+)",l)
        __,Cm,CD,CDf,CDp = b.groups() if b != None else [__,Cm,CD,CDf,CDp]
    #print(output_str)
    out = {
    "CL":float(CL),
    "CD":float(CD),
    "Cm":float(Cm),
    "CDf":float(CDf),
    "CDp":float(CDp),
    "out":output_str
    }
    return out

def run_xfoil_analysis(airfoil,reynolds,alpha,iterations=100,print_output = False):
    #print("running xfoil for %s,Re=%s,alpha=%s" % (airfoil,reynolds,alpha))
    with Popen(os.path.abspath(xfoil_path), stdin=PIPE, stdout=PIPE,
               universal_newlines=True) as process:

        def call(_str,process=process):
            #print(_str,file=process.stdin)
            process.stdin.write(_str+"\n")

        call("plop")
        call("G F")
        call("")

        if "naca" in airfoil.lower():
            #treat as standard naca airfoil
            call(airfoil)
        else:
            #open .dat file
            call("load %s" % airfoil)
        call("oper")
        call("re")
        #rint("settings reynolds")
        call("%s" % reynolds)
        call("v")
        call("iteration")
        call("%s" % iterations)
        call("alfa")
        call("%s" % alpha)
        call("")
        call("quit")
        #print(process.stdout.read())
        timer = Timer(0.5, process.kill)
        try:
            timer.start()
            output = process.communicate()[0]
            if print_output:
                for l in output.splitlines():
                    print(l)
                #print(output)
        finally:
            timer.cancel()
            #print("")
            #return False

    return get_coefficients_from_output(output)

#print(run_xfoil_analysis("airfoil_name.dat",26820,0.82))