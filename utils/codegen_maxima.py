#!/usr/bin/env python2
#----------------------------------
# generate dft(gga) code for dirac
# Andreas Hesselmann <andreas@theochem.uni-duesseldorf.de>, 2003
#----------------------------------
# needs: maxima 
# usage: codegen.py 'file.max'
# 'file.max' is a  maxima input file of the following general form:
#      ....some variables....
#      ....some functions...
#      K(rhoa,ghroa,rhob,gradb):= 'expression of rhoa,grhoa,xa
#                                                rhob,grhob,xb'
# where xa=sqrt(grada*grada)/rhoa^(4/3)
# output: fun-'file'.c

import sys, string, os, re;

funcfil=sys.argv[1];
funcname=string.split(funcfil,'.')[0]

#header
f=open(funcfil,'r')
func=f.readlines()
f.close()
f=open("fun-"+funcname+".c",'w')
f.write("/* Automatically generated functional code: "+funcname+"\n")
f.write("   Maxima input:\n")
for line in func:
    line=string.replace(line,"/*","//")
    line=string.replace(line,"*/","")
    f.write(4*" "+">> "+line)
f.write("*/\n\n")
f.write("// add \"extern Functional "+funcname+"Functional;\" to 'functionals.h'\n")
f.write("// add \"&"+funcname+"Functional,\" to 'functionals.c'\n\n")
f.write("#include <math.h>\n")
f.write("#include <stddef.h>\n\n")
f.write("#define __CVERSION__\n\n")
f.write("#include \"functionals.h\"\n")
f.write("#include \"dcbham.h\"\n\n")
f.write("#define LOG log\n")
f.write("#define ABS fabs\n")
f.write("#define ASINH asinh\n")
f.write("#define SQRT sqrt\n\n")
f.write("/* INTERFACE PART */\n")
f.write("static int "+funcname+"_read(const char* conf_line);\n")
f.write("static real "+funcname+"_energy(real rhoa, real rhob, real grada, real gradb);\n")
f.write("static void "+funcname+"_first(FirstDrv *ds, real factor, \n\
                        real rhoa, real rhob, real ngrada, real ngradb);\n")
f.write("static void "+funcname+"_second(SecondFuncDrv *ds, real factor,\n\
                        real rhoa, real rhob, real ngrada, real ngradb);\n")
f.write("static void "+funcname+"_third(ThirdFuncDrv *ds, real factor,\n\
                        real rhoa, real rhob, real ngrada, real ngradb);\n\n")
f.write("Functional "+funcname+"Functional = {\n")
f.write("  \""+funcname+"\",\n")
f.write("  fun_true,\n")
f.write("  "+funcname+"_read,\n")
f.write("  NULL,\n")
f.write("  "+funcname+"_energy,\n")
f.write("  "+funcname+"_first,\n")
f.write("  "+funcname+"_second,\n")
f.write("  "+funcname+"_third\n};\n\n")
f.write("/* IMPLEMENTATION PART */\n")
f.write("static int\n"+funcname+"_read(const char* conf_line)\n")
f.write("{\n#ifndef NO_BACKWARD_COMP\n  dcrham_.hfxfac = 0;\n")
f.write("#endif\n  return 1;\n}\n\n")
f.close()

#maxima commands
com1="string(float(subst(pow,\"^\",optimize("
com2="))));\n"
maxima_commands={
    "en": com1+"zk=K(rhoa,grada,rhob,gradb)"+com2,
    "first": com1+"[dfdra=diff(K(rhoa,grada,rhob,gradb),rhoa),"\
                 +" dfdga=diff(K(rhoa,grada,rhob,gradb),grada)]"+com2,
    "second":com1+"[dfdra=diff(K(rhoa,grada,rhob,gradb),rhoa),"\
                 +" dfdga=diff(K(rhoa,grada,rhob,gradb),grada),"\
                 +" d2fdraga=diff(K(rhoa,grada,rhob,gradb),rhoa,1,grada,1),"\
                 +" d2fdrara=diff(K(rhoa,grada,rhob,gradb),rhoa,2),"\
                 +" d2fdgaga=diff(K(rhoa,grada,rhob,gradb),grada,2)]"+com2,
    "third": com1+"[dfdra=diff(K(rhoa,grada,rhob,gradb),rhoa),"\
                 +" dfdga=diff(K(rhoa,grada,rhob,gradb),grada),"\
                 +" d2fdraga=diff(K(rhoa,grada,rhob,gradb),rhoa,1,grada,1),"\
                 +" d2fdrara=diff(K(rhoa,grada,rhob,gradb),rhoa,2),"\
                 +" d2fdgaga=diff(K(rhoa,grada,rhob,gradb),grada,2),"\
                 +" d3fdraraga=diff(K(rhoa,grada,rhob,gradb),rhoa,2,grada,1),"\
                 +" d3fdragaga=diff(K(rhoa,grada,rhob,gradb),rhoa,1,grada,2),"\
                 +" d3fdrarara=diff(K(rhoa,grada,rhob,gradb),rhoa,3),"\
                 +" d3fdgagaga=diff(K(rhoa,grada,rhob,gradb),grada,3)]"+com2
    }
                   

def make_maximafile(funcfil,com):
    f=open(funcfil,'r')
    func=f.readlines()
    f.close()
    f=open('maxima.in','w')
    f.write("xa:sqrt(grada*grada)/rhoa^(4/3);\n");
    f.write("xb:sqrt(gradb*gradb)/rhob^(4/3);\n");
    for line in func:
        f.write(line)
    f.write(maxima_commands[com])
    f.close()

def run_maxima(outfil):
    if os.system("maxima < maxima.in >"+outfil)!=0:
        sys.exit("Error in run_maxima!")

def vars_one2field(s):
    f=0
    s=list(s)
    for i in range(len(s)):
        if s[i]=="%":
            s[i]="t["
            f=1
        elif f==1 and not s[i] in ["0","1","2","3","4","5","6","7","8","9"]:
            s[i]="]"+s[i]
            f=0
    if f==1:
        s.append("]")
    string=""
    for c in s:
        string=string+c
    return string

def f2c(expr):
    s=string.split(expr,",")
    for i in range(len(s)):
        if s[i][0]=="%" and s[i][-1]=="]":
            nvars=int(s[i][1:-1])
            s=s[i+1:]
            break
        if s[i][-4:]=="[%1]":
            nvars=1
            s=s[i+1:]
            break
    print nvars
    s2=[]
    i=-1
    for s1 in s:
        if ":" in s1 or "=" in s1:
            s2.append(s1)
            i=i+1
        else:
            s2[i]=s2[i]+","+s1
    s=s2
    var={}
    for i in range(nvars):
        var["%"+str(i+1)]="t["+str(i+1)+"]"
    for i in range(len(s)):
        s[i]=vars_one2field(s[i]) #replace %... by t[...]
        if s[i][0:2]=="zk":
            s[i]=s[i][:-1]
        if s[i][0]=="[":
            s[i]=s[i][1:]
        s[i]=string.replace(s[i],":"," = ")+";"
    if s[-1][0:2] =="zk":
        s[-1]=s[-1][:-1]+";"
    else:
        s[-1]=s[-1][:-3]+";"
    return nvars,s
        
def read_maximaout(outfil):
    f=open(outfil,'r')
    block=0
    for line in f.readlines():
        data=string.split(line)
        if len(data)==0: continue
        if data[0]=="\n": continue
        if len(data)>=2:
            if len(data[1])>=5:
                if data[1][0:5]=="BLOCK":
                    block=1
                    expr=""
        if block==1:
            if data[0][0:2]!="(C":
                line=string.strip(line)
                line=string.replace(line,"#","")
                line=string.replace(line,"\n","")
                expr=expr+line
            else:
                block=0
    nvars,cexpr=f2c(expr)
    return nvars,cexpr


def make_energycode(nvars,code):
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic real\n")
    f.write(funcname+"_energy(real rhoa, real rhob, real grada, real gradb)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(nvars+1)+"],zk;\n")
    for i in range(len(code)):
        f.write(4*" "+code[i]+"\n")
    f.write(4*" "+"return zk;\n")
    f.write("}\n")
    f.close()

def make_firstcode(nvars,code):
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic void\n")
    f.write(funcname+"_first(FirstDrv *ds, real factor, real rhoa, real rhob, real grada, real gradb)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(nvars+1)+"];\n")
    f.write(4*" "+"real dfdra, dfdga;\n")
    for i in range(len(code)):
        f.write(4*" "+code[i]+"\n")
    f.write(4*" "+"ds->df1000 += factor*dfdra;\n")
    f.write(4*" "+"ds->df0010 += factor*dfdga;\n")
    f.write("}\n")
    f.close()    

def make_secondcode(nvars,code):
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic void\n")
    f.write(funcname+"_second(SecondFuncDrv *ds, real factor, real rhoa, real rhob, real grada, real gradb)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(nvars+1)+"];\n")
    f.write(4*" "+"real dfdra, dfdga;\n")
    f.write(4*" "+"real d2fdraga, d2fdrara, d2fdgaga;\n")    
    for i in range(len(code)):
        f.write(4*" "+code[i]+"\n")
    f.write(4*" "+"ds->df1000 += factor*dfdra;\n")
    f.write(4*" "+"ds->df0010 += factor*dfdga;\n")
    f.write(4*" "+"ds->df1010 += factor*d2fdraga;\n")
    f.write(4*" "+"ds->df2000 += factor*d2fdrara;\n")
    f.write(4*" "+"ds->df0020 += factor*d2fdgaga;\n")
    f.write("}\n")
    f.close()    

def make_thirdcode(nvars,code):
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic void\n")
    f.write(funcname+"_third(ThirdFuncDrv *ds, real factor, real rhoa, real rhob, real grada, real gradb)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(nvars+1)+"];\n")
    f.write(4*" "+"real dfdra, dfdga;\n")
    f.write(4*" "+"real d2fdraga, d2fdrara, d2fdgaga;\n")
    f.write(4*" "+"real d3fdraraga, d3fdragaga, d3fdrarara, d3fdgagaga;\n")  
    for i in range(len(code)):
        f.write(4*" "+code[i]+"\n")
    f.write(4*" "+"ds->df1000 += factor*dfdra;\n")
    f.write(4*" "+"ds->df0010 += factor*dfdga;\n")
    f.write(4*" "+"ds->df1010 += factor*d2fdraga;\n")
    f.write(4*" "+"ds->df2000 += factor*d2fdrara;\n")
    f.write(4*" "+"ds->df0020 += factor*d2fdgaga;\n")
    f.write(4*" "+"ds->df2010 += factor*d3fdraraga;\n")
    f.write(4*" "+"ds->df1020 += factor*d3fdragaga;\n")
    f.write(4*" "+"ds->df3000 += factor*d3fdrarara;\n")
    f.write(4*" "+"ds->df0030 += factor*d3fdgagaga;\n")
    f.write("}\n")
    f.close()    
            
            

#---------------------------------------------------------------------
parts=['en','first','second','third']
for p in parts:
    make_maximafile(funcfil,p)
    run_maxima(p+'.out')
    nvars,code=read_maximaout(p+'.out')
    if p=='en':
        make_energycode(nvars,code)
    elif p=='first':
        make_firstcode(nvars,code)
    elif p=='second':
        make_secondcode(nvars,code)
    elif p=='third':
        make_thirdcode(nvars,code)
    os.remove(p+'.out')

