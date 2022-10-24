import math
import cmath
import numpy as np
from tkinter import *
from tkinter import messagebox
import statistics
import matplotlib.pyplot as plt


ground=Tk()
#images
ground.geometry("200x200")
im_mean=PhotoImage(file='mean.png')
im_variance=PhotoImage(file='variance.png')
im=PhotoImage(file='statr.png')
im_logo=PhotoImage(file='logo.png')
im_standard=PhotoImage(file='standard.png')
im_n=PhotoImage(file='n.png')
form=PhotoImage(file="form.png")
back3=PhotoImage(file='back2.png')
back=PhotoImage(file='back3.png')
border=PhotoImage(file='border.png')

#********************************************************************************
#solver Equation
def openSolver():
    #desgin
    top=Toplevel()
    top.geometry('770x575')
    top.title('Equation solver by Doha Emad')
    label_back3=Label(top, image=back3).place(x=0,y=0)
    label_form=Label(top, image=form).place(x=100,y=100)
    Button(top,text='Home',font=30,command=top.destroy,fg='navy',bg='white',bd=4).place(x=670,y=520)
    title=Label(top, text="Equation solver",font="Times 30 italic bold",bg='white',fg='navy').place(x=250, y=15)
    a_label=Label(top, text="a",font="Times 25 italic bold",bg='white',fg='navy').place(x=50, y=190)
    a=Entry(top,width=10,font=("courier New",20,'bold'),bg='white',bd=2)
    a.place(x =80, y=195)
    b_label= Label(top, text="b",font="Times 25 italic bold",bg='white',fg='navy').place(x=50, y=240)
    b1=Entry(top,width=10,font=("courier New",20,'bold'),bg='white',bd=2)
    b1.place(x=80, y=250)
    c_label= Label(top, text="c",font="Times 25 italic bold",bg='white',fg='navy').place(x=50, y=300)
    c1=Entry(top,width=10,font=("courier New",20,'bold'),bg='white',bd=2)
    c1.place(x=80, y=305)
    d_label=Label(top, text="d",font="Times 25 italic bold",bg='white',fg='navy').place(x=50, y=350)
    d1=Entry(top,width=10,font=("courier New",20,'bold'),bg='white',bd=2)
    d1.place(x=80, y=355)
    e_label=Label(top, text="e",font="Times 25 italic bold",bg='white',fg='navy').place(x=50, y=400)
    e1=Entry(top,width=10,font=("courier New",20,'bold'),bg='white',bd=2)
    e1.place(x=80, y=405)
    label_sol=Label(top,image=border)
    label_sol.place(x=300, y=200)
    label_output=Label(top, text="Solution",font="Times 18 italic bold",bg='white',fg='navy')
    label_output.place(x=460, y=173)
    global full
    full=0
    #clear label
    def clear_label():
        if(a4==0 and a3==0 and a2==0 and a1==0 and a0==0):
            print("run")
            label_output1.place_forget()
        elif(a4==0 and a3==0 and a2==0 and a1==0):
            label_output1.place_forget()
        elif(a4==0 and a3==0 and a2==0):
            label_output1.place_forget()
        elif(a4==0 and a3==0):
            label_output1.place_forget()
            label_output2.place_forget()
        elif(a4==0):
            label_output1.place_forget()
            label_output2.place_forget()
            label_output3.place_forget()
            label_output4.place_forget()
        else:
            label_output1.place_forget()
            label_output2.place_forget()
            label_output3.place_forget()
            label_output4.place_forget()

    def solve():
        global full
        try:
            if(full==1):
                clear_label()
            #coefficient
            global a4
            global a3
            global a2
            global a1
            global a0
            global label_output1
            global label_output2
            global label_output3
            global label_output4
            a4=float(a.get())
            a3=float(b1.get())
            a2=float(c1.get())
            a1=float(d1.get())
            a0=float(e1.get())
            
            #no degree 
            if(a4==0 and a3==0 and a2==0 and a1==0 and a0==0):
                label_output1=Label(top, text="Trivial Equation",font="Times 20 italic bold",bg='white',fg='navy')
                label_output1.place(x=400, y=250)

            #degree 0
            elif(a4==0 and a3==0 and a2==0 and a1==0 ):
                label_output1=Label(top, text="No solution",font="Times 20 italic bold",bg='white',fg='navy')
                label_output1.place(x=400, y=250)

            #degree 1    
            elif(a4==0 and a3==0 and a2==0  ):
                label_output1=Label(top, text=" x1 =  "+str(round((-1*(a0)/a1),6)),font="Times 20 italic bold",
                bg='white',fg='navy')
                label_output1.place(x=400, y=250)

            #degree 2    
            elif(a4==0 and a3==0 ):
                def equation( a, b, c):
                    global label_output1
                    global label_output2
                    global label_output3
                    global label_output4 
                    dis = b * b - 4 * a * c 
                    sqrt_val = math.sqrt(abs(dis)) 
                    if dis > 0:
                        label_output1=Label(top, text="x1 =  "+str(round((-b + sqrt_val)/(2 * a), 8)),
                        font="Times 20 italic bold",bg='white',fg='navy')
                        label_output1.place(x=350, y=250)
                        label_output2=Label(top, text="x2 =  "+str(round((-b - sqrt_val)/(2 * a), 8)),
                        font="Times 20 italic bold",bg='white',fg='navy')
                        label_output2.place(x=350, y=300)

                    elif dis == 0: 
                        label_output1=Label(top, text="x1 =  "+str(round(-b / (2 * a), 8)),font="Times 20 italic bold",
                        bg='white',fg='navy')
                        label_output1.place(x=350, y=250)
                        label_output2=Label(top, text="x2 =  "+str(round(-b / (2 * a), 8)),font="Times 20 italic bold",
                        bg='white',fg='navy')
                        label_output2.place(x=350, y=300)
                
                    else:
                        label_output1=Label(top, text="x1 =  "+str(round(- b / (2 * a), 8))+' + i'+
                        str(round(sqrt_val/2*a, 8)),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output1.place(x=350, y=250)
                        label_output2=Label(top, text="x2 =  "+str(round(- b / (2 * a), 8))+' - i'+
                        str(round(sqrt_val/2*a, 8)),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output2.place(x=350, y=300) 
                equation(a2,a1,a0)
            

            #degree 3
            else:
                #Function of p coefficient ..... p *x^2
                def P(a3, a2):
                    p = a2/a3
                    return p

                #Functon of q coefficient ..... q*x
                def Q(a3, a1):
                    q = a1/a3
                    return q

                #Function of r coefficient ..... r
                def R(a3, a0):
                    r = a0/a3
                    return r

                #Function of a coefficient
                def a_cof(a3,a2,a1):
                    q = Q(a3, a1)
                    p = P(a3, a2)
                    a = q-(pow(p,2)/3.0)
                    return a

                #Function of b coefficient
                def b_cof(a3,a2, a1, a0):
                    q = Q(a3, a1)
                    p = P(a3, a2)
                    r = R(a3, a0)
                    b = r + (2/27)* pow(p, 3) - (1/3) * p * q
                    return b

                #Function of A.
                def A_function(a3,a2,a1,a0):
                    a = a_cof(a3, a2, a1)
                    b = b_cof(a3, a2, a1, a0)
                    A = np.cbrt(-b / 2.0 + math.sqrt(pow(b, 2) / 4.0 + pow(a, 3) / 27.0))
                    return A

                #Function of B.
                def B_function (a3,a2,a1,a0):
                    a = a_cof(a3, a2, a1)
                    b = b_cof(a3, a2, a1, a0)
                    B = np.cbrt(-b / 2.0 - (math.sqrt(pow(b, 2)/4.0 + pow(a, 3)/27.0)))
                    return B

                #The complex part of the solution.
                def complex_part(a3,a2,a1,a0):
                    A = A_function(a3, a2, a1, a0)
                    B = B_function(a3, a2, a1, a0)
                    complex = (math.sqrt(3) / 2) * (A - B)
                    return complex

                #The real part of the solution.
                def real_part(a3, a2, a1, a0):
                    A = A_function(a3, a2, a1, a0)
                    B = B_function(a3, a2, a1, a0)
                    p = P(a3, a2)
                    real= (-0.5 * (A + B)  - p / 3.0)
                    return real

                #The first solution.
                def solution_x1(a3,a2,a1,a0):
                    A = A_function(a3, a2, a1, a0)
                    B =  B_function(a3, a2, a1, a0)
                    p = P(a3, a2)
                    x1 = A + B - p / 3.0
                    return x1

                # ---CUBIC EQIATION---
                #Viete's Algorithm
                def result_2(a3,a2,a0,a1):
                    b = b_cof(a3,a2,a1,a0)
                    a = a_cof(a3,a2,a1)
                    p = P(a3, a2)
                    r = -(b / 2.0)
                    q = (a / 3.0)

                    if ((r**2)+(q**3))<= 0.0:
                        if q==0:
                            theta = 0
                        if q<0:
                            theta = cmath.acos(r/(-q**(3.0/2.0)))

                    phi1 = theta / 3.0
                    phi2 = phi1 - ((2*cmath.pi) / 3.0)
                    phi3 = phi1 + ((2*cmath.pi) / 3.0)
                    global label_output1
                    global label_output2
                    global label_output3
                    global label_output4
                    phi1 = theta / 3.0
                    phi2 = phi1 - ((2*cmath.pi) / 3.0)
                    phi3 = phi1 + ((2*cmath.pi) / 3.0)
                    label_output1=Label(top, text="x1 =  "+str("{0.real:.5f}".format(2*math.sqrt(-q)*cmath.cos(phi1)
                    -p/3.0)),font="Times 20 italic bold",bg='white',fg='navy')
                    label_output1.place(x=350, y=250)
                    label_output2=Label(top, text="x2 =  "+str("{0.real:.5f}".format(2*math.sqrt(-q)*cmath.cos(phi2)-
                    p/3.0)),font="Times 20 italic bold",bg='white',fg='navy')
                    label_output2.place(x=350, y=300)
                    label_output3=Label(top, text="x3 =  "+str("{0.real:.5f}".format(2*math.sqrt(-q)*cmath.cos(phi3)
                    -p/3.0)),font="Times 20 italic bold",bg='white',fg='navy')
                    label_output3.place(x=350, y=350)
                    label_output4=Label(top,text="Viete's Algorithm",font="Times 18 italic bold",bg='navy',fg='white')
                    label_output4.place(x=390, y=400)
                #Cardan's Algorithm
                def result_1(a3,a2,a1,a0):
                    b = b_cof (a3,a2, a1,a0)
                    a = a_cof (a3,a2, a1)
                    if((b**2) /4.0+(a**3 )/ 27.0 > 0.0):
                        x1 = solution_x1(a3,a2, a1, a0)
                        global label_output1
                        global label_output2
                        global label_output3
                        global label_output4
                        x1 = solution_x1(a3,a2, a1, a0)
                        label_output1=Label(top, text="x1 =  "+str(round(float(format(x1)),4)),
                        font="Times 20 italic bold",bg='white',fg='navy')
                        label_output1.place(x=350, y=250)
                        label_output2=Label(top, text="x2 =  "+str("{:.5f}".format(real_part(a3,a2, a1, a0)))+' + i'+str
                        ("{:.5f}".format(complex_part(a3,a2, a1, a0))),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output2.place(x=350, y=300)
                        label_output3=Label(top, text="x3 =  "+str("{:.5f}".format(real_part(a3,a2, a1, a0)))+' - i'+str
                        ("{:.5f}".format(complex_part(a3,a2, a1, a0))),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output3.place(x=350, y=350)
                        label_output4=Label(top,text="Cardan's Algorithm",font="Times 18 italic bold",bg='navy',
                        fg='white')
                        label_output4.place(x=390, y=400)
                    else:
                        result_2(a3,a2,a0,a1) #CASE 2

                #-------------------------
                #Coefficient of x^3
                def A_3(a4,a3):
                    A3 = a3 / a4
                    return A3

                #coefficient of x^2
                def A_2(a4,a2):
                    A2 = a2 / a4
                    return A2

                #Coefficient of x^1
                def A_1(a4,a1):
                    A1 = a1/ a4
                    return A1

                #Coefficient of x^0
                def A_0 (a4,a0):
                    A0 = a0/ a4
                    return A0

                def c(a4,a3):
                    A3 = A_3(a4,a3)
                    C = A3/4.0
                    return C

                def b_0(a4,a3,a2,a1,a0):
                    A0 = A_0(a4,a0)
                    A1 = A_1(a4,a1)
                    A2 = A_2 (a4,a2)
                    C = c(a4,a3)
                    b0 = A0-A1*C+A2*pow(C,2.0)-3.0*pow(C,4.0)
                    return b0

                def b_1(a4,a1,a3,a2):
                    A1 = A_1(a4,a1)
                    A2 = A_2 (a4,a2)
                    C = c(a4,a3)
                    b1 = A1-(2.0*A2*C)+(8.0*pow (C,3.0))
                    return b1

                def b_2(a4,a3,a2):
                    A2 = A_2 (a4,a2)
                    C = c(a4,a3)
                    b2 = A2-6.0*pow(C,2.0)
                    return b2

                #function pow(b, 2)/4.0 + pow(a, 3)/27.0 <= 0
                def max_m(a3,a2,a1,a0):
                    b = b_cof(a3,a2,a1,a0)
                    a = a_cof(a3,a2,a1)
                    p = P(a3, a2)
                    r = -(b / 2.0)
                    q = (a / 3.0)
                    if q==0:
                        theta = 0
                    if q<0:
                        theta = math.acos(r / (pow(-q, 3.0 / 2.0)))

                    phi1 = theta / 3.0
                    phi2 = phi1 - ((2 *math.pi) / 3.0)
                    phi3 = phi1 + ((2 *math.pi) / 3.0)
                    x1 = 2 * math.sqrt(-q) * math.cos(phi1)- p / 3.0
                    x2 = 2 * math.sqrt(-q) * math.cos(phi2)- p / 3.0
                    x3 = 2 * math.sqrt(-q) * math.cos(phi3)- p / 3.0
                    if (x1>x2) and (x1>x3):
                        m=x1
                    elif (x2>x1) and (x2>x3):
                        m=x2
                    else:
                        m=x3
                    return m

                #Function to get m
                def M(a4,a3,a2,a1,a0):
                    b0 = b_0(a4,a3,a2,a1,a0)
                    b1 = b_1 (a4,a1,a3,a2)
                    b2 = b_2(a4,a3,a2)
                    f = pow(b2,2)/4.0- b0
                    g = -(b1*b1)/8.0
                    b = b_cof(1.0, b2,f, g)
                    a = a_cof(1.0, b2,f)

                    if (pow(b, 2)/4.0 + pow(a, 3)/27.0) > 0.0:
                        m = solution_x1(1.0, b2,f, g)
                    else:
                        m = max_m(1.0, b2,f, g)

                    if m > 0:
                        m=m
                    else:
                        m=0

                    return m

                # ---QUADRATIC EQUATION---
                def result(a4,a3,a2,a1,a0):
                    b0 = b_0(a4,a3,a2,a1,a0)
                    b1 = b_1 (a4,a1,a3,a2)
                    b2 = b_2(a4,a3,a2)
                    C = c(a4,a3)
                    m = M (a4,a3,a2,a1,a0)

                    if b1>0:
                        Sum = 1
                    else:
                        Sum = -1

                    R = Sum*math.sqrt(pow(m,2)+b2*m+pow(b2,2)/4.0-b0)

                    if ((-m/2.0)-(b2/2.0)-R  < 0) and ((-m/2.0)-(b2/2.0)+ R < 0):
                        global label_output1
                        global label_output2
                        global label_output3
                        global label_output4
                        z1= math.sqrt(-((-m/2.0)-(b2/2.0)-R))
                        z2= math.sqrt(-((-m/2.0)-(b2/2.0)+R))
                        label_output1=Label(top, text="x1 =  "+str("{0.real:.5f}".format((math.sqrt(m/2.0)- C)))+' + i'
                        +str("{:.5f}".format(z1)),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output1.place(x=350, y=250)
                        label_output2=Label(top, text="x2 =  "+str("{0.real:.5f}".format((math.sqrt(m/2.0)- C)))+' - i'
                        +str("{:.5f}".format(z1)),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output2.place(x=350, y=300)
                        label_output3=Label(top, text="x3 =  "+str("{0.real:.5f}".format((-math.sqrt(m/2.0)- C)))+' + i'
                        +str("{:.5f}".format(z1)),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output3.place(x=350, y=350)
                        label_output4=Label(top, text="x4 =  "+str("{0.real:.5f}".format((-math.sqrt(m/2.0)- C)))+' - i'
                        +str("{:.5f}".format(z1)),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output4.place(x=350, y=400)

                    elif (-m/2.0)-(b2/2.0)-R < 0:
                        z1= math.sqrt(-((-m/2.0)-(b2/2.0)-R))
                        label_output1=Label(top, text="x1 =  "+str("{0.real:.5f}".format((math.sqrt(m/2.0)- C)))+' + i'
                        +str("{:.5f}".format(z1)),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output1.place(x=350, y=250)
                        label_output2=Label(top, text="x2 =  "+str("{0.real:.5f}".format((math.sqrt(m/2.0)- C)))+' - i'
                        +str("{:.5f}".format(z1)),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output2.place(x=350, y=300)
                        label_output3=Label(top, text="x3 =  "+str("{:.5f}".format((-math.sqrt(m/2.0)- C + 
                        math.sqrt((-m/2.0)-(b2/2.0)+R)))),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output3.place(x=350, y=350)
                        label_output4=Label(top, text="x4 =  "+str("{:.5f}".format((-math.sqrt(m/2.0)- C - 
                        math.sqrt((-m/2.0)-(b2/2.0)+R)))),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output4.place(x=350, y=400)

                    elif (-m/2.0)-(b2/2.0)+ R < 0:
                        z2 = math.sqrt(-((-m/2.0)-(b2/2.0)+R))
                        label_output1=Label(top, text="x1 =  "+str("{:.5f}".format((math.sqrt(m/2.0)- C 
                        +math.sqrt((-m/2.0)-(b2/2.0)-R)))),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output1.place(x=350, y=250)
                        label_output2=Label(top, text="x2 =  "+str("{:.5f}".format((math.sqrt(m/2.0)- C 
                        -math.sqrt((-m/2.0)-(b2/2.0)-R)))),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output2.place(x=350, y=300)
                        label_output3=Label(top, text="x3 =  "+str("{0.real:.5f}".format((-math.sqrt(m/2.0)- C)))
                        +' + i'+str("{:.5f}".format(z2)),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output3.place(x=350, y=350)
                        label_output4=Label(top, text="x4 =  "+str("{0.real:.5f}".format((-math.sqrt(m/2.0)- C)))
                        +' - i'+str("{:.5f}".format(z2)),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output4.place(x=350, y=400)

                    else:
                        label_output1=Label(top, text="x1 =  "+str(round(((math.sqrt(m/2.0)- C + math.sqrt(-m/2.0-(b2/2.0)-R))),3)),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output1.place(x=350, y=250)
                        label_output2=Label(top, text="x2 =  "+str("{:.3f}".format((math.sqrt(m/2.0)- C -math.sqrt(-m/2.0-(b2/2.0)-R)))),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output2.place(x=350, y=300)
                        label_output3=Label(top, text="x3 =  "+str("{:.3f}".format((-math.sqrt(m/2.0)- C + math.sqrt(-m/2.0-(b2/2.0)+R)))),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output3.place(x=350, y=350)
                        label_output4=Label(top, text="x4 =  "+str("{:.3f}".format((-math.sqrt(m/2.0)- C - math.sqrt(-m/2.0-(b2/2.0)+R)))),font="Times 20 italic bold",bg='white',fg='navy')
                        label_output4.place(x=350, y=400)
                #_______________________________________________________
                if (a4==0):
                    print("\n ----- OUTPUT -----\n")
                    result_1(a3, a2, a1, a0)

                else:
                    print("\n ----- OUTPUT -----\n")
                    result(a4, a3, a2, a1, a0)
            full=1
             
        except ValueError:
            messagebox.showinfo("Error","invalid value enter number")


    #clear coff       
    def clear2():
        #no degree 
        a.delete(0, END)
        b1.delete(0, END)
        c1.delete(0, END)
        d1.delete(0, END)
        e1.delete(0, END)
        clear_label()

    #button for solve equation
    Button(top,text='calculate',font=30,command=solve,fg='navy',bg='white',bd=4).place(x=70,y=460)
    Button(top,text='clear',font=30,command=clear2,fg='navy',bg='white',bd=4).place(x=170,y=460)
#**********************************************************************************************

#______________________________________________________________________________________________
#Stat Tool
def openStat():
    root= Toplevel()
    root.title("Stat Tool by Doha Emad")
    root.geometry("600x600")
    root.configure(bg='black')
    Button(root,text='Home',font=30,command=root.destroy,fg='black',bg='white',bd=4).place(x=520,y=550)
    global go
    go=0
    global clear
    clear=0
    value_count=0
    #logo
    logo_label=Label(root,image=im_logo,bg='black').place(x=500,y=8)

    #text
    head=Label(root,text="Stat Tool", font="Verdana 22 bold",bg='black',fg='deep sky blue')
    data=Label(root,text="Enter Data Set", font=("courier New",15,'bold'),bg='black',fg='white')
    head.place(x=200,y=10)
    data.place(x=40,y=60)
    d=Entry(root, width=40,font=("courier New",17,'bold'))
    d.place(x=20,y=100,height=50)
    d.focus()

    #label
    def label(d):
        global value
        global clear
        clear=1
        try:
            value=list(map(float,d.get().strip().split()))
        except ValueError:
            messagebox.showinfo("Error","invalid value Enter number")
        global go
        if(go==0):
            global label_output
            label_output=Label(root,text='',bg='black')
            label_output.place(x=20,y=500)
            go=1
        else:
            label_output.place_forget()
    

    #clear
    def clear():
        if(clear==1):
            d.delete(0, END)
            value.clear()
            label_output.place_forget()
        else:
            messagebox.showinfo("Error","There is no data to clear")
    want=Label(root,text="What do you want?", font=("courier New",14,'bold'),fg='white',bg='black')
    want.place(x=200,y=182)

    #histogram
    def Histogram():
        if(go==1):
            global label_output
            label_output.place_forget()
        plt.hist(value, bins =10)
        plt.show()
    histogram_b=b_mean=Button(root,text='Histogram',command=Histogram,bg='deep sky blue',fg='black',bd=4,width=14,height=2)
    histogram_b.place(x=230,y=400)

    #mean
    def mean_fun():
        if(go==1):
            global label_output
            label_output.place_forget()
            mean_l=statistics.mean(value)
            label_output=Label(root,text="Mean = "+ str(mean_l),font=("courier New",17,'bold'))
            label_output.place(x=20,y=500)
        else:
            messagebox.showinfo("Error","Click on Enter button first")
    b_mean=Button(root,image=im_mean,command=mean_fun,bg='deep sky blue',fg='black',bd=6,font=12,width=40,height=30 )
    b_mean.place(x=20,y=250)

    #variance
    def variance_fun():
        if(go==1):
            global label_output
            label_output.place_forget()
            variance_l=statistics.variance(value)
            label_output=Label(root,text="Variance = "+str(variance_l),font=("courier New",17,'bold'))
            label_output.place(x=20,y=500)
        else:
            messagebox.showinfo("Error","Click on Enter button first")

    b_variance=Button(root,image=im_variance,command=variance_fun,bg='deep sky blue',fg='black',bd=6,font=12,width=40,height=30)
    b_variance.place(x=140,y=250)

    #standard
    def standard_fun():
        if(go==1):
            global label_output
            label_output.place_forget()
            standard_l=statistics.stdev(value)
            label_output=Label(root,text="Standard Deviation  = "+str(standard_l),font=("courier New",14,'bold'))
            label_output.place(x=20,y=500)
        else:
            messagebox.showinfo("Error","Click on Enter button first")
    b_standard=Button(root,image=im_standard,command=standard_fun,bg='deep sky blue',fg='black',bd=6,font=12,width=40,height=30)
    b_standard.place(x=260,y=250)

    #lenght
    def n_fun():
        if(go==1):
            global label_output
            label_output.place_forget()
            label_output=Label(root,text="n = "+str(len(value)),font=("courier New",17,'bold'))
            label_output.place(x=20,y=500)
        else:
            messagebox.showinfo("Error","Click on Enter button first")
    b_n=Button(root,image=im_n,command=n_fun,bg='deep sky blue',fg='black',bd=6,font=12,width=40,height=30)
    b_n.place(x=380,y=250)

    #median
    def median_fun():
        if(go==1):
            global label_output
            label_output.place_forget()
            median_l=statistics.median(value)
            label_output=Label(root,text="Median  = "+str(median_l),font=("courier New",17,'bold'))
            label_output.place(x=20,y=500)
        else:
            messagebox.showinfo("Error","Click on Enter button first")
    b_median=Button(root,text='Median',command=median_fun,bg='deep sky blue',fg='black',bd=4,width=6,height=2)
    b_median.place(x=20,y=320)

    #mode
    def mode_fun():
        if(go==1):
            global label_output
            label_output.place_forget()
            count=0
            for i in value:
                tempnum=i
                tempcount=0
                for j in value:
                    if(j==tempnum):
                        tempcount+=1
                        if(tempcount>count):
                            number=tempnum
                            count=tempcount
            if(count==1):
                label_output=Label(root,text="No Mode",font=("courier New",17,'bold'))
            else:
                label_output=Label(root,text="Mode = "+str(number)+' is repeated '+str(count)+"times",font=("courier New",17,'bold'))
            label_output.place(x=20,y=500)
        else:
            messagebox.showinfo("Error","Click on Enter button first")
    b_mode=Button(root,text='Mode',command=mode_fun,bg='deep sky blue',fg='black',bd=4,width=6,height=2)
    b_mode.place(x=140,y=320)

    #max
    def max_fun():
        if(go==1):
            global label_output
            label_output.place_forget()
            label_output=Label(root,text="Max = "+str(max(value)),font=("courier New",17,'bold'))
            label_output.place(x=20,y=500)
        else:
            messagebox.showinfo("Error","Click on Enter button first")
    b_max=Button(root,text='Max',command=max_fun,bg='deep sky blue',fg='black',bd=4,width=6,height=2)
    b_max.place(x=260,y=320)

    #mini
    def mini_fun():
        if(go==1):
            global label_output
            label_output.place_forget()
            label_output=Label(root,text="Mini = "+str(min(value)),font=("courier New",17,'bold'))
            label_output.place(x=20,y=500)
        else:
            messagebox.showinfo("Error","Click on Enter button first")
    b_max=Button(root,text='Mini',command=mini_fun,bg='deep sky blue',fg='black',bd=4,width=6,height=2)
    b_max.place(x=380,y=320)

    #sort
    def sort_fun():
        if(go==1):
            global label_output
            label_output.place_forget()
            label_output=Label(root,text="Descending order = "+str(sorted(value,reverse=True)),font=("courier New",13,'bold'))
            label_output.place(x=20,y=500)
        else:
            messagebox.showinfo("Error","Click on Enter button first")
    b_sort=Button(root,text='Descending order',command=sort_fun,bg='deep sky blue',fg='black',bd=4,width=14,height=2)
    b_sort.place(x=470,y=250)

    #range
    def range_fun():
        if(go==1):
            global label_output
            label_output.place_forget()
            sort=sorted(value)
            range=sort[len(sort)-1]-sort[0]
            label_output=Label(root,text="Range = "+str(range),font=("courier New",17,'bold'))
            label_output.place(x=20,y=500)
        else:
            messagebox.showinfo("Error","Click on Enter button first")
    b_sort=Button(root,text='Range',command=range_fun,bg='deep sky blue',fg='black',bd=4,width=14,height=2)
    b_sort.place(x=470,y=320)

    #clear button
    b2=Button(root,text='Clear',command=clear,bg='white',fg='black',bd=6,width=4,height=0 )
    b2.place(x=530,y=152)
    b1=Button(root,image=im,command=lambda: label(d), bg='deep sky blue',fg='black',bd=6,font=12,width=20,height=20 )
    b1.place(x=20,y=155)
#________________________________________________________________________________________________

ground.geometry('770x575')
ground.title('Main Window')
# place a button on the root window
back_l=Label(ground,image=back).place(x=0,y=0)
label= Label(ground, text="What do you want?",font="Times 25 italic bold")
label.place(x=40, y=150)
l=Label(ground,text='Welcome',font="Times 35 italic bold",bg='white').place(x=300,y=40)
Button(ground,text='Equation Solver',font=30,command=openSolver,bg='white',bd=10).place(x=90,y=240)
Button(ground,text='Stat Tool',font=25,command=openStat,bg='white',bd=10).place(x=500,y=240)

# mainloop, runs infinitely
mainloop()