import numpy as np
import matplotlib.pyplot as plt




def func(x):
    y = x
    return y

def RoofThickness(x):
    """
    This function computes the roof thickness as function of the angle x
    """
    # Public Function f6(x_) As Double
    # 'For comments look ELEC5b.FORTRAN code
    #     t = tr / Sqr(1 - (td * Sin(x_) / ri) ^ 2)
    #     f6 = a1t * Exp(-t / lamat) * Sin(x_)
    #     f6 = f6 + b1t * Exp(-t / lambt) * Sin(x_)
    #     f6 = f6 + c1t * Exp(-t / lamct) * Sin(x_)
    # End Function

    c = 0.000225
    b = 0.00225
    a = 0.009545
    lambda_a = 30
    lambda_b = 55
    lambda_c= 120

    thickness = tr  / np.sqrt(1 - np.power(td*np.sin(x) / ri, 2))

    roof_thickness = np.sin(x) * (
                        a * np.exp(- thickness / lambda_a) + 
                        b * np.exp(- thickness / lambda_b) + 
                        c * np.exp(- thickness / lambda_c)
                        )

    return roof_thickness    

def Trap(a, b, s, j):
    if (j == 1):
        s = 0.5 * (b-a) * (RoofThickness(a) + RoofThickness(b))
        
    else:
        IT = np.power(2, j-2)
        tnm = IT
        DEL = (b - a) / tnm
        x = a + 0.5 * DEL
        SUMM = 0

        for j in range(1, IT):
            SUMM = SUMM + RoofThickness(x)
            x = x + DEL
        s = 0.5 * (s +  (b - a)*SUMM / tnm)
    
    return s

    # Public Sub TRAPZD(a_, b_, s_, n_)
    # Dim DEL, x, SUMM As Double   '?
    # Debug.Print ("a_ :" & a_ & ", b_ : " & b_)
    # Static IT, tnm As Integer
    #     Select Case n_
    #         Case 1
    #             s_ = 0.5 * (b_ - a_) * (f6(a_) + f6(b_))
    #             IT = 1
    #         Case Else
                
    #             tnm = IT
    #             DEL = (b_ - a_) / tnm
    #             x = a_ + 0.5 * DEL
    #             SUMM = 0
    #                 For j = 1 To IT
    #                 SUMM = SUMM + f6(x)
    #                 x = x + DEL
    #                 Next j
    #             s_ = 0.5 * (s_ + (b_ - a_) * SUMM / tnm)
    #              IT = IT * 2
    #     End Select
    # End Sub

def OTrap(a, b, EPS = 0, JMAX = 2, OLDS = -1e30):
    """
    functin to call trapezoidal integration over the range a to be
    convergence is compared to EPS provide the iterations are below JMAX
    """

    for i in range(1,JMAX):
        if (i == 1):
            s = Trap(a, b, 0, i)
        else:
            s = Trap(a, b, s, i)
        if (np.abs(s - OLDS) < EPS * np.abs(OLDS)):
            break
        OLDS = s
    
    if (i == JMAX):
        print("End of iterations reached, using last value")
    return s
# Public Sub OTRAP(a_, b_, s_)
# 'For comments look ELEC5b.FORTRAN code
# Dim OLDS As Single
#     Const EPS As Single = 0.0001      'parameter
#     Const JMAX As Integer = 50      'parameter
#     OLDS = -1E+30
#         For j = 1 To JMAX
#         TRAPZD a_, b_, s_, j
#         If Abs(s_ - OLDS) < EPS * Abs(OLDS) Then GoTo 12
#         OLDS = s_
#         Next j
#     MsgBox "TOO MANY STEPS (function OTRAP)"
# 12 End Sub


tr = 286.5
td = 36.5
ri = 56.4

if __name__ == "__main__":
    b = 1.107
    integ = (OTrap(0,b) *2 * np.pi / (1.4539473684210527e-14) )


    x = np.arange(1,3,0.001)
    y = func(x)
    print(np.trapz(y, dx=0.001))

    print(integ / 5552961931.032397)

    x = np.arange(0,b,0.001)
    y = RoofThickness(x)
    print(np.trapz(y,dx = 0.001))