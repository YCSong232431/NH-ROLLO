# Parameters of RQC and ROLLO-III for the 2-IRSD Problem
# Algebric method ---- MM modeling without Grobner basis  
# The MM modeling for the two Blockwise Rank Syndrome Decoding (2-RSD) Problem

from itertools import combinations, combinations_with_replacement
import math

def Equations_p(m,n,k,r1,r2,p):  # the maximal  number of equations   
    return m*binomial(n-k-p-1,r1+r2)

def Unkonwns_a_p(m,n,k,r1,r2,a1,a2,p): # the number of unkonwns 
    result = binomial(n/2-a1,r1)* binomial(n/2-a2,r2) - 1
    return result

def Equations_Rate_Unkonwns(m,n,k,r1,r2,a1,a2,p):  # If < 1, equations < unkonwns，underdetermined case
    result = Equations_p(m,n,k,r1,r2,p)/Unkonwns_a_p(m,n,k,r1,r2,a1,a2,p) 
    return result

def p_values(m,n,k,r1,r2): 
    delta = Equations_Rate_Unkonwns(m,n,k,r1,r2,0,0,0)
    list_p = [0] 
    for i in range(1,n/2-r2-1): 
        v = Equations_Rate_Unkonwns(m,n,k,r1,r2,0,0,i)
        if (n-i-k-1 >= r1+r2) and (n/2-i >= r3) and (1.000 <= v < delta): 
            list_p.append(i) 
    return list_p

def a_tuples(m,n,k,r1,r2): 
    list_a = [] 
    for comb in combinations_with_replacement(list(range(30)),2): 
        v = Equations_Rate_Unkonwns(m,n,k,r1,r2,comb[0],comb[1],0)
        if (n/2-comb[0] >= r1) and (n/2-comb[1] >= r2)  and (1.000 <= v <= 1.300):
            list_a.append(comb) 
    return list_a


def WF1(m,n,k,r1,r2,p):  # The complexity of overdetemined case if one knowns p
    U = Unkonwns_a_p(m,n,k,r1,r2,0,0,p)
    E = Equations_p(m,n,k,r1,r2,p) 
    WF = math.log(E,2) + (w-1) * math.log(U,2) 
    return WF

def WF2(m,n,k,r1,r2):  # The complexity of overdetemined case 
    list_p = p_values(m,n,k,r1,r2)
    Rate_list = []
    p_WF = { }
    for i in range(len(list_p)):
        p = list_p[i]
        U = Unkonwns_a_p(m,n,k,r1,r2,0,0,p)
        E = Equations_p(m,n,k,r1,r2,p) 
        WF = math.log(E,2) + (w-1) * math.log(U,2)
        p_WF[p] = float(WF)
        v = Equations_Rate_Unkonwns(m,n,k,r1,r2,0,0,p)
        Rate_list.append(round(v,4))
        
    p_WF_Rate = dict(zip(list(p_WF.items()),Rate_list))    
    Rate = min(Rate_list)
    p_WF = [key for key, value in p_WF_Rate.items() if value == Rate]
    return p_WF[0][0],p_WF[0][1], p_WF_Rate

def WF3(m,n,k,r1,r2,a1,a2): # The complexity of underdetermined case if one knowns a1 and a2
    U = Unkonwns_a_p(m,n,k,r1,r2,a1,a2,0)
    E = Equations_p(m,n,k,r1,r2,0) 
    WF = math.log(E,2) + (w-1)*math.log(U,2) + (a1*r1 + a2*r2)
    return WF

def WF4(m,n,k,r1,r2): # The complexity of underdetermined case
    E = Equations_p(m,n,k,r1,r2,0)
    a_WF = {}
    Rate_list = []
    list_a = a_tuples(m,n,k,r1,r2)
    for i in range(len(list_a)):
        a1 = list_a[i][0]
        a2 = list_a[i][1]
        U = Unkonwns_a_p(m,n,k,r1,r2,a1,a2,0)
        WF = math.log(E,2) + (w-1)*math.log(U,2) + (a1*r1 + a2*r2 )
        a_WF[list_a[i]] = round(WF,4)
        v = Equations_Rate_Unkonwns(m,n,k,r1,r2,a1,a2,0)
        Rate_list.append(round(v,4))
    
    a_WF_Rate = dict(zip(list(a_WF.items()),Rate_list))
    Rate = min(Rate_list)
    a_tuple_WF = [key for key, value in a_WF_Rate.items() if value == Rate]
    return a_tuple_WF[0][0], a_tuple_WF[0][1], a_WF_Rate

# m: extension degree; n: code_length, k: code_dimension; 
# w : the exponent of matrix multiplication 2 <= w <= 3 and a practical value is 2.81

# Our RQC; n1 = n2 = n/2, k = n/2
#(q,m,n,k,r1,r2,w)= (2,83,2*79,79,4,4,2.81)  #  2-IRSD(2n) 128    a1 = 1， a2 = 2,   WF = 127
#(q,m,n,k,r1,r2,w)= (2,101,2*97,97,5,4,2.81)  #  2-IRSD(2n)  128   a1 = 0, a2 = 12  more reliable in security level  180
# (q,m,n,k,r1,r2,w)= (2,127,2*113,113,5,5,2.81)  #  2-IRSD(2n)  192  a1 = 5，a2 = 16,   WF = 253
#(q,m,n,k,r1,r2,w)= (2,139,2*137,137,5,5,2.81)  #  2-IRSD(2n)  256    a1 = 5, a2 = 17,   WF = 267


# Our Rollo-III Ourobors; DFR is around 2**(-30); n1 = n2 = n/2, k = n/2
#(q,m,n,k,r1,r2,w)= (2,59,2*79,79,4,4,2.81) # 2-IRSD(2n)  # 128     a1 = 0, a2 = 9,   WF = 149
#(q,m,n,k,r1,r2,w)= (2,59,2*73,73,4,4,2.81) # 2-IRSD(2n)  # 128     a1 = 4, a2 = 5,   WF = 147
# (q,m,n,k,r1,r2,w)= (2,89,2*101,101,4,5,2.81) # 2-IRSD(2n)  # 192   a1 = 3 = a2 = 10,   WF = 195
# (q,m,n,k,r1,r2,w)= (2,97,2*103,103,5,5,2.81) # 2-IRSD(2n)  # 256   a1 = 3, a2 = 21,   WF =  263
#(q,m,n,k,r1,r2,w)= (2,97,2*109,109,5,5,2.81) # 2-IRSD(2n)  # 256   a1 = 3, a2 = 22,   WF = 270


print("Observe the Rate of Equations and Unkonwns:")
Rate = Equations_Rate_Unkonwns(m,n,k,r1,r2,0,0,0)
print("Rate：%f, If Rate >= 1：overdetermined case, otherwise, underdetermined case."% Rate)

if Rate >= 1: 
    print("Overdetermined case,  needing to set the value of p") 
    print("Workfactor when p = 0：%f"% WF1(m,n,k,r1,r2,0))
    print("Possible p ：%s"% p_values(m,n,k,r1,r2)) 
    p, WF, p_WF_Rate = WF2(m,n,k,r1,r2)
    Rate = Equations_Rate_Unkonwns(m,n,k,r1,r2,0,0,p)
    print("Updated the Minimal Rate of Equations and Unkonwns：%f"% Rate)
    print("\033[91m Optimal p: %s; Workfactor：%f \033[0m"% (p, WF)) 
    print("If need, print (p value, Workfactor, Updated Rate).")
    print("- p value,   Workfactor,  Updated Rate --------")
    for key, value in p_WF_Rate.items():
        print("     %s,      %f,  %f" %(key[0], key[1], value))
else: 
    print("Underdetermined case, needing to set the values of a1 and a2.") 
    print("Workfactor when a1 = a2 = 0：%f"% WF3(m,n,k,r1,r2,0,0))
    print("Possible a = (a1, a2)：%s"% a_tuples(m,n,k,r1,r2)) 
    a_tuple, WF, a_WF_Rate = WF4(m,n,k,r1,r2)
    a1 = a_tuple[0]; a2 = a_tuple[1]
    Rate = Equations_Rate_Unkonwns(m,n,k,r1,r2,a1,a2,0)
    print("Updated the Minimal Rate of Equations and Unkonwns：%f"% Rate)
    print("\033[91m Optimal a tuple: %s; Workfactor：%f \033[0m"% (a_tuple, WF)) 
    print("If need, print (tuple of a, Workfactor, Updated Rate).")
    print("- tuple of a,  Workfactor,  Updated Rate --------")
    for key, value in a_WF_Rate.items():
        print("  %s,   %f,  %f" %(key[0], key[1], value))
