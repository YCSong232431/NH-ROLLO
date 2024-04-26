# Complexity of 2-IRSD Problem for RQC and ROLLO
# Algebric method ---- MM modeling without Grobner basis  
# The MM modeling for the 2-RSD Problem

from itertools import combinations, combinations_with_replacement
import math

def Equations_p(m,n,k,r1,r2,p):  # the maximal  number of equations   
    return m*binomial(n-p-k-1,r1+r2)

def Unkonwns_a_p(m,n,k,n1,n2,r1,r2,a1,a2,p): # the number of unkonwns 
    result = binomial(n1-a1,r1) * binomial(n2-a2-p,r2) - 1
    return result

def Equations_Rate_Unkonwns(m,n,k,n1,n2,r1,r2,a1,a2,p):  # If < 1, equations < unkonwns，underdetermined case
    result = Equations_p(m,n,k,r1,r2,p)/Unkonwns_a_p(m,n,k,n1,n2,r1,r2,a1,a2,p) 
    return result

def p_values(m,n,k,n1,n2,r1,r2): 
    delta = Equations_Rate_Unkonwns(m,n,k,n1,n2,r1,r2,0,0,0)
    list_p = [0] 
    for i in range(1,n2-r2-1): 
        v = Equations_Rate_Unkonwns(m,n,k,n1,n2,r1,r2,0,0,i)
        if (n-i-k-1 >= r1+r2) and (n2-i >= r2) and (1.000 <= v < delta): 
            list_p.append(i) 
    return list_p

def a_tuples(m,n,k,n1,n2,r1,r2): 
    list_a = [] 
    for comb in combinations_with_replacement(list(range(40)),2): 
        v = Equations_Rate_Unkonwns(m,n,k,n1,n2,r1,r2,comb[0],comb[1],0)
        if (n1-comb[0] >= r1) and (n2-comb[1] >= r2) and (1.000 <= v <= 1.300):
            list_a.append(comb) 
    return list_a


def WF1(m,n,k,n1,n2,r1,r2,p):  # The complexity of overdetermined case if one knowns p
    U = Unkonwns_a_p(m,n,k,n1,n2,r1,r2,0,0,p)
    E = Equations_p(m,n,k,r1,r2,p) 
    WF = math.log(E,2) + (w-1)* math.log(U,2) 
    return WF

def WF2(m,n,k,n1,n2,r1,r2):  # The complexity of overdetermined case 
    list_p = p_values(m,n,k,n1,n2,r1,r2)
    Rate_list = []
    p_WF = { }
    for i in range(len(list_p)):
        p = list_p[i]
        U = Unkonwns_a_p(m,n,k,n1,n2,r1,r2,0,0,p)
        E = Equations_p(m,n,k,r1,r2,p) 
        WF = math.log(E,2) + (w-1)* math.log(U,2)
        p_WF[p] = float(WF)
        v = Equations_Rate_Unkonwns(m,n,k,n1,n2,r1,r2,0,0,p)
        Rate_list.append(round(v,4))
        
    p_WF_Rate = dict(zip(list(p_WF.items()),Rate_list))    
    Rate = min(Rate_list)
    p_WF = [key for key, value in p_WF_Rate.items() if value == Rate]
    return p_WF[0][0],p_WF[0][1], p_WF_Rate

def WF3(m,n,k,n1,n2,r1,r2,a1,a2): # The complexity of underdetermined case if one knowns a1, a2
    U = Unkonwns_a_p(m,n,k,n1,n2,r1,r2,a1,a2,0)
    E = Equations_p(m,n,k,r1,r2,0) 
    WF = math.log(E,2) + (w-1)* math.log(U,2) + (a1*r1 + a2*r2)
    return WF

def WF4(m,n,k,n1,n2,r1,r2): # The complexity of underdetermined case
    E = Equations_p(m,n,k,r1,r2,0)
    a_WF = {}
    Rate_list = []
    list_a = a_tuples(m,n,k,n1,n2,r1,r2)
    for i in range(len(list_a)):
        a1 = list_a[i][0]
        a2 = list_a[i][1]
        U = Unkonwns_a_p(m,n,k,n1,n2,r1,r2,a1,a2,0)
        WF = math.log(E,2) + (w-1) * math.log(U,2) + (a1*r1 + a2*r2)
        a_WF[list_a[i]] = round(WF,4)
        v = Equations_Rate_Unkonwns(m,n,k,n1,n2,r1,r2,a1,a2,0)
        Rate_list.append(round(v,4))
    
    a_WF_Rate = dict(zip(list(a_WF.items()),Rate_list))
    Rate = min(Rate_list)
    a_tuple_WF = [key for key, value in a_WF_Rate.items() if value == Rate]
    return a_tuple_WF[0][0], a_tuple_WF[0][1], a_WF_Rate


# m: extension degree; n: code_length, k: code_dimension; 
# w : the exponent of matrix multiplication 2 <= w <= 3 and a practical value is 2.81


# Our RQC, 2-IRSD: n1 = n2 = n/2, k = n/2
#(q,m,n,k,n1,n2,r1,r2,w)= (2,83,79*2,79,79,79,4,4,2.81) # 128    a1 = 1， a2 = 2,   WF = 127
#(q,m,n,k,n1,n2,r1,r2,w)= (2,127,113*2,113,113,113,5,5,2.81)  # 192   a1 = 5，a2 = 16,   WF = 253
#(q,m,n,k,n1,n2,r1,r2,w)= (2,139,137*2,137,137,137,5,5,2.81)  # 256   a1 = 5, a2 = 17,   WF = 267

# Our ROLLO-I  Lake, DFR is around 2**(-30);
# Message attack: c = e1 + h.e2; n1 = n2 = n/2, k = n/2
#(q,m,n,k,n1,n2,r1,r2,w)= (2,59,2*67,67,67,67,4,4,2.81) # 128    a1 = 4, a2 = 5,     WF = 144
#(q,m,n,k,n1,n2,r1,r2,w)= (2,79,2*89,89,89,89,5,5,2.81) # 192    a1 = 7,  a2 = 18,  WF = 261
#(q,m,n,k,n1,n2,r1,r2,w)= (2,79,2*97,97,97,97,5,5,2.81) # 256    a1 = 10, a2 = 17,  WF = 275


# Our Updated ROLLO-II  Locker;  DFR is around 2**(-130); 
# Message attack: c = e1 + h.e2; n1 = n2 = n/2, k = n/2
#(q,m,n,k,n1,n2,r1,r2,w)= (2,83,173*2,173,173,173,4,4,2.81) # 128   p = 1,     WF = 141
#(q,m,n,k,n1,n2,r1,r2,w)= (2,97,179*2,179,179,179,4,5,2.81)  # 192   a1 = 6, a2 = 10,     WF = 229
#(q,m,n,k,n1,n2,r1,r2,w)= (2,101,191*2,191,191,191,5,5,2.81)  # 256   a1 = 3, a2 = 35,     WF = 359


# Our Rollo-III Ouroboros; DFR is around 2**(-30); n1 = n2 = n/2, k = n/2
#(q,m,n,k,n1,n2,r1,r2,w)= (2,59,73*2,73,73,73,4,4,2.81) # 128    a1 = 4, a2 = 5,   WF = 147
#(q,m,n,k,n1,n2,r1,r2,w)= (2,89,101*2,101,101,101,4,5,2.81)  # 192   a1 = 3 = a2 = 10,   WF = 195
#(q,m,n,k,n1,n2,r1,r2,w)= (2,97,109*2,109,109,109,5,5,2.81)  # 256   a1 = 3, a2 = 22,   WF = 270


print("Observe the Rate of Equations and Unkonwns:")
Rate = Equations_Rate_Unkonwns(m,n,k,n1,n2,r1,r2,0,0,0)
print("Rate：%f, If Rate >= 1：overdetermined case, otherwise, underdetermined case."% Rate)

if Rate >= 1: 
    print("Overdetermined case, needing to set the value of p") 
    print("Workfactor when p = 0：%f"% WF1(m,n,k,n1,n2,r1,r2,0))
    print("Possible p ：%s"% p_values(m,n,k,n1,n2,r1,r2)) 
    p, WF, p_WF_Rate = WF2(m,n,k,n1,n2,r1,r2)
    Rate = Equations_Rate_Unkonwns(m,n,k,n1,n2,r1,r2,0,0,p)
    print("Updated the Minimal Rate of Equations and Unkonwns：%f"% Rate)
    print("\033[91m Optimal p: %s; Workfactor：%f \033[0m"% (p, WF)) 
    print("If need, print (p value, Workfactor, Updated Rate).")
    print("- p value,   Workfactor,  Updated Rate --------")
    for key, value in p_WF_Rate.items():
        print("     %s,      %f,  %f" %(key[0], key[1], value))
else: 
    print("Underdetermined case, needing to set the values of a1 and a2.") 
    print("Workfactor when a1 = a2 = 0：%f"% WF3(m,n,k,n1,n2,r1,r2,0,0))
    print("Possible a = (a1, a2)：%s"% a_tuples(m,n,k,n1,n2,r1,r2)) 
    a_tuple, WF, a_WF_Rate = WF4(m,n,k,n1,n2,r1,r2)
    a1 = a_tuple[0]; a2 = a_tuple[1]
    Rate = Equations_Rate_Unkonwns(m,n,k,n1,n2,r1,r2,a1,a2,0)
    print("Updated the Minimal Rate of Equations and Unkonwns：%f"% Rate)
    print("\033[91m Optimal a tuple: %s; Workfactor：%f \033[0m"% (a_tuple, WF))
    print("If need, print (tuple of a, Workfactor, Updated Rate).")
    print("- tuple of a,  Workfactor,  Updated Rate --------")
    for key, value in a_WF_Rate.items():
        print("  %s,     %f,  %f" %(key[0], key[1], value))
