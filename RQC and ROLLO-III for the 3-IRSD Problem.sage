# Parameters of RQC and ROLLO-III for the 3-IRSD Problem
# Algebric method ---- MM modeling without Grobner basis  
# The MM modeling for the two Blockwise Rank Syndrome Decoding (3-RSD) Problem

from itertools import combinations, combinations_with_replacement

def Equations_p(m,n,k,r1,r2,r3,p):  # the maximal  number of equations   
    return m*binomial(n-k-p-1,r1+r2+r3)

def Unkonwns_a_p(m,n,k,r1,r2,r3,a1,a2,a3,p): # the number of unkonwns 
    result = binomial(n/3-a1,r1)* binomial(n/3-a2,r2) * binomial(n/3-a3-p,r3)-1
    return result

def Equations_Rate_Unkonwns(m,n,k,r1,r2,r3,a1,a2,a3,p):  # If < 1,  equations < unkonwns，underdetermined case
    result = Equations_p(m,n,k,r1,r2,r3,p)/Unkonwns_a_p(m,n,k,r1,r2,r3,a1,a2,a3,p) 
    return result

def p_values(m,n,k,r1,r2,r3): 
    delta = Equations_Rate_Unkonwns(m,n,k,r1,r2,r3,0,0,0,0)
    list_p = [0] 
    for i in range(1,n/3-r3-1): 
        v = Equations_Rate_Unkonwns(m,n,k,r1,r2,r3,0,0,0,i)
        if (n-i-k-1 >= r1+r2+r3) and (n/3-i >= r3) and (1.000 <= v < delta): 
            list_p.append(i) 
    return list_p

def a_tuples(m,n,k,r1,r2,r3): 
    list_a = [] 
    for comb in combinations_with_replacement(list(range(40)),3): 
        v = Equations_Rate_Unkonwns(m,n,k,r1,r2,r3,comb[0],comb[1],comb[2],0)
        if (n/3-comb[0] >= r1) and (n/3-comb[1] >= r2) and (n/3-comb[2] >= r3) and (1.000 <= v <= 1.300):
            list_a.append(comb) 
    return list_a


def WF1(m,n,k,r1,r2,r3,p):  # The complexity of overdetemined case if one knowns p
    U = Unkonwns_a_p(m,n,k,r1,r2,r3,0,0,0,p)
    E = Equations_p(m,n,k,r1,r2,r3,p) 
    WF = log(E,2) + (w-1)*log(U,2) 
    return WF

def WF2(m,n,k,r1,r2,r3):  # The complexity of overdetemined case 
    list_p = p_values(m,n,k,r1,r2,r3)
    Rate_list = []
    p_WF = { }
    for i in range(len(list_p)):
        p = list_p[i]
        U = Unkonwns_a_p(m,n,k,r1,r2,r3,0,0,0,p)
        E = Equations_p(m,n,k,r1,r2,r3,p) 
        WF = log(E,2) + (w-1)*log(U,2)
        p_WF[p] = float(WF)
        v = Equations_Rate_Unkonwns(m,n,k,r1,r2,r3,0,0,0,p)
        Rate_list.append(round(v,4))
        
    p_WF_Rate = dict(zip(list(p_WF.items()),Rate_list))    
    Rate = min(Rate_list)
    p_WF = [key for key, value in p_WF_Rate.items() if value == Rate]
    return p_WF[0][0],p_WF[0][1], p_WF_Rate

def WF3(m,n,k,r1,r2,r3,a1,a2,a3): # The complexity of underdetermined case if one knowns a1, a2, a3 
    U = Unkonwns_a_p(m,n,k,r1,r2,r3,a1,a2,a3,0)
    E = Equations_p(m,n,k,r1,r2,r3,0) 
    WF = log(E,2) + (w-1)*log(U,2) + (a1*r1 + a2*r2 + a3*r3)
    return WF

def WF4(m,n,k,r1,r2,r3): # The complexity of underdetermined case
    E = Equations_p(m,n,k,r1,r2,r3,0)
    a_WF = {}
    Rate_list = []
    list_a = a_tuples(m,n,k,r1,r2,r3)
    for i in range(len(list_a)):
        a1 = list_a[i][0]
        a2 = list_a[i][1]
        a3 = list_a[i][2]
        U = Unkonwns_a_p(m,n,k,r1,r2,r3,a1,a2,a3,0)
        WF = log(E,2) + (w-1)*log(U,2) + (a1*r1 + a2*r2 + a3*r3)
        a_WF[list_a[i]] = round(WF,4)
        v = Equations_Rate_Unkonwns(m,n,k,r1,r2,r3,a1,a2,a3,0)
        Rate_list.append(round(v,4))
    
    a_WF_Rate = dict(zip(list(a_WF.items()),Rate_list))
    Rate = min(Rate_list)
    a_tuple_WF = [key for key, value in a_WF_Rate.items() if value == Rate]
    return a_tuple_WF[0][0], a_tuple_WF[0][1], a_WF_Rate

# m: extension degree of the field F_q to F_{q**m}; n: attacked code_length, k: attacked code_dimension; 
# w : the exponent of matrix multiplication 2 <= w <= 3 and a practical value is 2.81

# Our RQC; n1 = n2 = n3 = n/3, k = n/3
(q,m,n,k,r1,r2,r3,w)= (2,83,3*79,79,4,4,4,2.81) #  3-IRSD(3n)  128,    p= 40,    WF = 163 
#(q,m,n,k,r1,r2,r3,w)= (2,101,3*97,97,4,5,4,2.81)  #  3-IRSD(3n)  128,    p= 55 ,    WF = 182   
#(q,m,n,k,r1,r2,r3,w)= (2,127,3*113,113,5,5,5,2.81) #  3-IRSD(3n)  192,    p = 57 ,    WF = 214
#(q,m,n,k,r1,r2,r3,w)= (2,139,3*137,137,6,6,6,2.81) #  3-IRSD(3n)  256,    p = 19,    WF = 274
# p = 99 is a strange value 

# Our Rollo-III Ourobors; DFR is around 2**(-30); n1 = n2 = n3 = n/3, k = n/3
#(q,m,n,k,r1,r2,r3,w)= (2,59,3*79,79,4,4,5,2.81) #  3-IRSD(3n)  # 128,     p = 31,   WF = 175
#(q,m,n,k,r1,r2,r3,w)= (2,89,3*101,101,6,6,6,2.81) #  3-IRSD(3n)  # 192,   a1 = a2 = 0; a3 = 2,    WF = 266
#(q,m,n,k,r1,r2,r3,w)= (2,97,3*103,103,6,6,7,2.81) #  3-IRSD(3n)  # 256,   a1 = a2 = a3 = 2,    WF = 304
#(q,m,n,k,r1,r2,r3,w)= (2,97,3*109,109,6,6,7,2.81) #  3-IRSD(3n)  # 256,    a1 = a2 = 0; a3 = 5,  WF = 305



print("Observe the Rate of Equations and Unkonwns:")
Rate = Equations_Rate_Unkonwns(m,n,k,r1,r2,r3,0,0,0,0)
print("Rate：%f, If Rate >= 1：overdetermined case, otherwise, underdetermined case."% Rate)

if Rate >= 1: 
    print("Overdetermined case,  needing to set the value of p") 
    print("Workfactor when p = 0：%f"% WF1(m,n,k,r1,r2,r3,0))
    print("Possible p ：%s"% p_values(m,n,k,r1,r2,r3)) 
    p, WF, p_WF_Rate = WF2(m,n,k,r1,r2,r3)
    Rate = Equations_Rate_Unkonwns(m,n,k,r1,r2,r3,0,0,0,p)
    print("Updated the Minimal Rate of Equations and Unkonwns：%f"% Rate)
    print("\033[91m Optimal p: %s; Workfactor：%f \033[0m"% (p, WF)) 
    print("If need, print (p value, Workfactor, Updated Rate).")
    print("- p value,   Workfactor,  Updated Rate --------")
    for key, value in p_WF_Rate.items():
        print("     %s,      %f,  %f" %(key[0], key[1], value))
else: 
    print("Underdetermined case, needing to set the values of a1, a2 and a3.") 
    print("Workfactor when a1 = a2 = a3 = 0：%f"% WF3(m,n,k,r1,r2,r3,0,0,0))
    print("Possible a = (a1, a2, a3)：%s"% a_tuples(m,n,k,r1,r2,r3)) 
    a_tuple, WF, a_WF_Rate = WF4(m,n,k,r1,r2,r3)
    a1 = a_tuple[0]; a2 = a_tuple[1]; a3 = a_tuple[2]
    Rate = Equations_Rate_Unkonwns(m,n,k,r1,r2,r3,a1,a2,a3,0)
    print("Updated the Minimal Rate of Equations and Unkonwns：%f"% Rate)
    print("\033[91m Optimal a tuple: %s; Workfactor：%f \033[0m"% (a_tuple, WF)) 
    print("If need, print (tuple of a, Workfactor, Updated Rate).")
    print("- tuple of a,  Workfactor,  Updated Rate --------")
    for key, value in a_WF_Rate.items():
        print("  %s,   %f,  %f" %(key[0], key[1], value))