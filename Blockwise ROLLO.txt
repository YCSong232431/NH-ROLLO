
def random_small_vec_gen(n,t):
    B = matrix(Fqm.base_ring(),t,m,0)
    while B.rank() != t:
        B = matrix(Fqm.base_ring(),[vector(Fqm.random_element()) for i in range(t)])
    C = matrix(Fqm.base_ring(),n,m,0)
    while C.rank() != t:
        C = matrix(Fqm.base_ring(),n,t,[Fqm.base_ring().random_element() for _ in range(n*t)]) * B
    return vector(Fqm,[C[i] for i in range(n)])

def vector_matrix(small_vector): 
    length = len(list(small_vector))
    return matrix(Fqm.base_ring(),length,m,[vector(small_vector[j]) for j in range(length)])

def R_to_Matrix(z):
    return matrix(Fqm.base_ring(),n,m,[vector(z[i]) for i in range(n)])

def rank_R(z):
    return R_to_Matrix(z).rank()

def ideal_matrix(v): 
    list1 = []
    for i in range(n/2):
        c1 = list(X**i * v)
        list1.append(c1)
    return matrix(n/2,n/2,list1)

def H_to_basis_list(B):
    b = vector(B)
    c = B.nrows() * B.ncols()
    BB =  matrix(Fqm.base_ring(),c,m,[vector(b[j]) for j in range(c)]).row_space().basis_matrix()
    return [Fqm(x) for x in BB.rows()]

def vector_matrix(small_vector): 
    length = len(list(small_vector))
    return matrix(Fqm.base_ring(),length,m,[vector(small_vector[j]) for j in range(length)])

def vector_space(small_vector): 
    return vector_matrix(small_vector).row_space()

def vector_to_basis_list(small_vector):
    return [Fqm(x) for x in vector_space(small_vector).basis_matrix().rows()]

def basis_list_to_basis_matrix(b):
    return vector_space(b).basis_matrix() 

def basis_matrix_to_list(B):
    return [Fqm(x) for x in B.rows()]

def basis_product_of_two_space(basis_list1,basis_list2):
    product_list = []
    for i in range(len(basis_list1)): 
        for j in range(len(basis_list2)):
            product_list.append(basis_list1[i] * basis_list2[j])
    return product_list

def find_E(S_list, F_list):
    tmp = matrix(Fqm.base_ring(),[vector(Fqm.gen()**i) for i in range(m)]).row_space()
    for Fi in F_list:
        FiS_list = [(1/Fi)*i for i in S_list]
        SSBi = matrix(Fqm.base_ring(),[vector(i) for i in FiS_list])
        tmp = tmp.intersection(SSBi.row_space())
    return tmp

----------------------------------------------------------------------------------------------------------------


# Blockwise_ROLLO_I
#(q,m,n,n1,n2,r1,r2,d1,d2) = (2,61,134,67,67,4,4,5,4) # Blockwise_ROLLO_I - 128
#(q,m,n,n1,n2,r1,r2,d1,d2) = (2,71,158,79,79,5,5,5,5) # Blockwise_ROLLO_I - 192
#(q,m,n,n1,n2,r1,r2,d1,d2) = (2,79,178,89,89,5,5,6,5) # Blockwise_ROLLO_I -256

# Updated Blockwise_ROLLO_I
#(q,m,n,n1,n2,r1,r2,d1,d2) = (2,61,2*67,67,67,4,4,5,5) # Blockwise_ROLLO_I - 128
#(q,m,n,n1,n2,r1,r2,d1,d2) = (2,79,2*89,89,89,5,5,5,6) # Blockwise_ROLLO_I - 192
#(q,m,n,n1,n2,r1,r2,d1,d2) = (2,79,2*97,97,97,5,5,6,6) # Blockwise_ROLLO_I -256

import hashlib
Fqm = GF(q**m)
P1  = GF(q)['XX'].irreducible_element(n1, algorithm="minimal_weight")
P.<XX> = Fqm[]
R.<X> = P.quotient(P1)

# Key Generation
def Blockwise_ROLLO_I_KGen(q,m,n,n1,n2,r1,r2,d1,d2):
    x = R(list(random_small_vec_gen(n1,d1)))
    y = R(list(random_small_vec_gen(n1,d2)))
    h = x**(-1) * y
    return x,y,h

%time Key_List = Blockwise_ROLLO_I_KGen(q,m,n,n1,n2,r1,r2,d1,d2)

# Encapsulation
def Blockwise_ROLLO_I_Encap(Public_Key):  
    e1 = R(list(random_small_vec_gen(n1,r1)))
    e2 = R(list(random_small_vec_gen(n1,r2)))
    E = vector_space(list(e1)+list(e2))
    E_string = "".join(map(str,list(E.basis_matrix()))) # list(E.basis_matrix()) can be replace by E.basis_matrix().list()
    c = e1 + Public_Key * e2
    return c,hashlib.sha256(E_string.encode()).hexdigest()

Public_Key = Key_List[2]
%time Encapsulation_Key_List = Blockwise_ROLLO_I_Encap(Public_Key)

# Decapsulation 
def Blockwise_ROLLO_I_Decap(Private_Key,Ciphertext):  
    xlist = list(Private_Key[0])
    ylist = list(Private_Key[1])
    F1_list = vector_to_basis_list(xlist)
    F2_list = vector_to_basis_list(ylist)
    slist = list(Private_Key[0]*Ciphertext)
    S_list = vector_to_basis_list(slist)
    E_test = find_E(S_list,F1_list) + find_E(S_list,F2_list)
    E_test_string = "".join(map(str,list(E_test.basis_matrix())))
    return hashlib.sha256(E_test_string.encode()).hexdigest()

Private_Key = Key_List[0:2]; Ciphertext = Encapsulation_Key_List[0]
%time Key_test = Blockwise_ROLLO_I_Decap(Private_Key,Ciphertext)

# Check Correctness
Key_test == Encapsulation_Key_List[1]

--------------------------------------------------------------------------------------------------------------

# Blockwise_ROLLO_II

import hashlib

# This function turns Hexadecimal into Binary
def Hexadecimal_to_Binary(Hex_str): 
    bin_str = ""
    for i in Hex_str:
        bin_str += bin(int(i,16))[2:].zfill(4)
    return bin_str

digest = hashlib.sha256().hexdigest()
lamda = len(Hexadecimal_to_Binary(digest)) # the output length of  hashlib.sha256()

# Blockwise_ROLLO_II
#(q,m,n,n1,n2,r1,r2,d1,d2) = (2,89,2*163,163,163,4,4,4,4) # Blockwise_ROLLO_II-128
#(q,m,n,n1,n2,r1,r2,d1,d2) = (2,97,2*179,179,179,4,5,5,5) # Blockwise_ROLLO_II-192
#(q,m,n,n1,n2,r1,r2,d1,d2) = (2,101,2*181,181,181,5,5,5,5) # Blockwise_ROLLO_II-256

# Updated Blockwise_ROLLO_II
#(q,m,n,n1,n2,r1,r2,d1,d2) = (2,89,2*163,163,163,4,4,5,5) # Blockwise_ROLLO_II - 128
#(q,m,n,n1,n2,r1,r2,d1,d2) = (2,97,2*179,179,179,4,5,5,6) # Blockwise_ROLLO_II - 192
(q,m,n,n1,n2,r1,r2,d1,d2) = (2,101,2*191,191,191,5,5,5,6) # Blockwise_ROLLO_II -256

Fqm = GF(q^m)
P1  = GF(q)['XX'].irreducible_element(n1, algorithm="minimal_weight")
P.<XX> = Fqm[]
R.<X> = P.quotient(P1)

# Key Generation
def Blockwise_ROLLO_II_KGen(q,m,n,n1,n2,r1,r2,d1,d2):
    x = R(list(random_small_vec_gen(n1,d1)))
    y = R(list(random_small_vec_gen(n1,d2)))
    h = x**(-1) * y
    return x,y,h

%time Key_List = Blockwise_ROLLO_II_KGen(q,m,n,n1,n2,r1,r2,d1,d2)

# Encryption
def Blockwise_ROLLO_II_Enc(Public_Key,Message):  
    e1 = R(list(random_small_vec_gen(n1,r1)))
    e2 = R(list(random_small_vec_gen(n1,r2)))
    E = vector_space(list(e1)+list(e2))
    E_string = "".join(map(str,list(E.basis_matrix())))
    hash1 = hashlib.sha256(E_string.encode()).hexdigest()
    b1 = Hexadecimal_to_Binary(hash1)
    #b1 = bin(int(hash1,16))[2:]  # [2:]删去二进制串前面的0b,删去后是长度256比特的二进制串
    b1_list = [Fqm(b1[i]) for i in range(lamda)]
    c = e1 + Public_Key * e2
    return c,vector(Fqm.base_ring(),b1_list) + Message

Public_Key = Key_List[2]; Message = random_vector(Fqm.base_ring(),lamda)
%time Ciphertext_List = Blockwise_ROLLO_II_Enc(Public_Key,Message)

# Decryption
def Blockwise_ROLLO_II_Dec(Private_Key,Ciphertext):  
    #H1 = ideal_matrix(Private_key[0]).transpose()
    #H2 = ideal_matrix(Private_key[1]).transpose()
    #H = block_matrix(1,2,[H1,H2])
    xlist = list(Private_Key[0])
    ylist = list(Private_Key[1])
    F1_list = vector_to_basis_list(xlist)
    F2_list = vector_to_basis_list(ylist)
    slist = list(Private_Key[0]*Ciphertext[0])
    S_list = vector_to_basis_list(slist)
    E_test = find_E(S_list,F1_list) + find_E(S_list,F2_list)
    E_test_string = "".join(map(str,list(E_test.basis_matrix())))
    hash2 = hashlib.sha256(E_test_string.encode()).hexdigest()
    b2 = Hexadecimal_to_Binary(hash2)
    #b2 = bin(int(hash2,16))[2:]
    b2_list = [Fqm(b2[i]) for i in range(lamda)]
    return vector(Fqm.base_ring(),b2_list) + Ciphertext[1]

Private_Key = Key_List[0:2]; Ciphertext = Ciphertext_List
%time Message_test = Blockwise_ROLLO_II_Dec(Private_Key,Ciphertext)

# check correctness
Message_test == Message

--------------------------------------------------------------------------------


# Blockwise_ROLLO_III

#(q,m,n,n1,n2,n3,r1,r2,r3,d1,d2) = (2,53,3*73,73,73,73,4,4,5,4,4) # Blockwise_ROLLO_III-128
#(q,m,n,n1,n2,n3,r1,r2,r3,d1,d2) = (2,89,3*101,101,101,101,6,6,6,4,5) # Blockwise_ROLLO_III-192
(q,m,n,n1,n2,n3,r1,r2,r3,d1,d2) = (2,97,3*103,103,103,103,6,6,7,5,5) # Blockwise_ROLLO_III-256

import hashlib
Fqm = GF(q^m)
P1  = GF(q)['XX'].irreducible_element(n/3, algorithm="minimal_weight")
P.<XX> = Fqm[]
R.<X> = P.quotient(P1)

# Key Generation
def Blockwise_ROLLO_III_KGen(q,m,n,n1,n2,n3,r1,r2,r3,d1,d2):
    f1 = R.random_element()
    h0 = R(list(random_small_vec_gen(n1,d1)))
    h1 = R(list(random_small_vec_gen(n1,d2)))
    f0 = h1 + h0*f1
    return h0,h1,f0,f1

%time Key_List = Blockwise_ROLLO_III_KGen(q,m,n,n1,n2,n3,r1,r2,r3,d1,d2)

# Encapsulation 
def Blockwise_ROLLO_III_Encap(Public_Key):
    e0 = R(list(random_small_vec_gen(n1,r1)))
    e1 = R(list(random_small_vec_gen(n1,r2)))
    e = R(list(random_small_vec_gen(n1,r3)))
    E = vector_space(list(e0)+list(e1))
    E_string = "".join(map(str,list(E.basis_matrix())))
    c0 = e + Public_Key[0] * e1
    c1 = e0 + Public_Key[1] * e1
    return c0,c1,hashlib.sha256(E_string.encode()).hexdigest()

Public_Key = Key_List[2:4]
%time Encapsulation_Key_List = Blockwise_ROLLO_III_Encap(Public_Key)

# Decapsulation 
def Blockwise_ROLLO_III_Decap(Private_Key,Ciphertext):  
    h0_list = list(Private_Key[0])
    h1_list = list(Private_Key[1])
    F1_list = vector_to_basis_list(h0_list)
    F2_list = vector_to_basis_list(h1_list)
    slist = list(Ciphertext[0] - Private_Key[0]*Ciphertext[1])
    S_list = vector_to_basis_list(slist)
    E_test = find_E(S_list,F1_list) + find_E(S_list,F2_list)
    E_test_string = "".join(map(str,list(E_test.basis_matrix())))
    return hashlib.sha256(E_test_string.encode()).hexdigest()

Private_Key = Key_List[0:2]; Ciphertext = Encapsulation_Key_List[0:2]
%time Key_test = Blockwise_ROLLO_III_Decap(Private_Key,Ciphertext)

# Check correctness
Key_test == Encapsulation_Key_List[2]


