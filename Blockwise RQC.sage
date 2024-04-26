def random_small_vec_gen(n,t):
    B = matrix(Fqm.base_ring(),t,m,0)
    while B.rank() != t:
        B = matrix(Fqm.base_ring(),[vector(Fqm.random_element()) for i in range(t)])
    C = matrix(Fqm.base_ring(),n,t,0)
    while C.rank() != t:
        C = matrix(Fqm.base_ring(),n,t,[Fqm.base_ring().random_element() for _ in range(n*t)]) 
    CC = C * B
    return vector(Fqm,[CC[i] for i in range(n)])

def rank_R(z):
    return matrix(Fqm.base_ring(),n,m,[vector(z[i]) for i in range(n)]).rank()

def Encoding_Gabidulin(Message, Gabidulin_Support):
    f = S(Message.list())  # The message polynomial 
    return vector(f.multi_point_evaluation(Gabidulin_Support))

def Decoding_Gabidulin(Noisy_Word, Gabidulin_Support,r): 
    g_monomials = [Gabidulin_Support[i]**(q**j) for i in range(n) for j in range(k+r)] 
    SC1 = matrix(Fqm,n,k+r,g_monomials) 
    y_monomials = [Noisy_Word[i]**(q**j) for i in range(n) for j in range(r+1)] 
    SC2 = matrix(Fqm,n,r+1,y_monomials) 
    SC = block_matrix(Fqm,1,2,[SC1,SC2])
    Solution = SC.right_kernel_matrix().list()  
    N = S(Solution[0:k+r])
    V_vector = vector(Solution[k+r:k+2*r+1])
    V = S(list(-V_vector))
    ff,re = N.left_quo_rem(V)
    return vector(ff.list())
    
# Key Generation
def Blockwise_RQC_KGen(q,m,n,k,w_x,w_y,w_r1,w_r2,w_e):
    h = R.random_element()
    x = R(list(random_small_vec_gen(n,w_x)))
    y = R(list(random_small_vec_gen(n,w_y)))
    s = x + h*y
    pk = [h,s]; sk = [x,y]
    return pk,sk

# Encryption
def Blockwise_RQC_Enc(Public_Key,Message,Gabidulin_Support):  
    r1 = R(list(random_small_vec_gen(n,w_r1)))
    r2 = R(list(random_small_vec_gen(n,w_r2)))
    e = R(list(random_small_vec_gen(n,w_e)))
    u = r1 + Public_Key[0]*r2
    v = Encoding_Gabidulin(Message,Gabidulin_Support)+  vector(e + Public_Key[1]*r2)
    ct = [vector(u),v]
    return ct

# Decryption
def Blockwise_RQC_Dec(Private_Key,Ciphertext,Gabidulin_Support,r):  
    u = R(list(Ciphertext[0]))
    Noisy_Word = Ciphertext[1] - vector(u*Private_Key[1])
    return Decoding_Gabidulin(Noisy_Word, Gabidulin_Support,r)


# Blockwise_RQC
(q,m,n,k,w_x,w_y,w_r1,w_r2,w_e) = (2,83,79,7,4,4,4,4,4) # Blockwise_RQC - 128
#(q,m,n,k,w_x,w_y,w_r1,w_r2,w_e) = (2,127,113,3,5,5,5,5,5) # Blockwise_RQC - 192
#(q,m,n,k,w_x,w_y,w_r1,w_r2,w_e) = (2,139,137,3,5,5,6,6,7) # Blockwise_RQC - 256


Fqm = GF(q**m)
P1  = GF(q)['XX'].irreducible_element(n, algorithm = "minimal_weight")
P.<XX> = Fqm[]
R.<X> = P.quotient(P1)
Frob = Fqm.frobenius_endomorphism()
S = OrePolynomialRing(Fqm, Frob, 'x')

%time Public_Key, Private_Key = Blockwise_RQC_KGen(q,m,n,k,w_x,w_y,w_r1,w_r2,w_e)

Message = random_vector(Fqm,k); g = random_small_vec_gen(n,n)
%time Ciphertext = Blockwise_RQC_Enc(Public_Key,Message,g)

r = w_x*w_r2 + w_y*w_r1 + w_e
%time Message_test = Blockwise_RQC_Dec(Private_Key,Ciphertext,g,r)

# check correctness
Message_test == Message
