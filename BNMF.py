import numpy as np
import numpy.matlib as mtl
import scipy.special as spc
import time 


def gnmf_solvebynewton(c, a0, cond):
    
    (M ,N) = np.shape(a0)
    c = np.array(c)
    a0 = np.array(a0)
    if c.ndim > 1:
        (Mc ,Nc) = np.shape(c)
    elif cond == 3:
        Mc = np.size(c)
        Nc = 1
    elif cond == 2:
        Mc = 1
        Nc = np.size(c)
    elif cond == 4:
        Mc = 1
        Nc = 1
    
            
    
    if M == Mc and N == Nc: #no tie
        a1 = a0
#        cond = 1
        
    elif Mc == 1 and Nc > 1: # tie rows
        a1 = a0[0]
#        cond = 2
        
    elif Mc > 1 and Nc == 1 : #tie columns
        a1 = a0.T[0]
#        cond = 3
        
    elif Mc == 1 and Nc == 1:  #tie all
        a1 = a0[0,0]
#        cond = 4
    a1 = np.array(a1)
    for i in range(10):
        a2 = a1 - (np.log(a1) - spc.digamma(a1) + 1 - c )/(1/a1 - spc.polygamma(1,a1))
        idx = np.where(a2<0)
        if np.array(idx).size > 0:
            if isinstance(a2, float) :
                a2 = a1/2
            else  :
                a2[idx] = a1[idx]/2
            
        a1 = a2
        
    if cond == 2:
        a1 = mtl.repmat(a1, M, 1)
    elif cond == 3:
        a1 = mtl.repmat(np.reshape(a1,(M,1)),1,N)
    elif cond == 4:
        a1 = np.ones((M,N))*a1
        
    return a1
    
    
        
    
    


def gnmf_vb_poisson_mult_fast(x, M, a_tm, b_tm, a_ve, b_ve, EPOCH=2000, METHOD="vb", UPDATE= 10000000, tie_a_ve = "clamp", tie_b_ve = "clamp",  tie_a_tm = "clamp", tie_b_tm = "clamp", print_period = 500, t_initializer = "gamma", v_initializer = "gamma") :
    
    if t_initializer == "gamma":
        t_init = np.random.gamma(a_tm, b_tm/a_tm)
    if v_initializer == "gamma":
        v_init = np.random.gamma(a_ve, b_ve/a_ve)
        
    W = x.shape[0]
    
    K = x.shape[1]
    I = b_tm.shape[1]
    
    X = x
    
    L_t = t_init
    L_v = v_init
    E_t = t_init
    E_v = v_init
    Sig_t = t_init
    Sig_v = v_init
    
    B = np.zeros((1,EPOCH))
    
    rec = 0
    g_all = []
    
    gammalnX = spc.gammaln(X + 1)
    
    tic = time.time()
    ticall = time.time()
    
    for e in range(EPOCH) :
        
        LtLv = np.dot(L_t,L_v)
        tmp = X/(LtLv)
        Sig_t = L_t*np.dot(tmp,L_v.T)
        Sig_v = L_v*np.dot(L_t.T, tmp)
        
        alpha_tm = a_tm + Sig_t
        beta_tm = 1/(a_tm/b_tm + np.dot(M,E_v.T) )
        E_t = alpha_tm*beta_tm
        
        alpha_ve = a_ve + Sig_v
        beta_ve = 1/(a_ve/b_ve + np.dot(E_t.T,M))
        E_v = alpha_ve*beta_ve
        
        #COMPUTE THE BOUND
        
        if e % 100 == 1 :
            print("*")
        if e % print_period == 1 or e == EPOCH:
            g = {'E_T' : E_t}
            g['E_logT'] = np.log(L_t)
            g['E_V'] = E_v
            g['E_logV'] = np.log(L_v)
            
            g['Bound'] = -np.sum(M*np.dot(g['E_T'],g['E_V']) + gammalnX) 
            + np.sum(-X*  ( (np.dot(L_t*g['E_logT'],L_v)+ np.dot(L_t,L_v*g['E_logV'])) / LtLv - np.log(LtLv)))
            +np.sum( -a_tm/b_tm*g['E_T'] -  spc.gammaln(a_tm) +  a_tm*np.log(a_tm/b_tm)  )
            +np.sum( -a_ve/b_ve*g['E_V'] -  spc.gammaln(a_ve) +  a_ve*np.log(a_ve/b_ve)  )
            +np.sum( spc.gammaln(alpha_tm) + alpha_tm*(np.log(beta_tm) + 1)      )
            +np.sum( spc.gammaln(alpha_ve) + alpha_ve*(np.log(beta_ve) + 1)      )
            
            g['a_ve'] = a_ve
            g['b_ve'] = b_ve
            g['a_tm'] = a_tm
            g['b_tm'] = b_tm
            
            g_all = np.append(g_all,g)
            print('\nBound {} \t a_ve = {} \t b_ve = {} \t a_tm = {}\t b_tm = {}\n'.format(g['Bound'], g['a_ve'][0,0], g['b_ve'][0,0],g['a_tm'][0,0] ,g['b_tm'][0,0]  ))
            
        if e == EPOCH:
            break
        
        L_t = np.exp(spc.psi(alpha_tm))*beta_tm
        L_v = np.exp(spc.psi(alpha_ve))*beta_ve
        
        if e>UPDATE:
            if ~(tie_a_tm == 'clamp') :
                Z = E_t/b_tm - (np.log(L_t) - np.log(b_tm))
                if tie_a_tm == 'free' :
                    a_tm = gnmf_solvebynewton(Z, a_tm, 1)
                elif tie_a_tm == 'rows' :
                    a_tm = gnmf_solvebynewton(np.sum(Z,0)/W, a_tm, 2)
                elif tie_a_tm == 'cols' :
                    a_tm = gnmf_solvebynewton(np.sum(Z,1)/I, a_tm, 3)
                elif tie_a_tm == 'tie_all' :
                    a_tm = gnmf_solvebynewton(np.sum(Z)/(W*I), a_tm, 4)
                    
            if tie_b_tm == 'free' :
                b_tm = E_t
            elif tie_b_tm == 'rows' :
                b_tm = mtl.repmat(np.sum(a_tm*E_t, 0)/np.sum(a_tm,0), W, 1)
            elif tie_b_tm == 'cols' :
                b_tm = mtl.repmat(np.reshape(np.sum(a_tm*E_t, 1)/np.sum(a_tm,1),(W,1)), 1, I)
            elif tie_b_tm == 'tie_all' :
                b_tm = np.sum(a_tm*E_t)/np.sum(a_tm)*np.ones((W, I))
                
            if ~(tie_a_ve == 'clamp') :
                Z = E_v/b_ve - (np.log(L_v) - np.log(b_ve))
                if tie_a_ve == 'free' :
                    a_ve = gnmf_solvebynewton(Z, a_ve, 1)
                elif tie_a_ve == 'rows' :
                    a_ve = gnmf_solvebynewton(np.sum(Z,0)/I, a_ve, 2)
                elif tie_a_ve == 'cols' :
                    a_ve = gnmf_solvebynewton(np.sum(Z,1)/K, a_ve, 3)
                elif tie_a_ve == 'tie_all' :
                    a_ve = gnmf_solvebynewton(np.sum(Z)/(I*K), a_ve, 4)
                    
            if tie_b_ve == 'free' :
                b_ve = E_v
            elif tie_b_ve == 'rows' :
                b_ve = mtl.repmat(np.sum(a_ve*E_v, 0)/np.sum(a_ve,0), I, 1)
            elif tie_b_ve == 'cols' :
                b_ve = mtl.repmat(np.reshape(np.sum(a_ve*E_v, 1)/np.sum(a_ve,1),(I,1)), 1, K)
            elif tie_b_ve == 'tie_all' :
                b_ve = np.sum(a_ve*E_v)/np.sum(a_ve)*np.ones((I, K))
        if e % 100 == 0:
            toc = time.time()
            print('number of epoch',e,'  time elapsed',toc-tic)
            tic = time.time()
    tocall = time.time()   
    print('total time elapsed',tocall-ticall)
    return g

x = [[1,0,0,0],[0,1,1,0],[0,0,1,0],[1,1,0,0],[0,0,0,1]]
x = np.array(x)
M = [[0.9,0.1,0.1,0.1],[0.1,0.9,0.9,0.1],[0.1,0.1,0.9,0.1],[0.9,0.9,0.1,0.1],[0.1,0.1,0.1,0.9]]
M = np.array(M)
a_tm = np.ones((5,3))
b_tm = np.ones((5,3))*10
b_ve = np.ones((3,4))*10
a_ve = np.ones((3,4))*10

g = gnmf_vb_poisson_mult_fast(x, M, a_tm, b_tm, a_ve, b_ve,tie_a_tm='free' ,tie_b_tm='free' ,tie_a_ve='free' ,tie_b_ve='free' , UPDATE = 490)
 
#S = [ 4.88370603,  4.88370603,  4.88370603]
#gnmf_solvebynewton(S,a_tm, 2)

    