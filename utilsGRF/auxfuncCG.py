import numpy as np
import itertools

def get_omega_i_S(c,N,site,S,intrinsiccoop=False):
    numerator=''    
    rho_S_den=''
    Kai0=''
    
    Zconf='1'
    for cn in range(1,c+1):
        if cn>1:
            Zconf+='+l%d'%(cn)
    if intrinsiccoop:
        bg='_0'
    else:
        bg=''
    for cn in range(1,c+1):
        Kai0+='K_%d_%d%s'%(cn,site,bg)
        if cn>1:
            Kai0+='*l%d+'%(cn)
        else:
            Kai0+='+'
    Kai0=Kai0.strip('+')
    
    for cn in range(1,c+1):
        
        mu_S=''
        nu0=len(S)
        
        if intrinsiccoop:
            numerator+='K_%d_%d_%s*'%(cn,site,''.join(map(str,S)))
            if nu0>1:#go "upwards" along the graph with kiS,i<S
                for site2 in S:
                    #print(site2)
                    uppersites=np.arange(site2+1,N+1) 
                    nu=len(uppersites)
                    if nu>0:
                        #print(nu)
                        for i in range(1,nu+1):
                            combis=itertools.combinations(np.arange(site2+1,N+1),i)
                            for combi in combis:
                                S2=''.join(map(str,combi))
                                name='K_%d_%d_%s*'%(cn,site2,S2)
                                mu_S+=name
                    else:
                        name='K_%d_%d_0*'%(cn,site2)
                        mu_S+=name
                    
            
            else:
                site2=S[0]
                name='K_%d_%d_0*'%(cn,site2)
                mu_S+=name
        else:
            numerator+='K_%d_%d*%s'%(cn,site,mu_S)
            for site2 in S:
                mu_S+='K_%d_%d*'%(cn,site2)
            
        if cn>1:    
            mu_S+='l%d'%(cn)
        else:
            mu_S=mu_S.strip('*')
            
        numerator+=mu_S+'+'
        
        rho_S_den+=mu_S+'+'
    
    numerator='((%s)*(%s))'%(numerator.strip('+'),Zconf)
    rho_S_den=rho_S_den.strip('+').strip('*')
    
    denominator='((%s)*(%s))'%(rho_S_den,Kai0)
    return [numerator, denominator]

def write_effective_Ks0(c,N,intrinsiccoop=False,pybind=False):
    Knames=[]
    Niter=N+1
    lines=[]
    #effective affinities with nothing else bound
    
    if intrinsiccoop:
        bg='_0'
    else:
        bg=''
    lines=[]
    for i in range(N):
        #numerator
        expr='('
        for cn in range(1,c+1):
            expr+='K_%d_%d%s'%(cn,i+1,bg)
            if cn>1:
                expr+='*l%d+'%(cn)
            else:
                expr+='+'
        expr=expr.strip('+')
        expr+=')/'
        
        #denominator
        expr+='('
        for cn in range(1,c+1):
            if cn>1:
                expr+='l%d+'%(cn)
            else:
                expr+='1+'
        expr=expr.strip('+')
        expr+=')'
        Kn='K%d'%(i+1)
        if pybind:
            line='    T %s=%s;\n'%(Kn,expr)
            #lines.append(line)
            #fh.write(line)
        else:
            line='    %s=%s\n'%(Kn,expr)
            #fh.write(line)
        lines.append(line)    
        Knames.append(Kn)
    return [Knames,lines]

def write_effective_w(c,N,intrinsiccoop=False,pybind=False):
    # effective cooperativities
    wnames=[]
    lines=[]
    for site in range(1,N): #w_i_S
        uppersites=np.arange(site+1,N+1)
        nu=len(uppersites)
       
        #print(nu)
        if nu>0:
            for i in range(1,nu+1):
                combis=itertools.combinations(np.arange(site+1,N+1),i)
                for S in combis:
                    num,den=get_omega_i_S(c,N,site,S,intrinsiccoop=intrinsiccoop)
                    wn='w%d%s'%(site,''.join(map(str,S)))
                    #print(context)
                    if pybind:
                        line='    T %s=%s/%s;\n'%(wn,num,den)
                        #fh.write()
                        #lines.append(line)
                    else:
                        line='    %s=%s/%s\n'%(wn,num,den)
                        #fh.write()
                    lines.append(line)
                    wnames.append(wn)
    return [wnames,lines]

def get_parslist_atconf(c,N,intrinsiccoop=False,samesites=False):
    outstr=''
    if samesites:
        for cn in range(1,c+1): #for each conformation
            outstr+='K_%d,'%(cn)
                        #print('K_%d_%d_%s'%(cn,site,S),end=',')

        for cn in range(1,c+1): #I'm looping again so that these transitions will be printed at the end
                if cn>1:
                    outstr+='l%d,'%(cn)
                    #print('l%d'%(cn),end=',') #l2 is transition from conf1 to conf2, l3 is from 1 to 3 and so on
    else:
        for cn in range(1,c+1): #for each conformation
            for site in range(1,N+1): #for each site
                #combination of contexts, considering K_i,S with i<S
                if intrinsiccoop is True:
                    uppersites=np.arange(site+1,N+1) #S
                    nu=len(uppersites)
                    #print(nu)
                    if nu>0: #at least one site bound
                        for i in range(nu+1):
                            combis=itertools.combinations(np.arange(site+1,N+1),i)
                            for combi in combis:
                                if len(combi)==0:
                                    S='0'
                                else:
                                    S=''.join(map(str,combi))
                                #print('K_%d_%d_%s'%(cn,site,S),end=',')
                                outstr+='K_%d_%d_%s,'%(cn,site,S)
                    else:
                        S='0'
                        outstr+='K_%d_%d_%s,'%(cn,site,S)
                        #print('K_%d_%d_%s'%(cn,site,S),end=',')
                else:
                    outstr+='K_%d_%d,'%(cn,site)
                    #print('K_%d_%d_x'%(cn,site),end=',')
            

        for cn in range(1,c+1): #I'm looping again so that these transitions will be printed at the end
                if cn>1:
                    outstr+='l%d,'%(cn)
                    #print('l%d'%(cn),end=',') #l2 is transition from conf1 to conf2, l3 is from 1 to 3 and so on
    return outstr.strip(',')

        