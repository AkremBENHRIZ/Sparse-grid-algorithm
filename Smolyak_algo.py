#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
import numpy as np
import matplotlib.pyplot as plt 
import time
from scipy.stats import norm


# In[2]:


def permutliste(seq, er=False):
    p = [seq]
    n = len(seq)
    for k in range(0,n-1):
        for i in range(0, len(p)):
            z = p[i][:]
            for c in range(0,n-k-1):
                z.append(z.pop(k))
                if er==False or (z not in p):
                    p.append(z[:])
    return p


# In[3]:


def possibilite(mu,d): #OK

    
#PART 1  #OK
    if(mu==d):
        return [[1]*d] 
    else:
#PART 2  #OK
        p=(mu-d)
        one=[1]*(d-p) 
       
            
        j=0
        if(d-p>0):
            
            l=[[1]*(p-1)+[p+1]]
        else:
            
            l=[[1]*(d-1)+[p+1]]
        
        while(j<len(l[0])-1):
            for i in range(l[-1][-1-j]-2):
                l1=l[-1].copy()

                l1[-1-j]=l1[-1-j]-1

                l1[-1-(j+1)]=l1[-1-(j+1)]+1
                l.append(l1)
            j+=1

            
        final=[]
        
        j=0
        i=0
        div=l[0][-1]%2
        
#Part 4 #OK

        if(len(l)==2)or(len(l)==1):
            final=l
        else:
            
            while(i<len(l)):
                if(div==1):
                    k=int(l[i][-1-j]/2)
                    if(k>=2):
                        for cpt in range(k+1):
                            
                            final.append(l[i+cpt])
                    else:

                        final.append(l[i])

                    
                   
                    i=l[i][-1-j]+i-1
                    j=j+1

                else:
                    
                    k=int(l[i][-1-j]/2)
                    if(k>=2):
                        for cpt in range(k):
                            
                            final.append(l[i+cpt])
                    else:
                        final.append(l[i])
                    
                    i=l[i][-1-j]+i-1
                    j=j+1

#Part 5 #OK
       
        possibilities=[]
        
        for l1 in final:
            if(d-p>0):
                element=one+l1
            else:
                element=l1
           
            possibilities.append(permutliste(element, True))
          
            
        return possibilities
   


# In[5]:


possibilite(13,9))


# In[ ]:


def Chebyshev(n,x): #OK
    return np.cos((n-1)*np.arccos(x))


# In[ ]:


def basis_function(mu,x): #OK
    m=2**(mu-1)+1
    basis=[]
    for i in range(m):
        basis.append(Chebyshev(n,x))


# In[ ]:


def between(vec): #OK
    
    l=[vec[i] for i in range(1,len(vec),2)]
    
    return l        


# In[ ]:


def smolyak_points(mu): #OK
    points=[]
    points.append([0])
    points.append([-1,1])
    
    for i in range(3,mu+2):
        
        points.append(between([-np.cos(np.pi*(j-1)/2**(i-1)) for j in range(1,2**(i-1)+2)]))
    return points


# In[ ]:


def produit(args, repeat=1): #OK
    iterables = args*repeat  
    result=[()] 
  
    for it in iterables: 
        result = produit_cartesien(result,it) 
  
    return result 
  
def produit_cartesien(iterable1,iterable2): 
    for e1 in iterable1:  
        for e2 in iterable2: 
            yield (e1 + (e2,)) 


# In[ ]:


def Tchebyshev_point(list_of_n):
    point=[]
    X=[]
    for n in list_of_n:
        point.append([np.cos(k *np.pi/n) for k in range(n+1)])
    iterateur=produit(point)
                    
    for it in iterateur:
        X.append(list(it))
    return X


# In[ ]:


def smolyak_grid(mu,d): #OK
    
    SP=smolyak_points(mu)
    
    X=[]
    
    
    for i in range(d,d+mu+1):
           
            every=possibilite(i,d)
            
           
            if(i==d):
               
                X.append([0]*d)
            else:
                for ensemble in every:
                    
                    for t in ensemble:
                        please=[]
                        
                        l=list(t)
                       
                        for j in range(len(l)):

                            please.append(SP[l[j]-1])
                        
                        iterateur=produit(please)
                    
                        for it in iterateur:
                            X.append(list(it))
    
    return X


# In[ ]:


def compute_m(n):
    if(n==0):
        return 0
    elif(n==1):
        return 1
    else:
        return 2**(n-1)+1


# In[ ]:


def Smolyak(d,mu,X,b): #OK
    
    final=[]
    
    for x in X:
        
       
        val=0
        entree=0
        previous=0
        start=0
        mat=[0]*len(b)
        for i in range(d,d+mu+1):
           

            every=possibilite(i,d)
           
            for l in every:

                if(entree==0):
                   
                    big=[]
                    for j in range(len(l)):
                       
                        big.append(np.array([Chebyshev(inc,x[j]) for inc in range(compute_m(l[j]-1)+1,compute_m(l[j])+1)]))
                   
                    res=big[0]
                    
                    for j in range(1,len(big)):
                       
                        res=np.dot(res.reshape((-1,1)),big[j])
                   
                    for j in range(previous,previous+len(res.reshape((-1,1)))):
                                #print(j)

                                mat[j]=res.reshape((-1,1))[j-previous][0]
                    previous=previous+len(res.reshape((-1,1)))
                    entree=1

                else:
                   

                    for l1 in l :
                        big=[]
                        for j in range(len(l1)):
                            big.append(np.array([Chebyshev(inc,x[j]) for inc in range(compute_m(l1[j]-1)+1,compute_m(l1[j])+1)]))
                    
                        res=big[0]
                        for j in range(1,len(big)):
                          
                            res=np.dot(res.reshape((-1,1)),big[j].reshape((1,-1)))
                           
                        for j in range(previous,previous+len(res.reshape((-1,1)))):
                               
                                mat[j]=res.reshape((-1,1))[j-previous][0]
                        previous=previous+len(res.reshape((-1,1)))

        final.append((mat*b).sum())
       
    return final


# In[ ]:


def calibration(true_V,X,mu): #OK
    
    d=len(X[0])
    M=len(X)
    mat=np.zeros((M,M))
    k=-1
    
    
    for x in X:
        
        entree=0
        previous=0
        k=k+1
       
        for i in range(d,d+mu+1):
            
            every=possibilite(i,d)
            
            for l in every:
                if(entree==0):
                    big=[]
                    for j in range(len(l)):
                        big.append(np.array([Chebyshev(inc,x[j]) for inc in range(compute_m(l[j]-1)+1,compute_m(l[j])+1)]))
                    res=big[0]
                    for j in range(1,len(big)):
                        res=np.dot(res.reshape((-1,1)),big[j].reshape((1,-1)))
                    for j in range(previous,previous+len(res.reshape((-1,1)))):
                        mat[k,j]=res[j]
                    previous=len(res.reshape((-1,1)))
                    entree=1
                else:
                    
                    
                    for l1 in l :
                        big=[]
                        for j in range(len(l1)):
                            big.append(np.array([Chebyshev(inc,x[j]) for inc in range(compute_m(l1[j]-1)+1,compute_m(l1[j])+1)]))
                   
                        res=big[0]
                        for j in range(1,len(big)):
                            res=np.dot(res.reshape((-1,1)),big[j].reshape((1,-1)))
                       
                        for j in range(previous,previous+len(res.reshape((-1,1)))):
                         
                            mat[k,j]=res.reshape((-1,1))[j-previous]
                        previous=previous+len(res.reshape((-1,1)))
                    
   
    return np.dot(np.linalg.inv(mat),true_V)

