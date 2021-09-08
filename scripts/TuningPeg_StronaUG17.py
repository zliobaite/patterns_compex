import copy
from numpy import array,zeros,mean,std,arange,delete,where,argsort,ones,sort
from random import randrange,random,sample,shuffle,randint
from os import listdir
from copy import deepcopy
from math import isnan
from numpy.random import binomial,normal,choice

import pdb

#Function to compute C-score
def CH(mch):    #the higher the more checkerboardness
    ch=[]
    for i in range(len(mch)):
        for j in range(len(mch)):
            Ri=mch[i].sum()
            Rj=mch[j].sum()
            Sh=(mch[i]*mch[j]).sum()
            if i<j:
                #ch.append(((Ri - Sh)*(Rj - Sh))/float(Ri*Rj))    #use this line for normalized C-score
                ch.append((Ri - Sh)*(Rj - Sh))
    return mean(ch)




###NODF functions
def p_func(l):
    if len(l[0])<=len(l[1]):
        return 0
    else:
        return len(l[0]&l[1])/float(len(l[1]))




def pack(m):
    col_ord=argsort(m.sum(0))[::-1]
    row_ord=argsort(m.sum(1))[::-1]
    pm=m[:,col_ord]
    pm=pm[row_ord,:]
    pm=delete(pm,where(pm.sum(0)==0),1)
    pm=delete(pm,where(pm.sum(1)==0),0)
    return pm





def NODF(a):
    aaa=pack(array(a))
    aa=[set(where (i==1)[0]) for i in aaa]
    b=[[aa[i] for j in range(len(aa)-i-1)] for i in range(len(aa))]
    c=[aa[i+1:] for i in range(len(aa)-1)]
    bb=[item for sublist in b for item in sublist] 
    cc=[item for sublist in c for item in sublist] 
    d=[[bb[i],cc[i]] for i in range(len(bb))]
    mr=map(p_func,d)
    aaa=aaa.transpose()
    aa=[set(where (i==1)[0]) for i in aaa]
    b=[[aa[i] for j in range(len(aa)-i-1)] for i in range(len(aa))]
    c=[aa[i+1:] for i in range(len(aa)-1)]
    bb=[item for sublist in b for item in sublist] 
    cc=[item for sublist in c for item in sublist] 
    d=[[bb[i],cc[i]] for i in range(len(bb))]
    mc=map(p_func,d)
    return mean(mr+mc)



#function to compute row/column discrepancy
def ndisc(exp,obs,sort_='no'):
    if sort_=='no':
        return mean(abs(obs-exp)/exp)
    else:
        obs=sort(obs)
        exp=sort(exp)
        return mean(abs(obs-exp)/exp)



#function to identify max(M_r) and max(M_c)
def randomness(a,sort_='no'):
    exp_c=a.sum(0)
    exp_r=a.sum(1)
    rand_r,rand_c=[],[]
    for rep in range(1000):
        ra=EE(a)
        obs_c=ra.sum(0)
        obs_r=ra.sum(1)
        rand_c.append(ndisc(exp_c,obs_c,sort_))
        rand_r.append(ndisc(exp_r,obs_r,sort_))
    return mean(rand_r),mean(rand_c)








####classic null models
def EE(matrix):
    matrix=array(matrix)
    Ctot=matrix.sum(0)
    Rtot=matrix.sum(1)
    R,C=matrix.shape
    occs=matrix.sum()
    fill=occs/float(R*C)
    rm=zeros([R,C])
    while rm.sum()<occs:
        rr,rc=randrange(R),randrange(C)
        if random()<=fill:
            rm[rr][rc]=1
    return rm





###function needed by UG
def bin(n,N):
    if n==N:
        return n
    p=0.5
    rep=n-1
    att=rep*p
    adj=n-att
    if rep>0:
        rval=binomial(rep,p)
    else:
        return 1
    return rval+adj




###unbiased proportional null model by Ulrich and Gotelli (2012)
def UG(M):
    R,C=M.shape
    ct_,rt_=sorted(M.sum(0),reverse=True),sorted(M.sum(1),reverse=True)
    ct=[bin(i,R) for i in ct_]
    rt=[bin(i,C) for i in rt_]
    while sum(ct)<sum(ct_):
        try:
            bin_c=binomial(C-1,0.5)
            if ct[bin_c]<R:
                ct[bin_c]+=1
        except:
            pass
    while sum(ct)>sum(ct_):
        try:
            bin_c=binomial(C-1,0.5)
            if ct[bin_c]>1:
                ct[bin_c]-=1
        except:
            pass
    while sum(rt)<sum(rt_):
        try:
            bin_r=binomial(R-1,0.5)
            if rt[bin_r]<C:
                rt[bin_r]+=1
        except:
            pass
    while sum(rt)>sum(rt_):
        try:
            bin_r=binomial(R-1,0.5)
            if rt[bin_r]>1:
                rt[bin_r]-=1
        except:
            pass
    occs=float(M.sum())
    rm=zeros([R,C],dtype=int)
    RR=[]
    for i in range(len(rt)):
        for j in range(int(rt[i])):
            RR.append(i)
    RC=[]
    for i in range(len(ct)):
        for j in range(int(ct[i])):
            RC.append(i)
    while rm.sum()<occs:
        rr,rc=sample(RR,1)[0],sample(RC,1)[0]
        rm[rr][rc]+=1
    sc=0
    more=where(rm>1)
    # while len(more[0])>0 and sc<100000:
    while len(more[0])>1 and sc<100000:
        a,b=sample(range(len(more[0])),2)
        rr1,rc1=more[0][a],more[1][a]
        rr2,rc2=more[0][b],more[1][b]
        sc+=1
        if rm[rr1][rc2]<(rm[rr1][rc1]) and rm[rr1][rc2]<(rm[rr2][rc2]) and (rm[rr2][rc1])<(rm[rr1][rc1]) and rm[rr2][rc1]<(rm[rr2][rc2]):
            rm[rr1][rc1]-=1
            rm[rr1][rc2]+=1
            rm[rr2][rc1]+=1
            rm[rr2][rc2]-=1
            more=where(rm>1)
    if sc==100000:
        more=where(rm>1)
        while len(more[0])>0:
            z=where(rm==rm.min())
            rr1,rc1=more[0][0],more[1][0]
            a=randrange(len(z[0]))
            rr2,rc2=z[0][a],z[1][a]
            if random()<=(float(ct[rc2])*float(rt[rr2]))/float(occs**2):
                rm[rr1][rc1]-=1
                rm[rr2][rc2]+=1
                more=where(rm>1)
    while 0 in list(rm.sum(0))+list(rm.sum(1)):
        try:
            zc=sample(list(where(rm.sum(0)==0)[0]),1)[0]
        except:
            zc=sample(range(C),1)[0]
        try:
            zr=sample(list(where(rm.sum(1)==0)[0]),1)[0]
        except:
            zr=sample(range(R),1)[0]

        rc=sample(list(where(rm.sum(0)>1)[0]),1)[0]
        rr=sample(list(where(rm.sum(1)>1)[0]),1)[0]
        rm[rr][rc]-=1
        rm[zr][zc]+=1
    return CB(rm)

#####################################################
#####################################################
###unbiased proportional null model by Ulrich and Gotelli (2012)
### cutting off if not getting anywhere
def UG_mod(M, sc_max=100000, max_budget_loop=-1):
    log = False
    R,C=M.shape
    ct_,rt_=sorted(M.sum(0),reverse=True),sorted(M.sum(1),reverse=True)
    pmc = sum(ct_ > mean(ct_))/len(ct_)
    pmr = sum(rt_ > mean(rt_))/len(rt_)
    ##### LOG MATRIX
    if log:
        nb0c = len([c for c in ct_ if c==0])
        nb0r = len([r for r in rt_ if r==0])
        nb1c = len([c for c in ct_ if c==1])
        nb1r = len([r for r in rt_ if r==1])
        multir = [r for r in rt_ if r>1]
        multic = [c for c in ct_ if c>1]
        print("-- UG on matrix %d %.4f %.4f %.4f\tC %d/%d/%d\tR %d/%d/%d" % (M.sum(), sum(multir)/(len(multir)*len(multic)), pmr, pmc, C-nb0c-nb1c, nb1c, nb0c, R-nb0r-nb1r, nb1r, nb0r))        
        # print("\t%s\n\t%s" % (multir, multic))
    # #####
    ct= array([bin(i,R) for i in ct_], dtype=float)
    rt= array([bin(i,C) for i in rt_], dtype=float)
    # if log: print("I) Assignment of margin totals")
    b=0
    # succes = 0
    while sum(ct)<sum(ct_):
        b+=1
        if max_budget_loop > 0 and b > max_budget_loop: return None
        try:
            # bin_c=binomial(C-1,0.5)
            bin_c=binomial(C-1, pmc)
            if ct[bin_c]<R:
                ct[bin_c]+=1
                # succes += 1
        except:
            pass
    # if log: print("IIa) Adjustment of column totals, up", b, succes)
    b=0
    # succes = 0
    while sum(ct)>sum(ct_):
        b+=1
        if max_budget_loop > 0 and b > max_budget_loop: return None
        try:
            # bin_c=binomial(C-1,0.5)
            bin_c=binomial(C-1, pmc)
            if ct[bin_c]>1:
                ct[bin_c]-=1
                # succes += 1
        except:
            pass
    # if log: print("IIb) Adjustment of column totals, down", b, succes)
    b=0
    # succes = 0
    while sum(rt)<sum(rt_):
        b+=1
        if max_budget_loop > 0 and b > max_budget_loop: return None
        try:
            # bin_r=binomial(R-1,0.5)
            bin_r=binomial(R-1, pmr)
            if rt[bin_r]<C:
                rt[bin_r]+=1
                # succes += 1
        except:
            pass
    # if log: print("IIc) Adjustment of row totals, up", b, succes)
    b=0
    # succes = 0
    while sum(rt)>sum(rt_):
        b+=1
        if max_budget_loop > 0 and b > max_budget_loop: return None
        try:
            # bin_r=binomial(R-1,0.5)
            bin_r=binomial(R-1, pmr)
            if rt[bin_r]>1:
                rt[bin_r]-=1
                # succes += 1
        except:
            pass
    # if log: print("IId) Adjustment of row totals, down", b, succes)
    # b=0
    occs=float(M.sum())
    rm=zeros([R,C],dtype=int)
    RR=[]
    for i in range(len(rt)):
        for j in range(int(rt[i])):
            RR.append(i)
    RC=[]
    for i in range(len(ct)):
        for j in range(int(ct[i])):
            RC.append(i)
    while rm.sum()<occs:
        # b+=1
        if max_budget_loop > 0 and b > max_budget_loop: return None
        rr,rc=sample(RR,1)[0],sample(RC,1)[0]
        rm[rr][rc]+=1
    # if log: print("III) Placement of cell occurrences", b)
    sc=0
    more=where(rm>1)
    # while len(more[0])>0 and sc<100000:
    while len(more[0])>1 and sc<sc_max:
        a,b=sample(range(len(more[0])),2)
        rr1,rc1=more[0][a],more[1][a]
        rr2,rc2=more[0][b],more[1][b]
        sc+=1
        if rm[rr1][rc2]<(rm[rr1][rc1]) and rm[rr1][rc2]<(rm[rr2][rc2]) and (rm[rr2][rc1])<(rm[rr1][rc1]) and rm[rr2][rc1]<(rm[rr2][rc2]):
            rm[rr1][rc1]-=1
            rm[rr1][rc2]+=1
            rm[rr2][rc1]+=1
            rm[rr2][rc2]-=1
            more=where(rm>1)
    # if log: print("IVa) SSR reduction", sc)
    b=0
    if sc==sc_max:
        more=where(rm>1)
        if log: print("IVb) Nb irreducibles", len(more[0]))
        while len(more[0])>0:
            b+=1
            if max_budget_loop > 0 and b > max_budget_loop:
                if log: print("IVb) Stop irreducibles, left", len(more[0]))
                return None
            z=where(rm==rm.min())
            rr1,rc1=more[0][0],more[1][0]
            # a=randrange(len(z[0]))
            # rr2,rc2=z[0][a],z[1][a]
            # if random()<=(float(ct[rc2])*float(rt[rr2]))/float(occs**2):
            #     rm[rr1][rc1]-=1
            #     rm[rr2][rc2]+=1
            #     more=where(rm>1)
            ### MODIFIED
            #######################
            ps = rt[z[0]]*ct[z[1]]/occs**2
            a = choice(len(ps), p=ps/sum(ps))
            rr2,rc2=z[0][a],z[1][a]
            rm[rr1][rc1]-=1
            rm[rr2][rc2]+=1
            more=where(rm>1)

    # if log: print("IVb) Reassign irreducible", b)
    # b=0
    while 0 in list(rm.sum(0))+list(rm.sum(1)):
        # b+=1
        # if max_budget_loop > 0 and b > max_budget_loop: return None
        try:
            zc=sample(list(where(rm.sum(0)==0)[0]),1)[0]
        except:
            zc=sample(range(C),1)[0]
        try:
            zr=sample(list(where(rm.sum(1)==0)[0]),1)[0]
        except:
            zr=sample(range(R),1)[0]
            
        if random() > 0.5:
            rc=sample(list(where(rm.sum(0)>1)[0]),1)[0]
            rr=sample(list(where(rm[:,rc]>0)[0]),1)[0]
        else:
            rr=sample(list(where(rm.sum(1)>1)[0]),1)[0]
            rc=sample(list(where(rm[rr,:]>0)[0]),1)[0]

        rm[rr][rc]-=1
        rm[zr][zc]+=1
    # if log: print("IVc) Fill zeros", b)
    ## reorder columns and rows
    ros = argsort(argsort(-M.sum(1)))
    cos = argsort(argsort(-M.sum(0)))
    return CB(rm[ros,:][:,cos])
    # return CB(rm)
#####################################################
#####################################################


def PP(matrix):
    matrix=array(matrix)
    Ctot=matrix.sum(0)
    Rtot=matrix.sum(1)
    R,C=matrix.shape
    occs=matrix.sum()
    rm=zeros([R,C])
    while rm.sum()<occs or (0 in list(rm.sum(0))+list(rm.sum(1))):
        rr,rc=randrange(R),randrange(C)
        if random()<=(float(Rtot[rr])*float(Ctot[rc]))/float(occs**2):
            rm[rr][rc]=1
    while rm.sum()>occs:
        rr,rc=randrange(R),randrange(C)
        if rm.sum(0)[rc]>1 and rm.sum(1)[rr]>1:
            rm[rr][rc]=0
    return rm




def PE(matrix):
    matrix=array(matrix)
    Ctot=matrix.sum(0)
    Rtot=matrix.sum(1)
    R,C=matrix.shape
    occs=matrix.sum()
    rm=zeros([R,C])
    while rm.sum()<occs:
        rr,rc=randrange(R),randrange(C)
        if random()<=float(Rtot[rr])/float(occs*C):
            rm[rr][rc]=1
    return rm




def PF(matrix):
    matrix=array(matrix)
    Ctot=matrix.sum(0)
    Rtot=matrix.sum(1)
    R,C=matrix.shape
    occs=matrix.sum()
    rm=zeros([R,C])
    for c in range(C):
        while rm.sum(0)[c]<Ctot[c]:
            rr=randrange(R)
            if random()<=float(Rtot[rr])/float(occs):
                rm[rr][c]=1
    return rm





def EP(matrix):
    matrix=array(matrix)
    Ctot=matrix.sum(0)
    Rtot=matrix.sum(1)
    R,C=matrix.shape
    occs=matrix.sum()
    rm=zeros([R,C])
    while rm.sum()<occs:
        rr,rc=randrange(R),randrange(C)
        if random()<=float(Ctot[rc])/float(occs*R):
            rm[rr][rc]=1
    return rm




def FP(matrix):
    matrix=array(matrix)
    Ctot=matrix.sum(0)
    Rtot=matrix.sum(1)
    R,C=matrix.shape
    occs=matrix.sum()
    rm=zeros([R,C])
    for r in range(R):
        while rm.sum(1)[r]<Rtot[r]:
            rc=randrange(C)
            if random()<=float(Ctot[rc])/float(occs):
                rm[r][rc]=1
    return rm







def EF(matrix):
    matrix=array(matrix)
    matrix=matrix.transpose()
    tsM=[]
    for i in matrix:
        row=list(i)
        shuffle(row)
        tsM.append(list(row))
    return array(tsM).transpose() 




def FE(matrix):
    sM=[]
    for i in matrix:
        sr=list(i)
        shuffle(sr)
        sM.append(list(sr))
    return array(sM)




###FF null model, Curveball (Strona et al. 2014)

def find_presences(input_matrix):
    num_rows, num_cols = len(input_matrix),len(input_matrix[0])
    hp = []
    iters = num_rows if num_cols >= num_rows else num_cols
    if num_cols >= num_rows:
        input_matrix_b = input_matrix
    else:
        input_matrix_b = input_matrix.T
    for r in range(iters):
        hp.append(where(array(input_matrix_b[r]) == 1)[0])
    return hp



def CB(m,num_iterations=-1):
    r_hp=find_presences(m)
    num_rows, num_cols = len(m),len(m[0])
    l = range(len(r_hp))
    num_iters = 5 * min(num_rows, num_cols)
    if num_iterations > num_iters:
        num_iters = num_iterations
    # num_iters = 5 * min(num_rows, num_cols) if num_iterations == -1 else 
    for rep in range(num_iters):
        ab = sample(l, 2)
        a = ab[0]
        b = ab[1]
        ab = set(r_hp[a]) & set(r_hp[b])
        a_ba = set(r_hp[a]) - ab
        if len(a_ba) != 0:
            b_aa = set(r_hp[b]) - ab
            if len(b_aa) != 0:
                ab = list(ab)
                a_ba = list(a_ba)
                b_aa = list(b_aa)
                shuffle(a_ba)
                shuffle(b_aa)
                swap_extent = randint(1, min(len(a_ba), len(b_aa)))
                r_hp[a] = ab+a_ba[:-swap_extent]+b_aa[-swap_extent:]
                r_hp[b] = ab+b_aa[:-swap_extent]+a_ba[-swap_extent:]
    out_mat = zeros([num_rows, num_cols], dtype='int8') if num_cols >= num_rows else zeros([num_cols,num_rows], dtype='int8')
    for r in range(min(num_rows, num_cols)):
        out_mat[r, r_hp[r]] = 1
    result = out_mat if num_cols >= num_rows else out_mat.T
    return result
