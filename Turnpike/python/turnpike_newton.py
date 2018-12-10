#Kevin Wasiluk, 2018
#performs one iteration of babylonian method, after finding least squares soln to v*v2 = D
import numpy
from numpy.linalg import pinv 
from random import randrange
from numpy import sqrt 
from collections import Counter as mset
from _bisect import bisect_left


def turnpike_newton_iterate(v, D, req, possible):
    v2 = least_squares(v,D)
    ret = []
    for i in range(0,len(v)):
        k = float(v[i][0] + v2[i][0])
        k = k/2
        k = round(k,6) #Not one hundered percent sure how to handle numerical issues.
        if(v[i][1] in req):
            k = 1
        if(v[i][1] not in possible):
            k = 0
        ret.append((k, v[i][1]))
    return ret


#P.dot(v)
#find least squares soln to v*v2 = D, should return a vector...
#v and D must be sorted.
def least_squares(v,D):
    A = autoc_matrix(v,D)#np.zeros((len(v),len(v)))
    P = pinv(A)
    dm = mvec(D)
    v2 = P.dot(dm)
    return [(v2[i],v[i][1]) for i in range(0,len(v))]

#|D| == |v|
def autoc_matrix(v,D):
    A = numpy.zeros((len(v),len(v)))
    Dl = list(set(D))
    for i in range(0,len(v)):
        for j in range(0,len(v)):
            z = abs(v[j][1] - v[i][1])
            k = bisect_left(Dl,z)
            if(k < len(v) and v[k][1] == z):
                A[k][j] += v[i][0]
    return A

#Important D remains sorted
def m_element(D):
    c = mset(D)
    return [(c[a],a) for a in c]


#Gives the counts of elements in order.
def mvec(D):
    c = mset(D)
    ret = []
    for a in c:
        ret.append(c[a])
    return numpy.array(ret)

#n = num elements
#calculate generating element := 2*D + n
def aug_distances(D):
    z = sqrt(8*len(D) + 1)
    dn = 0.5
    n = (1+z)*(dn)
    n = int(n)
    return sorted(2*D + [0]*n)


def correct_solution(v,D):
    return mult_magma(v,v) == m_element(D)

#Not a very efficient method.
def mult_magma(v1,v2):
    interm = []
    ret = []
    for e1 in v1:
        for e2 in v2:
            interm.append((e1[0]*e2[0], abs(e1[1] - e2[1])))
    interm = sorted(interm,key=lambda x: x[1])
    #Gather like terms...
    curr = -1
    for i in range(0,len(interm)):
        if(curr > i):
            continue
        j = i
        z = 0
        while(j < len(interm) and interm[i][1] == interm[j][1]):
            z += interm[j][0]
            j+=1
        ret.append((z,interm[i][1]))
        curr = j
    return ret

def reversed(P):
    return [P[len(P)-1] - P[i] for i in xrange(-1,-len(P)-1, -1)]


def ini_placements(W,req,D):
    if(len(W) <= (1+sqrt(len(D)*8 + 1))/2):
        return W,req
    D1 = mset(D)
    D2 = mset(distances(W))
    D3 = mset(distances(req))
    UD = unused_distances(D3,D)
    k = max(UD)
    if(index(W,k) > -1 and index(W,D[len(D)-1]-k) > -1 ):
        return W,req
    elif(index(W,k) > -1):
        req.append(k)
        req.sort()
        return ini_placements(possible_placements(req,D),req,D)
    elif(index(W,D[len(D)-1] - k) > -1):
        req.append(D[len(D)-1] - k)
        req.sort()
        return ini_placements(possible_placements(req,D),req,D)
    return W,req


    for j in range(0,len(D)):
        if(D[j] not in UD or (j > 0 and D[j] == D[j-1])):
            continue
        i = D[j]
        cd2 = D2[i]
        if(cd2 == D1[i] and D3[i] != cd2):
            V = place_points_with_d(i,W,req)
            Wn = possible_placements(V,D)
            if(Wn == W):
                return W,req
            return ini_placements(Wn,V,D)

def count(c,D):
    ret = 0
    i = index(D,c)
    if(i == -1):
        return 0
    while(i < len(D)):
        if c == D[i]:
            ret+=1
            i+=1
        else:
            break
    return ret

def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    else:
        return -1



def initial_placements(D):
    req = [0,D[len(D)-2],D[len(D)-1]]
    return ini_placements(first_try(D),req,D)

def first_try(D):
    Wp = [0,D[len(D)-2],D[len(D)-1]]
    Wp = possible_placements(Wp,D)
    return Wp

def possible_placements(P,D):
    ret = []
    for i in P:
        ret.append(i)
    for i in D:
        if(can_place_point(i,P,D)):
            if(i not in ret):
                ret.append(i)
    ret.sort()
    return ret
def unused_distances(D_p,D):
    ret = mset(D) - mset(D_p)
    return ret

def can_place_point(point,W,D):
    checkW = W[:]
    checkW.insert(bisect_left(W,point),point)
    CWD = distances(checkW)
    UCWD = unique_elements(CWD)
    for i in UCWD:
        if(count(i,CWD) > count(i,D)):
            return False
    return True
#Ok, should be able to just do list(set(X)), not sure why I did this... but there must be a reason.
def unique_elements(D):
    ret = []
    ret.append(D[0])
    for i in range(1,len(D)):
        if(D[i] != D[i-1]):
            ret.append(D[i])
    return ret

def distances(v):
    return sorted([i-j for i in v for j in v if(i-j > 0)])


def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    else:
        return -1



def get_random_arr(num,max_length):
    ret = []
    while(len(ret) < num):
        ri = randrange(0,max_length)
        if(ri not in ret):
            ret.append(ri)
    ret.sort()
    if(ret[0] != 0):
        sub = ret[0]
        for i in range(0,len(ret)):
            ret[i] = ret[i]-sub
    return ret

def mag_element(v,D):
    #return [(0,a) for a in D if a not in v else (1,a)]
    #Problem, D doesn't have 0
    Dl = list(set(D))
    ret = [(1,0)]
    for a in Dl:
        if(a in v):
            ret.append((1,a))
        else:
            ret.append((0,a))
    return ret


def reconstruct(D, maxiter=100):
    z = sqrt(8*len(D) + 1)
    dn = 0.5
    n = (1+z)*(dn)
    n = int(n)
    Da = aug_distances(D)
    v = initial_placements(D)
    if(len(v[0]) == n):
        return v[0]
    else:
        vp = mag_element(v[0],D)
        i = 0
        while(i < maxiter):
            #print(vp)
            vp = turnpike_newton_iterate(vp,Da,v[1],v[0])
            i+=1
        #Warning, this will not work for multiple points at the same location
       # ret = [v[i][1] for i in range(0,len(v)) if round(v[i][0],0) >= 1]
        ret = []
        for i in range(0,len(vp)):
            if(round(vp[i][0],0) >= 1):
                #L = [[vp[i][1]]*int(round(vp[i][0],0))]
                #ret += L
                ret.append(vp[i][1])

        if(distances(ret) == D):
            print("returning here")
            return ret
        else:
            #Fail
            print("Fail")
            return vp,ret


P = get_random_arr(10,100)
D = distances(P)
print("Desired reconstruction:")
print(P)
print("Actual reconstruction:")
print(reconstruct(D))



'''

def get_zero_vector(D):
    ret = []
    for i in D:
        ret.append((0,i))
    return ret






v = [(0,0), (1,1), (1,2), (2,3)]
v2 = [(1,0), (2,1), (1,2), (0,3)]
ret = []
for i in range(0,len(v)):
    k = float(v[i][0] + v2[i][0])
    k = k/2
    ret.append((k, v[i][1]))
print(ret)
print(mult_magma(v,v2))
'''
