#Written by Kevin Wasiluk, some written pre-2018
#call function acaaf_turnpike for the always correct, almost always fast algorithm.
from random import randrange
from collections import Counter as mset
from _bisect import *

def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    else:
        return -1

def first_try(D):
	Wp = [0,D[len(D)-2],D[len(D)-1]]
	Wp = possible_placements(Wp,D)
	return Wp

def gen_rand_get_ft(num,length):
	A = get_random_arr(num,length)
	D = distances(A)
	W = first_try(D)
	return D,W
#Could be written as a one liner, but oh well.
def distances(A):
	ret = []
	for i in range(len(A)):
		for j in range(i+1,len(A)):
			ret.append(A[j]-A[i])
	ret.sort()
	return ret

def unused_distances(D_p,D):
	ret = mset(D) - mset(D_p)
	return ret

def satisfies_distances(W,D):
	D1 = distances(W)
	D = list(D)
	D.sort()
	return D1 == D


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


def can_place_point(point,W,D):
	checkW = W[:]
	checkW.insert(bisect_left(W,point),point)
	CWD = distances(checkW)
	UCWD = unique_elements(CWD)
	for i in UCWD:
		if(count(i,CWD) > count(i,D)):
			return False
	return True


def generating_polynomial(D):
	R.<x> = PolynomialRing(ZZ)
	ret = (1 + sqrt((8*len(D)) + 1))/2
	for i in D:
		ret+= x^i + x^(-i)
	return ret

def generating_polynomial_numerator(D):
	ret = generating_polynomial(D)
	return ret.numerator()

def pol_no_inv(A):
	R.<x> = PolynomialRing(ZZ)
	ret = 0
	for i in A:
		ret += x^i
	return ret


def factored_distances(D):
	return factor(generating_polynomial(D))


def reversed_points(W):
	ret = []
	for i in W:
		ret.append(W[len(W)-1] - i)
	ret.sort()
	return ret


def initial_placements(D):
	req = [0,D[len(D)-2],D[len(D)-1]]
	return ini_placements(first_try(D),req,D)

#appears to benchmark very similar to 'unoptimized' version

def b_initial(D):
	req = [0,D[len(D)-2], D[len(D)-1]]
	W = first_try(D)
	D1 = mset(D)
	D2 = mset(distances(W))
	D3 = mset(distances(req))
	UD = unused_distances(D3,D)
	n = (1+sqrt(len(D)*8 + 1))/2
	m = D[len(D)-1]
	return ini_placements_m(req,UD,n,m)

def can_place_point_UD(k,req,UD):
	DD = []
	for i in req:
		if i > k:
			DD.append(i-k)
		else:
			DD.append(k-i)
	DD.sort()
	for i in DD:
		if count(i,DD) > UD[i]:
			return false
	return true

#Different implementation of sparse heuristic.
def ini_placements_m(req,UD,n,m):
	if(len(req) == n):
		return req
	k = max(k for k, v in UD.iteritems() if v != 0) # This is a lot slower than it ought to be
	dk = m - k
	a1 = can_place_point_UD(k,req,UD)
	a2 = can_place_point_UD(dk,req,UD)
	a3 = a1 and a2
	if(a3):
		if(UD[k] == 2):
			if(k == dk):
				req.insert(bisect_left(req,k),k)
				return req
			for i in req:
				if(i > k):
					UD[i-k] = UD[i-k] - 1
				else:
					UD[k-i] = UD[k-i] - 1
				if(i > dk):
					UD[i-dk] = UD[i-dk] - 1
				else:
					UD[dk-i] = UD[dk-i] - 1
			req.insert(bisect_left(req,k),k)
			req.insert(bisect_left(req,dk),dk)
			return ini_placements_m(req,UD,n,m)
		else:
			return req
	elif(a1):
		for i in req:
			if(i > k):
				UD[i-k] = UD[i-k] - 1
			else:
				UD[k-i] = UD[k-i] - 1
		req.insert(bisect_left(req,k),k)
		return ini_placements_m(req,UD,n,m)
	elif(a2):
		for i in req:
			if(i > dk):
				UD[i-dk] = UD[i-dk] - 1
			else:
				UD[dk-i] = UD[dk-i] - 1
		req.insert(bisect_left(req,dk),dk)
		return ini_placements_m(req,UD,n,m)
	return req

def unique_elements(D):
    ret = []
    ret.append(D[0])
    for i in range(1,len(D)):
        if(D[i] != D[i-1]):
            ret.append(D[i])
    return ret

#Main method for the sparse heuristic.
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

	return W,req

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

def place_points_with_d(c,W,req):
	for i in range(0,len(W)):
		j = index(W,W[i] + c)
		if(j != -1):
			if(W[j] not in req):
				req.insert(bisect_left(req,W[j]),W[j])
				req.sort()
			if(W[i] not in req):
				req.insert(bisect_left(req,W[i]),W[i])
				req.sort()
	return req


def possible_placements_min(W,D):
	ret = []
	PP = possible_placements(W,D)
	for i in PP:
		if i not in W:
			ret.append(i)
	ret.sort()
	return ret

#Rosenblatt-Seymour algorithm, first get the generating polynomial.
def reconstructions(D):
	gp = generating_polynomial(D)
	gf = gp.factor() 
	r = []
	gf_list = list(gf) # .factor() and this change p(x)q(x) into such a list: [(p(x), 1), (q(x), 1)]
	collec = []
	collecount = []
	for i in range(1,len(gf_list)):
		collec.append(gf_list[i][0])
	for i in range(1,len(gf_list)):
		collecount.append(gf_list[i][1])
	g = 1
	find_all_reconstructions(collec,collecount,r,g,0)
	ret = []
	for i in range(1,len(r)): # doesn't work for g = 1 ? 
		bb = r[i].exponents()
		if satisfies_distances(bb,D):
			ret.append(bb)
	newlist = [i for n,i in enumerate(ret) if i not in ret[:n] and reversed_points(i) not in ret[:n]]
	return newlist


def reconstructions_pol(gp):
	gf = gp.factor() 
	r = []
	gf_list = list(gf) # .factor() and this change p(x)q(x) into such a list: [(p(x), 1), (q(x), 1)]
	collec = []
	collecount = []
	for i in range(1,len(gf_list)):
		collec.append(gf_list[i][0])
	for i in range(1,len(gf_list)):
		collecount.append(gf_list[i][1])
	g = 1
	find_all_reconstructions(collec,collecount,r,g,0)
	ret = []
	for i in range(1,len(r)): # doesn't work for g = 1 ? 
		bb = r[i].exponents()
		ret.append(bb)
	newlist = [i for n,i in enumerate(ret) if i not in ret[:n] and reversed_points(i) not in ret[:n]]
	return newlist

#Tests for homometric sets, not necessary for anything here.
def hom_pairs(L):
	ret = []
	for A in L:
		for B in L:
			if(A!=B and A!=reversed_points(B) and distances(A) == distances(B)):
				ret.append(A)
				ret.append(B)
	return ret

#Tries every multiplication of factors.
def find_all_reconstructions(collec,collecount,ret,g,index):
	ret.append(g)
	if(index >= len(collec)):
		return
	if(sum(collecount) < 2): #base case, we've got it down to p(x)q(x) when sum of collecount = 2
		return
	else:
		for i in range(index,len(collec)):
			if(collecount[i] > 0):
				p = g*collec[i] #new polynomial for recursive case
				newcollecount = collecount[:] # just a copy
				newcollecount[i] = newcollecount[i] - 1 # decrement
				find_all_reconstructions(collec,newcollecount,ret,p,i+1)
	return

#Algorithm that is always correct and almost always fast
def acaaf_turnpike(D):
    I = initial_placements(D)
    if(len(I[0])*(len(I[0])-1)/2 == len(D)):
        return I[0]
    else:
        return reconstructions(D)
#For testing
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
