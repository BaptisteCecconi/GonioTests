; ===========================================================================
;  Quaternion.pro
; ---------------------------------------------------------------------------
; Set of Functions that enables quaternions calculations basics : 
; Q_PROD : product of 2 quaternions or vector of quaternions
; Q_CONJ : conjugate of a quaternion or a vector of quaternions
; Q_MAKE : builds a quaternion or a vector of quaternons with angle(s) and
;          vector(s)
; Q_VMAKE : builds a pure quaternion or vector of quaterion with vector(s)
; Q_ROT : rotation by qrot of qvect
; ===========================================================================
; USE : to compile the whole set of function, just type 'quaternion' :
;
; IDL> quaternion
; The quaternion.pro procedure set has been successfully compiled and loaded !
; (c) BC, Jun 11 2002
; IDL>
;
; ---------------------------------------------------------------------------
  FUNCTION Q_prod,q1in,q2in
; ---------------------------------------------------------------------------
s1 = size(q1in)
s2 = size(q2in)

n1 = (s1(0) eq 1) + (s1(0) eq 2)*s1(2)
n2 = (s2(0) eq 1) + (s2(0) eq 2)*s2(2)

err = bytarr(3)
err(0) = s1(1) ne 4 or n1 eq 0
err(1) = s2(1) ne 4 or n2 eq 0
err(2) = not (n1 eq 1 or n2 eq 1 or (n1 eq n2 and n1 ne 1 and n2 ne 1) )
err = err mod 2

if (where(err ne 0))(0) ne -1 then begin
    print,"Q_PROD : error, wrong argument size..."
    print," Use : Q_PROD,q1in,q2in"
    if err(0) then print," > q1in must be a quaternion (4) or a vector of quaternions (4,n)"
    if err(1) then print," > q2in must be a quaternion (4) or a vector of quaternions (4,n)"
    if err(2) then print," > in case of vector of quaternions (4,n1) and (4,n2), n1=n2 is required"
    return,-1
endif

n3 = max([n1,n2])
q3 = dblarr(4,n3)
q1 = rebin(reform(double(q1in),4,n1),4,n3)
q2 = rebin(reform(double(q2in),4,n2),4,n3)

q3(0,*) = q1(0,*)*q2(0,*) - q1(1,*)*q2(1,*) - q1(2,*)*q2(2,*) - q1(3,*)*q2(3,*)
q3(1,*) = q1(0,*)*q2(1,*) + q2(0,*)*q1(1,*) + q1(2,*)*q2(3,*) - q1(3,*)*q2(2,*)
q3(2,*) = q1(0,*)*q2(2,*) + q2(0,*)*q1(2,*) + q1(3,*)*q2(1,*) - q1(1,*)*q2(3,*)
q3(3,*) = q1(0,*)*q2(3,*) + q2(0,*)*q1(3,*) + q1(1,*)*q2(2,*) - q1(2,*)*q2(1,*)

return,reform(q3)
end

; ---------------------------------------------------------------------------
  FUNCTION Q_conj,q1in
; ---------------------------------------------------------------------------
s1 = size(q1in)

n1 = (s1(0) eq 1) + (s1(0) eq 2)*s1(2) 

if s1(1) ne 4 or n1 eq 0 then begin
    print,"Q_CONJ : error, wrong argument size..."
    print," Use : Q_CONJ,q1in"
    print," > q1in must be a quaternion (4) or a vector of quaternions (4,n)"
    return,-1
endif

q1 = reform(double(q1in),4,n1)
q2 = -q1
q2(0,*) = q1(0,*)

return,reform(q2)
end

; ---------------------------------------------------------------------------
  FUNCTION Q_make,thin,vin,deg=deg
; ---------------------------------------------------------------------------

th = [thin]
sth = size(th)
sv = size(vin)
nth = sth(1)
nv = (sv(0) eq 1) + (sv(0) eq 2)*sv(2)

err = bytarr(3)
err(0) = sth(0) ne 1
err(1) = sv(1) ne 3 or nv eq 0
err(2) = not (nth eq 1 or nv eq 1 or (nth eq nv and nth ne 1 and nv ne 1) )
err = err mod 2

if (where(err))(0) ne -1 then begin
    print,"Q_MAKE : error, wrong argument size..."
    print," Use : Q_MAKE,thin,vin,[/deg]"
    if err(0) then print," > angle (thin) must be a scalar (1) or a vector (n)"
    if err(1) then print," > vector (vin) must be a vector (3) or a vector of vector (3,n)"
    if err(2) then print," > in case of vector of angles thin(4,nth) and vectors vin(4,nv), nth=nv is required"
    return,-1
endif

nq = max([nth,nv])

th = rebin(reform(double(th),1,nth),1,nq)
v = rebin(reform(double(vin),3,nv),3,nq)

if keyword_set(deg) then th = th/!radeg
Q = dblarr(4,nq)
Q(0,*) = cos(th/2.)
Q(1,*) = v(0,*) * sin(th/2.)
Q(2,*) = v(1,*) * sin(th/2.)
Q(3,*) = v(2,*) * sin(th/2.)

return,reform(Q)
end

; ---------------------------------------------------------------------------
  FUNCTION Q_vmake,vin
; ---------------------------------------------------------------------------

sv = size(vin)
nv = (sv(0) eq 1) + (sv(0) eq 2)*sv(2)

if sv(1) ne 3 or nv eq 0 then begin
    print,"Q_VMAKE : error, wrong argument size..."
    print," Use : Q_VMAKE,vin"
    print," > vector (vin) must be a vector (3) or a vector of vectors (3,n)"
    return,-1
endif

Q = dblarr(4,nv)
Q(0,*) = 0.d0
Q(1:3,*) = reform(double(vin),3,nv) 

return,reform(Q)
end

; ---------------------------------------------------------------------------
  FUNCTION Q_rot,qrot,qvec
; ---------------------------------------------------------------------------

sr = size(qrot)
sv = size(qvec)

nr = (sr(0) eq 1) + (sr(0) eq 2)*sr(2)
nv = (sv(0) eq 1) + (sv(0) eq 2)*sv(2)

err = bytarr(3)
err(0) = sr(1) ne 4 or nr eq 0
err(1) = sv(1) ne 4 or nv eq 0
err(2) = not (nr eq 1 or nv eq 1 or (nr eq nv and nr ne 1 and nv ne 1) )
err = err mod 2

if (where(err ne 0))(0) ne -1 then begin
    print,"Q_ROT : error, wrong argument size..."
    print," Use : Q_ROT,qrot,qvec"
    if err(0) then print," > qrot must be a quaternion (4) or a vector of quaternions (4,n)"
    if err(1) then print," > qvec must be a quaternion (4) or a vector of quaternions (4,n)"
    if err(2) then print," > in case of vector of quaternions (4,n1) and (4,n2), n1=n2 is required"
    return,-1
endif

qres = Q_prod(Q_prod(qrot,qvec),Q_conj(Qrot))


return,qres
end

; ----------------------------------------------------------------------------
 PRO Quaternion, quiet=quiet
; ----------------------------------------------------------------------------
if not keyword_set(quiet) then begin
    print,"The quaternion.pro procedure set has been successfully compiled and loaded ! "
    print," (c) BC, Jun 11 2002"
endif

end
