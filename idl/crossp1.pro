; ===========================================================================
;  crossp1.pro
; ---------------------------------------------------------------------------
; Replacement function for CROSSP function that enables cross-product
; of arrays of vectors 
; ===========================================================================
; ---------------------------------------------------------------------------
  FUNCTION crossp1,v1in,v2in,norm=norm
; ---------------------------------------------------------------------------
; V1in : (3) or (3,n) array
; V2in : (3) or (3,n) array
; /norm = normalize the result
; ---------------------------------------------------------------------------

; testing consistency of inputs

s1 = size(v1in)
s2 = size(v2in)

n1 = (s1(0) eq 1) + (s1(0) eq 2)*s1(2)
n2 = (s2(0) eq 1) + (s2(0) eq 2)*s2(2)

err = bytarr(3)
err(0) = s1(1) ne 3 or n1 eq 0
err(1) = s2(1) ne 3 or n2 eq 0
err(2) = not (n1 eq 1 or n2 eq 1 or (n1 eq n2 and n1 ne 1 and n2 ne 1) )
err = err mod 2

if (where(err ne 0))(0) ne -1 then begin
    print,"CROSSP1 : error, wrong argument size..."
    print," Use : CROSSP1,v1,v2"
    if err(0) then print," > v1 must be a vector (3) or a 2-D array (3,n)"
    if err(1) then print," > v2 must be a vector (3) or a 2-D array (3,n)"
    if err(2) then print," > in case of vector of vectors (3,n1) and (3,n2), n1=n2 is required"
    return,-1
endif

; preparation of data

n3 = max([n1,n2])
v3 = dblarr(3,n3)
v1 = rebin(reform(double(v1in),3,n1),3,n3)
v2 = rebin(reform(double(v2in),3,n2),3,n3)

; cross product

for i=0l,n3-1l do v3(*,i) = crossp(v1(*,i),v2(*,i))

; normlization (if required)

if keyword_set(norm) then v3 = v3 / rebin(reform(sqrt(total(v3^2.,1)),1,n3),3,n3)

return,reform(v3)
end
