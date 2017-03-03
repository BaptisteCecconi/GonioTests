; ----------------------------------------------------------------------------
  FUNCTION angular_distance, v1_in,v2_in 
; ----------------------------------------------------------------------------
; computes the angular distance (arccos of scalar product) between 2 vectors
; ----------------------------------------------------------------------------
; INPUTS : vectors or array of vectors
;   V1_IN = (3) or (3,n)
;   V2_IN = (3) or (3,n)
; NB : if V1_IN is (3,n1) and V2_IN is (3,n2), then n1 and n2 must be equal


; testing consistency of inputs

sv1 = size(v1_in)
sv2 = size(v2_in)

nv1 = (sv1(0) eq 1) + (sv1(0) eq 2)*sv1(2)
nv2 = (sv2(0) eq 1) + (sv2(0) eq 2)*sv2(2)

err = bytarr(3)
err(0) = sv1(1) ne 3 or nv1 eq 0
err(1) = sv2(1) ne 3 or nv2 eq 0
err(2) = ((nv1-nv2)*(nv1-1)*(nv2-1)) ne 0.

if (where(err))(0) ne -1 then begin
    message,'Angular_distance : error, wrong argument size ...',/info
    message,' Use : Result = Angular_Distance(V1,V2)',/info
    if err(0) then message,'> V1 must be a vector (3) or an array (3,n)',/info
    if err(1) then message,'> V2 must be a vector (3) or an array (3,n)',/info
    if err(2) then message,'> in case of vector arrays (3,n1) and (3,n2), n1=n2 is required',/info
    return,-1
endif

; normalization

if nv1 eq 1 then v1 = reform(v1_in,3,1) else v1 = v1_in
if nv2 eq 1 then v2 = reform(v2_in,3,1) else v2 = v2_in

nv = max([nv1,nv2])
v1 = rebin(v1,3,nv)
v2 = rebin(v2,3,nv)

norm1 = sqrt(total(v1^2.,1))
norm2 = sqrt(total(v2^2.,1))
v1 = v1/reform(transpose(rebin(reform([norm1],nv,1),nv,3)),3,nv)
v2 = v2/reform(transpose(rebin(reform([norm2],nv,1),nv,3)),3,nv)

; distance

distance = acos( v1(0,*)*v2(0,*) + v1(1,*)*v2(1,*) + v1(2,*)*v2(2,*) <1.)

distance = reform(distance)

return,distance
end
