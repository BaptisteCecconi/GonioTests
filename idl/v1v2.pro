; -----------------------------------------------------------------------------
  function V1V2, S,h1,h2,al1i,al2i,be1i,be2i,Q,U,thi,phi,deg=deg
; -----------------------------------------------------------------------------
; autocorrelation & Re (cross-correlation)

if keyword_set(deg) then begin
    al1 = al1i/!radeg
    al2 = al2i/!radeg
    be1 = be1i/!radeg
    be2 = be2i/!radeg
    th = thi/!radeg
    ph = phi/!radeg
endif else begin
    al1 = al1i
    al2 = al2i
    be1 = be1i
    be2 = be2i
    th = thi
    ph = phi
endelse


V12=0.5*S*h1*h2*( (1.+Q)* (cos(al1)*sin(th)-sin(al1)*cos(th)*cos(ph-be1))* $
                  (cos(al2)*sin(th)-sin(al2)*cos(th)*cos(ph-be2))- $
                  U*sin(al2)*sin(ph-be2)* (cos(al1)*sin(th)-sin(al1)*cos(th)*cos(ph-be1))- $
                  U*sin(al1)*sin(ph-be1)* (cos(al2)*sin(th)-sin(al2)*cos(th)*cos(ph-be2))+ $
                  (1.-Q)* sin(al1)*sin(al2)*sin(ph-be1)*sin(ph-be2) )

 return, V12
end
