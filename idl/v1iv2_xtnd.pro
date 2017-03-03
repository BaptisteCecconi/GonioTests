; -------------------
; integrant functions
; -------------------

function F_cosine,xx
common param, gam,kk
const = 4./!pi-8./!pi^2
return, cos(!pi/2.*tan(xx)/tan(gam))*sin(kk*xx)/const
end

function F_cosine_small_source,xx
common param, gam,kk
const = 4./!pi-8./!pi^2
return, cos(!pi/2.*xx/gam)*sin(kk*xx)/const
end

function F_spherical,xx
common param, gam,kk
const =  2./3.
return, sqrt(1.-(tan(xx)/tan(gam))^2.)*sin(kk*xx)/const
end

function F_spherical_small_source,xx
common param, gam,kk
const = 2./3.
return, sqrt(1.-xx^2./gam^2.)*sin(kk*xx)/const
end

function F_gaussian_kk1,xx
common param, gam,kk
const = 1./alog(2.)
return,exp(-alog(2)*(xx/tan(gam))^2.)*(1.+xx)^(-1.5)/const
end

function F_gaussian_kk2,xx
common param, gam,kk
const = 1./alog(2.)
return,exp(-alog(2)*(xx/tan(gam))^2.)*(1.+xx)^(-2.5)/const
end

function F_gaussian_kk3,xx
common param, gam,kk
const = 1./alog(2.)
return,exp(-alog(2)*(xx/tan(gam))^2.)*(1.+xx)^(-2.5)*(3.-xx)/2./const
end

function F_gaussian,xx
common param, gam,kk
const = 1./alog(2.)
return,exp(-alog(2)*(tan(xx)/tan(gam))^2.)*sin(kk*xx)/const
end

function F_gaussian_small_source,xx
common param, gam,kk
const = 1./alog(2.)
return,exp(-alog(2)*(xx/gam)^2.)*sin(kk*xx)/const
end

; -----------------------------------------------------------------------------
  function V1iV2_xtnd, S,h1,h2,al1i,al2i,be1i,be2i,V,thi,phi,gami, $
  	deg=deg,src_model=src_model
; -----------------------------------------------------------------------------
; Im (cross-correlation) 
common param,gam,kk

if keyword_set(deg) then begin
    al1 = al1i/!radeg
    al2 = al2i/!radeg
    be1 = be1i/!radeg
    be2 = be2i/!radeg
    th = thi/!radeg
    ph = phi/!radeg
    gam = gami/!radeg
endif else begin
    al1 = al1i
    al2 = al2i
    be1 = be1i
    be2 = be2i
    th = thi
    ph = phi
    gam = gami
endelse

if not keyword_set(src_model) then src_model = 'uniform'

A1 = -sin(al1)*cos(th)*cos(ph-be1) + cos(al1)*sin(th)
A2 = -sin(al2)*cos(th)*cos(ph-be2) + cos(al2)*sin(th)

B1 = -sin(al1)*sin(ph-be1)
B2 = -sin(al2)*sin(ph-be2)


case src_model of 
  'uniform' : begin
    G2 = 1. + cos(gam)
  end
  'spherical' : begin
    kk = 2 & G2 = qromb('f_spherical',0.0,gam)
    if gam le 0.001 then begin
      G2 = G2/(gam^2./2)
    endif else begin
      G2 = G2/(1-cos(gam))
    endelse
  end
  'spherical_small_source' : begin
    kk = 2 & G2 = qromb('f_spherical_small_source',0.0,gam)
    if gam le 0.001 then begin
      G2 = G2/(gam^2./2)
    endif else begin
      G2 = G2/(1-cos(gam))
    endelse
  end
  'gaussian' : begin
    kk = 2 & G2 = qromb('f_gaussian',0.0,!pi/2.)
    if gam le 0.001 then begin
      G2 = G2/(gam^2./2)
    endif else begin
      G2 = G2/(1-cos(gam))
    endelse
  end
  'gaussian_small_source' : begin
    kk = 2 & G2 = qromb('f_gaussian_small_source',0.0,!pi/2.)
    if gam le 0.001 then begin
      G2 = G2/(gam^2./2)
    endif else begin
      G2 = G2/(1-cos(gam))
    endelse
  end
endcase

V1i2=0.5*S*h1*h2*( V*(-A1*B2 + A2*B1)*G2/2. )

 return, V1i2
end

