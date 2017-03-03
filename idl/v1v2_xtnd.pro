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
  function V1V2_xtnd, S,h1,h2,al1i,al2i,be1i,be2i,Q,U,thi,phi,gami, $
  	deg=deg,src_model=src_model,n_int_steps=n_int_steps
; -----------------------------------------------------------------------------
; autocorrelation & Re (cross-correlation)
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

C1 =  sin(al1)*sin(th)*cos(ph-be1) + cos(al1)*cos(th)
C2 =  sin(al2)*sin(th)*cos(ph-be2) + cos(al2)*cos(th)

if not keyword_set(n_int_steps) then n_int_steps=101

case src_model of 
  'uniform' : begin
    G1 = 1.
    G2 = 1. + cos(gam)
    G3 = 4./3.*(1.+cos(gam)+cos(gam)^2.) - 1.
  end
  'spherical' : begin
    kk = 1 & G1 = qromb('f_spherical',0.0,gam)
    kk = 2 & G2 = qromb('f_spherical',0.0,gam)
    kk = 3 & G3 = qromb('f_spherical',0.0,gam)     
    if gam le 0.001 then begin
      G1 = G1/(gam^2./2)
      G2 = G2/(gam^2./2)
      G3 = G3/(gam^2./2)
    endif else begin
      G1 = G1/(1-cos(gam))
      G2 = G2/(1-cos(gam))
      G3 = G3/(1-cos(gam))
    endelse
  end
  'spherical_small_source' : begin
    kk = 1 & G1 = qromb('f_spherical_small_source',0.0,gam)
    kk = 2 & G2 = qromb('f_spherical_small_source',0.0,gam)
    kk = 3 & G3 = qromb('f_spherical_small_source',0.0,gam)     
    if gam le 0.001 then begin
      G1 = G1/(gam^2./2)
      G2 = G2/(gam^2./2)
      G3 = G3/(gam^2./2)
    endif else begin
      G1 = G1/(1-cos(gam))
      G2 = G2/(1-cos(gam))
      G3 = G3/(1-cos(gam))
    endelse
  end
  'gaussian' : begin
    kk = 1 & G1 = qromb('f_gaussian',0.0,!pi/2.)
    kk = 2 & G2 = qromb('f_gaussian',0.0,!pi/2.)
    kk = 3 & G3 = qromb('f_gaussian',0.0,!pi/2.)
    if gam le 0.001 then begin
      G1 = G1/(gam^2./2)
      G2 = G2/(gam^2./2)
      G3 = G3/(gam^2./2)
    endif else begin
      G1 = G1/(1-cos(gam))
      G2 = G2/(1-cos(gam))
      G3 = G3/(1-cos(gam))
    endelse
  end
  'gaussian_small_source' : begin
    kk = 1 & G1 = qromb('f_gaussian_small_source',0.0,!pi/2.)
    kk = 2 & G2 = qromb('f_gaussian_small_source',0.0,!pi/2.)
    kk = 3 & G3 = qromb('f_gaussian_small_source',0.0,!pi/2.)
    if gam le 0.001 then begin
      G1 = G1/(gam^2./2)
      G2 = G2/(gam^2./2)
      G3 = G3/(gam^2./2)
    endif else begin
      G1 = G1/(1-cos(gam))
      G2 = G2/(1-cos(gam))
      G3 = G3/(1-cos(gam))
    endelse
  end
endcase

V12=0.5*S*h1*h2*( (1.+Q) * ( A1*A2 * G2/2.                   $
                           + C1*C2 * (G1-G2/2.)              $
                           )                                 $
                + U * ( A1*B2 + A2*B1) * G2/2.               $
                + (1.-Q) * ( A1*A2 * (G1-G2+(G3+G1)/4.)/2.   $
                           + B1*B2 * (G1+(G3+G1)/4.)/2.      $
                           + C1*C2 * (G2/2.-(G3+G1)/4.)      $ 
                           )                                 $
                )

 return, V12
end
