; ====================================================================
PRO dfb_main,df_result,data_p,data_m,ephem,antenna_file,sn, $
             VERBOSE=VERBOSE,TEST=TEST
; ====================================================================
; all angles are in radian
; ====================================================================
; v1.0 [BC, 2008-jan-28] : first version
; v1.1 [BC, 2009-jun-02] : best solution selection modified 
; --------------------------------------------------------------------

if antenna_file eq 'calDec04' then begin 
  read_antenna_set,antenna_HF,'calDec04_H12',/rad
  read_antenna_set,antenna_LF,'calDec04_ABC',/rad
endif else begin
  read_antenna_set,antenna_HF,antenna_file,/rad
  antenna_LF = antenna_HF
endelse

f_antenna_lim = 320.

ndim = n_elements(df_result)

; --------------------------------------------------------------------
;                    antenna parameters initialization
; --------------------------------------------------------------------

al_tmp  = [antenna_LF.Z.al,antenna_HF.Z.al]
be_tmp  = [antenna_LF.Z.be,antenna_HF.Z.be]
h__tmp  = [antenna_LF.Z.h, antenna_HF.Z.h ]
alZ = al_tmp( data_p.f ge f_antenna_lim )
beZ = be_tmp( data_p.f ge f_antenna_lim )
hZ  = h__tmp( data_p.f ge f_antenna_lim )

al_tmp = [antenna_LF.Xp.al,antenna_HF.Xp.al]
be_tmp = [antenna_LF.Xp.be,antenna_HF.Xp.be]
h__tmp = [antenna_LF.Xp.h, antenna_HF.Xp.h ]
alXp = al_tmp( data_p.f ge f_antenna_lim )
beXp = be_tmp( data_p.f ge f_antenna_lim )
hXp  = h__tmp( data_p.f ge f_antenna_lim )

al_tmp = [antenna_LF.Xm.al,antenna_HF.Xm.al]
be_tmp = [antenna_LF.Xm.be,antenna_HF.Xm.be]
h__tmp = [antenna_LF.Xm.h, antenna_HF.Xm.h ]
alXm = al_tmp( data_p.f ge f_antenna_lim )
beXm = be_tmp( data_p.f ge f_antenna_lim )
hXm  = h__tmp( data_p.f ge f_antenna_lim )

; --------------------------------------------------------------------
; reference frame rotation : S/C to Antenna
; --------------------------------------------------------------------

vhZ  = make_vect_sph(1.,alZ,beZ)
vhXp = make_vect_sph(1.,alXp,beXp)
vhXm = make_vect_sph(1.,alXm,beXm)

e3 = [0.,0.,1.] 

vr1 = crossp1(vhZ,e3,/norm)
qr1 = q_make(alZ,vr1)

qhZ = q_vmake(vhZ)
qhXp = q_vmake(vhXp)
qhXm = q_vmake(vhXm)

qhZ1 = q_rot(qr1,qhZ)
qhXp1 = q_rot(qr1,qhXp)
qhXm1 = q_rot(qr1,qhXm)

beXp1 = reform(atan(qhXp1(2,*),qhXp1(1,*)))     ; in radian
beXm1 = reform(atan(qhXm1(2,*),qhXm1(1,*)))     ; in radian

beta = ( !pi - (beXm1 - beXp1) )/2. ; in radian

qr2 = q_make(-beXp1+beta,e3)

qhZ2 = q_rot(qr2,qhZ1)
qhXp2 = q_rot(qr2,qhXp1)
qhXm2 = q_rot(qr2,qhXm1)

alXp2 = reform(acos(qhXp2(3,*)))                ; in radian
alXm2 = reform(acos(qhXm2(3,*)))                ; in radian

beXp2 = reform(atan(qhXp2(2,*),qhXp2(1,*)))     ; in radian
beXm2 = reform(atan(qhXm2(2,*),qhXm2(1,*)))     ; in radian

;if keyword_set(verbose) then begin
;  print,'Antenna Frame : '
;  print,'Alpha (X+) = ',alXp2(0)*!radeg
;  print,'Beta (X+)  = ',beXp2(0)*!radeg
;  print,'Alpha (X-) = ',alXm2(0)*!radeg
;  print,'Beta (X-)  = ',beXm2(0)*!radeg
;endif

; Computing ephemeris in antenna frame
; ------------------------------------

veph = make_vect_sph(1.,ephem.th,ephem.ph)
qeph = q_vmake(veph)

qeph1 = q_rot(qr1,qeph)
qeph2 = q_rot(qr2,qeph1)

ph_eph2 = (atan(qeph2(2,*),qeph2(1,*)) + 2*!pi) mod (2.*!pi)
th_eph2 = acos(qeph2(3,*))

;if keyword_set(test) then begin
  ph_eph2i = ph_eph2
  th_eph2i = th_eph2
;endif else begin
;  ph_eph2i = interpol(ph_eph2,ephem.time,ti_t97(data.ti))
;  th_eph2i = interpol(th_eph2,ephem.time,ti_t97(data.ti))
;endelse

; --------------------------------------------------------------------------
; Determination of Phi in the Antenna frame
; --------------------------------------------------------------------------

ph2 = (atan((hXp*sin(alXp2)*data_m.CrossI - $
             hXm*sin(alXm2)*data_p.CrossI)*tan(beta), $
            (hXp*sin(alXp2)*data_m.CrossI +  $
             hXm*sin(alXm2)*data_p.CrossI) ) + 2.*!pi) mod (2.*!pi)

ph2_a = reform(ph2,ndim) 
ph2_b = reform(ph2,ndim) + !pi

; ---------------------------------------------------------------------------
; Determination of Theta in the antenna frame 
; ---------------------------------------------------------------------------


vvz = (data_p.AutoZ+data_m.AutoZ)/2.

Th2_a = atan(vvz*hXp*hXm*sin(alXp2)*sin(alXm2)*sin(2.*beta), $
           ((hXp*vvz*cos(alXp2)-hZ*data_p.crossR)*hXm*sin(alXm2)* $
            sin(ph2_a+beta)+(hXm*vvz*cos(alXm2)-hZ*data_m.crossR)* $
            hXp*sin(alXp2)*sin(ph2_a-beta)))

Th2_a =( Th2_a + !pi ) mod !pi

Th2_b = atan(vvz*hXp*hXm*sin(alXp2)*sin(alXm2)*sin(2.*beta), $
           ((hXp*vvz*cos(alXp2)-hZ*data_p.crossR)*hXm*sin(alXm2)* $
            sin(ph2_b+beta)+(hXm*vvz*cos(alXm2)-hZ*data_m.crossR)* $
            hXp*sin(alXp2)*sin(ph2_b-beta)))

Th2_b =( Th2_b + !pi ) mod !pi

; ---------------------------------------------------------------------------
; Selecting source with input guessed ephemeris
; ---------------------------------------------------------------------------

dtha = angular_distance(make_vect_sph(1,th2_a,ph2_a),qeph2(1:3,*))
dthb = angular_distance(make_vect_sph(1,th2_b,ph2_b),qeph2(1:3,*))

dth = reform([dtha,dthb],ndim,2)
dthmin = min(dth,idthmin,dimension=2)

th2 = ([th2_a,th2_b])(idthmin)
ph2 = ([ph2_a,ph2_b])(idthmin)


; ---------------------------------------------------------------------------
; Preparing Useful Quantities...
; ---------------------------------------------------------------------------

sxp = hXp*(cos(alXp2)*sin(th2) - sin(alXp2)*cos(th2)*cos(ph2-beta))
sxm = hXm*(cos(alXm2)*sin(th2) + sin(alXm2)*cos(th2)*cos(ph2+beta))

syp = -hXp*sin(alXp2)*sin(ph2-beta)
sym = +hXm*sin(alXm2)*sin(ph2+beta)

sz  = hZ*sin(th2)

; ---------------------------------------------------------------------------
; Determination of Flux (S)
; ---------------------------------------------------------------------------

Sp = (data_p.autoX*sz^2. - 2.*data_p.crossR*sxp*sz + $
      data_p.autoZ*(sxp^2.+syp^2.)) / (syp^2.*sz^2.)
Sm = (data_m.autoX*sz^2. - 2.*data_m.crossR*sxm*sz + $
      data_m.autoZ*(sxm^2.+sym^2.)) / (sym^2.*sz^2.)

; ---------------------------------------------------------------------------
; Determination of Stokes Parameter V (in antenna frame)
; ---------------------------------------------------------------------------
;
;Vp = -2.*data_p.crossI/(Sp*sz*syp)
;Vm = -2.*data_m.crossI/(Sm*sz*sym)
;
; ===== NB =====
; le signe moins vient du fait qu'on calcul th,ph en tant que direction 
; de la source et non direction de k
; ==============
;
; => V is now computed in the Wave frame
;
; ---------------------------------------------------------------------------
; Theta and Phi in S/C frame
; ---------------------------------------------------------------------------

vsrc2 = make_vect_sph(1.,th2,ph2)
qsrc2 = Q_vmake(vsrc2)

qsrc1 = Q_rot(Q_conj(qr2),qsrc2)
qsrc = Q_rot(Q_conj(qr1),qsrc1)

th = acos(qsrc(3,*))
ph = (atan(qsrc(2,*),qsrc(1,*)) + 2*!pi) mod (2.*!pi)

ith = where(abs(th - !pi/2.) - !pi/2. ge -0.01)
if ith(0) ne -1 then ph(ith) = 0.

; ---------------------------------------------------------------------------
; Determination of Stokes Parameters Q,U (linear polarization) in S/C frame
; ---------------------------------------------------------------------------

sxz = hZ*(cos(alZ)*sin(th)-sin(alZ)*cos(th)*cos(ph-beZ))
sxp = hXp*(cos(alXp)*sin(th)-sin(alXp)*cos(th)*cos(ph-beXp))
sxm = hXm*(cos(alXm)*sin(th)-sin(alXm)*cos(th)*cos(ph-beXm))

syz = -hZ*sin(alZ)*sin(ph-beZ)
syp = -hXp*sin(alXp)*sin(ph-beXp)
sym = -hXm*sin(alXm)*sin(ph-beXm)

Up = (data_p.autoZ*(sxp^2.-syp^2.) - data_p.autoX*(sxz^2.-syz^2.))/ $
     (Sp*(syz*sxp-syp*sxz)*(sxz*sxp+syz*syp)) - (syz*sxp+syp*sxz)/ $
     (sxz*sxp+syz*syp)
Um = (data_m.autoZ*(sxm^2.-sym^2.) - data_m.autoX*(sxz^2.-syz^2.))/ $
     (Sm*(syz*sxm-sym*sxz)*(sxz*sxm+syz*sym)) - (syz*sxm+sym*sxz)/ $
     (sxz*sxm+syz*sym)


Qp =  2.*(data_p.autoX*sxz*syz - data_p.autoZ*sxp*syp)/ $
      (Sp*(syz*syp+sxz*sxp)*(syz*sxp-syp*sxz)) + $
      (syz*syp-sxz*sxp)/(syz*syp+sxz*sxp)
Qm =  2.*(data_m.autoX*sxz*syz - data_m.autoZ*sxm*sym)/ $
      (Sm*(syz*sym+sxz*sxm)*(syz*sxm-sym*sxz)) + $
      (syz*sym-sxz*sxm)/(syz*sym+sxz*sxm)

Vp = -2.*data_p.crossI/(Sp*(syz*sxp-syp*sxz))
Vm = -2.*data_m.crossI/(Sm*(syz*sxm-sym*sxz))

df_result.ydh = data_p.ydh
df_result.num = transpose(reform([data_p.num,data_m.num],ndim,2))
df_result.s  = transpose(reform([sp,sm],ndim,2))
df_result.q  = transpose(reform([qp,qm],ndim,2))
df_result.u  = transpose(reform([up,um],ndim,2))
df_result.v  = transpose(reform([vp,vm],ndim,2))
df_result.th = reform(th,ndim)
df_result.ph = reform(ph,ndim)
df_result.zr = reform((data_p.autoZ-data_m.autoZ)/(data_p.autoZ+data_m.autoZ),ndim)
df_result.sn = sn

return
end
