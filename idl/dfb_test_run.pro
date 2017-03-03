; ====================================================================
  PRO dfb_test_run,     q0=q0, u0=u0, v0=v0, sdb0=sdb0,           $
                        nth=nth, nph=nph,                         $
                        dhz=dhz,   dalz=dalz,   dbez=dbez,        $
                        dhxp=dhxp, dalxp=dalxp, dbexp=dbexp,      $
                        dhxm=dhxm, dalxm=dalxm, dbexm=dbexm,      $
                        agcNbit = agcNbit,                        $
                        snr=snr, Bsnr=Bsnr, Tsnr=Tsnr, Nsnr=Nsnr, $
                        dazz=dazz, dbg=dbg, bglev=bglev,          $
                        gam0=gam0, src_model=src_model,           $
                        cross_phase=cross_phase,				  $
                        dipole_null_level=dipole_null_level,      $
                        antenna_file=antenna_file,                $
                        file_ext=file_ext, output_path=output_path
; ====================================================================
; VERSION 3.3 [03-mar-2017]
; --------------------------------------------------------------------
; --------------------------------------------------------------------
; PURPOSE : Tests the dfb procedure set with a simulated grid of data
; with various alterations.
; --------------------------------------------------------------------
; KEYWORDS :
; [definition of the grid]
;  - Q0 : 1D array containing the Q values for the simulated grid
;  - U0 : 1D array containing the U values for the simulated grid
;  - V0 : 1D array containing the V values for the simulated grid
;  # NB: - default values if not set :
;  #       U,Q = [-1.,-.75,-5,-.25,0.,.25,.5,.75,1.]
;  #       V   = [-1.,-.75,-5,-.25,   .25,.5,.75,1.]
;  #       NB: V = 0.0 is excluded for dfb inversion
;  #     - automatic deletion of (U^2+Q^2+V^2)>1 points
;  - sbd0 : 1D array containig the S [dB] values for the simul. grid
;           default = [-145,-150]
;  - gam0 : gamma source angular extension (in deg) (float not array)
;           default = 0.
;  - nth : number of colatitude steps (th) 
;          [default = 36 <-> 5 deg steps]
;  - nph : number of azimuth steps (ph)
;          [default = 72 <-> 5 deg steps]
; [errors on antenna lengths]
;  - dhz : error on hZ parameter (in percent)
;  - dhxp : error on hX+ parameter (in percent)
;  - dhxm : error on hX- parameter (in percent)
; [errors on antenna colatitudes]
;  - dalz : error on alZ parameter (in deg)
;  - dalxp : error on alX+ parameter (in deg)
;  - dalxm : error on alX- parameter (in deg)
; [error on antenna azimuths]
;  - dbez : error on beZ parameter (in deg)
;  - dbexp : error on beX+ parameter (in deg)
;  - dbexm : error on beX- parameter (in deg)
; [AGC noise : digitization error introduced by the AGC]
;  - agcNbit : set to 8 or 12
;    introduces - 8bit AGC digitization noise as on Cassini 
;               - 12bit AGC digitization noise as on Stereo 
; [SNR]
;  - /snr : activate SNR simulation
;  - Bsnr : Band width of measurements [default = 25.e3 Hz]
;  - Tsnr : Real integration time [default = 16.e-3 sec]
;  - Nsnr : size of the array of gaussian noise which is added to
;    the replicated simul. grid [default = 100]
; [DAzz : changing in S during the antenna switching]
;  - dazz : value in percent
; [Dbg : error on background determination (constant offset on
;  autocorrelations with threshold at bg level)]
;  - dbg : value on the offset
;  - bglev : value of bg level [default = -160 dB]
; [error on the cross-correlation phase calibration]
;  - cross_phase in deg. [default=0.]
; [level of the gain null in the antenna direction]
;  - dipole_null_level [default=0. (perfect dipole)]
;    in units of antenna gain: 1. = perfect dipole max gain
;    antenna gain is re-normalized to 8!pi/3 sr.
; [antenna_file]
;  - path of the antenna file (containing antenna parameters) 
; [output file]
;  - file_ext : (string) suffix for ouput file :
;    dfb_test_[file_ext].bin 
; [output_path]
;  - specify the output_path for files. 
;    [default=$DATA_RPWS/test/]
; --------------------------------------------------------------------


; compiling quaternion routines
; --------------------------------------------------------------------

  quaternion

; ====================================================================
; ====================   Beginning initialisation  ===================
; ====================================================================

; ====================================================================
;   building data : 
; ====================================================================


;     Polarization parameters :
;
;---------------------------------------------------------------------

 
if not(keyword_set(q0)) then q00 = [-1.,-.75,-.5,-.25,0.,.25,.5,.75,1.] else q00 = [q0]
nq = n_elements(q00)

if not(keyword_set(u0)) then u00 = [-1.,-.75,-.5,-.25,0.,.25,.5,.75,1.] else u00 = [u0]
nu = n_elements(u00)

if not(keyword_set(v0)) then v00 = [-1.,-.75,-.5,-.25,.25,.5,.75,1.] else v00 = [v0]
nv = n_elements(v00)

q0 = rebin(reform(q00,nq,1,1),nq,nu,nv)
u0 = rebin(reform(u00,1,nu,1),nq,nu,nv)
v0 = rebin(reform(v00,1,1,nv),nq,nu,nv)

w_quv = where(sqrt(q0^2.+u0^2.+v0^2.) le 1.)
npol = n_elements(w_quv)

q0 = q0(w_quv)
u0 = u0(w_quv)
v0 = v0(w_quv)

;     Angular parameters :
;
; ---------------------------------------------------------------------

; colatitude : 5 deg to 175 deg (5 deg steps) + 0 deg + 180 deg
;             (5 deg = 36 total steps)
; azimuths   : 0 deg to 355 deg (5 deg steps) + 0 deg +   0 deg
;             (5 deg = 72 total steps)


if not(keyword_set(nth))   then nth =    36l
if not(keyword_set(nph))   then nph =    72l

th00 = ((findgen(nth-1.)+1.)/(nth))*180.
ph00 = findgen(nph)/nph*360.

th0 = [0.,(rebin(reform(th00,nth-1,1),nth-1,nph))(*),180.] /!radeg
ph0 = [0.,(rebin(reform(ph00,1,nph),nth-1,nph))(*),0.]     /!radeg

;     Source Extension parameter :
;
; ---------------------------------------------------------------------

if not(keyword_set(gam0)) then gam00 = 0. else gam00 = gam0

gam = gam00/!radeg

;     Flux intensity parameter : 
;
; ---------------------------------------------------------------------

if not(keyword_set(sdb0)) then sdb00 = [-145,-150] $
else sdb00 = [sdb0]

ns = n_elements(sdb00)
s0 = 10.^(sdb00/10.)

;     Phase calibration Error on cross-correlation :
;
; --------------------------------------------------------------------

if not(keyword_set(cross_phase)) then cross_phase00 = 0. $
else cross_phase00=cross_phase

cross_phase0 = cross_phase00/!radeg

; ====================================================================
;   Merging :
; ====================================================================

ndim = (2+(nth-1)*nph)*npol*ns

s = (rebin(reform(s0,1,1,ns),npol,2+(nth-1)*nph,ns))(*)
q = (rebin(reform(q0,npol,1,1),npol,2+(nth-1)*nph,ns))(*)
u = (rebin(reform(u0,npol,1,1),npol,2+(nth-1)*nph,ns))(*)
v = (rebin(reform(v0,npol,1,1),npol,2+(nth-1)*nph,ns))(*)
th = (rebin(reform(th0,1,2+(nth-1)*nph,1),npol,2+(nth-1)*nph,ns))(*)
ph = (rebin(reform(ph0,1,2+(nth-1)*nph,1),npol,2+(nth-1)*nph,ns))(*)


print,'Polar : ',npol,' -- Angles : ',2+(nth-1)*nph,' -- Flux : ',ns
print,' -> Total : ',ndim

; --------------------------------------------------------------------
; SNR preparation :
; --------------------------------------------------------------------
if keyword_set(snr) then begin
  if not keyword_set(Nsnr) then Nsnr = 100l
  Nsnr = long(Nsnr)
  print,' SNR simulation : x',strtrim(string(Nsnr),2)
  print,' -> Total :',ndim*Nsnr
  s = reform(rebin(reform(s,1,ndim),Nsnr,ndim),ndim*Nsnr)
  q = reform(rebin(reform(q,1,ndim),Nsnr,ndim),ndim*Nsnr)
  u = reform(rebin(reform(u,1,ndim),Nsnr,ndim),ndim*Nsnr)
  v = reform(rebin(reform(v,1,ndim),Nsnr,ndim),ndim*Nsnr)
  th = reform(rebin(reform(th,1,ndim),Nsnr,ndim),ndim*Nsnr)
  ph = reform(rebin(reform(ph,1,ndim),Nsnr,ndim),ndim*Nsnr)
  ndim = ndim*Nsnr
endif

; ====================================================================
; Filling data structures :  
; ====================================================================

; --------------------------------------------------------------------
; loading input antenna parameters 

;if keyword_set(antenna_file) then begin
  read_antenna_set,ant,antenna_file,/rad,path=output_path+'temp/'
;endif else begin
;  antenna_file='dfb_test'
;  read_antenna_set,ant,'dfb_test',/rad
;endelse

;  antenna parameter error : 
;
; NB: this simulates a wrong antenna calibration
; ---------------------------------------------------------------------

if keyword_set(dhz) then ant.Z.h = ant.Z.h * (1+dhz/100.)   
if keyword_set(dalz) then ant.Z.al = ant.Z.al + dalz/!radeg
if keyword_set(dbez) then ant.Z.be = ant.Z.be + dbez/!radeg

if keyword_set(dhxp) then ant.Xp.h = ant.Xp.h * (1+dhxp/100.)
if keyword_set(dalxp) then ant.Xp.al = ant.Xp.al + dalxp/!radeg
if keyword_set(dbexp) then ant.Xp.be = ant.Xp.be + dbexp/!radeg

if keyword_set(dhxm) then ant.Xm.h = ant.Xm.h * (1+dhxm/100.)
if keyword_set(dalxm) then ant.Xm.al = ant.Xm.al + dalxm/!radeg
if keyword_set(dbexm) then ant.Xm.be = ant.Xm.be + dbexm/!radeg

; --------------------------------------------------------------------

datap = replicate({data_n2},ndim)
datam = replicate({data_n2},ndim)
df = replicate({data_n3b},ndim)
angle = replicate({data_ephem},ndim)

time = dindgen(ndim)/86400.d0+1.d0

df.num = transpose([[lindgen(ndim)*2l],[lindgen(ndim)*2l+1l]])
df.ydh = 0l
df.s = rebin(reform(s,1,ndim),2,ndim)
df.q = rebin(reform(q,1,ndim),2,ndim)
df.u = rebin(reform(u,1,ndim),2,ndim)
df.v = rebin(reform(v,1,ndim),2,ndim)
df.th = th  
df.ph = ph


datap.num = df.num(0)
datap.ydh = df.ydh
datam.num = df.num(1)
datam.ydh = df.ydh
datap.ant = 11b
datam.ant = 12b

datap.autoX = v1v2_xtnd(s,ant.Xp.h,ant.Xp.h,ant.Xp.al,ant.Xp.al,ant.Xp.be, $
                 ant.Xp.be,Q,U,th,ph,gam,src_model=src_model)
datap.autoZ = v1v2_xtnd(s,ant.Z.h,ant.Z.h,ant.Z.al,ant.Z.al,ant.Z.be, $
                 ant.Z.be,Q,U,th,ph,gam,src_model=src_model)
datam.autoX = v1v2_xtnd(s,ant.Xm.h,ant.Xm.h,ant.Xm.al,ant.Xm.al,ant.Xm.be, $
                 ant.Xm.be,Q,U,th,ph,gam,src_model=src_model)
datam.autoZ = datap.autoZ
datap.crossR = v1v2_xtnd(s,ant.Xp.h,ant.Z.h,ant.Xp.al,ant.Z.al,ant.Xp.be, $
                 ant.Z.be,Q,U,th,ph,gam,src_model=src_model)
datap.crossI = v1iv2_xtnd(s,ant.Xp.h,ant.Z.h,ant.Xp.al,ant.Z.al,ant.Xp.be, $
                 ant.Z.be,v,th,ph,gam,src_model=src_model)
datam.crossR = v1v2_xtnd(s,ant.Xm.h,ant.Z.h,ant.Xm.al,ant.Z.al,ant.Xm.be, $
                 ant.Z.be,Q,U,th,ph,gam,src_model=src_model)
datam.crossI = v1iv2_xtnd(s,ant.Xm.h,ant.Z.h,ant.Xm.al,ant.Z.al,ant.Xm.be, $
                 ant.Z.be,v,th,ph,gam,src_model=src_model)

angle.time = time
angle.th = th
angle.ph = ph

; =====================================================================
; Adding Errors : 
; =====================================================================


;   non-zero gain in dipole null :
; ---------------------------------------------------------------------

if keyword_set(dipole_null_level) then begin

   dipole_null_normalization_coeff = 3.*dipole_null_level/2. + 1
   dthXp = angular_distance(make_vect_sph(1.,ant.Xp.al,ant.Xp.be),make_vect_sph(1.,th,ph)) 
   dthXm = angular_distance(make_vect_sph(1.,ant.Xm.al,ant.Xm.be),make_vect_sph(1.,th,ph)) 
   dthZ  = angular_distance(make_vect_sph(1.,ant.Z.al,ant.Z.be),make_vect_sph(1.,th,ph))

   coefXp = (dipole_null_level + sin(dthXp)^2)/sin(dthXp)^2 * dipole_null_normalization_coeff
   coefXm = (dipole_null_level + sin(dthXm)^2)/sin(dthXm)^2 * dipole_null_normalization_coeff
   coefZ  = (dipole_null_level + sin(dthZ)^2)/sin(dthZ)^2 * dipole_null_normalization_coeff

   datap.autoX = datap.autoX * coefXp
   datap.autoZ = datap.autoZ * coefZ
   datam.autoX = datam.autoX * coefXm
   datam.autoZ = datam.autoZ * coefZ
   datap.crossR = datap.crossR * sqrt(coefXp*coefZ)
   datap.crossI = datap.crossI * sqrt(coefXp*coefZ)
   datam.crossR = datam.crossR * sqrt(coefXm*coefZ)
   datam.crossI = datam.crossI * sqrt(coefXm*coefZ)

endif

;   phase calibration error on cross-correlation :
; ---------------------------------------------------------------------

if keyword_set(cross_phase) then begin

    ; +XZ antenna pair :
    ; ------------------

    cross_rr = sqrt(datap.crossR^2. + datap.crossI^2.)
    cross_ph = atan(datap.crossI,datap.crossR)
    
    datap.crossR = cross_rr * cos(cross_ph+cross_phase0)
    datap.crossI = cross_rr * sin(cross_ph+cross_phase0)
    
    ; -XZ antenna pair :
    ; ------------------

    cross_rr = sqrt(datam.crossR^2. + datam.crossI^2.)
    cross_ph = atan(datam.crossI,datam.crossR)
    
    datam.crossR = cross_rr * cos(cross_ph+cross_phase0)
    datam.crossI = cross_rr * sin(cross_ph+cross_phase0)
    
    
endif

;   agc digitization error :
; ---------------------------------------------------------------------


if keyword_set(agcNbit) then begin
    if agcNbit ne 8 and agcNbit ne 12 then stop,'dfb_test - FATAL ERROR : AgcNbit must 8 or 12 ...'
    print,'AGC/LogCompression on ',string(format='(I2)',agcNbit),' bit words ...'

    ; +XZ antenna pair :
    ; ------------------

    axx  = datap.autoX
    azz  = datap.autoZ
    Crxz = datap.crossR
    Cixz = datap.crossI

    agc_error,Axx,Azz,Crxz,Cixz,nbits=agcNbit

    datap.autoX = axx 
    datap.autoZ = azz 
    datap.crossR = Crxz
    datap.crossI = Cixz

    ; -XZ antenna pair :
    ; ------------------

    axx  = datam.autoX
    azz  = datam.autoZ
    Crxz = datam.crossR
    Cixz = datam.crossI

    agc_error,Axx,Azz,Crxz,Cixz,nbits=agcNbit

    datam.autoX = axx 
    datam.autoZ = azz 
    datam.crossR = Crxz
    datam.crossI = Cixz

endif


;   SNR error :
; ---------------------------------------------------------------------

if keyword_set(snr) then begin
  if not keyword_set(bglev) then bglev = -160.
  if not keyword_set(Bsnr) then Bsnr = 25.e3 
  if not keyword_set(Tsnr) then Tsnr = 16.e-3

  noise_lev = 10.^(bglev/10.)/sqrt(Bsnr*Tsnr)
  noise = reform(rebin(reform(randomn(seed,4,Nsnr),4,Nsnr,1),4,Nsnr,Ndim/Nsnr),4,Ndim)*noise_lev
  df.sn(0) = 10.*alog10(datap.autoX/noise_lev)
  df.sn(1) = 10.*alog10(datap.autoZ/noise_lev)
  df.sn(2) = 10.*alog10(datam.autoX/noise_lev)
  df.sn(3) = 10.*alog10(datam.autoZ/noise_lev)
  
  datap.autoX = ( datap.autoX  + reform(noise(0,*),Ndim) ) > 0.
  datap.autoZ = ( datap.autoZ  + reform(noise(1,*),Ndim) ) > 0.
  datam.autoX = ( datam.autoX  + reform(noise(2,*),Ndim) ) > 0.
  datam.autoZ = ( datam.autoZ  + reform(noise(3,*),Ndim) ) > 0.
endif

;  Delta Azz :
; ---------------------------------------------------------------------

if keyword_set(dazz) then begin
  datam.autoX  = datam.autoX  * (1+dazz/100.)
  datam.autoZ  = datam.autoZ  * (1+dazz/100.)
  datam.crossR = datam.crossR * (1+dazz/100.)
  datam.crossI = datam.crossI * (1+dazz/100.)
endif

;   Background error :
; ---------------------------------------------------------------------

if keyword_set(dbg) then begin
  if not keyword_set(bglev) then bglev = -155.
  datap.autoX = (datap.autoX - 10^(bglev/10.) + 10^((bglev+dbg)/10.))>10^((bglev+dbg)/10.)
  datap.autoZ = (datap.autoZ - 10^(bglev/10.) + 10^((bglev+dbg)/10.))>10^((bglev+dbg)/10.)
  datam.autoX = (datam.autoX - 10^(bglev/10.) + 10^((bglev+dbg)/10.))>10^((bglev+dbg)/10.)
  datam.autoZ = (datam.autoZ - 10^(bglev/10.) + 10^((bglev+dbg)/10.))>10^((bglev+dbg)/10.)
endif

; ====================================================================
; Internal stuff : 
; ====================================================================
; NaN exceptions (cos > 1, for instance)
; --------------------------------------
nan_exception = bytarr(ndim)


if keyword_set(output_path) then path = output_path else path='$DATA_RPWS/test/'

if strmid(path,0,/rev) ne '/' then path = path+'/'



; writing files

if not keyword_set(file_ext) then begin
  file_ext=''
  print,'Please specify the file name : '
  read,file_ext
endif

data = replicate({data_n2},ndim*2l)
data(df.num(0)) = datap
data(df.num(0)).ant = 11
data(df.num(1)) = datam
data(df.num(1)).ant = 12


write_data_binary,path+'n2/Pdfb_test_'+file_ext+'.00',data
write_data_binary,path+'ephem/dfb_test_'+file_ext+'.ephem',angle
write_data_binary,path+'n3b/N3b_Ixx_test_'+file_ext+'.00',df
write_antenna_set,file='dfb_test_'+file_ext,ant_set=ant,/rad,path=path+'temp'

; ====================================================================
; ====================  End of initialisation  =======================
; ====================================================================


; ====================================================================
; ====================  Call dfb_main   ==============================
; ====================================================================

dfb_main,df,datap,datam,angle,antenna_file,df.sn,path=output_path
write_data_binary,path+'n3b/N3b_Oxx_test_'+file_ext+'.00',df


end
