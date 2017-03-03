;============================================================
PRO AGC_ERROR,Axx,Azz,Crxz,Cixz, $
              nbits = nbits
;------------------------------------------------------------
; computes AGC+LogCompression digitization error as in HFR
;------------------------------------------------------------
; INPUT/OUTPUT :
; Axx    (I/O) : auto(XX) (V^2/Hz)
; Azz    (I/O) : auto(ZZ) (V^2/Hz)
; Crxz   (I/O) : crossR(XZ) (V^2/Hz)
; Cixz   (I/O) : crossI(XZ) (V^2/Hz)
;------------------------------------------------------------
; LOCAL :
; a1,a2,a3 : AGC coefficient (see HFR doc)
; cal      : dBcal (see HFR doc)
; dbmin    : zero equivalent level in dB 
; agcX,   agcZ   : hfr agc ouput value
; autoX,  autoZ  : hfr logCompressed output values for Autos
; CrossR, CrossI : hfr logCompressed output values for Cross
; AX, AZ : cross computed from DSP input voltages (after AGC)
; CR, CI : cross computed from DSP input voltages (after AGC)
;============================================================

if keyword_set(nbits) then bb = nbits else bb=8
if bb ne 8 and bb ne 12 then stop,'Error ! nbits must be 8 or 12 !'

case bb of
    8 : begin
        a1 = 176.35
        a2 = 94.79
        a3 = 0.
        cal = 68.23  
    end                                              
    12 :begin
        a1 = 176.35
        a2 = 1528.21
        a3 = 0.
        cal = 68.335 
    end  
endcase

; set dbmin to a value right up -a1 :
; -----------------------------------

dbmin = -176.                                                     ; dB V^2/Hz [real]

;**********************************************************************
; -> NB on dBmin :
;   agc_db_invert is undefined for signal <= -a1 [dB]
;   -176.35 dB[V^2/Hz] <=> 2.32 10^-18 V^2/Hz
;   -163.10 dB[V^2/Hz] <=> 49.0 10^-18 V^2/Hz 
;                     [<=> 7 nV.Hz^-1/2 (noise receiver level)]
;**********************************************************************

Azz0 = Azz
Axx0 = Axx
Crxz0 = Crxz
Cixz0 = Cixz

; put auto measurements in dB and keep sign of Cross
; --------------------------------------------------

Axx = reform(Axx)>(10.^(dbmin/10.))
Azz = reform(Azz)>(10.^(dbmin/10.))

Axx_db = 10.*alog10(Axx)              ; dB V^2/Hz [real]
Azz_db = 10.*alog10(Azz)              ; dB V^2/Hz [real]
sign_Crxz = reform((Crxz ge 0.)*2. - 1.)
sign_Cixz = reform((Cixz ge 0.)*2. - 1.)

; compute agcX and agcZ
; ---------------------

agcX = agc_db_invert(Axx_db,a1,a2,a3)>0
agcZ = agc_db_invert(Azz_db,a1,a2,a3)>0

; compute AutoX and AutoZ
; -----------------------

autoX = auto_db_invert(Axx_db-agc_db(agcX,a1,a2,a3)+cal,nbits=bb) ; pseudoLog Nbits
autoZ = auto_db_invert(Azz_db-agc_db(agcZ,a1,a2,a3)+cal,nbits=bb) ; pseudoLog Nbits
AX_db = auto_db(autoX,nbits=bb)                                   ; dB V^2/Hz [DSP]  
AZ_db = auto_db(autoZ,nbits=bb)                                   ; dB V^2/Hz [DSP] 
AX = 10.^(AX_db/10.)                                              ; V^2/Hz [DSP]  
AZ = 10.^(AZ_db/10.)                                              ; V^2/Hz [DSP]

; compute CrossR and CrossI
; -------------------------

CR = (Crxz*sign_Crxz)/sqrt(Axx*Azz)*sqrt(AX*AZ)                   ; V^2/Hz [DSP] 
CI = (Cixz*sign_Cixz)/sqrt(Axx*Azz)*sqrt(AX*AZ)                   ; V^2/Hz [DSP] 
CR_db = 10.*alog10(CR>2^(bb-5))                                   ; dB V^2/Hz [DSP]
CI_db = 10.*alog10(CI>2^(bb-5))                                   ; dB V^2/Hz [DSP]
CrossR = auto_db_invert(CR_db,nbits=bb)                           ; pseudoLog Nbits
CrossI = auto_db_invert(CI_db,nbits=bb)                           ; pseudoLog Nbits

; Reconstruct Axx and Azz
; -----------------------

Axx_db = agc_db(agcX,a1,a2,a3)+AX_db-cal                          ; dB V^2/Hz [real]
Azz_db = agc_db(agcZ,a1,a2,a3)+AZ_db-cal                          ; dB V^2/Hz [real]
Axx = 10.^(Axx_db/10.)                                            ; V^2/Hz [real]
Azz = 10.^(Azz_db/10.)                                            ; V^2/Hz [real]

; Reconstruct Crxz and Cixz
; -------------------------

Crxz = sign_Crxz*10.^(auto_db(CrossR,nbits=bb)/10.)/sqrt(AX*AZ)*sqrt(Axx*Azz) ; V^2/Hz [real]
Cixz = sign_Cixz*10.^(auto_db(CrossI,nbits=bb)/10.)/sqrt(AX*AZ)*sqrt(Axx*Azz) ; V^2/Hz [real]

stop

return
end
