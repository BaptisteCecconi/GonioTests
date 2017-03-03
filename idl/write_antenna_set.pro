PRO write_antenna_set,file=file,ant_set=ant_set,help=help,rad=rad,path=path
;===========================================================================
; help : 
if keyword_set(help) then begin
    print,'% PRO write_antenna_set[,file=''file''][,ant_set=ant_set][,/rad][,PATH=path][,/help]'
    print,'% Description : '
    print,'%    This procedure writes *.ant  binary files containing antenna '
    print,'%    sets parameters for the HFR-RPWS radio reciever.'
    print,'% Keywords : '
    print,'%    file = file name, must be a string with neither path nor extension'
    print,'%    ant_set = antenna set you want to save, '
    print,'%              must be defined as an {antenna_set} structure'
    print,'%              angles must be in degrees except if /RAD is set'
    print,'%    /rad  = set keyword to set input in radian'
    print,'%    Path  = set specific path for output default is $ROOT_RPWS/ant/'
    print,'%    /help = displays this help'
    print,'% NOTE : if the keywords are not set, you are asked to enter '
    print,'% the parameters manually.'
    print,'% Convention : angles in degrees unless specified in input.'
    print,'%              angles are in degrees in the output file.'
    
endif else begin 
    antenna = {antenna_set}
    
    if keyword_set(ant_set) then antenna = ant_set else begin
        print,"Please prepare your antenna parameters ... "
        if keyword_set(rad) then print,'(in rad.)' else print,'(in deg.)'
        print,""
        
        tmp=0.
        read,"X+ : h  ? ",tmp
        antenna.Xp.h = tmp
        read,"X+ : al ? ",tmp
        antenna.Xp.al = tmp
        read,"X+ : be ? ",tmp
        antenna.Xp.be = tmp
        read,"X- : h  ? ",tmp
        antenna.Xm.h = tmp
        read,"X- : al ? ",tmp
        antenna.Xm.al = tmp
        read,"X- : be ? ",tmp
        antenna.Xm.be = tmp
        read,"Z  : h  ? ",tmp
        antenna.Z.h = tmp
        read,"Z  : al ? ",tmp
        antenna.Z.al = tmp
        read,"Z  : be ? ",tmp
        antenna.Z.be = tmp
        read,"Dip: h  ? ",tmp
        antenna.Dip.h = tmp
        read,"Dip: al ? ",tmp
        antenna.Dip.al = tmp
        read,"Dip: be ? ",tmp
        antenna.Dip.be = tmp
        
    endelse
    
    if keyword_set(rad) then begin
        antenna.Xp.al = antenna.Xp.al*!radeg
        antenna.Xp.be = antenna.Xp.be*!radeg
        antenna.Xm.al = antenna.Xm.al*!radeg
        antenna.Xm.be = antenna.Xm.be*!radeg
        antenna.Z.al = antenna.Z.al*!radeg
        antenna.Z.be = antenna.Z.be*!radeg
        antenna.Dip.al = antenna.Dip.al*!radeg
        antenna.Dip.be = antenna.Dip.be*!radeg
    endif

    message,"You entered the followind parameters : ",/info
    message,string(format='("X+ antenna : h=",f4.2," al=",f5.1," be=",f5.1)',antenna.Xp.h,antenna.Xp.al,antenna.Xp.be),/info
    message,string(format='("X- antenna : h=",f4.2," al=",f5.1," be=",f5.1)',antenna.Xm.h,antenna.Xm.al,antenna.Xm.be),/info
    message,string(format='("Z  antenna : h=",f4.2," al=",f5.1," be=",f5.1)',antenna.Z.h,antenna.Z.al,antenna.Z.be),/info
    message,string(format='("Dipole     : h=",f4.2," al=",f5.1," be=",f5.1)',antenna.Dip.h,antenna.Dip.al,antenna.Dip.be),/info
    
    if keyword_Set(file) then filant = file else begin
        filant=""
        read,"File name for this antenna set ? [no path, no extension] ",filant
    endelse
    
    if not keyword_set(path) then path = getenv('ROOT_RPWS')+'/ant'
    file = path+'/'+filant+'.ant'
    
    openw,lun,file,/get_lun,/swap_if_big_endian
    writeu,lun,antenna
    close,lun
    free_lun,lun
    message,"Antenna set parameters writen in "+file,/info

endelse    
end
