; ================================================================
  PRO read_data_binary,file,data,valid,ndata=ndata,level=level
; ================================================================
; PARAMETERS
; - file [IN]   : fullpath data file name
; - data [OUT]  : output data. Either a array of data structure
;                 or a long integer (see ndata keyword)   
; - valid [OUT] : 1b if file is ok
;                 0b if file is 0 byte long
; KEYWORDS
; /ndata        : return the number of records in data parameter 
;                 if set.
; level=level   : specify the type of binary data to read. No 
;                 default value. This keyword MUST be set.
; ================================================================

if not keyword_set(level) then begin
    message,'>>>   ERROR : you MUST select a data level !!!'
    valid=0b
    return
endif else begin
  case level of
    'n1'    : data = {data_n1}
    'n2'    : data = {data_n2}
    'n3a'   : data = {data_n3a}
    'n3b'   : data = {data_n3b}
    'n3c'   : data = {data_n3c}
    'n3d'   : data = {data_n3d}
    'n3e'   : data = {data_n3e}
    'n3f'   : data = {data_n3f}
    'n3g'   : data = {data_n3g}
    'n3'    : data = {data_n3}
    'ephem' : data = {data_ephem}  ; obsolete...
    'qeph'  : data = {data_qeph}
    'veph'  : data = {data_veph}
    'bgold' : data = {data_bg_old}
    'bg320' : data = {data_bg_320}
    'bg640' : data = {data_bg_640}
    'bg800' : data = {data_bg_800}
    'bg1600': data = {data_bg_1600}
    'bg'    : data = {data_bg_1600}
    'bg3200': data = {data_bg_3200}
    'sed'   : data = {data_sed}
    'loc'	: data = {data_loc}
    'index' : data = {data_index}
    else    : begin
      message,'>>>   ERROR : This type of data "'+level+'" does not exist !!! <<<'
      valid = 0b
      return
    end
  endcase
endelse

; size of 1 record of the structure :
nbytesdata = n_tags(data, /data_length) 
; works only on IDL 5.6 and up
; doesn't work in GDL 0.9.6... !! 

; >>>> QUICK AND DIRTY FIX FOR GDL >>>>
case level of 
	'n2': nbytesdata = 45l
	'n3b': nbytesdata = 72
	'ephem': nbytesdata = 44l
endcase
; <<<< QUICK AND DIRTY FIX FOR GDL <<<<

; opening binary data file :
; (swapping byte order if CPU is big endian)

openr,lun,file,/get_lun,/swap_if_big_endian

; number of records in the file :
s = fstat(lun)
nd = s.size/nbytesdata
; checking if file length is consistent with single record length
if long(nd)*long(nbytesdata) ne s.size then begin
  message,'>>>   WARNING : you may be using a wrong binary data format !!!',/info
  message,'>>>   file='+file,/info  
endif

if keyword_set(ndata) then begin
  data = nd
endif else begin
  if nd gt 0 then begin 
    valid = 1b
;   final data structure definition : 
    data = replicate(data,nd)
;   reading data : 
    readu,lun,data
  endif else valid=0l
endelse

; closing
close,lun
free_lun,lun

return
end
