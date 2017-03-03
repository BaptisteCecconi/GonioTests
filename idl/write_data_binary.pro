; ===========================================
  PRO write_data_binary,file,data,null=null
; ===========================================
; file : output file name
; data : data structure (built with data_XXXX__define)
; /null : if set, writes a 0b file length

; opening binary data file :
; (swapping byte order if CPU is big endian)

openw,lun,file,/get_lun,/swap_if_big_endian

if not keyword_set(null) then writeu,lun,data

; closing

close,lun
free_lun,lun

return
end





