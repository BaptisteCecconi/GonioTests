PRO read_antenna_set,ant,file,quiet = quiet,rad=rad,help=help,path=path

if keyword_set(help) then begin
  print, '--------------------------------------------------------------------'
  print, 'READ_ANTENNA_SET, antenna, file_name[, /quiet][, /rad][, PATH=path]'
  print, '--------------------------------------------------------------------'
  print, '  antenna   : output structure containing antenna parameters'
  print, '  file_name : string of the antenna set file name (no path, no ext.)'
  print, '  /quiet    : no verbose output'
  print, '  /rad      : convert angles into radian (default is degree)'
  print, '  PATH      : specify Path for reading ant file.'
  print, '--------------------------------------------------------------------'
  return
endif

if keyword_set(quiet) then notquiet = 0b else notquiet = 1b 

ant = {antenna_set}
if not keyword_set(path) then path = getenv('ROOT_RPWS')+'/ant'
filelong = path+'/'+file+".ant"

if notquiet then begin 
  message,"Reading Antenna Parameters from : "+filelong,/info
endif
  
openr,lun,filelong,/get_lun,/swap_if_big_endian
  readu,lun,ant
close,lun
free_lun,lun

if keyword_set(rad) then begin
  if notquiet then message," Converting angles into radian ...",/info
  ant.Xp.al  = ant.Xp.al/!radeg
  ant.Xm.al  = ant.Xm.al/!radeg
  ant.Z.al   = ant.Z.al/!radeg
  ant.Dip.al = ant.Dip.al/!radeg
  ant.Xp.be  = ant.Xp.be/!radeg
  ant.Xm.be  = ant.Xm.be/!radeg
  ant.Z.be   = ant.Z.be/!radeg
  ant.Dip.be = ant.Dip.be/!radeg
  
endif

if notquiet then message,"Antenna Parameters loaded.",/info

return
end
