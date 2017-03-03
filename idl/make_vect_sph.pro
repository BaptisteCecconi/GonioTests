  FUNCTION make_vect_sph,h_in,th_in,ph_in,deg=deg
  
  h = double([h_in])
  th = double([th_in])    
  ph = double([ph_in])

  sh = size(h)
  sth = size(th)
  sph = size(ph)
  
  nh = sh(1)
  nth = sth(1)
  nph = sph(1)

  err = bytarr(4)
  err(0) = sh(0) ne 1
  err(1) = sth(0) ne 1
  err(2) = sph(0) ne 1
  err(3) = ((nh-nth)*(nh-1)*(nth-1) + (nh-nph)*(nh-1)*(nph-1) + $
    (nth-nph)*(nth-1)*(nph-1)) ne 0

  if (where(err))(0) ne -1 then begin
    print,"MAKE_VECT_SPH: error, wrong argument size..."
    if err(0) then print," > length must a scalar (1) or a vector (n)"
    if err(1) then print," > colat. must a scalar (1) or a vector (n)"
    if err(2) then print," > azimuth must a scalar (1) or a vector (n)"
    if err(3) then print," > in case of vectors, argument must have same dimensions (n or 1)"
    return,-1
  endif

  if keyword_set(deg) then begin
    th=th/!radeg & ph=ph/!radeg 
  endif
  nv = max([nh,nth,nph])
  vect = dblarr(3,nv)
  h = rebin(reform(h,1,nh),3,nv)
  th = rebin(reform(th,1,nth),1,nv)
  ph = rebin(reform(ph,1,nph),1,nv)

  vect = h*[cos(ph)*sin(th),sin(th)*sin(ph),cos(th)]
  return,vect
  end


