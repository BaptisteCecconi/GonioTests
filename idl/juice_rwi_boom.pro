; IDL Version 7.0, Mac OS X (darwin i386 m32)
; Journal File for baptiste@MacBookBC.local
; Working directory: /Users/baptiste/Projets/JUICE-Laplace/Etude Antennes/error_analysis
; Date: Tue Jul 16 09:42:49 2013

data_path = '../data/'

ant_boom={antenna_set}
ant_boom.xp.h=1.
ant_boom.xm.h=1.
ant_boom.z.h=1.
ant_boom.xp.al=36
ant_boom.xp.be=213
ant_boom.xm.al=117
ant_boom.xm.be=249
ant_boom.z.al=69
ant_boom.z.be=345
write_antenna_set,file='juice_boom1',ant_set=ant_boom,path=data_path+'temp'

ant_boom={antenna_set}
ant_boom.xp.h=1.
ant_boom.xm.h=1.
ant_boom.z.h=1.
ant_boom.z.al=36
ant_boom.z.be=213
ant_boom.xp.al=117
ant_boom.xp.be=249
ant_boom.xm.al=69
ant_boom.xm.be=345
write_antenna_set,file='juice_boom2',ant_set=ant_boom,path=data_path+'temp'

ant_boom={antenna_set}
ant_boom.xp.h=1.
ant_boom.xm.h=1.
ant_boom.z.h=1.
ant_boom.xm.al=36
ant_boom.xm.be=213
ant_boom.z.al=117
ant_boom.z.be=249
ant_boom.xp.al=69
ant_boom.xp.be=345
write_antenna_set,file='juice_boom3',ant_set=ant_boom,path=data_path+'temp'

nsnr = 50l
dfb_test_run,antenna_file='juice_boom1',output_path=data_path,file_ext='juice_boom1_snr',/snr,sdb0=[-150],u0=[0.],q0=[0.],v0=[-1.,1.],Nsnr=nsnr,Bsnr=25.e3,Tsnr=1.e-3,nth=36,nph=72
dfb_test_run,antenna_file='juice_boom1',output_path=data_path,file_ext='juice_boom1_snr0',sdb0=[-150],u0=[0.],q0=[0.],v0=[-1,1.],nth=36,nph=72
read_data_binary,data_path+'ephem/dfb_test_juice_boom1_snr0.ephem',eph,level='ephem'
read_data_binary,data_path+'n2/Pdfb_test_juice_boom1_snr0.00',n2,level='n2'
read_data_binary,data_path+'n3b/N3b_Ixx_test_juice_boom1_snr0.00',n3i,level='n3b'
read_data_binary,data_path+'n3b/N3b_Oxx_test_juice_boom1_snr.00',n3o1,level='n3b'

dfb_test_run,antenna_file='juice_boom2',output_path=data_path,file_ext='juice_boom2_snr',/snr,sdb0=[-150],u0=[0.],q0=[0.],v0=[-1.,1.],Nsnr=nsnr,Bsnr=25.e3,Tsnr=1.e-3,nth=36,nph=72
dfb_test_run,antenna_file='juice_boom2',output_path=data_path,file_ext='juice_boom2_snr0',sdb0=[-150],u0=[0.],q0=[0.],v0=[-1,1.],nth=36,nph=72
read_data_binary,data_path+'ephem/dfb_test_juice_boom2_snr0.ephem',eph,level='ephem'
read_data_binary,data_path+'n2/Pdfb_test_juice_boom2_snr0.00',n2,level='n2'
read_data_binary,data_path+'n3b/N3b_Ixx_test_juice_boom2_snr0.00',n3i,level='n3b'
read_data_binary,data_path+'n3b/N3b_Oxx_test_juice_boom2_snr.00',n3o2,level='n3b'

dfb_test_run,antenna_file='juice_boom3',output_path=data_path,file_ext='juice_boom3_snr',/snr,sdb0=[-150],u0=[0.],q0=[0.],v0=[-1.,1.],Nsnr=nsnr,Bsnr=25.e3,Tsnr=1.e-3,nth=36,nph=72
dfb_test_run,antenna_file='juice_boom3',output_path=data_path,file_ext='juice_boom3_snr0',sdb0=[-150],u0=[0.],q0=[0.],v0=[-1,1.],nth=36,nph=72
read_data_binary,data_path+'ephem/dfb_test_juice_boom3_snr0.ephem',eph,level='ephem'
read_data_binary,data_path+'n2/Pdfb_test_juice_boom3_snr0.00',n2,level='n2'
read_data_binary,data_path+'n3b/N3b_Ixx_test_juice_boom3_snr0.00',n3i,level='n3b'
read_data_binary,data_path+'n3b/N3b_Oxx_test_juice_boom3_snr.00',n3o3,level='n3b'


read_antenna_set,ant1,'juice_boom1',/rad,path=data_path+'temp'
read_antenna_set,ant2,'juice_boom2',/rad,path=data_path+'temp'
read_antenna_set,ant3,'juice_boom3',/rad,path=data_path+'temp'
ant = ant1

ndim=n_elements(n3i)
n3o1=reform(n3o1,nsnr,ndim)
n3o2=reform(n3o2,nsnr,ndim)
n3o3=reform(n3o3,nsnr,ndim)

n3o_med = n3i
n3o_med.s = [median(n3o1.s,dim=2)+median(n3o2.s,dim=2)+median(n3o3.s,dim=2)]/3.
n3o_med.q = [median(n3o1.q,dim=2)+median(n3o2.q,dim=2)+median(n3o3.q,dim=2)]/3.
n3o_med.u = [median(n3o1.u,dim=2)+median(n3o2.u,dim=2)+median(n3o3.u,dim=2)]/3.
n3o_med.v = [median(n3o1.v,dim=2)+median(n3o2.v,dim=2)+median(n3o3.v,dim=2)]/3.
n3o_med.th = [median(n3o1.th,dim=1)+median(n3o2.th,dim=1)+median(n3o3.th,dim=1)]/3.
n3o_med.ph = [median(n3o1.ph,dim=1)+median(n3o2.ph,dim=1)+median(n3o3.ph,dim=1)]/3.

n3o_ave = n3i
n3o_ave.s = [total(n3o1.s,2)+total(n3o2.s,2)+total(n3o3.s,2)]/3./nsnr
n3o_ave.q = [total(n3o1.q,2)+total(n3o2.q,2)+total(n3o3.q,2)]/3./nsnr
n3o_ave.u = [total(n3o1.u,2)+total(n3o2.u,2)+total(n3o3.u,2)]/3./nsnr
n3o_ave.v = [total(n3o1.v,2)+total(n3o2.v,2)+total(n3o3.v,2)]/3./nsnr
n3o_ave.th = [total(n3o1.th,1)+total(n3o2.th,1)+total(n3o3.th,1)]/3./nsnr
n3o_ave.ph = [total(n3o1.ph,1)+total(n3o2.ph,1)+total(n3o3.ph,1)]/3./nsnr

n3o_max = n3i
n3o_max.s = max([[[max(n3o1.s,dim=2)]],[[max(n3o2.s,dim=2)]],[[max(n3o3.s,dim=2)]]],dim=3)
n3o_max.q = max([[[max(n3o1.q,dim=2)]],[[max(n3o2.q,dim=2)]],[[max(n3o3.q,dim=2)]]],dim=3)
n3o_max.u = max([[[max(n3o1.u,dim=2)]],[[max(n3o2.u,dim=2)]],[[max(n3o3.u,dim=2)]]],dim=3)
n3o_max.v = max([[[max(n3o1.v,dim=2)]],[[max(n3o2.v,dim=2)]],[[max(n3o3.v,dim=2)]]],dim=3)
n3o_max.th = max([[max(n3o1.th,dim=1)],[max(n3o2.th,dim=1)],[max(n3o3.th,dim=1)]],dim=2)
n3o_max.ph = max([[max(n3o1.ph,dim=1)],[max(n3o2.ph,dim=1)],[max(n3o3.ph,dim=1)]],dim=2)

dth_med01 = angular_distance(make_vect_sph(1.,n3i.th,n3i.ph),median(reform(make_vect_sph(1.,reform(n3o1.th,nsnr*ndim),reform(n3o1.ph,nsnr*ndim)),3,nsnr,ndim),dim=2))
dth_max01 = angular_distance(make_vect_sph(1.,n3i.th,n3i.ph),max(reform(make_vect_sph(1.,reform(n3o1.th,nsnr*ndim),reform(n3o1.ph,nsnr*ndim)),3,nsnr,ndim),dim=2))
dth_ave01 = angular_distance(make_vect_sph(1.,n3i.th,n3i.ph),total(reform(make_vect_sph(1.,reform(n3o1.th,nsnr*ndim),reform(n3o1.ph,nsnr*ndim)),3,nsnr,ndim),2)/nsnr)

dth_med02 = angular_distance(make_vect_sph(1.,n3i.th,n3i.ph),median(reform(make_vect_sph(1.,reform(n3o2.th,nsnr*ndim),reform(n3o2.ph,nsnr*ndim)),3,nsnr,ndim),dim=2))
dth_max02 = angular_distance(make_vect_sph(1.,n3i.th,n3i.ph),max(reform(make_vect_sph(1.,reform(n3o2.th,nsnr*ndim),reform(n3o2.ph,nsnr*ndim)),3,nsnr,ndim),dim=2))
dth_ave02 = angular_distance(make_vect_sph(1.,n3i.th,n3i.ph),total(reform(make_vect_sph(1.,reform(n3o2.th,nsnr*ndim),reform(n3o2.ph,nsnr*ndim)),3,nsnr,ndim),2)/nsnr)

dth_med03 = angular_distance(make_vect_sph(1.,n3i.th,n3i.ph),median(reform(make_vect_sph(1.,reform(n3o3.th,nsnr*ndim),reform(n3o3.ph,nsnr*ndim)),3,nsnr,ndim),dim=2))
dth_max03 = angular_distance(make_vect_sph(1.,n3i.th,n3i.ph),max(reform(make_vect_sph(1.,reform(n3o3.th,nsnr*ndim),reform(n3o3.ph,nsnr*ndim)),3,nsnr,ndim),dim=2))
dth_ave03 = angular_distance(make_vect_sph(1.,n3i.th,n3i.ph),total(reform(make_vect_sph(1.,reform(n3o3.th,nsnr*ndim),reform(n3o3.ph,nsnr*ndim)),3,nsnr,ndim),2)/nsnr)


dth_med1 = angular_distance(make_vect_sph(1.,n3i.th,n3i.ph),make_vect_sph(1.,n3o_med.th,n3o_med.ph))
dth_max1 = angular_distance(make_vect_sph(1.,n3i.th,n3i.ph),make_vect_sph(1.,n3o_max.th,n3o_max.ph))
dth_ave1 = angular_distance(make_vect_sph(1.,n3i.th,n3i.ph),make_vect_sph(1.,n3o_ave.th,n3o_ave.ph)) 

average_ant = (make_vect_sph(1.,ant1.xp.al,ant1.xp.be)+make_vect_sph(1.,ant1.xm.al,ant1.xm.be)+make_vect_sph(1.,ant1.z.al,ant1.z.be))/3.
ant_dist = angular_distance(average_ant,make_vect_sph(1.,n3i.th,n3i.ph))

set_plot,'ps'
device,file=data_path+'plot/juice_rwi_boom_error.ps'
contour,dth_med01*!radeg,n3i.ph*!radeg,n3i.th*!radeg,/irr,lev=[1,2,5,10],xs=1,ys=1,xtit='azimuth [deg]',ytit='colatitude [deg]'; ,c_annotation=['1','2','5','10']
contour,dth_med02*!radeg,n3i.ph*!radeg,n3i.th*!radeg,/irr,lev=[1,2,5,10],/over,xs=1,ys=1; ,c_annotation=['1','2','5','10']
contour,dth_med03*!radeg,n3i.ph*!radeg,n3i.th*!radeg,/irr,lev=[1,2,5,10],/over,xs=1,ys=1; ,c_annotation=['1','2','5','10']
contour,ant_dist*!radeg,n3i.ph*!radeg,n3i.th*!radeg,/irr,lev=findgen(9)*10.+10.,/over,xs=1,ys=1,c_line=fltarr(9)+2; ,c_annotation=string(format='(17I2)',findgen(9)*10.+10.)
contour,180-ant_dist*!radeg,n3i.ph*!radeg,n3i.th*!radeg,/irr,lev=findgen(8)*10.+10.,/over,xs=1,ys=1,c_line=fltarr(8)+2; ,c_annotation=string(format='(17I2)',findgen(8)*10.+10.)
oplot,[ant.xp.be,ant.xm.be,ant.z.be]*!radeg,[ant.xp.al,ant.xm.al,ant.z.al]*!radeg,psym=4,thick=5
oplot,(!pi+[ant.xp.be,ant.xm.be,ant.z.be])*!radeg,(!pi-[ant.xp.al,ant.xm.al,ant.z.al])*!radeg,psym=4,thick=5
oplot,(-!pi+[ant.xp.be,ant.xm.be,ant.z.be])*!radeg,(!pi-[ant.xp.al,ant.xm.al,ant.z.al])*!radeg,psym=4,thick=5
device,/close
set_plot,'x'
set_plot,'ps'
device,file=data_path+'plot/juice_rwi_boom_plot_v.ps'
contour,abs(n3o_med.v(0)),n3i.ph*!radeg,n3i.th*!radeg,/irr,lev=[0.7,0.8,0.9,1,1.1],/fill,xs=1,ys=1,xtit='azimuth [deg]',ytit='colatitude [deg]'
contour,abs(n3o_med.v(0)),n3i.ph*!radeg,n3i.th*!radeg,/over,/irr,lev=[0.7,0.8,0.9,1,1.1],xs=1,ys=1;,c_annotation=['0.7','0.8','0.9','1.0','1.1']
oplot,[ant.xp.be,ant.xm.be,ant.z.be]*!radeg,[ant.xp.al,ant.xm.al,ant.z.al]*!radeg,psym=4,thick=5
oplot,(!pi+[ant.xp.be,ant.xm.be,ant.z.be])*!radeg,(!pi-[ant.xp.al,ant.xm.al,ant.z.al])*!radeg,psym=4,thick=5
oplot,(-!pi+[ant.xp.be,ant.xm.be,ant.z.be])*!radeg,(!pi-[ant.xp.al,ant.xm.al,ant.z.al])*!radeg,psym=4,thick=5
device,/close
set_plot,'x'

end
