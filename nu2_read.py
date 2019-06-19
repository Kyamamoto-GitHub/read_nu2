import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import tkinter as tk
import tkinter.messagebox, tkinter.filedialog
import os


Tintegr = 5 #integration time (s)
Ncycle = 40 #number of cycles
Nhead = 60 #number of header lines
Ndetect = 22 #number of detectors
Zn68 = 6 #ordinal number of detector for 68Zn
Zn67 = 8 #ordinal number of detector for 67Zn
Zn66 = 10 #ordinal number of detector for 66Zn
Cu65 = 12 #ordinal number of detector for 65Cu
Zn64 = 14 #ordinal number of detector for 64Zn
Cu63 = 16 #ordinal number of detector for 63Cu
m68 = 67.92484 #atomic mass of 68Zn
m67 = 66.92713 #atomic mass of 68Zn
m66 = 65.92603 #atomic mass of 66Zn
m65 = 64.92779 #atomic mass of 65Cu
m64 = 63.92914 #atomic mass of 64Zn
m63 = 62.92960 #atomic mass of 63Cu
m_68_66 = np.log(m65/m63)/np.log(m68/m66)
m_66_64 = np.log(m65/m63)/np.log(m66/m64)
abundance_70 = 0.0061
abundance_68 = 0.1845
abundance_67 = 0.0404
abundance_66 = 0.2773
abundance_65 = 0.3085
abundance_64 = 0.4917
abundance_63 = 0.6915
#m = np.log(m65/m63)/np.log(m66/m64)
R_ref_68_66 = abundance_68/abundance_66 #reference isotope ratio for external correction
R_ref_66_64 = abundance_66/abundance_64 #reference isotope ratio for external correction
dif_value = 0. #standard d65Cu relative to primary standard
err_value = 0. #standard deviation of standard d65Cu


'''
dif = input('Please enter the standard d65Cu relative to NIST SRM 976.\t')
err = input('Please enter the standard deviation of standard d65Cu.\t')
dif_value = float(dif)
err_value = float(err)
'''

def NP2_read(dirname,Ndetect=22,Nhead=60):
	files = glob.glob(dirname+'/Data_*.csv')
	sheets = [[line.rstrip().split(',') for line in open(files[filenumber],'r')] for filenumber in range(len(files))]
	voltages = np.reshape([float(sheet[cycle+Nhead][column+1]) for sheet in sheets for cycle in range(len(sheet)-Nhead) for column in range(Ndetect*2)],(len(sheets),len(sheets[0])-Nhead,Ndetect*2))
	signals = voltages[:,:,Ndetect:]-voltages[:,:,:Ndetect]
	signals_mean = {'Mean':np.average(signals,axis=1),'SE':np.std(signals,axis=1,ddof=1)/np.sqrt(len(signals[0]))}
	return signals, signals_mean

def calc_IR(signals,signals_mean,dividend,divisor):
	R = (signals[1:len(signals)-1,:,dividend-1]-np.average([signals_mean['Mean'][0,dividend-1],signals_mean['Mean'][len(signals)-1,dividend-1]]))/(signals[1:len(signals)-1,:,divisor-1]-np.average([signals_mean['Mean'][0,divisor-1],signals_mean['Mean'][len(signals)-1,divisor-1]]))
	R_mean = {'Mean':np.average(R,axis=1),'SD':np.std(R,axis=1,ddof=1),'SE':np.std(R,axis=1,ddof=1)/np.sqrt(len(R[0]))}
	R_within2sigma = [[r for r in R[filenumber] if R_mean['Mean'][filenumber]-R_mean['SD'][filenumber]*2<=r<=R_mean['Mean'][filenumber]+R_mean['SD'][filenumber]*2] for filenumber in range(len(R))]
	R_mean_within2sigma = {'Mean':[np.average(r_2s) for r_2s in R_within2sigma],'SE':[np.std(r_2s,ddof=1)/np.sqrt(len(r_2s)) for r_2s in R_within2sigma]}
	return R, R_mean_within2sigma

def IR_corr(R,R_ex,R_ref,m_exp):
	R_corr = R*(R_ref/R_ex)**m_exp
	R_mean_corr = {'Mean':np.average(R_corr,axis=1),'SD':np.std(R_corr,axis=1,ddof=1),'SE':np.std(R_corr,axis=1,ddof=1)/np.sqrt(len(R[0]))}
	R_corr_within2sigma = [[r_corr for r_corr in R_corr[filenumber] if R_mean_corr['Mean'][filenumber]-R_mean_corr['SD'][filenumber]*2<=r_corr<=R_mean_corr['Mean'][filenumber]+R_mean_corr['SD'][filenumber]*2] for filenumber in range(len(R_corr))]
	R_mean_corr_within2sigma = {'Mean':[np.average(rc_2s) for rc_2s in R_corr_within2sigma],'SE':[np.std(rc_2s,ddof=1)/np.sqrt(len(rc_2s)) for rc_2s in R_corr_within2sigma]}
	return R_corr, R_mean_corr_within2sigma

def calc_delta(R_mean,dif_value=0,err_value=0):
	times = int((len(R_mean['Mean'])-1)/2)
	M1 = [R_mean['Mean'][time*2+1] for time in range(times)]; E1 = [R_mean['SE'][time*2+1] for time in range(times)]
	m2 = 1+dif_value/1000; e2 = err_value/1000
	M3 = [R_mean['Mean'][time*2] for time in range(times)]; E3 = [R_mean['SE'][time*2] for time in range(times)]
	M4 = [R_mean['Mean'][time*2+2] for time in range(times)]; E4 = [R_mean['SE'][time*2+2] for time in range(times)]
	delta_SE = [np.sqrt((m2/(m3+m4)*e1)**2+(m1/(m3+m4)*e2)**2+(m1*m2/(m3+m4)**2*e3)**2+(m1*m2/(m3+m4)**2*e4)**2)*2*1000 for m1,e1,m3,e3,m4,e4 in zip(M1,E1,M3,E3,M4,E4)]
	delta_value = {'Mean':[(m1*m2/(m3+m4)*2-1)*1000 for m1,m3,m4 in zip(M1,M3,M4)],'SE':delta_SE,'weight':[1/d_SE**2 for d_SE in delta_SE]}
	return delta_value


win = tk.Tk()
win.withdraw()
tkinter.messagebox.showinfo('Information','Please select a folder containing data files.')
dirname = tkinter.filedialog.askdirectory(initialdir=os.getcwd(),title='Select a folder containing data files')
print(dirname)

try:signals, signals_mean = NP2_read(dirname,Ndetect,Nhead)
except:print('Please select a folder containing data files.')

R_Cu, R_Cu_mean = calc_IR(signals,signals_mean,Cu65,Cu63)
#R_68_67, R_68_67_mean = calc_IR(signals,signals_mean,Zn68,Zn67)
R_68_66, R_68_66_mean = calc_IR(signals,signals_mean,Zn68,Zn66)
R_68_64, R_68_64_mean = calc_IR(signals,signals_mean,Zn68,Zn64)
R_67_66, R_67_66_mean = calc_IR(signals,signals_mean,Zn67,Zn66)
R_67_64, R_67_64_mean = calc_IR(signals,signals_mean,Zn67,Zn64)
R_66_64, R_66_64_mean = calc_IR(signals,signals_mean,Zn66,Zn64)

R_Cu_68_66, R_Cu_mean_68_66 = IR_corr(R_Cu,R_68_66,R_ref_68_66,m_68_66)
R_Cu_66_64, R_Cu_mean_66_64 = IR_corr(R_Cu,R_66_64,R_ref_66_64,m_66_64)

d65Cu = calc_delta(R_Cu_mean,0,0)
d65Cu_68_66 = calc_delta(R_Cu_mean_68_66,0,0)
d65Cu_66_64 = calc_delta(R_Cu_mean_66_64,0,0)
d_68_66 = calc_delta(R_68_66_mean,0,0)
d_68_64 = calc_delta(R_68_64_mean,0,0)
d_66_64 = calc_delta(R_66_64_mean,0,0)

ans = np.empty((3,2),float)
ans[0,0] = np.average(d65Cu['Mean'],weights=d65Cu['weight'])
ans[0,1] = np.sqrt(1/np.sum(d65Cu['weight']))
ans[1,0] = np.average(d65Cu_68_66['Mean'],weights=d65Cu_68_66['weight'])
ans[1,1] = np.sqrt(1/np.sum(d65Cu_68_66['weight']))
ans[2,0] = np.average(d65Cu_66_64['Mean'],weights=d65Cu_66_64['weight'])
ans[2,1] = np.sqrt(1/np.sum(d65Cu_66_64['weight']))
print('d65Cu = %.2f +- %.2f (2SE)\nd65Cu_68_66 = %.2f +- %.2f (2SE)\nd65Cu_66_64 = %.2f +- %.2f (2SE)'%(ans[0,0],ans[0,1]*2,ans[1,0],ans[1,1]*2,ans[2,0],ans[2,1]*2))
print('d65Cu = %.2f +- %.2f (2SE)\nd65Cu_68_66 = %.2f +- %.2f (2SE)\nd65Cu_66_64 = %.2f +- %.2f (2SE)'%(ans[0,0],ans[0,1]*2,ans[1,0],ans[1,1]*2,ans[2,0],ans[2,1]*2),file=open(dirname+'_delta-Cu.txt','w'))

plt.rcParams['font.size'] = 12
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['axes.grid'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig0 = plt.figure(figsize=(12,9))
std_meas = [meas for meas in range(1,len(R_Cu_mean['Mean'])+1,2)]
smpl_meas = [meas for meas in range(2,len(R_Cu_mean['Mean'])+1,2)]
std_Cu = [R_Cu_mean['Mean'][time] for time in range(0,len(R_Cu_mean['Mean']),2)]
std_Cu_err = [R_Cu_mean['SE'][time] for time in range(0,len(R_Cu_mean['SE']),2)]
smpl_Cu = [R_Cu_mean['Mean'][time] for time in range(1,len(R_Cu_mean['Mean']),2)]
smpl_Cu_err = [R_Cu_mean['SE'][time] for time in range(1,len(R_Cu_mean['SE']),2)]
std_68_66 = [R_68_66_mean['Mean'][time] for time in range(0,len(R_68_66_mean['Mean']),2)]
std_68_66_err = [R_68_66_mean['SE'][time] for time in range(0,len(R_68_66_mean['Mean']),2)]
smpl_68_66 = [R_68_66_mean['Mean'][time] for time in range(1,len(R_68_66_mean['Mean']),2)]
smpl_68_66_err = [R_68_66_mean['SE'][time] for time in range(1,len(R_68_66_mean['Mean']),2)]
std_68_64 = [R_68_64_mean['Mean'][time] for time in range(0,len(R_68_64_mean['Mean']),2)]
std_68_64_err = [R_68_64_mean['SE'][time] for time in range(0,len(R_68_64_mean['Mean']),2)]
smpl_68_64 = [R_68_64_mean['Mean'][time] for time in range(1,len(R_68_64_mean['Mean']),2)]
smpl_68_64_err = [R_68_64_mean['SE'][time] for time in range(1,len(R_68_64_mean['Mean']),2)]
std_66_64 = [R_66_64_mean['Mean'][time] for time in range(0,len(R_68_64_mean['Mean']),2)]
std_66_64_err = [R_66_64_mean['SE'][time] for time in range(0,len(R_68_64_mean['Mean']),2)]
smpl_66_64 = [R_66_64_mean['Mean'][time] for time in range(1,len(R_68_64_mean['Mean']),2)]
smpl_66_64_err = [R_66_64_mean['SE'][time] for time in range(1,len(R_68_64_mean['Mean']),2)]
ax1 = fig0.add_subplot(2,2,1)
ax2 = fig0.add_subplot(2,2,2)
ax3 = fig0.add_subplot(2,2,3)
ax4 = fig0.add_subplot(2,2,4)
ax1.errorbar(std_meas,std_Cu,yerr=std_Cu_err,marker='o',mfc='blue',mec='blue',linestyle='',markersize=8,ecolor='k',capsize=4,linewidth=0.8)
ax1.errorbar(smpl_meas,smpl_Cu,yerr=smpl_Cu_err,marker='o',mfc='orange',mec='orange',linestyle='',markersize=8,ecolor='k',capsize=4,linewidth=0.8)
ax2.errorbar(std_meas,std_68_66,yerr=std_68_66_err,marker='o',mfc='blue',mec='blue',linestyle='',markersize=8,ecolor='k',capsize=4,linewidth=0.8)
ax2.errorbar(smpl_meas,smpl_68_66,yerr=smpl_68_66_err,marker='o',mfc='orange',mec='orange',linestyle='',markersize=8,ecolor='k',capsize=4,linewidth=0.8)
ax3.errorbar(std_meas,std_68_64,yerr=std_68_64_err,marker='o',mfc='blue',mec='blue',linestyle='',markersize=8,ecolor='k',capsize=4,linewidth=0.8)
ax3.errorbar(smpl_meas,smpl_68_64,yerr=smpl_68_64_err,marker='o',mfc='orange',mec='orange',linestyle='',markersize=8,ecolor='k',capsize=4,linewidth=0.8)
ax4.errorbar(std_meas,std_66_64,yerr=std_66_64_err,marker='o',mfc='blue',mec='blue',linestyle='',markersize=8,ecolor='k',capsize=4,linewidth=0.8)
ax4.errorbar(smpl_meas,smpl_66_64,yerr=smpl_66_64_err,marker='o',mfc='orange',mec='orange',linestyle='',markersize=8,ecolor='k',capsize=4,linewidth=0.8)
ax1.set_xlabel('Measurement')
ax1.set_ylabel('$\mathrm{^{65}Cu/^{63}Cu}$')
ax1.get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
ax2.set_xlabel('Measurement')
ax2.set_ylabel('$\mathrm{^{68}Zn/^{66}Zn}$')
ax2.get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
ax3.set_xlabel('Measurement')
ax3.set_ylabel('$\mathrm{^{68}Zn/^{64}Zn}$')
ax3.get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
ax4.set_xlabel('Measurement')
ax4.set_ylabel('$\mathrm{^{66}Zn/^{64}Zn}$')
ax4.get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
plt.tight_layout()
ax1.set_yticklabels(['{:g}'.format(x) for x in ax1.get_yticks()])
ax2.set_yticklabels(['{:g}'.format(x) for x in ax2.get_yticks()])
ax3.set_yticklabels(['{:g}'.format(x) for x in ax3.get_yticks()])
ax4.set_yticklabels(['{:g}'.format(x) for x in ax4.get_yticks()])
plt.savefig(dirname+'/R_2SE.png')
plt.clf(); plt.close()

fig1 = plt.figure(figsize=(12,9))
ax1 = fig1.add_subplot(5,1,1); ax6 = ax1.twinx()
ax2 = fig1.add_subplot(5,1,2); ax7 = ax2.twinx()
ax3 = fig1.add_subplot(5,1,3); ax8 = ax3.twinx()
ax4 = fig1.add_subplot(5,1,4); ax9 = ax4.twinx()
ax5 = fig1.add_subplot(5,1,5); ax10 = ax5.twinx()
meas_t = [(cycle+filenumber*len(signals[filenumber]))*Tintegr for filenumber in range(1,len(signals)-1) for cycle in range(len(signals[filenumber]))]
ax1.scatter(meas_t,[signals[filenumber,cycle,Cu63-1] for filenumber in range(1,len(signals)-1) for cycle in range(len(signals[filenumber]))],s=10,color='red')
ax2.scatter(meas_t,[signals[filenumber,cycle,Cu65-1] for filenumber in range(1,len(signals)-1) for cycle in range(len(signals[filenumber]))],s=10,color='blue')
ax3.scatter(meas_t,[signals[filenumber,cycle,Zn64-1] for filenumber in range(1,len(signals)-1) for cycle in range(len(signals[filenumber]))],s=10,color='brown')
ax4.scatter(meas_t,[signals[filenumber,cycle,Zn66-1] for filenumber in range(1,len(signals)-1) for cycle in range(len(signals[filenumber]))],s=10,color='orange')
ax5.scatter(meas_t,[signals[filenumber,cycle,Zn68-1] for filenumber in range(1,len(signals)-1) for cycle in range(len(signals[filenumber]))],s=10,color='purple')
#ax6.scatter(meas_t,[signals[filenumber,cycle,Cu63-1]/1.602*10**8 for filenumber in range(1,len(signals)-1) for cycle in range(len(signals[filenumber]))],s=10,color='red')
#ax7.scatter(meas_t,[signals[filenumber,cycle,Cu65-1]/1.602*10**8 for filenumber in range(1,len(signals)-1) for cycle in range(len(signals[filenumber]))],s=10,color='blue')
#ax8.scatter(meas_t,[signals[filenumber,cycle,Zn64-1]/1.602*10**8 for filenumber in range(1,len(signals)-1) for cycle in range(len(signals[filenumber]))],s=10,color='brown')
#ax9.scatter(meas_t,[signals[filenumber,cycle,Zn66-1]/1.602*10**8 for filenumber in range(1,len(signals)-1) for cycle in range(len(signals[filenumber]))],s=10,color='orange')
#ax10.scatter(meas_t,[signals[filenumber,cycle,Zn68-1]/1.602*10**8 for filenumber in range(1,len(signals)-1) for cycle in range(len(signals[filenumber]))],s=10,color='purple')
ax1.tick_params(labelbottom=False)
ax2.tick_params(labelbottom=False)
ax3.tick_params(labelbottom=False)
ax4.tick_params(labelbottom=False)
ax1.get_xaxis().set_major_locator(ticker.MultipleLocator(len(signals[0])*Tintegr))
ax2.get_xaxis().set_major_locator(ticker.MultipleLocator(len(signals[0])*Tintegr))
ax3.get_xaxis().set_major_locator(ticker.MultipleLocator(len(signals[0])*Tintegr))
ax4.get_xaxis().set_major_locator(ticker.MultipleLocator(len(signals[0])*Tintegr))
ax5.get_xaxis().set_major_locator(ticker.MultipleLocator(len(signals[0])*Tintegr))
ax6.set_ylim([y/1.602*10**8 for y in ax1.get_ylim()])
ax7.set_ylim([y/1.602*10**8 for y in ax2.get_ylim()])
ax8.set_ylim([y/1.602*10**8 for y in ax3.get_ylim()])
ax9.set_ylim([y/1.602*10**8 for y in ax4.get_ylim()])
ax10.set_ylim([y/1.602*10**8 for y in ax5.get_ylim()])
ax5.set_xlabel('Time /s')
ax1.set_ylabel('$\mathrm{^{63}Cu}$ /V'); ax6.set_ylabel('$\mathrm{^{63}Cu}$ /cps')
ax2.set_ylabel('$\mathrm{^{65}Cu}$ /V'); ax7.set_ylabel('$\mathrm{^{65}Cu}$ /cps')
ax3.set_ylabel('$\mathrm{^{64}Zn}$ /V'); ax8.set_ylabel('$\mathrm{^{64}Zn}$ /cps')
ax4.set_ylabel('$\mathrm{^{66}Zn}$ /V'); ax9.set_ylabel('$\mathrm{^{66}Zn}$ /cps')
ax5.set_ylabel('$\mathrm{^{68}Zn}$ /V'); ax10.set_ylabel('$\mathrm{^{68}Zn}$ /cps')
#plt.grid(False,axis='y')
'''
ax1.set_yticklabels(['{:g}'.format(x) for x in ax1.get_yticks()])
ax2.set_yticklabels(['{:g}'.format(x) for x in ax2.get_yticks()])
ax3.set_yticklabels(['{:g}'.format(x) for x in ax3.get_yticks()])
ax4.set_yticklabels(['{:g}'.format(x) for x in ax4.get_yticks()])
ax5.set_yticklabels(['{:g}'.format(x) for x in ax5.get_yticks()])
ax6.set_yticklabels(['{:g}'.format(x) for x in ax6.get_yticks()])
ax7.set_yticklabels(['{:g}'.format(x) for x in ax7.get_yticks()])
ax8.set_yticklabels(['{:g}'.format(x) for x in ax8.get_yticks()])
ax9.set_yticklabels(['{:g}'.format(x) for x in ax9.get_yticks()])
ax10.set_yticklabels(['{:g}'.format(x) for x in ax10.get_yticks()])
'''
plt.tight_layout()
plt.savefig(dirname+'/temporal_signal.png')
plt.clf(); plt.close()

fig2 = plt.figure(figsize=(12,9))
ax1 = fig2.add_subplot(4,1,1); ax5 = ax1.twinx()
ax2 = fig2.add_subplot(4,1,2); ax6 = ax2.twinx()
ax3 = fig2.add_subplot(4,1,3); ax7 = ax3.twinx()
ax4 = fig2.add_subplot(4,1,4); ax8 = ax4.twinx()
ratio_temp_Cu = [R_Cu[filenumber,cycle] for filenumber in range(len(R_Cu)) for cycle in range(len(R_Cu[filenumber]))]
ratio_temp_68_66 = [R_68_66[filenumber,cycle] for filenumber in range(len(R_68_66)) for cycle in range(len(R_68_66[filenumber]))]
ratio_temp_68_64 = [R_68_64[filenumber,cycle] for filenumber in range(len(R_68_64)) for cycle in range(len(R_68_64[filenumber]))]
ratio_temp_66_64 = [R_66_64[filenumber,cycle] for filenumber in range(len(R_66_64)) for cycle in range(len(R_66_64[filenumber]))]
ax1.scatter(meas_t,ratio_temp_Cu,s=10,color='red',label='$\mathrm{^{65}Cu/^{63}Cu}$')
ax2.scatter(meas_t,ratio_temp_68_66,s=10,color='blue',label=('$\mathrm{^{68}Zn/^{66}Zn}$'))
ax3.scatter(meas_t,ratio_temp_68_64,s=10,color='orange',label=('$\mathrm{^{68}Zn/^{64}Zn}$'))
ax4.scatter(meas_t,ratio_temp_66_64,s=10,color='purple',label=('$\mathrm{^{66}Zn/^{64}Zn}$'))
ax1.tick_params(labelbottom=False)
ax2.tick_params(labelbottom=False)
ax3.tick_params(labelbottom=False)
ax1.get_xaxis().set_major_locator(ticker.MultipleLocator(len(signals[0])*Tintegr))
ax2.get_xaxis().set_major_locator(ticker.MultipleLocator(len(signals[0])*Tintegr))
ax3.get_xaxis().set_major_locator(ticker.MultipleLocator(len(signals[0])*Tintegr))
ax4.get_xaxis().set_major_locator(ticker.MultipleLocator(len(signals[0])*Tintegr))
ax1.set_ylim(min(ratio_temp_Cu)*2-np.average(ax1.get_ylim()),max(ratio_temp_Cu)*2-np.average(ax1.get_ylim()))
ax2.set_ylim(min(ratio_temp_68_66)*2-np.average(ax2.get_ylim()),max(ratio_temp_68_66)*2-np.average(ax2.get_ylim()))
ax3.set_ylim(min(ratio_temp_68_64)*2-np.average(ax3.get_ylim()),max(ratio_temp_68_64)*2-np.average(ax3.get_ylim()))
ax4.set_ylim(min(ratio_temp_66_64)*2-np.average(ax4.get_ylim()),max(ratio_temp_66_64)*2-np.average(ax4.get_ylim()))
ax5.set_ylim((ax1.get_ylim()[0]/np.average(ax1.get_ylim())-1)*1000,(ax1.get_ylim()[1]/np.average(ax1.get_ylim())-1)*1000)
ax6.set_ylim((ax2.get_ylim()[0]/np.average(ax2.get_ylim())-1)*1000,(ax2.get_ylim()[1]/np.average(ax2.get_ylim())-1)*1000)
ax7.set_ylim((ax3.get_ylim()[0]/np.average(ax3.get_ylim())-1)*1000,(ax3.get_ylim()[1]/np.average(ax3.get_ylim())-1)*1000)
ax8.set_ylim((ax4.get_ylim()[0]/np.average(ax4.get_ylim())-1)*1000,(ax4.get_ylim()[1]/np.average(ax4.get_ylim())-1)*1000)
ax4.set_xlabel('Time /s')
ax1.set_ylabel('$\mathrm{^{65}Cu/^{63}Cu}$'); ax5.set_ylabel('Deviation (‰)')
ax2.set_ylabel('$\mathrm{^{68}Zn/^{66}Zn}$'); ax6.set_ylabel('Deviation (‰)')
ax3.set_ylabel('$\mathrm{^{68}Zn/^{64}Zn}$'); ax7.set_ylabel('Deviation (‰)')
ax4.set_ylabel('$\mathrm{^{66}Zn/^{64}Zn}$'); ax8.set_ylabel('Deviation (‰)')
plt.tight_layout()
ax1.set_yticklabels(['{:g}'.format(x) for x in ax1.get_yticks()])
ax2.set_yticklabels(['{:g}'.format(x) for x in ax2.get_yticks()])
ax3.set_yticklabels(['{:g}'.format(x) for x in ax3.get_yticks()])
ax4.set_yticklabels(['{:g}'.format(x) for x in ax4.get_yticks()])
plt.savefig(dirname+'/temporal_ratio.png')
plt.clf(); plt.close()

fig3 = plt.figure()
plt.errorbar(d_68_64['Mean'],d_66_64['Mean'],xerr=d_68_64['SE'],yerr=d_66_64['SE'],marker='o',linestyle='',markersize=8,color='k',capsize=4,linewidth=0.8)
#plt.scatter(d_68_64[:,0],d_66_64[:,0],s=10,color='black')
x = np.array(plt.xlim())
y = x*((66-64)/66*64)/((68-64)/68*64)
plt.plot(x,y,color='black')
plt.xlabel('$\mathrm{\delta^{68}Zn/^{64}Zn}$ (‰)')
plt.ylabel('$\mathrm{\delta^{66}Zn/^{64}Zn}$ (‰)')
plt.tight_layout()
plt.gca().set_xticklabels(['{:g}'.format(x) for x in plt.gca().get_xticks()])
plt.gca().set_yticklabels(['{:g}'.format(x) for x in plt.gca().get_yticks()])
plt.savefig(dirname+'/Zn_three_isotope.png')
plt.clf(); plt.close()

'''
#import tkinter as tk
#import tkinter.messagebox, tkinter.filedialog
import sys

win = tk.Tk()
scr_w, scr_h = win.winfo_screenwidth(), win.winfo_screenheight()
#win.geometry('%dx%d+%d+%d'%(scr_w/2,scr_h/2,scr_w/4,scr_h/4))
win.title(u'Menu') #Japanese can be used with unicode.

def callback(Press_OK):
	win.withdraw();	win.quit() #Error occurs replacing this line with win.destroy().
def close(Press_esc):
	sys.exit(u'Processing was discontinued.') #sys.exit() can be replaced by quit() or exit() without exit status.
def delete():
	sys.exit(u'Processing was discontinued.') #sys.exit() can be replaced by quit() or exit() without exit status.

f0 = tk.LabelFrame(win,text='Reference',font=('',12,'bold'))
f1 = tk.LabelFrame(win,text='Number of decimal places',font=('',12,'bold'))

radio = tk.IntVar()
radio.set(0)
tk.Radiobutton(f0,text=u'Select a batch result file',font=('',12),variable=radio,value=0).grid(row=0,column=0)
tk.Radiobutton(f0,text=u'Select a folder containing data files',font=('',12),variable=radio,value=1).grid(row=0,column=1,sticky=tk.W+tk.E)

digit = tk.IntVar()
digit.set(2)
tk.Radiobutton(f1,text=u'1',font=('',12),variable=digit,value=1).grid(row=0,column=0)
tk.Radiobutton(f1,text=u'2',font=('',12),variable=digit,value=2).grid(row=0,column=1)
tk.Radiobutton(f1,text=u'3',font=('',12),variable=digit,value=3).grid(row=0,column=2)
tk.Radiobutton(f1,text=u'4',font=('',12),variable=digit,value=4).grid(row=0,column=3)
tk.Radiobutton(f1,text=u'5',font=('',12),variable=digit,value=5).grid(row=0,column=4)

f0.columnconfigure(0,weight=1); f0.columnconfigure(1,weight=1)
f1.columnconfigure(0,weight=1); f1.columnconfigure(1,weight=1); f1.columnconfigure(2,weight=1); f1.columnconfigure(3,weight=1); f1.columnconfigure(4,weight=1)
f0.pack(padx=scr_w/80,pady=scr_h/80,fill='both')
f1.pack(padx=scr_w/80,pady=scr_h/80,fill='both')

check = tk.BooleanVar()
check.set(True)
tk.Checkbutton(text=u'External correction',font=('',12),variable=check).pack(padx=scr_w/80,pady=scr_h/80)

button = tk.Button(win,text=u' OK ')
button.bind('<ButtonPress-1>',callback)
button.bind('<KeyPress-Return>',callback)
button.bind('<KeyPress-Escape>',close)
button.pack(padx=scr_w/80,pady=scr_h/80)
button.focus_set()

win.protocol('WM_DELETE_WINDOW',delete) #overriding Tkinter "X" button control
win.mainloop()
print(radio.get(),digit.get(),check.get())

tkinter.messagebox.showinfo('Information','Please select a batch result file.')
filename = tkinter.filedialog.askopenfilename(filetypes=[('CSV File','*.csv')],initialdir=os.getcwd(),title='Open a batch result file')
tkinter.messagebox.showinfo('Information','Please select a folder containing data files.')
dirname = tkinter.filedialog.askdirectory(initialdir=os.getcwd(),title='Select a folder containing data files')
print('%s\n%s'%(filename,dirname))
'''
