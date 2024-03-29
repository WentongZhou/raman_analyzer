import os,collections,glob,bisect
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import matplotlib.font_manager as fm
import rampy as rp
import scipy
from scipy import signal
import pandas as pd
import lmfit
import warnings
warnings.filterwarnings("ignore")


class raman_analyzer:
    """
    Author: Wentong Zhou
    This class is used to analyze the Raman spectrum.

    parameters:

    name: the name of the spectrum file
    min: the minimum wavenumber of the spectrum
    max: the maximum wavenumber of the spectrum
    filter: the filter used to identify the peak
    noise_min: the baseline noise minimum
    noise_max: the baseline noise maximum
    
    """
    def __init__(self,name,min,max,filter,noise_min=1700,noise_max=2400):
        self.name = name
        self.min = min
        self.max = max
        self.noise_min = noise_min
        self.noise_max = noise_max
        self.filter = filter
        self.peak_pos = []
        if self.name[-2] == 'p':
            self.spectrum = np.genfromtxt(self.name,delimiter=',')
        else:
            self.spectrum = np.genfromtxt(self.name)
        self.raman()
        self.peak_finder()
        self.peak_organizer()
        self.peak_fwhm()
        self.noise_filter()
    
    def raman(self):
        """
        This function is used to load the Raman spectrum raw data and perform baseline correction, smoothing and normalization.
        For graphene-based materials, the baseline correction is performed using the drPLS method, 
        and the smoothing is performed using the Whittaker method.

        """
        x_new = np.arange(self.min, self.max, 0.5)
        y_new = rp.resample(self.spectrum[:,0], self.spectrum[:,1], x_new)
        self.spectrum_resample = np.vstack((x_new,y_new)).T
        y_smo_10 = rp.smooth(self.spectrum_resample[:,0],self.spectrum_resample[:,1],method="whittaker",Lambda=3000,window_length=7)
        self.spectrum_resample = np.vstack((x_new,y_smo_10)).T
        d = scipy.signal.argrelextrema(self.spectrum_resample, np.less)
        points = len(d[0])
        bir=np.zeros((points,2))
        for i in range(len(d[0])):
            wavenumber=self.spectrum_resample[d[0][i]][0]
            bir[i][0]=wavenumber
            bir[i][1]=wavenumber+5
        y_corr, self.y_base = rp.baseline(self.spectrum_resample[:,0],self.spectrum_resample[:,1],bir,'arPLS', lam=10**7, ratio=0.001)
        x = self.spectrum_resample[:,0]
        x_fit = x[np.where((x > self.min)&(x < self.max))]
        y_fit = y_corr[np.where((x > self.min)&(x < self.max))]
        self.spectrum_corr = np.column_stack((x_fit,y_fit))
        self.ese0 = np.sqrt(abs(y_fit[:,0]))/abs(y_fit[:,0]) # the relative errors after baseline subtraction
        y_fit[:,0] = y_fit[:,0]/np.amax(y_fit[:,0])*10 # normalise spectra to maximum intensity, easier to handle
        self.spectrum_fit = np.column_stack((x_fit,y_fit))
        self.sigma = abs(self.ese0*y_fit[:,0])

    def peak_finder(self):
        """
        use scipy.signal.find_peaks to find the peaks in the spectrum
        
        """
        # find_peaks adn argrelmin are exporting the indices of the peaks and saddles
        peaks_1 = scipy.signal.find_peaks(self.spectrum_fit.T[1])
        saddles = scipy.signal.argrelmin(self.spectrum_fit.T[1])
        peaks = np.insert(peaks_1[0],0,saddles[0])
        peaks = np.sort(peaks)
        self.peak_pos = peaks
        self.peak_wavenumber,self.peak_signal,self.peak_abs_signal = [],[],[]
        for peak in peaks:
            self.peak_wavenumber.append(self.spectrum_fit.T[0][peak])
            self.peak_signal.append(self.spectrum_fit.T[1][peak])
            self.peak_abs_signal.append(self.spectrum_corr.T[1][peak])
    
    def peak_organizer(self):
        """
        Get the peak position and peak intensity of the spectrum based on a filter value.
        
        """
        self.peak_pos_filtered,self.peak_filtered,self.signal_filtered,self.abs_signal_filtered = [],[],[],[]
        for i in range(len(self.peak_wavenumber)):
            if self.peak_signal[i] > self.filter:
                self.peak_pos_filtered.append(self.peak_pos[i])
                self.peak_filtered.append(self.peak_wavenumber[i])
                self.signal_filtered.append(self.peak_signal[i])
                self.abs_signal_filtered.append(self.peak_abs_signal[i])


    def peak_fwhm(self):
        """
        Calculate the FWHM of the peaks in the spectrum.
        """
        pw = signal.peak_widths(self.spectrum_fit[:,1],self.peak_pos_filtered,rel_height=0.495)
        pw_1 = signal.peak_widths(self.spectrum_fit[:,1],self.peak_pos_filtered,rel_height=0.97)
        f = scipy.interpolate.interp1d(range(len(self.spectrum_fit[:,0])),self.spectrum_fit[:,0])
        x_left = f(pw[2])
        x_right = f(pw[3])
        self.fwhm = x_right - x_left
        self.fwhm_hline = np.vstack((pw[1],x_left,x_right))

    def noise_filter(self):
        """
        Calculate the baseline noise of the spectrum.
        """
        noise = self.spectrum_fit[:, 1][np.where((self.spectrum_fit[:, 0] > self.noise_min) & (self.spectrum_fit[:, 0] < self.noise_max))]
        self.noise_up = np.max(noise)
        self.noise_down = np.min(noise)
        self.peak_threshold = (np.max(noise) - np.min(noise))




def raman_batch(dirs,min,max,filter,show=False,normalized=True,export=False,noise_test=False,noise_min=1700,noise_max=2400):
    """
    This function is used to analyze the Raman spectrum in batch.
    dirs: the directory of the spectrum files
    min: the minimum wavenumber of the spectrum
    max: the maximum wavenumber of the spectrum
    filter: the filter used to identify the peak
    show: whether to show the spectrum
    export: whether to export the baseline corrected spectrum data
    noise_test: whether to show the baseline noise
    noise_min: the baseline noise wavenumber minimum
    noise_max: the baseline noise wavenumber maximum

    return values:

    peaks: the peak position of the spectrum
    signals: the peak intensity of the spectrum
    fwhm: the FWHM of the peaks in the spectrum

    """
    g = os.walk(dirs)
    if export == True:
        os.system(f'mkdir -p {dirs}_trans')
    peaks = []
    signals = []
    fwhm = []
    fwhm_test = []
    name = []
    peak_threshold = []
    for path,dir_list,file_list in g:
        for file_name in file_list:
            a = os.path.join(path,file_name)
            raman = raman_analyzer(a,min,max,filter,noise_min,noise_max)
            if export == False:
                peaks.append(raman.peak_filtered)
                if normalized == True:
                    signals.append(raman.signal_filtered)
                else:
                    signals.append(raman.abs_signal_filtered)
            peak_threshold.append(raman.peak_threshold)
            name.append(file_name)
            fwhm.append(raman.fwhm)
            fwhm_test.append(raman.fwhm_hline[0])
            if export == True:
                np.savetxt(f'./{dirs}_trans/{file_name}.txt', raman.spectrum_fit,fmt='%10.5f')
                peaks.append(raman.spectrum_fit[:,0])
                if normalized == True:
                    signals.append(raman.spectrum_fit[:,1])
                else:
                    signals.append(raman.spectrum_corr[:,1])
            if show == True:
                raman_plot(raman,min,max)
            if noise_test == True:
                noise_plot(raman,min,max)
    return peaks,signals,fwhm,fwhm_test,name,peak_threshold

def raman_plot(a,min,max):
    """
    This function is used to plot the spectrum with raw data, baseline corrected normalized data.
    """
    plt.figure(figsize=(10,5))
    plt.subplot(1,2,1)
    plt.plot(a.spectrum_resample[:,0],a.spectrum_resample[:,1],'ko',markersize=1,label='spectrum_resample')
    plt.plot(a.spectrum_corr[:,0],a.spectrum_corr[:,1],'ro',markersize=1,label='spectrum_corr')
    plt.plot(a.spectrum_resample[:,0],a.y_base,'g--',label='baseline',linewidth=2)
    plt.xlim(min,max)
    plt.xlabel("Raman shift, cm$^{-1}$", fontsize = 12)
    plt.ylabel("Normalized intensity, a. u.", fontsize = 12)
    plt.legend(fontsize=8)
    plt.subplot(1,2,2)  
    plt.plot(a.spectrum_fit[:,0],a.spectrum_fit[:,1],'#F5420a',marker='.',markersize=3,label='spectrum_fit')
    plt.scatter(a.peak_filtered,a.signal_filtered,s=100,c='#0af57c',marker='o')
    plt.hlines(*a.fwhm_hline,colors='#2596be').set_linewidth(2)
    plt.ylim(-1,11)
    plt.xlim(min,max)
    plt.title(a.name)
    plt.xlabel("Raman shift, cm$^{-1}$", fontsize = 12)
    plt.show()


def noise_plot(a,min,max):
    """
    This function is used to plot the spectrum with baseline noise
    """
    plt.figure(figsize=(10,5)) 
    plt.plot(a.spectrum_fit[:,0],a.spectrum_fit[:,1],'#F5420a',markersize=3,label='spectrum_fit')
    plt.scatter(a.peak_filtered,a.signal_filtered,s=100,c='#0af57c',marker='o')
    plt.hlines(*a.fwhm_hline,colors='#2596be').set_linewidth(2)
    plt.axhline(y=a.noise_up,linestyle='--',color='#2596be')
    plt.axhline(y=a.noise_down,linestyle='--',color='#2596be')
    plt.ylim(-1,11)
    plt.xlim(min,max)
    plt.title(a.name)
    plt.ylabel("Normalized intensity, a. u.", fontsize = 12)
    plt.xlabel("Raman shift, cm$^{-1}$", fontsize = 12)
    plt.show()



def ratio_calculator(directory,ref = "g",normalized = True,**peak_param):
    """
    This function is used to calculate the ratio of the peaks in batch
    parameters:

    directory: the directory of the spectrum files
    seive: the seive of the spectrum files
    peak_info: the peak information of the spectrum files
    
    """
    peak_data = collections.defaultdict(dict)

    for peak_name in ['total',ref]:
        if peak_name == 'total':
            peaks,signals,fwhm,fwhm_test,names,peak_threshold = raman_batch(f'{directory}',peak_param[peak_name][0],\
                                                                            peak_param[peak_name][1],peak_param[peak_name][2],\
                                                                            normalized=normalized,noise_min=peak_param[peak_name][3],\
                                                                            noise_max=peak_param[peak_name][4],export=True)
        else:
            peaks,signals,fwhm,fwhm_test,names,peak_threshold = raman_batch(f'{directory}_trans',peak_param[peak_name][0],\
                                                                            peak_param[peak_name][1],peak_param[peak_name][2],\
                                                                            normalized=normalized,noise_min=peak_param[peak_name][0],\
                                                                            noise_max=peak_param[peak_name][1],export=False)
        # store the peak information in a dictionary
        for idx, name in enumerate(names):
            if peak_name == 'total':
                peak_data[f'{name}.txt'][peak_name] = {"peak":peaks[idx],"signal":signals[idx],"fwhm":fwhm[idx]}
            else:
                peak_value = peaks[idx][0]
                peak_signal = peak_data[f'{name}']['total']['signal'][int((peak_value - peak_param['total'][0])*2)]
                peak_data[name][peak_name] = {"peak":peak_value,"signal":peak_signal,'fwhm':fwhm[idx][0],f'{peak_name}_{ref}':1}

    for peak_name, peak_info in peak_param.items():
        if len(peak_info) == 3 and peak_info[2] != "auto" and peak_name != ref:
            peaks,signals,fwhm,fwhm_test,names,peak_threshold = raman_batch(f'{directory}_trans',peak_info[0],\
                                                                            peak_info[1],peak_info[2],noise_min=peak_info[0],\
                                                                            normalized=normalized,noise_max=peak_info[1],export=False)
            for idx, name in enumerate(names):
                peak_value = peaks[idx][0]
                peak_signal = peak_data[f'{name}']['total']['signal'][int((peak_value - peak_param['total'][0])*2)]
                peak_ratio = peak_signal/peak_data[name][ref]['signal']
                peak_data[name][peak_name] = {"peak":peak_value,"signal":peak_signal,'fwhm':fwhm[idx][0],f'{peak_name}_{ref}':peak_ratio}

        elif len(peak_info) == 3 and peak_info[2] == "auto" and peak_name != ref:
            for file in glob.glob(f"{directory}/*.dpt"):
                name = file.split('/')[-1] + '.txt'
                test = raman_analyzer(file,peak_param['total'][0],peak_param['total'][1],0.2)
                if normalized == True:
                    fit_spectrum = pd.DataFrame(test.spectrum_fit)
                else:
                    fit_spectrum = pd.DataFrame(test.spectrum_corr)
                # fit_spectrum.plot(x=0,y=1,title=file_name)
                d_peak = fit_spectrum[(fit_spectrum.iloc[:,0] <= peak_info[1]) & (fit_spectrum.iloc[:,0] >= peak_info[0])]
                # d_peak.plot(x=0,y=1,title=file_name)
                d_peak_max = d_peak.iloc[:,1].max()
                peak_value = float(d_peak[fit_spectrum.iloc[:,1] == d_peak_max].iloc[:,0].tolist()[0])
                peak_signal = peak_data[f'{name}']['total']['signal'][int((peak_value - peak_param['total'][0])*2)]
                peak_signal = peak_signal if peak_signal > test.peak_threshold else 0
                peak_ratio = peak_signal/peak_data[name][ref]['signal']
                peak_indice = bisect.bisect_left(test.peak_filtered,peak_value)
                fwhm = test.fwhm[peak_indice] if peak_signal != 0 else 0
                peak_data[name][peak_name] = {"peak":peak_value,"signal":abs(peak_signal),'fwhm':fwhm,f'{peak_name}_{ref}':peak_ratio}

    file_name = [peak for peak in peak_data.keys()]
    all_peaks = []
    columns = ['name']
    for peak_name in peak_param:
        if peak_name != 'total':
            try:
                all_peaks.append([peak_data[name][peak_name]["peak"] for name in file_name])
                all_peaks.append([peak_data[name][peak_name]["signal"] for name in file_name])
                all_peaks.append([peak_data[name][peak_name]["fwhm"] for name in file_name])
                all_peaks.append([peak_data[name][peak_name][f'{peak_name}_{ref}'] for name in file_name])
                columns.append(f'{peak_name}_peak')
                columns.append(f'{peak_name}_signal')
                columns.append(f'{peak_name}_fwhm')
                columns.append(f'I{peak_name}_I{ref}')
            except:
                print(f'Error: {name} of {peak_name} is not in the spectrum')

    all_peaks = pd.DataFrame(np.array(all_peaks).T)
    # add two more rows about average and std
    all_peaks.loc['average'] = all_peaks.mean()
    all_peaks.loc['std'] = all_peaks.std()
    # add a new column about the file name
    file_name += ['average','std']
    all_peaks.insert(0,'name',file_name)
    all_peaks.columns = columns
    # export the data as csv file
    all_peaks.to_csv(f'{directory}_ratio.csv',index=False)
      
    return peak_data, all_peaks


class raman_fitting(raman_analyzer):

    def __init__(self,name,fit_min,fit_max,kw_fn={},\
                 fit_algo='nelder',fit_type='gaussian',\
                 peak_num=5,peak_init=[],min=1100,max=3000,filter=5):
        self.fit_min = fit_min
        self.fit_max = fit_max
        self.fit_algo = fit_algo
        self.kw_fn = kw_fn
        self.peak_num = peak_num
        self.peak_init = peak_init
        self.fit_type = fit_type
        super().__init__(name, min, max, filter, min, max)
        self.fit_gen()
        self.fit_minimize()
        self.fit_export()

    def fit_preplot(self):
        """
        This function is used to plot the spectrum with raw data, baseline corrected normalized data.
        """
        raman_plot(self, self.fit_min, self.fit_max)

    def fit_gen(self,amplitude:float=3,wide:float=15):
        """
        This function is generate inital guess for the peak fitting.
        """
        # get self.spectrum_fit within range of fit_min and fit_max
        self.fit_range = self.spectrum_fit[(self.spectrum_fit[:,0] >= self.fit_min) & (self.spectrum_fit[:,0] <= self.fit_max)]
        self.fit_peaks = np.array(self.peak_filtered)
        self.fit_peaks = self.fit_peaks[(self.fit_peaks >= self.fit_min) & (self.fit_peaks <= self.fit_max)]
        self.fit_params = lmfit.Parameters()
        if self.peak_init == []:
            n = len(self.fit_peaks)
            if n > self.peak_num:
                self.fit_peaks = self.fit_peaks[:self.peak_num]
            elif n < self.peak_num:
                range_gap = self.peak_num - len(self.fit_peaks)
                self.fit_peaks = np.append(self.fit_peaks, np.linspace(self.fit_min + 20, self.fit_max - 20, range_gap))
        else:
            self.fit_peaks = self.peak_init
            n = len(self.peak_init)
            self.peak_num = len(self.peak_init)
        # build the initial guess for the fitting and check lmfit documentation for more details
        for i in range(self.peak_num):
            if i < n and self.peak_init[i] != 1590:
                self.fit_params.add_many((f'a{i}',amplitude,True,0,None,None),
                                         (f'f{i}',self.fit_peaks[i],True,self.fit_peaks[i]-5,self.fit_peaks[i]+5,None),
                                         (f'l{i}',wide,True,0,30,None))
            elif i < n and self.peak_init[i] == 1590:
                self.fit_params.add_many((f'a{i}',amplitude,True,0,None,None),
                                         (f'f{i}',self.fit_peaks[i],True,self.fit_peaks[i]-5,self.fit_peaks[i]+5,None),
                                         (f'l{i}',wide,True,0,30,None))

            else:
                self.fit_params.add_many((f'a{i}',amplitude,True,0,None,None),
                                         (f'f{i}',self.fit_peaks[i],True,self.fit_peaks[i]-55,self.fit_peaks[i]+55,None),
                                         (f'l{i}',wide,True,0,50,None))
            
    def fit_residue(self,par,x,data=None,eps=None):
        """
        This function is used to calculate the residue of the fitting.
        for more details, check lmfit documentation.
        """
        model = np.zeros(len(self.fit_range[:,0]))
        fit_y = {}
        for i in par:
            locals()[i] = par[i].value
        for i in range(self.peak_num):
            if self.fit_type == 'gaussian':
                locals()[f'peak{i}'] = rp.gaussian(x,locals()[f'a{i}'],locals()[f'f{i}'],locals()[f'l{i}'])
                fit_y[f'peak{i}'] = locals()[f'peak{i}']
            elif self.fit_type == 'lorentzian':
                locals()[f'peak{i}'] = rp.lorentzian(x,locals()[f'a{i}'],locals()[f'f{i}'],locals()[f'l{i}'])
                fit_y[f'peak{i}'] = locals()[f'peak{i}']
            model += locals()[f'peak{i}']
        # if we don't have data, the function only returns the direct calculation
        if data is None: 
            return model, fit_y
        # without errors, no ponderation
        if eps is None: 
            return (model - data)
        # with errors, the difference is ponderated
        return (model - data) / eps

    def fit_minimize(self):
        """
        This function is used to minimize the residue of the fitting.
        """
        for i in range(self.peak_num):
            self.fit_params[f'f{i}'].vary = False
        self.fit_result = lmfit.minimize(self.fit_residue,self.fit_params,method = self.fit_algo,\
                                         args=(self.fit_range[:,0],self.fit_range[:,1]),**self.kw_fn)
        for i in range(self.peak_num):
            self.fit_params[f'f{i}'].vary = True
        self.fit_result1 = lmfit.minimize(self.fit_residue,self.fit_params,method = self.fit_algo,\
                                         args=(self.fit_range[:,0],self.fit_range[:,1]),**self.kw_fn)

    def fit_plot(self, save=False):
        """
        This function is used to plot the fitting result.
        """
        # get the script path
        colors = ["#C4421A", "#1AC441", "#1A42C4", "#C4C41A", "#C41AC4", "#1AC4C4"]
        script_path = os.path.dirname(os.path.abspath(__file__))
        # Add the font files to Matplotlib's font manager.
        font_path = script_path + '/Times_New_Roman.ttf'
        bold_font_path = script_path + '/Times_New_Roman_Bold.ttf'
        mpl.rcParams['axes.labelweight'] = 'bold'
        # Add the font to Matplotlib's font manager.
        fm.fontManager.addfont(font_path)
        fm.fontManager.addfont(bold_font_path)
        plt.rcParams['font.family'] = 'Times New Roman'
        plt.rcParams['mathtext.fontset'] = 'stix'
        fig, ax = plt.subplots(figsize=(6, 4))
        plt.gca().spines['top'].set_linewidth(1.5)
        plt.gca().spines['right'].set_linewidth(1.5)
        plt.gca().spines['bottom'].set_linewidth(1.5)
        plt.gca().spines['left'].set_linewidth(1.5)
        plt.tick_params(which='major',direction='out',width=1,length=4.5,labelsize=12)

        plt.tick_params(which='minor',direction='out',width=1,length=3,labelsize=12)
        # set the minor ticks
        ax.xaxis.set_major_locator(MultipleLocator((self.fit_max - self.fit_min)//5))
        minor_locator_x = AutoMinorLocator(2)
        minor_locator_y = AutoMinorLocator(2)
        ax.xaxis.set_minor_locator(minor_locator_x)
        ax.yaxis.set_minor_locator(minor_locator_y)
        x_fit = self.fit_range[:,0]
        y_fit = self.fit_range[:,1]
        self.model = lmfit.fit_report(self.fit_result1.params)
        yout,fit_y = self.fit_residue(self.fit_result1.params,x_fit)
        plt.scatter(x_fit[::6],y_fit[::6],s=16,marker='o',color='k',label='experimental data')
        plt.plot(x_fit,yout,color='#E6750A',linewidth=2,label='fitting data')
        for i in range(self.peak_num):
            plt.plot(x_fit,fit_y[f'peak{i}'],color=colors[i%len(colors)],linewidth=2)
        plt.ylim(-0.5,10.5)
        plt.xlim(self.fit_min,self.fit_max)
        plt.xlabel("Raman shift (cm$^{\mathbf{-1}}$)", fontsize = 14,fontweight='bold')
        plt.ylabel("Normalized intensity (a. u.)", fontsize = 14,fontweight='bold')
        # add border of legend
        plt.legend(frameon=True,edgecolor='black',facecolor='white',fontsize='small')
        if save == True:
            plt.savefig(f'{self.name.split("/")[-1]}_fit.jpg',dpi=1000,bbox_inches='tight')
        plt.title(self.name.split('/')[-1],fontsize=10)
        plt.show()

    def fit_export(self):
        """
        This function is used to export the fitting result.
        """
        self.fit_report = []
        for i in range(self.peak_num):
            if i > 0:
                single_report = []
                single_report.append(self.fit_result1.params[f'f{i}'].value)
                single_report.append(self.fit_result1.params[f'l{i}'].value * 2)
                single_report.append(self.fit_result1.params[f'a{i}'].value)
                single_report.append(rp.peakarea(self.fit_type,\
                                                 pos = self.fit_result1.params[f'f{i}'].value,\
                                                 amp=self.fit_result1.params[f'a{i}'].value,\
                                                 HWHM=self.fit_result1.params[f'l{i}'].value))
                single_report.append(self.fit_result1.params[f'a{i}']/self.fit_report[0][2])
                single_report.append(single_report[3]/self.fit_report[0][3])
            else:
                single_report = []
                single_report.append(self.fit_result1.params[f'f{i}'].value)
                single_report.append(self.fit_result1.params[f'l{i}'].value * 2)
                single_report.append(self.fit_result1.params[f'a{i}'].value)
                single_report.append(rp.peakarea(self.fit_type,\
                                                 pos=self.fit_result1.params[f'f{i}'].value,\
                                                 amp=self.fit_result1.params[f'a{i}'].value,\
                                                 HWHM=self.fit_result1.params[f'l{i}'].value))
                single_report.append(1)
                single_report.append(1)
            self.fit_report.append(single_report)
        self.fit_report.sort()


def raman_fitting_batch(directory,fit_min,fit_max,min=1100,max=3000,filter=5,kw_fn={},peak_init=[],peak_num=5,export=False,save=False,fit_algo='powell',fit_type='gaussian'):
    """
    This function is used to perform peak fitting in batch and export the fitting result.
    params:
    directory: the directory of the spectrum files
    name: the name list of the spectrum files
    fit_min: the minimum wavenumber of the fitting
    fit_max: the maximum wavenumber of the fitting
    peak_num: the number of peaks to be fitted
    fit_algo: the fitting algorithm
    fit_type: the fitting type, gaussian or lorentzian
    """
    # create an empty array with 3 * peak_num columns
    fit_reports = []
    names = []
    redchi = []
    # loop all the files in the directory
    for file in glob.glob(f"{directory}/*.dpt"):
        result = raman_fitting(file,fit_min,fit_max,kw_fn,fit_algo,fit_type,peak_num,peak_init,min,max,filter)
        names.append(file.split('/')[-1])
        redchi.append(result.fit_result1.redchi)
        if export == True:
            result.fit_plot(save=save)
            if not os.path.exists(f'{directory}_fit'):
                os.mkdir(f'{directory}_fit')
            os.system(f'mv {file.split("/")[-1]}_fit.jpg {directory}_fit')
        # export the fitting result
        fit_report = np.array(result.fit_report).T
        fit_report = fit_report.reshape(1,-1)
        fit_reports.append(fit_report[0])
    columns = [f'peak{i}' for i in range(result.peak_num)] + \
              [f'fwhm{i}' for i in range(result.peak_num)] + \
              [f'intensity{i}' for i in range(result.peak_num)] + \
              [f'area{i}' for i in range(result.peak_num)] + \
              [f'I_ratio{i}' for i in range(result.peak_num)] + \
              [f'area_ratio{i}' for i in range(result.peak_num)]



    fit_reports = pd.DataFrame(fit_reports,columns=columns)
    # add the reduced chi square as the last column
    fit_reports.insert(len(fit_reports.columns),'redchi',redchi)

    # calculate the average and std of the fitting result
    fit_reports.loc['average'] = fit_reports.mean()
    fit_reports.loc['std'] = fit_reports.std()
    # add the names as the first column
    names.append('average')
    names.append('std')
    fit_reports.insert(0,'name',names)
    # export the fitting result as csv file
    fit_reports.to_csv(f'{directory}_fit.csv',index=False)
    return fit_reports


        
        
    







    

