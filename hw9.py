from __future__ import print_function
from __future__ import division
#get the sample rate
import matplotlib.pyplot as plt
import numpy as np
import csv

t = [[] for i in range(4)]
data = [[] for i in range(4)]
sr = [[] for i in range(4)]

##############################################
graphs_path = './Graphs/'

def make_fft(DAT, SR):
    N = len(DAT) # length of the signal
    k = np.arange(N) #frequencies
    T = N/SR
    frq = k/T # two sides frequency range
    frq = frq[range(int(N/2))] # one side frequency range
    Y = np.fft.fft(DAT)/N # fft computing and normalization
    Y = Y[range(int(N/2))]
    return [frq, Y]

for file_num in range(4):
    filename = 'sig{}.csv'.format(chr(file_num+65))
    with open(filename) as f:
        # open the csv file
        reader = csv.reader(f)
        for row in reader:
            t[file_num].append(float(row[0])) # leftmost column
            data[file_num].append(float(row[1])) # second column

    sample_rate = len(t[file_num]) / (t[file_num][-1] - t[file_num][0])
    sr[file_num] = sample_rate

############################################################################
# Functions for low pass filters
def low_pass_MAF(t, data, X):
    filtered = np.zeros([len(data)-X,1])
    for i in range(0, len(data)-X):
        filtered[i] = np.average(data[i:i+X-1])
    return [t[X:],filtered]

def low_pass_IIR(data, a, b):
    filtered = np.zeros([len(data),1])
    filtered[0] = b * data[0]
    for i in range(1, len(data)):
        filtered[i] = b*data[i] + a*filtered[i-1]
    return filtered

def apply_FIR(SR, cutoff, filt_len, DAT):
    #code for calculating filter from the FIR website
    fS = SR  # Sampling rate.
    fL = cutoff  # Cutoff frequency.
    N = filt_len  # Filter length, must be odd.
    h = np.sinc(2 * fL / fS * (np.arange(N) - (N - 1) / 2))
    h *= np.blackman(N)
    h /= np.sum(h)
    s = np.convolve(DAT,h)
    return s


#make some arrays for the filtered data, and number of points to average
tf_MAF = [[] for i in range(4)]
df_MAF = [[] for i in range(4)]
df_IIR = [[] for i in range(4)]
df_FIR = [[] for i in range(4)]

X = [1200,1000,100,100]
A = 0
B = 1
C = 2
D = 3
a = [0.999, 0.998, 0.02, 0.992]
b = [0]*4
for ii in range(4):
    b[ii] = 1 - a[ii]
    
#Cutoffs
#A  < 3Hz, SR=10kHz
#B: < 1Hz, SR=3.3kHz
#C: < 50Hz, SR=2.5kHz
#D: < 5Hz, SR=400Hz

#[[cutoff, filt_length]]
FIR_cutoffs= [[3, 1001],[0.5,4001],[75,11],[5,201]]



for file_num in range(4):
    [tf_MAF[file_num], df_MAF[file_num]] = low_pass_MAF(t[file_num], data[file_num], X[file_num])
    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True)
    plt.suptitle('MAF data for sig'+chr(file_num+65)+'  X = ' + str(X[file_num]))
    # plt.figure()

##########################################################################
    #MAF
    ax1.plot(t[file_num],data[file_num],'k',label='unfiltered')
    ax1.plot(tf_MAF[file_num],df_MAF[file_num],'r',label='filtered')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Amplitude')
    ax1.set_title('Unfiltered vs Filtered Data')
    ax1.legend()

    #make the ffts
    [ff,mf] = make_fft(df_MAF[file_num], sr[file_num]) #make fft of filtered data
    ax2.loglog(ff,abs(mf),'r', label='filtered') # plotting the fft
    [f,m] = make_fft(data[file_num], sr[file_num]) #make fft of original data
    ax2.loglog(f,abs(m),'b', label='unfiltered') # plotting the fft

    #plot the ffts
    ax2.set_xlabel('Freq (Hz)')
    ax2.set_ylabel('|Y(freq)|')
    ax2.set_title('Unfiltered vs filtered FFT')
    ax2.legend()
    plt.savefig(graphs_path + 'MAFfiltered_' + chr(file_num+65) +'.png')
    

###############################################################################
    #IIR
    df_IIR[file_num] = low_pass_IIR(data[file_num], a[file_num], b[file_num])
    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True)
    plt.suptitle('IIR on sig'+chr(file_num+65)+'  A = ' + str(a[file_num]) + ', B = ' + str(b[file_num]))

    ax1.plot(t[file_num],data[file_num],'k',label='unfiltered')
    ax1.plot(t[file_num],df_IIR[file_num],'r',label='filtered')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Amplitude')
    ax1.set_title('Unfiltered vs Filtered Data')
    ax1.legend()

    [ff,mf] = make_fft(df_IIR[file_num], sr[file_num]) #make fft of filtered data
    ax2.loglog(ff,abs(mf),'r', label='filtered') # plotting the fft
    ax2.loglog(f,abs(m),'b', label='unfiltered') # plotting the fft

    ax2.set_xlabel('Freq (Hz)')
    ax2.set_ylabel('|Y(freq)|')
    ax2.set_title('Unfiltered vs filtered FFT')
    ax2.legend()
    plt.savefig(graphs_path + 'IIRfiltered_' + chr(file_num+65) +'.png')

    #########################################################################
    # Use the FIR low pass

    df_FIR[file_num] = apply_FIR(sr[file_num], FIR_cutoffs[file_num][0],FIR_cutoffs[file_num][1], data[file_num])
    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True)
    plt.suptitle('FIR data for sig'+chr(file_num+65)+'  Cutoff = ' + str(FIR_cutoffs[file_num][0]) + ' Filter length = ' + str(FIR_cutoffs[file_num][1]))
    ax1.plot(t[file_num],data[file_num],'k',label='unfiltered')
    ax1.plot(t[file_num],df_FIR[file_num][:len(df_FIR[file_num])-FIR_cutoffs[file_num][1]+1],'r',label='filtered')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Amplitude')
    ax1.set_title('Unfiltered vs Filtered Data')
    ax1.legend()

    [ff,mf] = make_fft(df_FIR[file_num], sr[file_num]) #make fft of filtered data
    ax2.loglog(ff,abs(mf),'r', label='filtered') # plotting the fft
    # [f,m] = make_fft(data[file_num], sr[file_num]) #make fft of original data
    ax2.loglog(f,abs(m),'b', label='unfiltered') # plotting the fft
    ax2.set_xlabel('Freq (Hz)')
    ax2.set_ylabel('|Y(freq)|')
    ax2.set_title('Unfiltered vs filtered FFT')
    ax2.legend()
    plt.savefig(graphs_path + 'FIRfiltered_' + chr(file_num+65) +'.png')


    #fft for FIR data
