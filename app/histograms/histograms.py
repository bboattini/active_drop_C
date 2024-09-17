import os
import app.aux_func as af
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")

PARAM = {
    'var': 'B',
    'a': 0,
    'h': 0,
    'samp_frac': 10,
    'bin': 1000
}

def histogram_of(var=None, a=None, h=None):
    print("-----------------------------------Start----------------------------------\n")
    files_to_read = af.file_crawler()['measures']
    fo_values = np.sort(np.unique(np.array([float(af.fo_from_file(f)) for f in files_to_read])))
    a_values = np.sort(np.unique(np.array([af.a_from_file(f) for f in files_to_read])))
    h_values = np.sort(np.unique(np.array([af.h_from_file(f) for f in files_to_read])))
    #--------------------------------------------------------------------------------------
    # User input
    if var == None:
        print("Choose between 'B' or 'theta': ")
        PARAM['var'] = input()
    else:
        PARAM['var'] = var

    if a == None:
        print("Chose a value for 'a': ")
        for a in a_values:
            print(f'\t{a:.2f}')
        PARAM['a'] = input()
    else:
        PARAM['a'] = a

    if h == None:
        print("Chose a value for 'h': ")
        for h in h_values:
            print(f'\t{h}')
        PARAM['h'] = input()
    else:
        PARAM['h'] = h 
    #--------------------------------------------------------------------------------------
    # Selecting the files to read
    h = PARAM['h']
    a = PARAM['a']
    files_to_read = [f for f in files_to_read if f'h_{h}' in f]
    files_to_read = [f for f in files_to_read if f'a_{a}' in f]

    columns = len(fo_values)
    #--------------------------------------------------------------------------------------
    # Figure init
    fig = plt.figure(figsize=(7*columns,10))
    fig.suptitle('a = '+ str(a) + ' and h = ' + str(h), fontsize=25)
    plt.rcParams['font.size'] = '15'
    save_dir = str(os.path.abspath(__file__)) 
    save_dir = save_dir.replace('/'+save_dir.split('/')[-1], "")

    #--------------------------------------------------------------------------------------
    # Graphics generator loop
    for f in files_to_read:
        
        t, V, Vw, E, bxw, byw, rxw, ryw, txw, tyw, vbw, vrw, vpw, x_CM, y_CM, z_CM = np.loadtxt(f, unpack=True)
        
        if PARAM['var'] == 'theta':
            xTS = txw
            yTS = tyw
            x_limits = (80, 160)
            labeling_x = r"$\theta_x$"
            labeling_y = r"$\theta_y$"
        elif PARAM['var'] == 'B':
            xTS = bxw
            yTS = byw
            x_limits = (15, 70)
            labeling_x = r"$B_x$"
            labeling_y = r"$B_y$"
        elif PARAM['var'] == 'Vp':
            xTS = vbw
            yTS = []
            x_limits = (0, 0.2)
            labeling_x = r"$V_f$"
        else:
            print("Invalid 'var' parameter")
            return

        state = af.state_from_file(f)
        N = 1
        leg_loc = 'lower left'
        if state == 'WE':
            N = columns+1 
            leg_loc = 'upper left'
        
        # Find the list index of fo in fo_values
        fo = af.fo_from_file(f)
        fo_loc = np.where(fo_values == float(fo))[0][0]
        print(f"Plotting histogram a = {PARAM['a']}, f = {fo} and {state}")

        #================================HISTOGRAMS================================
        sample = len(t)//PARAM['samp_frac']
        plt.subplot(2,columns,N+fo_loc)

        xTS = xTS[sample:]
        if not yTS==[]:
            yTS = yTS[sample:]

        plt.yscale("log")
        xTS_freq = af.time_series_to_frequency_array(xTS, bin = PARAM['bin'])
        if not yTS==[]:
            yTS_freq = af.time_series_to_frequency_array(yTS, bin = PARAM['bin'])
        else:
            yTS_freq = np.array([[0, 0]])
        maximum = np.max(np.append(np.max(xTS_freq[:,1]), np.max(yTS_freq[:,1])))
        plt.bar(xTS_freq[:,0], xTS_freq[:,1]/maximum, color='b', width = 1/PARAM['bin'], align='center', label=labeling_x)
        if not yTS==[]:
            plt.bar(yTS_freq[:,0], yTS_freq[:,1]/maximum, color='g', width = 1/PARAM['bin'], align='center', label=labeling_y)

        plt.xlim(x_limits)
        
        plt.ylim(0, 1)
        plt.legend(loc="upper left")
        plt.title(PARAM['var'] + r' $\mu$='+ str(fo) +' state ' + state)
        if PARAM['var']=='Vp':
            means = np.array([np.mean(xTS)])
            stds = np.array([np.std(xTS)])
            with open(f"{save_dir}/{PARAM['var']}_means_h_{PARAM['h']}_a_{PARAM['a']}.txt", 'a') as out:
                out.write(f"{state}, {fo}, {means}, {stds}\n")

    #--------------------------------------------------------------------------------
    # Final adjustments
    fig.tight_layout()
    fig.savefig(f"{save_dir}/{PARAM['var']}_hist_h_{PARAM['h']}_a{PARAM['a']}.pdf", format='pdf')
    print("-----------------------------------End------------------------------------\n")

histogram_of(var='Vp', a=8, h=10)
