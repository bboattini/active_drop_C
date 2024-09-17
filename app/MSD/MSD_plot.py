import numba
import numpy as np
import pwlf
import matplotlib.pyplot as plt
import app.aux_func as af
import os
from typing import Dict, List, Optional, Union

PATH = str(os.path.abspath(__file__))
PATH = PATH.replace(PATH.split("/")[-1], "")

CI_DICT = {}
CI_DICT['WE'] = 'WE'
CI_DICT['CB'] = 'CB'

G_DICT = {}
G_DICT['5'] = 'a=5'
G_DICT['11'] = 'a=11'
G_DICT['8'] = 'a=8' 

def MSD_plotter(state=None, a=None, h=None, limit=None):
    print("-----------------------------------Start----------------------------------\n")
    measures = af.file_crawler()['measures']
    print(f"Files: {len(measures)}\n")
    # Check if there is a "msd_files" folder in this directory in PATH
    if 'msd_files' not in os.listdir(PATH):
        MSD_file_writer(measures)

    # Select one of the msd_files
    if a==None or h==None or state==None:
        msd_files = os.listdir(f'{PATH}msd_files')
        print("Arquivos:")
        for i in range(len(msd_files)):
            if ".txt" in msd_files[i]:
                print("\t"+str(i)+") "+ msd_files[i])
        msd_file = f'{PATH}msd_files/' + msd_files[int(input("Digite o índice do arquivo desejado: "))]
        # get h value from msd file
        h = msd_file.split("/")[-1].split("_")[4]
        a = msd_file.split("/")[-1].split("_")[2]
        state = msd_file.split("/")[-1].split("_")[5]
        state = state.replace(".txt", "")
    else:
        msd_file = f'{PATH}msd_files/msd_a_{a}_h_{h}_{state}.txt'

    plt.rcParams.update({'font.size': 13})
    fo_values = np.sort(np.unique(np.array([float(af.fo_from_file(f)) for f in measures])))
    print("\nResultado:"+str(fo_values))
    print("\nOpções:")
    if limit == None:
        for i in range(len(fo_values)):
            print('\t'+fo_values[i])
        limit = input(r"Digite os valores de mu separados por espaço: ")
    # Transform the values into a list of floats
    limit = [float(i) for i in limit.split()]
    limit = np.sort(np.array(limit))
    print("\nlimit: "+str(limit))
    limit_index = np.array([i for i in range(len(fo_values)) if fo_values[i] in limit]) +1
    print("limit_index: "+str(limit_index))

    # fit = input("Digite fit se deseja fazer a regressão linear por partes:") 
    data = np.loadtxt(msd_file, unpack=True)
    # filter the limit_index's in data
    dados = {}
    dados = {'dt': data[0,:]}
    dt = dados['dt']
    data = data[limit_index,:]

    '''# dictionary that relates fo_values with color values
    color_dict = {}
    for i in range(len(fo_values)):
        color_dict[str(fo_values[i])] = COLORS[i]'''

    # create a dictionary with the folowing indexes [t] + [fo_values]
    for i in range(len(limit)):
        dados[limit[i]] = data[i,:]

    for i in range(1, len(dados.items())):
        msd = np.array(list(dados.values())[i])
        # Create a logarithmically spaced list of indices
        indices = np.logspace(0, np.log10(len(dt)-1), 20).astype(int)
        # Define the color for the current m
        fo = str(list(dados.keys())[i])
        # Plot MSD over time with a label containing the value of a and the defined color
        plt.plot(dt, msd, 'o-', label=r'$\mu$'+f'={fo}', color=af.COLORS[i-1], markevery = indices)

    plt.yscale("log")
    plt.xscale("log")
    plt.xlim(0, dt[-1]/10)
    gide = np.arange(1, 10, 1)
    plt.plot(gide, 0.12*gide*gide, 'k--', label='slope = 2')
    gide = np.arange(50, 500, 1)
    plt.plot(gide, 0.05*gide, 'k:', label='slope = 1')

    # Add labels and legend
    plt.title(f'CM MSD for {G_DICT[a]} and {state}')
    plt.xlabel(r'$\Delta$ t')
    plt.ylabel(r'MSD')
    #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', prop={'size': 5})
    plt.legend(loc='upper left', prop={'size': 7})
    # Create some space on the right side of the plot for the legend
    plt.subplots_adjust(right=0.15)

    # Display the plot
    plt.tight_layout()
    plt.savefig(msd_file.replace(msd_file.split("/")[-1], f"plot_msd_a_{a}_h_{h}_{state}.png"), dpi=300)
    # clear plt
    plt.clf()

def MSD_fixed_plotter(var, values, limit=None):
    print("-----------------------------------Start----------------------------------\n")
    measures = af.file_crawler()['measures']
    fo_values = np.sort(np.unique(np.array([float(af.fo_from_file(f)) for f in measures])))
    a_values = np.sort(np.unique(np.array([af.a_from_file(f) for f in measures])))
    h_values = np.sort(np.unique(np.array([af.h_from_file(f) for f in measures])))
    CI_values = ['WE', 'CB']
    print(f"Files: {len(measures)}\n")

    # Check if there is a "msd_files" folder in this directory in PATH
    if 'msd_files' not in os.listdir(PATH):
        MSD_file_writer(measures)
    
    values_dict = {'a': a_values, 'h': h_values, 'CI': ['WE', 'CB']}
    assert var in ['a', 'h', 'CI'], "Invalid 'var' parameter"
    #>>>>>>>>>>>>>
    '''for v in np.array(fixed['values']):
        assert str(v) in values_dict[fixed['var']] , f"Invalid 'value' parameter: {v}"'''
    #>>>>>>>>>>>>>
    if not isinstance(values, list):
        values = list(values)

    #list of all files in {PATH}msd_files/
    msd_files = os.listdir(f'{PATH}msd_files')
    msd_files = [f for f in msd_files if '.txt' in f]
    search_dict = {'a': -4, 'h': -2, 'CI': -1}
    curve = {}
    print(msd_files)
    if var == 'CI':
        msd_files = [f for f in msd_files if f.split("_")[search_dict[var]].replace(".txt","") in values]
        curve['var'] = 'a'
        curve['values'] = list(np.sort(np.append(a_values,np.array([5, 11]))))
    else:
        msd_files = [f for f in msd_files if int(f.split("_")[search_dict[var]]) in values]
        curve['var'] = 'CI'
        curve['values'] = ['WE', 'CB']
    print(msd_files)

    plt.rcParams.update({'font.size': 13})
    # Create a sub plot with one line and len(fixed['values']) columns
    fig, axs = plt.subplots(1, len(values), figsize=(len(values)*5, 4))

    print("\nResultado:"+str(fo_values))
    print("\nOpções:")
    if limit == None:
        for i in range(len(fo_values)):
            print('\t'+fo_values[i])
        limit = input(r"Digite os valores de mu separados por espaço: ")
    # Transform the values into a list of floats
    limit = [float(i) for i in limit.split()]
    limit = np.sort(np.array(limit))
    print("\nlimit: "+str(limit))
    limit_index = np.array([i for i in range(len(fo_values)) if fo_values[i] in limit]) +1
    print("limit_index: "+str(limit_index))
    label_register = []
    for msd_file in msd_files:
        val = msd_file.split("_")[search_dict[var]].replace('.txt','')
        if var != 'CI':
            val = float(val)
        ax = axs[values.index(val)]
        data = np.loadtxt(PATH +"msd_files/"+ msd_file, unpack=True)
        # filter the limit_index's in data
        dados = {}
        dados = {'dt': data[0,:]}
        dt = dados['dt']
        data = data[limit_index,:]

        # create a dictionary with the folowing indexes [t] + [fo_values]
        for i in range(len(limit)):
            dados[limit[i]] = data[i,:]

        for i in range(1, len(dados.items())):
            msd = np.array(list(dados.values())[i])
            # Create a logarithmically spaced list of indices
            indices = np.logspace(0, np.log10(len(dt)-1), 20).astype(int)
            # Define the color for the current m
            fo = str(list(dados.keys())[i])
            cur = msd_file.split("_")[search_dict[curve['var']]].replace('.txt','')
            cur = float(cur) if curve['var'] != 'CI' else cur
            # Plot MSD over time with a label containing the value of a and the defined color
            ax.plot(dt, msd, af.MARKERS[curve['values'].index(cur)]+"-", label=r'$\mu$'+f"={fo} {curve['var']}={cur}", color=af.COLORS[i-1], markevery = indices, alpha=1/(1+curve['values'].index(cur)))
        
        if values.index(val) not in label_register:
            ax.set_xlim(1, dt[-1]/10)
            ax.set_yscale("log")
            ax.set_xscale("log")
            gide = np.arange(1, 10, 1)
            ax.plot(gide, 0.12*gide*gide, 'k--', label='slope = 2')
            gide = np.arange(50, 500, 1)
            ax.plot(gide, 0.05*gide, 'k:', label='slope = 1')

            # Add labels and legend
            ax.set_title(f"CM MSD for fixed {var}={val}")
            ax.set_xlabel(r'$\Delta$ t')
            ax.set_ylabel(r'MSD')
            #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', prop={'size': 5})
            ax.legend(loc='upper left', prop={'size': 6})
        label_register.append(values.index(val))
    # Create some space on the right side of the plot for the legend
    fig.subplots_adjust(right=0.15)

    # Display the plot
    fig.tight_layout()
    fig.savefig(f"{PATH}fixed_msd_{var}.png", dpi=300)
    # clear plt
    plt.clf()

def MSD_file_writer(measures):
    # Create a directory to store the files
    os.makedirs(f'{PATH}msd_files', exist_ok=True)
    h_values = np.sort(np.unique(np.array([af.h_from_file(f) for f in measures])))
    a_values = np.sort(np.unique(np.array([af.a_from_file(f) for f in measures])))
    fo_values = np.sort(np.unique(np.array([float(af.fo_from_file(f)) for f in measures])))

    for state in ['WE','CB']:
        files_to_read = [f for f in measures if f'{state}_' in f]
        for h_val in h_values:
            for a_val in a_values:
                # Initialize an empty dict for the dt values
                dt_values = {}
                # Initialize an empty dictionary to store the msd values for each fo value
                msd_dict = {}
                for fo_val in fo_values:
                    m = [f for f in files_to_read if f'a_{a_val}' in f and f'h_{h_val}' in f and f'fo_{float(fo_val):.2f}.' in f][0]
                    t, V, Vw, E, bxw, byw, rxw, ryw, txw, tyw, vbw, vrw, vpw, xm_CM, ym_CM, zm_CM = np.loadtxt(m, unpack=True)   
                    a = a_val
                    fo = fo_val
                    l = 240
                    frac_data = 1/10
                    print(f'Calculating MSD for a = {a}, h = {h_val}, fo = {fo} and {state}')
                    msd, dt = MSD_calculator(xm_CM, ym_CM, t, l, frac_data)
                    # Store the dt values
                    dt_values[fo] = dt
                    # Store the msd values for the current fo value
                    msd_dict[fo] = msd
                    
                # Find the length of the shortest msd array
                min_length = min(len(msd) for msd in msd_dict.values())

                # Find the shortest array in dt_values dictionary
                dt_values = [dt for dt in dt_values.values() if len(dt) == min_length][0]

                # Truncate the msd arrays to make them the same length
                for fo in msd_dict:
                    msd_dict[fo] = msd_dict[fo][:min_length]

                # Create a 2D array with dt in the first column and the msd values in the other columns
                data = np.column_stack((dt_values, *[msd_dict[fo] for fo in fo_values]))

                # Save the data to a file
                np.savetxt(f'{PATH}msd_files/msd_a_{a_val}_h_{h_val}_{state}.txt', data)

@numba.njit
def MSD_calculator(xm_CM, ym_CM, t, l, frac_data):
    xm_CM = xm_CM.copy()
    ym_CM = ym_CM.copy()

    for i in range(1,len(xm_CM)):
        if abs(xm_CM[i] - xm_CM[i-1]) > l/2:
            delta = xm_CM[i] - xm_CM[i-1]
            xm_CM[i:] += -(delta/abs(delta))*l
        if abs(ym_CM[i] - ym_CM[i-1]) > l/2:
            delta = ym_CM[i] - ym_CM[i-1]
            ym_CM[i:] += -(delta/abs(delta))*l

    xm_CM1 = xm_CM[int(len(xm_CM)*frac_data):]
    ym_CM1 = ym_CM[int(len(ym_CM)*frac_data):]
    xm_CM2 = xm_CM[int(len(xm_CM)*frac_data):]
    ym_CM2 = ym_CM[int(len(ym_CM)*frac_data):]
    
    d_wind = t[1] - t[0]
    # Calculate the length of the arrays
    length = int(len(t)*(1-frac_data)-1)

    # Preallocate the arrays
    msd = np.zeros(length)
    dt = np.zeros(length)

    for phase in range(1, length):
        # Calculate the CM MSD per window
        dx = xm_CM1[phase:] - xm_CM2[:-phase]
        dy = ym_CM1[phase:] - ym_CM2[:-phase]
        dr2 = dx**2 + dy**2
        msd[phase] = np.mean(dr2)
        dt[phase] = phase*d_wind
        
    return msd, dt

if __name__ == '__main__':
    #MSD_plotter('WE', '8', '10', '0 0.5 1 2 5 10')
    #MSD_plotter('CB', '8', '10', '0 0.5 1 2 5 10')
    MSD_fixed_plotter('a',  [5, 8, 11], '0 0.5 2 5 10')
    MSD_fixed_plotter('CI', ["WE", "CB"], '0 0.5 2 5 10')



