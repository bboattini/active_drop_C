import numba
import numpy as np
import pwlf
import matplotlib.pyplot as plt
import app.aux_func as af
import os

PATH = str(os.path.abspath(__file__))
PATH = PATH.replace(PATH.split("/")[-1], "")

CI_DICT = {}
CI_DICT['WE'] = 'WE'
CI_DICT['CB'] = 'CB'

G_DICT = {}
G_DICT['5'] = 'a=5'
G_DICT['11'] = 'a=11'
G_DICT['8'] = 'a=8' 

def MSD_linear_fit(state=None, a=None, h=None, limit=None, break_n=3, lim_frac=0.1):
    print("-----------------------------------Start----------------------------------\n")
    measures = af.file_crawler()['measures']
    print(f"Files: {len(measures)}\n")
    # Check if there is a "msd_files" folder in this directory in PATH
    if 'msd_files' not in os.listdir(PATH):
        print("Run MSD_plot first")
        return

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
    fo_values = np.sort(np.unique(np.array([af.fo_from_file(f) for f in measures])))
    print("\nResultado:"+str(fo_values))
    print("\nOpções:")
    if limit == None:
        for i in range(len(fo_values)):
            print('\t'+fo_values[i])
        limit = input(r"Digite os valores de mu separados por espaço: ")
    # Transform the values into a list of floats
    limit = [i for i in limit.split()]
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
        # Index of the neerest value to dt[-1]*lim_frac
        lim_index = np.where(dt < dt[-1]*lim_frac)[0][-1]

        # Fit the MSD data with a piecewise linear function
        my_pwlf = pwlf.PiecewiseLinFit(np.log(dt[1:lim_index]), np.log(msd[1:lim_index]))
        breaks = my_pwlf.fit(break_n)
        slopes = my_pwlf.calc_slopes()
        dt_fit = np.linspace(dt[1], dt[-1], num=10000)
        msd_fit = my_pwlf.predict(np.log(dt_fit))
        print(f"fo={fo}, slopes={slopes}")

        # Plot MSD over time with a label containing the value of a and the defined color
        plt.plot(dt_fit, np.exp(msd_fit), '--', label=f'slope={slopes[0]:.2f}', color=af.COLORS[i-1])
        plt.plot(dt, msd, 'o', label=r'$\mu$'+f'={fo}', color=af.COLORS[i-1], markevery = indices)
        print(breaks)
        plt.axvline(x=np.exp(breaks[1]), color=af.COLORS[i-1], linestyle='-')

    plt.yscale("log")
    plt.xscale("log")
    plt.xlim(1, dt[-1]*lim_frac)
    gide = np.arange(1, 10, 1)
    plt.plot(gide, 0.12*gide*gide, 'k--', label='slope = 2')
    gide = np.arange(50, 500, 1)
    plt.plot(gide, 0.05*gide, 'k:', label='slope = 1')

    # Add labels and legend
    plt.title(f'CM MSD fit for {G_DICT[a]} and {state}')
    plt.xlabel(r'$\Delta$ t')
    plt.ylabel(r'MSD')
    #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', prop={'size': 5})
    plt.legend(loc='upper left', prop={'size': 7})
    # Create some space on the right side of the plot for the legend
    plt.subplots_adjust(right=0.15)

    # Display the plot
    plt.tight_layout()
    plt.savefig(f"{PATH}fit_msd_a_{a}_h_{h}_{state}.png", dpi=300)
    # clear plt
    plt.clf()
PARAM = {'vert': 'a', 
         'hori': 'h',
         'multi': 'CI',
         'var': 'theta'}
def MSD_diagram(var = PARAM['var'], config=None, limit='all', break_n=3, lim_frac=0.1):
    print("-----------------------------------Start----------------------------------\n")
    measures = af.file_crawler()['measures']
    fo_values = np.sort(np.unique(np.array([float(af.fo_from_file(f)) for f in measures])))
    print("\nResultado:"+str(fo_values))
    print("\nOpções:")
    if limit == None:
        for i in range(len(fo_values)):
            print('\t'+fo_values[i])
        limit = input(r"Digite os valores de mu separados por espaço: ")
    if limit == 'all':
        limit = ''
        for i in fo_values:
            limit += str(i)+' '
    # Transform the values into a list of floats
    limit = [i for i in limit.split()]
    limit = np.sort(np.array(limit))
    print("\nlimit: "+str(limit))
    limit_index = np.array([i for i in range(len(fo_values)) if fo_values[i] in limit]) +1
    print("limit_index: "+str(limit_index))
    
    # Check if there is a "msd_files" folder in this directory in PATH
    if 'msd_files' not in os.listdir(PATH):
        print("Run MSD_plot first")
        return
    # --------------------------------------------------------------------------------------
    # User input
    if config == None:
        print(
'''Insert the variables(CI, a or h) to be used in vertical and horizontal axis, respectively
(separate by a space):''')
        # Convert the input into a tuple
        config = tuple(input().split())
    PARAM['vert'] = config[0]
    PARAM['hori'] = config[1]
    PARAM['multi'] = [l for l in ['a', 'h', 'CI'] if l not in [PARAM['vert'], PARAM['hori']]][0]

    breaks = []
    slopes = []
    dr = []
    var_dict = {'p_length': dr, 'p_time': breaks, 'slope': slopes}

    if var == 'p_length':
        y_limits = (0, 5000)
        labeling_y = r'$\Delta r_{p}$'
        leg_loc = "upper left"
    elif var == 'p_time':
        y_limits = (0, 3000)
        labeling_y = r"$\Delta t_{p}$"
        leg_loc = "upper left"
    elif var == 'slope':
        y_limits = (0, 2)
        labeling_y = r"$1^{st}$ slope"
        leg_loc = "upper right"
    else:
        print("Invalid 'var' parameter")
        return

    #list of all files in {PATH}msd_files/
    msd_files = os.listdir(f'{PATH}msd_files')
    msd_files = [f for f in msd_files if '.txt' in f]
    search_dict = {'a': -4, 'h': -2, 'CI': -1}
    a_values = np.sort(np.unique(np.array([float(f.split("/")[-1].split("_")[search_dict['a']]) for f in msd_files])))
    h_values = np.sort(np.unique(np.array([float(f.split("/")[-1].split("_")[search_dict['h']]) for f in msd_files])))
    CI_values = ['WE', 'CB']
    print(msd_files)

    config_dict = {'a': a_values, 'h': h_values, 'fo': fo_values, 'CI': np.array(CI_values)}
    vertical = config_dict[PARAM['vert']]
    horizontal = config_dict[PARAM['hori']]
    curves = config_dict[PARAM['multi']]

    #--------------------------------------------------------------------------------------
    # Create a subplot
    plt.rcParams.update({'font.size': 15})
    fig, axs = plt.subplots(len(vertical), len(horizontal), figsize=(int(len(horizontal)*5), int(len(vertical)*4)))
    legend_added = {}

    for msd_file in msd_files:
        vert_index = msd_file.split("/")[-1].split("_")[search_dict[PARAM['vert']]].replace(".txt", "")
        try:
            vert_index = float(vert_index)
        except:
            vert_index = vert_index
        hori_index = msd_file.split("/")[-1].split("_")[search_dict[PARAM['hori']]].replace(".txt", "")
        try:
            hori_index = float(hori_index)
        except:
            hori_index = hori_index
        curve_index = msd_file.split("/")[-1].split("_")[search_dict[PARAM['multi']]].replace(".txt", "")
        try:
            curve_index = float(curve_index)
        except:
            curve_index = curve_index

        # Find the corresponding subplot
        Row = vertical.tolist().index(vert_index)
        Col = horizontal.tolist().index(hori_index)
        cur = curves.tolist().index(curve_index)

        if len(vertical) == 1 and len(horizontal) == 1:
            ax = axs
        elif isinstance(axs, np.ndarray):
            if axs.ndim > 1:
                ax = axs[Row, Col]
            else:
                if len(vertical) == 1: # or axs[Col/Row], depending on which dimension is 1
                    ax = axs[Col]
                else:
                    ax = axs[Row]
        if ax not in legend_added:
            legend_added[ax] = []
            if var in ['p_length', 'p_time']:
                ax.set_yscale("log")
                ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xscale("log")
            ax.set_ylim(y_limits)
            ax.set_xlim(fo_values[1], fo_values[-2]+0.5)
            # Add labels and legend
            ax.set_title(f'{PARAM["vert"]}={vert_index} and {PARAM["hori"]}={hori_index}')
            ax.set_xlabel(r'$\mu$')
            ax.set_ylabel(labeling_y)

        # fit = input("Digite fit se deseja fazer a regressão linear por partes:") 
        data = np.loadtxt(PATH+'msd_files/'+msd_file, unpack=True)
        # filter the limit_index's in data
        dados = {}
        dados = {'dt': data[0,:]}
        dt = dados['dt']
        data = data[1:,:]
        
        # create a dictionary with the folowing indexes [t] + [fo_values]
        for i in range(len(limit)):
            dados[limit[i]] = data[i,:]

        var_dict = {}
        for i in range(1, len(dados.items())):
            msd = np.array(list(dados.values())[i])
            # Create a logarithmically spaced list of indices
            indices = np.logspace(0, np.log10(len(dt)-1), 20).astype(int)
            # Define the color for the current m
            fo = str(list(dados.keys())[i])
            # Index of the neerest value to dt[-1]*lim_frac
            lim_index = np.where(dt < dt[-1]*lim_frac)[0][-1]

            # Fit the MSD data with a piecewise linear function
            my_pwlf = pwlf.PiecewiseLinFit(np.log(dt[1:lim_index]), np.log(msd[1:lim_index]))
            if float(fo) > 7:
                breaks = my_pwlf.fit(2)
            else:
                breaks = my_pwlf.fit(break_n)
            slopes = my_pwlf.calc_slopes()
            # Find the coresponding msd point for the break[0] point
            break_index = np.where(dt < np.exp(breaks[1]))[0][-1]
            dr = msd[break_index]
            var_dict['p_length'] = np.sqrt(dr)
            var_dict['p_time'] = np.exp(breaks[1])
            var_dict['slope'] = slopes[0]

            print(f"fo={fo}, slopes={slopes}")
            if f"{af.COLORS[cur]}" not in legend_added[ax]:
                ax.plot(float(fo), var_dict[var], af.MARKERS[cur], label=f"{PARAM['multi']}={curve_index}", color=af.COLORS[cur])
                legend_added[ax].append(f"{af.COLORS[cur]}")
            else:
                ax.plot(float(fo), var_dict[var], af.MARKERS[cur], color=af.COLORS[cur])

    for ax in legend_added:
        ax.legend(loc=leg_loc, prop={'size': 7})
    
    # Create some space on the right side of the plot for the legend
    plt.subplots_adjust(right=0.15)
    # Display the plot
    plt.tight_layout()
    plt.savefig(f"{PATH}fit_{var}_diagram_{PARAM['vert']}-vs-{PARAM['hori']}.png", dpi=300)
    # clear plt
    plt.clf()


if __name__ == "__main__":
    #MSD_linear_fit('WE', '5', '10', '1 3 5 8 10')
    #MSD_linear_fit('CB', '5', '10', '1 3 5 8 10')
    #MSD_linear_fit('WE', '11', '10', '1 3 5 8 10')
    #MSD_linear_fit('CB', '11', '10', '1 3 5 8 10')
    MSD_diagram('p_length', ('h', 'a'))
    MSD_diagram('p_time', ('h', 'a'))
    MSD_diagram('slope', ('h', 'a'))
    MSD_diagram('p_length', ('h', 'CI'))
    MSD_diagram('p_time', ('h', 'CI'))
    MSD_diagram('slope', ('h', 'CI'))