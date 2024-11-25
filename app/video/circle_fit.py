import numpy as np
import matplotlib.pyplot as plt
from itertools import islice
import os
import app.aux_func as af
import imageio.v2 as imageio
from scipy.optimize import least_squares

PATH = str(os.path.abspath(__file__))
PATH = PATH.replace(PATH.split("/")[-1], "")

COLORS = {
    'light_blue': [173, 216, 230],
    'dark_grey': [169, 169, 169]
}

PARAMS = {
    'l': 240,
    'w': 5,
    'h_base': 3,
    'samp_frac': 10,
    'bin': 1000
}

MAX = 1000 #should be the mcs_MAX/d_wind

def video_maker(state, fo_values, interval=1, start_line=20):
    # Define the colors
    wather_color = COLORS['light_blue']  # RGB values for light bl
    pillar_color = COLORS['dark_grey']  # RGB values for dark grey

    # Define the size of the lattice
    l = PARAMS['l']  # adjust this value as needed
    l2 = PARAMS['l']*PARAMS['l']
    w = PARAMS['w']
    h_base = PARAMS['h_base']

    # Transform the values into a list of floats
    fo_values = [i for i in fo_values.split()]

    config_dict = af.file_crawler()['config']
    # filter all the config files that has the selected state
    files_to_read = [f for f in config_dict if f'{state}_' in f]
    # filter all the config files that has the selected fo in the list fo_values
    files_to_read = [f for f in files_to_read if any(f'fo_{fo}.' in f for fo in fo_values)]
    config_dict = files_to_read

    measures_dict = af.file_crawler()['measures']
    files_to_read = [f for f in measures_dict if f'{state}_' in f]
    files_to_read = [f for f in files_to_read if any(f'fo_{fo}.' in f for fo in fo_values)]
    measure_file = files_to_read

    # Create a plot with only one frame
    # fig, axs = plt.subplots(figsize=(16, 16), dpi = 20)
    fig, axs = plt.subplots(figsize=(16, 16))
    plt.rcParams.update({'font.size': 30})

    axs.set_title(f'Base')

    for f_index in range(len(config_dict)):
        t, V, Vw, E, bxw, byw, rxw, ryw, txw, tyw, vbw, vrw, vpw, xm_CM, ym_CM, zm_CM = np.loadtxt(measure_file[f_index], unpack=True)
        d_wind = t[1] - t[0]

        file = config_dict[f_index]

        # Get the values of a and h from the file
        a = af.a_from_file(file)
        h = af.h_from_file(file)
        fo = af.fo_from_file(file)
        state = af.state_from_file(file)
        plt.suptitle(f"h={h} a={a} " + r"$\mu$" + f"={fo} {state}")
        print(f"for h={h} a={a} " + r"$\mu$" + f"={fo} {state}")

        #Create movie folder inside the current working directory
        movie_path = f"{PATH}movie{state}_h{h}_a{a}_fo{fo}"
        if not os.path.exists(movie_path):
            os.makedirs(movie_path)

        # CReate a empty ellipse file
        with open(f"{movie_path}/ellipses.txt", "w") as ellipse_file:
            ellipse_file.write("")

        # Create an empty image
        base = 255*np.ones((l, l, 3), dtype=np.uint8)
        X_base = np.array([])
        Y_base = np.array([])
        X_CM = np.array([])
        Y_CM = np.array([])

        # Read the file using numpy
        data = np.loadtxt(file, skiprows=1, dtype=int)
        time = -1
        mcs = 0
        with open(file, 'r') as data:
            for row in islice(data, start_line, None):
                if row.startswith("# tempo"):
                    if time % interval == 0 and time > -1: 
                        for x in range(l):
                            for y in range(l):
                                if (x % (w + a) < w and y % (w + a) < w):
                                    base[x, y] = pillar_color

                        # Rotate the image 270 degrees
                        rotated = np.rot90(base, k=1)
                        # Invert the y axis
                        rotated = np.flipud(rotated)                        
                        # Display the rotated image in the corresponding subplot
                        axs.imshow(rotated)
                        axs.set_title(f'Base mcs = {row.split()[2]}')
                        axs.plot(X_base, Y_base, '.')
                        
                        # Fit the ellipse
                        X, Y, Z, x0, y0, a1, a2 = ellipse_fit(X_base, Y_base)
                        axs.contour(X, Y, Z, levels=[1], colors='r')
                        axs.set_xlim(0, l)
                        axs.set_ylim(0, l)
                        
                        # Write in an ellipse file inside the movie folder
                        with open(f"{movie_path}/ellipses.txt", "a") as ellipse_file:
                            ellipse_file.write(f"{time*d_wind} {x0} {y0} {a1} {a2} {a2/a1}\n")
                    
                        # Save the frame as a PNG file inside move folder
                        fig.savefig(f"{movie_path}/{time}_frame.png", format='png')
                        plt.cla()

                        # Create an empty image
                        base = 255*np.ones((l, l, 3), dtype=np.uint8)
                        X_base = np.array([])
                        Y_base = np.array([])
                        X_CM = np.array([])
                        Y_CM = np.array([])

                        mcs = int(int(row.split()[2])/d_wind)
                    time += 1

                elif row.strip() and mcs%interval == 0:
                    # Iterate over the data
                    site, spin = map(int, row.split())
                    x = site % l
                    y = (site // l) % l
                    z = site // (l * l)

                    if z == int(zm_CM[mcs]):  # adjust this value as needed
                        if spin == 1:
                            X_CM = np.append(X_CM, x)
                            Y_CM = np.append(Y_CM, y)
                    if z == h + h_base + 3:  # adjust this value as needed
                        if spin == 1:
                            base[x, y] = wather_color
                            X_base = np.append(X_base, x)
                            Y_base = np.append(Y_base, y)

                if time == MAX:
                    break

def ellipse_fit(x, y, l=240):
    # Correct periodic boundary conditions
    if np.abs(np.min(x) - np.max(x)) > l/2:
        x = np.where(x < l/2, x + l, x)
    if np.abs(np.min(y) - np.max(y)) > l/2:
        y = np.where(y < l/2, y + l, y)

    A = np.stack([x**2, x * y, y**2, x, y]).T
    b = np.ones_like(x)
    w = np.linalg.lstsq(A, b)[0].squeeze()

    xlin = np.linspace(np.min(x), np.max(x), 300)
    ylin = np.linspace(np.min(y), np.max(y), 300)
    X, Y = np.meshgrid(xlin, ylin)

    Z = w[0]*X**2 + w[1]*X*Y + w[2]*Y**2 + w[3]*X + w[4]*Y

    # Ellipse center
    x0 = (w[1]*w[3] - 2*w[2]*w[4])/(4*w[0]*w[2] - w[1]**2)
    y0 = (w[1]*w[4] - 2*w[0]*w[3])/(4*w[0]*w[2] - w[1]**2)
    
    # Calculate the terms inside the square roots
    term1 = 2 * (w[0]*x0**2 + w[2]*y0**2 + w[3]*x0 + w[4]*y0 - 1)
    term2 = (w[0] + w[2]) + np.sqrt(np.abs((w[0] - w[2])**2 + w[1]**2))
    term3 = (w[0] + w[2]) - np.sqrt(np.abs((w[0] - w[2])**2 + w[1]**2))

    # Ensure terms inside sqrt are non-negative
    term1 = np.abs(term1)
    term2 = np.abs(term2)
    term3 = np.abs(term3)

    # Ellipse semi-axes
    a1 = np.sqrt(term1 / term2)
    a2 = np.sqrt(term1 / term3)

    # Determine the semi-major and semi-minor axes
    if a1 > a2:
        semi_major_axis = a1
        semi_minor_axis = a2
    else:
        semi_major_axis = a2
        semi_minor_axis = a1

    return X, Y, Z, x0, y0, semi_major_axis, semi_minor_axis

def ellipse_statistics(var, a, h=10):
    # List all the folders inside the current working directory
    folders = [f for f in os.listdir(PATH) if f.startswith('movie')]
    # Find each a values in each folder name
    a_values = np.unique(np.array([int(folder.split('_')[2].replace("a","")) for folder in folders]))
    fo_values = np.unique(np.array([float(folder.split('_')[3].replace("fo","")) for folder in folders]))
    state_values = np.unique(np.array([folder.split('_')[0].replace("movie","") for folder in folders]))
    if int(a) in a_values:
        folders = [f for f in folders if f'a{a}' in f]

    columns = len(fo_values)
    # Figure init
    fig = plt.figure(figsize=(7*columns,5))
    fig.suptitle(f'a = {a} and h = {h}', fontsize=25)
    plt.rcParams['font.size'] = '15'

    # Iterate over the folders
    for folder in folders:
        t, x0, y0, a1, a2, e = np.loadtxt(f"{PATH}{folder}/ellipses.txt", unpack=True)

        vard_dict = {'axis': a1, 'e': e}
        xTS = vard_dict[var]
        yTS = []
        labeling = {'axis': 'major axis', 'e': 'eccentricity'}
        x_limits = {'axis': (np.min(a1), np.max(a2)), 'e': (0, 1)}

        state = folder.split('_')[0].replace("movie","")
        N = 1
        leg_loc = 'lower left'
        if state == 'WE':
            N = columns+1 
            leg_loc = 'upper left'

        # Find the list index of fo in fo_values
        fo = float(folder.split('_')[3].replace("fo",""))
        fo_loc = np.where(fo_values == float(fo))[0][0]
        print(f"Plotting histogram a = {a}, f = {fo} and {state}")

        #================================HISTOGRAMS================================
        sample = len(t)//PARAMS['samp_frac']
        plt.subplot(2,columns,N+fo_loc)

        xTS = xTS[sample:]
        if var == 'axis':
            yTS = yTS[sample:]

        plt.yscale("log")
        xTS_freq = af.time_series_to_frequency_array(xTS, bin = PARAMS['bin'])
        if not yTS==[]:
            yTS_freq = af.time_series_to_frequency_array(yTS, bin = PARAMS['bin'])
            maximum = np.max(np.append(np.max(xTS_freq[:,1]), np.max(yTS_freq[:,1])))
        else:
            maximum = np.max(xTS_freq[:,1])

        plt.bar(xTS_freq[:,0], xTS_freq[:,1]/maximum, color='b', width = 1/PARAMS['bin'], align='center', label=labeling[var])
        if not yTS==[]:
            plt.bar(yTS_freq[:,0], yTS_freq[:,1]/maximum, color='g', width = 1/PARAMS['bin'], align='center', label="minor axis")

        plt.xlim(x_limits[var])
        
        plt.ylim(0, 1)
        plt.legend(loc="upper left")
        plt.title(var + r' $\mu$='+ str(fo) +' state ' + state)

    # Final adjustments
    fig.tight_layout()
    fig.savefig(f"{PATH}{var}_hist_h_{h}_a{a}.pdf", format='pdf')
    print("-----------------------------------End------------------------------------\n")

if __name__ == '__main__':
    #video_maker('WE', '1', interval=1, start_line=20)
    #video_maker('WE', '2', interval=1, start_line=20)
    #video_maker('WE', '5', interval=1, start_line=20)
    #video_maker('WE', '10', interval=1, start_line=20)
    ellipse_statistics('e', 5)
    ellipse_statistics('e', 11)
