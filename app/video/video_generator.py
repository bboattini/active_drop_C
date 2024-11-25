import numpy as np
import matplotlib.pyplot as plt
from itertools import islice
import os
import app.aux_func as af
import imageio.v2 as imageio


PATH = str(os.path.abspath(__file__))
PATH = PATH.replace(PATH.split("/")[-1], "")

def video_maker():
    # Define the colors
    light_blue = [173, 216, 230]  # RGB values for light blue
    light_grey = [211, 211, 211]  # RGB values for light grey
    dark_grey = [169, 169, 169]  # RGB values for dark grey
    red = [255, 0, 0]  # RGB values for red

    # Define the size of the lattice
    l = 240  # adjust this value as needed
    l2 = l*l
    w = 5
    h_base = 3

    state = input("Digite o estado desejado: (WE/CB)")
    fo_values = input("Digite os valores de fo desejado: ")
    # Transform the values into a list of floats
    fo_values = [i for i in fo_values.split()]

    config_dict = af.file_crawler()['base']
    # filter all the config files that has the selected state
    files_to_read = [f for f in config_dict if f'{state}_' in f]
    # filter all the config files that has the selected fo in the list fo_values
    files_to_read = [f for f in files_to_read if any(f'fo_{fo}' in f for fo in fo_values)]
    config_dict = files_to_read

    measures_dict = af.file_crawler()['measures']
    files_to_read = [f for f in measures_dict if f'{state}_' in f]
    files_to_read = [f for f in files_to_read if any(f'fo_{fo}' in f for fo in fo_values)]
    measure_file = files_to_read

    # Create a subplot with 2 rows and 6 columns
    fig, axs = plt.subplots(1, 3, figsize=(48, 16), dpi = 20)
    plt.rcParams.update({'font.size': 30})

    axs[0].set_title(f'Base')
    axs[1].set_title(f'Fixed y')
    axs[2].set_title(f'Fixed x')

    '''for i in range(len(config_dict)):
        print(f"{i}) h={h_from_file(config_dict[i])} a={a_from_file(config_dict[i])} fo={fo_from_file(config_dict[i])}")

    file_index = input("Digite o número do arquivo que deseja visualizar: ")'''
    for f_index in range(len(config_dict)):
        t, V, Vw, E, bxw, byw, rxw, ryw, txw, tyw, vbw, vrw, vpw, xm_CM, ym_CM, zm_CM = np.loadtxt(measure_file[f_index], unpack=True)
        d_wind = t[1] - t[0]

        file = config_dict[f_index]

        # Get the values of a and h from the file
        a = af.a_from_file(file)
        h = af.h_from_file(file)
        fo = af.fo_from_file(file)
        state = af.state_from_file(file)
        plt.suptitle(f"h={h} a={a} "+r"$\mu$"+f"={fo} {state}")

        #Create movie folder inside the current working directory
        movie_path = f"{PATH}movie100{state}_h{h}_a{a}_fo{fo}"
        if not os.path.exists(movie_path):
            os.makedirs(movie_path)

        # Create an empty image
        base = np.zeros((l, l, 3), dtype=np.uint8)
        y_fix = np.zeros((l, l, 3), dtype=np.uint8)
        x_fix = np.zeros((l, l, 3), dtype=np.uint8)

        y = ym_CM[0]
        # Iterate over the range
        for z in range(h + h_base):
            for x in range(l):
                site = z*l2 + y*l + x
                #if z < h_base or (x % (w + a) < w and y % (w + a) < w):
                if z < h_base or (x % (w + a) < w):
                    y_fix[x, z] = light_grey  # add a gray pixel
                    x_fix[x, z] = light_grey

        # Read the file using numpy
        data = np.loadtxt(file, skiprows=1, dtype=int)
        interval = 1
        MAX = 1000 #should be the mcs_MAX/d_wind
        start_line = 20  # replace with the line number you want to start from
        time = -1
        mcs = 0
        with open(file, 'r') as data:
            for row in islice(data, start_line, None):
                if row.startswith("# tempo"):
                    if time % interval == 0 and time > -1:
                        # Iterate over the range
                        for z in range(h + h_base):
                            y = ym_CM[mcs]
                            for x in range(l):
                                site = z*l2 + y*l + x
                                if z < h_base or (x % (w + a) < w and y % (w + a) < w):
                                #if z < h_base or (x % (w + a) < w):
                                    y_fix[x, z] = light_grey  # add a gray pixel
                            x = xm_CM[mcs]
                            for y in range(l):
                                site = z*l2 + y*l + x
                                if z < h_base or (x % (w + a) < w and y % (w + a) < w):
                                #if z < h_base or (y % (w + a) < w):
                                    x_fix[y, z] = light_grey    
                        for x in range(l):
                            for y in range(l):
                                if (x % (w + a) < w and y % (w + a) < w):
                                    base[x, y] = dark_grey

                        # Rotate the image 90° anti-clockwise
                        rotated_image = np.rot90(base, k=1, axes=(0, 1))
                        # Display the rotated image in the corresponding subplot
                        axs[0].imshow(rotated_image)
                        axs[0].set_title(f'Base mcs = {row.split()[2]}')

                        # Rotate the image 90° anti-clockwise
                        rotated_image = np.rot90(y_fix, k=1, axes=(0, 1))
                        # Display the rotated image in the corresponding subplot
                        axs[1].imshow(rotated_image)
                        axs[1].set_title(f'Fixed y CM={int(ym_CM[mcs])}')

                        # Rotate the image 90° anti-clockwise
                        rotated_image = np.rot90(x_fix, k=1, axes=(0, 1))
                        # Display the rotated image in the corresponding subplot
                        axs[2].imshow(rotated_image)
                        axs[2].set_title(f'Fixed x CM={int(xm_CM[mcs])}')

                        # Save the frame as a PNG file inside move folder
                        fig.savefig(f"{movie_path}/{time}_frame.png", format='png')

                        # Create an empty image
                        base = np.zeros((l, l, 3), dtype=np.uint8)
                        y_fix = np.zeros((l, l, 3), dtype=np.uint8)
                        x_fix = np.zeros((l, l, 3), dtype=np.uint8)
                        mcs = int(int(row.split()[2])/d_wind)
                    time += 1

                elif row.strip() and mcs%interval == 0:
                    # Iterate over the data
                    site, spin = map(int, row.split())
                    x = site % l
                    y = (site // l) % l
                    z = site // (l * l)
                    if x == int(xm_CM[mcs]):  # adjust this value as needed
                        base[x, y] = red
                        if spin == 1:
                            x_fix[y, z] = light_blue
                    if y == int(ym_CM[mcs]):  # adjust this value as needed
                        base[x, y] = red
                        if spin == 1:
                            y_fix[x, z] = light_blue
                    if z == h + h_base + 3 :  # adjust this value as needed
                        if spin == 1:
                            base[x, y] = light_blue

                if time == MAX:
                    break

def video_render(movie_path=None):
    print("----------------- Start -----------------\n")
    print(f'Reading movie frames of: {movie_path.split("/")[-1]}')
    fileList = []
    movie_paths = [f for f in os.listdir(PATH) if 'movie' in f]
    if movie_path is None:
        for i, movie in enumerate(movie_paths):
            print(f"\t{i}) {movie.replace(f'{PATH}', '')}")
        movie_path = movie_paths[int(input("Digite o número do arquivo que deseja visualizar: "))]
    #sorting the frames in the order of frame number
    listdir = os.listdir(movie_path)
    #remove ellipses.txt from the list
    listdir = [f for f in listdir if 'ellipses.txt' not in f]
    listdir.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

    for file in listdir:
        if "frame" in file:
            complete_path = movie_path +"/"+ file
            fileList.append(complete_path)
    movie = movie_path.split("/")[-1]
    state = movie_path.split("/")[-1].split("_")[0].replace("movie", "")
    writer = imageio.get_writer(f"{PATH}{state}{movie}.mp4", fps=20, codec="libx264")

    for frame in fileList:
        #print(f"Reading frame{frame.split('/')[-1]}")
        img = imageio.imread(frame, pilmode='RGB')
        writer.append_data(img)
    writer.close()
    print(f"\n----------------- End -----------------\n")

if __name__ == '__main__':
    # DEACTIVATE CONDA TO RENDER THE VIDEO
    #video_render('/home/bernardo/workspace/ACTIVE_WETTING/greene_videos_janela10/app/video/movieCB_h10_a5_fo1')
    #video_render('/home/bernardo/workspace/ACTIVE_WETTING/greene_videos_janela10/app/video/movieCB_h10_a11_fo1')
    #video_render('/home/bernardo/workspace/ACTIVE_WETTING/greene_videos_janela10/app/video/movieCB_h10_a5_fo10')
    #video_render('/home/bernardo/workspace/ACTIVE_WETTING/greene_videos_janela10/app/video/movieCB_h10_a11_fo10')
    video_maker()
