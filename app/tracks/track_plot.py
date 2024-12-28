from mpl_toolkits.mplot3d.art3d import Line3DCollection
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
from matplotlib.patches import Circle
import app.aux_func as af
#ax.clear()

PATH = str(os.path.abspath(__file__))
PATH = PATH.replace(PATH.split("/")[-1], "")

def Track_plotter_2D():
    state = input("Digite o estado desejado: (WE/CB)")
    l = 240
    plt.rcParams['font.size'] = '15'

    measures_dict = af.file_crawler()['measures']
    measures_dict = [f for f in measures_dict if f'{state}_' in f]
    for i in range(len(measures_dict)):
        print(f"{i}) h={af.h_from_file(measures_dict[i])} a={af.a_from_file(measures_dict[i])} fo={af.fo_from_file(measures_dict[i])}")

    file_index = input("Digite o número do arquivo que deseja visualizar: ")

    measure_file = measures_dict[int(file_index)]
    t, V, Vw, E, bxw, byw, rxw, ryw, txw, tyw, vbw, nulo1, nulo2, xm_CM, ym_CM, z_CM = np.loadtxt(measure_file, unpack=True)

    B_init = (bxw[0] + byw[0])/2
    Bx_stat = np.mean(bxw[int(len(bxw)//2):])
    By_stat = np.mean(byw[int(len(byw)//2):])
    B_stat = (Bx_stat + By_stat)/2

    file = measures_dict[int(file_index)]

    # Get the values of a and h from the file
    a = af.a_from_file(file)
    h = af.h_from_file(file)
    fo = af.fo_from_file(file)
    state = af.state_from_file(file)

    for i in range(1,len(xm_CM)):
        if abs(xm_CM[i] - xm_CM[i-1]) > l/2:
            delta = xm_CM[i] - xm_CM[i-1]
            xm_CM[i:] += -(delta/abs(delta))*l
        if abs(ym_CM[i] - ym_CM[i-1]) > l/2:
            delta = ym_CM[i] - ym_CM[i-1]
            ym_CM[i:] += -(delta/abs(delta))*l

    R = np.column_stack((t, xm_CM, ym_CM))
    # Setting the plots with gradient lines
    segments = np.stack([R[:-1,1:3], R[1:,1:3]], axis=1)

    cmap = plt.cm.get_cmap("viridis")
    norm = plt.Normalize(0, t[-1])

    fig = plt.figure()
    ax = fig.add_subplot(111)

    xmin = np.min(xm_CM)
    xmax = np.max(xm_CM)
    ymin = np.min(ym_CM)
    ymax = np.max(ym_CM)

    lc = mcoll.LineCollection(segments, linewidths=0.5, colors=cmap(norm(R[:-1,0])))
    ax.add_collection(lc)

    # Set specific ticks for the x-axis and y-axis
    ax.set_xticks([(xmin//l)*l, ((xmin + xmax)/2*l)*l, (xmax//l)*l])
    ax.set_yticks([(ymin//l)*l, ((ymin + ymax)/2*l)*l, (ymax//l)*l])

    # Optionally, set the tick labels if you want custom labels
    #ax.set_xticklabels([f"{(xmin//l)*l:.0f}", f"{((xmin + xmax) / 2*l)*l:.2f}L", f"{(xmax//l)*l:.0f}L"])
    #ax.set_yticklabels([f"{(xmin//l)*l:.0f}", f"{((ymin + ymax) / 2*l)*l:.2f}L", f"{(ymax//l)*l:.0f}L"])

    ax.set_xlabel("x")
    if (xmax) < l and (xmin) > 0:
        ax.set_xlim(0, l)
        ax.set_xticks([0,l/2, l])
        ax.set_xticklabels([0,"0.5L","L"])
        flag = False
    else:
        ax.set_xlim(xmin , xmax)
        ax.set_xticks([int((xmin//l)*l), ((xmin + xmax)//(2*l))*l, (xmax//l)*l])
        ax.set_xticklabels([f"{int((xmin//l))}L",f"{int(((xmin + xmax)/(2*l)))}L",f"{int((xmax//l))}L"])
        flag = True
    ax.set_ylabel("y")
    if (ymax) < l and (ymin) > 0:
        ax.set_ylim(0, l+1)
        ax.set_yticks([0,l/2, l])
        ax.set_yticklabels([0,"0.5L","L"])
    else:
        ax.set_ylim(ymin , ymax)
        ax.set_yticks([int((ymin//l)*l), ((ymin + ymax)//(2*l))*l, (ymax//l)*l])
        ax.set_yticklabels([f"{int((ymin//l))}L",f"{int(((ymin + ymax)/(2*l)))}L",f"{int((ymax//l))}L"])
    G_dict = {}
    G_dict['5'] = '5'
    G_dict['11'] = '11' 

    CI_dict = {}
    CI_dict['WE'] = 'W'
    CI_dict['CB'] = 'D'

    fig.colorbar(lc, label=r't', ax=ax)
    ax.set_title(f"a={G_dict[str(a)]} CI={CI_dict[state]} "+r"$\mu$="+f"{fo}")

    # Create a circle
    circle_init = Circle((xm_CM[0], ym_CM[0]), B_init, color=cmap(norm(0)), fill=True)
    circle = Circle((xm_CM[-1], ym_CM[-1]), B_stat, color=cmap(norm(t[-1])), fill=True)

    # Add the circle to the axes
    ax.add_patch(circle_init)
    ax.add_patch(circle)
    fig.tight_layout()

    plt.savefig(f"{PATH}2D_{state}_CM_a-{a}_h-{h}_fo-{fo}.jpeg")

if __name__ == '__main__':
    Track_plotter_2D()