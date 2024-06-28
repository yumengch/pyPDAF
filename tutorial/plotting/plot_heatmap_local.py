
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw=None, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current Axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar

def rmse(filename1, filename2):
    field1 = np.loadtxt(filename1, delimiter=';')
    field2 = np.loadtxt(filename2)
    field1 = field1.reshape(18,36)
    field2 = field2.reshape(18,36)
    rmse = 0;
    for i in range(16):
       for j in range(36):
	       rmse = rmse + (field1[i,j]-field2[i,j])**2
    rmse = np.sqrt(1/(18*36)*(rmse));

    return rmse

if __name__ == "__main__":
    nsteps = ["2","4","6","8","10","12","14","16","18"]
    forgets = ["0.0","1.0","5.0","10.0","15.0","20.0","25.0","30.0","35.0","40.0"]

    locweight = 2

    rmse_arr = np.full((len(forgets)+1, len(nsteps)), np.nan)
    for i in range(len(forgets)):
        for j in range(len(nsteps)):
            rmse_arr[i+1,j] = rmse(f'out_ensB_N4_lw{locweight}_r{forgets[i]}/state_step{nsteps[j]}_ana.txt', f'inputs_online/true_step{nsteps[j]}.txt')
            #rmse_arr[i+1,j] = rmse(f'out_N4_lw{locweight}_r{forgets[i]}/state_step{nsteps[j]}_ana.txt', f'inputs_online/true_step{nsteps[j]}.txt')

#    minval = np.min(rmse_arr)
#    print ('minval ', minval)

    fig, ax1 = plt.subplots(1, 1, figsize=(6, 6),facecolor='.9')
    divnorm = colors.TwoSlopeNorm(vmin=0,vcenter=.5,vmax=1)
    # Replicate the above example with a different font size and colormap.
    coefficients = rmse_arr[1:].T
    im1, cbar1 = heatmap(coefficients, nsteps, forgets, norm=divnorm, ax=ax1,cmap='RdYlGn_r')
    ax1.set_xlabel('cut-off radius')
    ax1.set_ylabel('analysis time step')
    plt.gcf().set_size_inches(10, 4)
    plt.show()
