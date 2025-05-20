import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    '''for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            if float(valfmt(data[i, j]))>0.0000001:
                text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                print(text)
                texts.append(text)'''
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            # Get the raw numeric value first
            current_value = data[i, j]
        
            #   Check if it's a valid number and meets your threshold
            valid_number = (
                current_value is not None and
                not (isinstance(current_value, float) and np.isnan(current_value)) and
                isinstance(current_value, (int, float)) and
                current_value > 0.0000001
                )
        
            if valid_number:
                # Now it's safe to update colors and format the value
                kw.update(color=textcolors[int(im.norm(current_value) > threshold)])
                formatted_value = valfmt(current_value, None)
                
                # Add the text to the plot
                text = im.axes.text(j, i, formatted_value, **kw)
                print(text)
                texts.append(text)

    return texts

def plot_micelle_text(path,data,ion,resname,name,x_max,y_max):
    """
    

    Parameters
    ----------
    path : str
        the destiantion path to save the file
    data : array
        numpy array of the data to plot
    resname : str
        residue name of the ligand
    name : str
        name of the file to be saved
    x_max : int
        the maximum number of ions in the aggregate
    y_max : int
        the maximum number of ligands in the aggregate

    Returns
    -------
    None.

    """
    
    #custom colormap because I want to start with white color - https://jdherman.github.io/colormap/
    n=0
    
    distribution_AB=data
    distribution_AB_mod=np.where(np.array(distribution_AB) < 0.0005, np.nan, distribution_AB)
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams.update({'font.size': 25})
    plt.rc('legend',fontsize=22)
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)
    
    #print(len(distribution_AB))
    fig, ax1=plt.subplots(figsize=(15, 8))
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    majorLocator_y1 = MultipleLocator(1)
    minorLocator_y1 = MultipleLocator(0.5)
    #create custom colormap
    
    #cmap=[(0, 'white'), (0.5, 'green'), (1, 'darkgreen')]
    cmap=[(0, 'blue'),  (0.1, 'green'), (0.3, 'gold'), (0.6, 'red'),(0.8,'magenta'),(1,'darkviolet')]
    custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', cmap)
    
    majorLocator_y1 = MultipleLocator(5)
    minorLocator_y1 = MultipleLocator(2)
    
    im=ax1.imshow(distribution_AB_mod[0], cmap=custom_cmap, origin='lower',vmin=0, vmax=25, aspect='auto')
    #cbar=fig.colorbar(im,aspect=5, pad=0.01, fraction=0.15 * 2)
    cbar  = fig.colorbar(im, cax=cax)
    cbar.set_label(r'$N$', fontsize=22, labelpad = 10)
    #cbar.set_label(r'$N_{\mathrm{cluster}}$', fontsize=15, labelpad = 10)
    cbar.ax.minorticks_on()
    cbar.ax.yaxis.set_major_locator(majorLocator_y1)
    cbar.ax.yaxis.set_minor_locator(minorLocator_y1)
    
    cbar.ax.tick_params(which='both', width=1)
    cbar.ax.tick_params(which='major', length=10)
    cbar.ax.tick_params(which='minor', length=5)
    
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(22)
    ax1.set_xlabel('N {}'.format(ion), fontsize=25)
    ax1.set_ylabel('N {}'.format(resname),fontsize=25)
    ax1.set_xticks(np.arange(0,int(x_max),1))
    ax1.set_yticks(np.arange(0,int(y_max),1))
    ax1.tick_params(axis='x', labelsize=27)
    ax1.tick_params(axis='y', labelsize=27)
    text=annotate_heatmap(im, valfmt="{x:.3f}",threshold=25,fontsize=15)
    
    majorLocator_x = MultipleLocator(1)
    #minorLocator_x = MultipleLocator(0.5)

    majorLocator_y2 = MultipleLocator(2)
    minorLocator_y2 = MultipleLocator(1)
    
    
    
    ax1.xaxis.set_major_locator(majorLocator_x)
    #ax1.xaxis.set_minor_locator(minorLocator_x)

    ax1.yaxis.set_major_locator(majorLocator_y2)
    ax1.yaxis.set_minor_locator(minorLocator_y2)

    ax1.tick_params(which='both', width=2)
    ax1.tick_params(which='major', length=14)
    ax1.tick_params(which='minor', length=7)
    ax1.tick_params(axis = 'y')

    ax1.spines["bottom"].set_linewidth(2)
    ax1.spines["left"].set_linewidth(2)
    ax1.spines["top"].set_linewidth(2)
    ax1.spines["right"].set_linewidth(2)
    
    ax1.spines["bottom"].set_linewidth(2)
    ax1.spines["left"].set_linewidth(2)
    ax1.spines["top"].set_linewidth(2)
    ax1.spines["right"].set_linewidth(2)
    
    ax1.set_xlim(-1,x_max)
    ax1.set_ylim(0,y_max)
    filename=name+resname+'_text_parallel.pdf'
    plt.savefig(path / filename,bbox_inches='tight')
    #plt.show()
    n+=1
