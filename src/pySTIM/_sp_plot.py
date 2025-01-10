import pandas as pd
import numpy as np
import matplotlib
import seaborn as sns
from anndata import AnnData
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.font_manager as fm
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
import random
from descartes import PolygonPatch
import scanpy as sc
from ._utils import to_hex, crop, get_color_map, split_field
from typing import Any, List, Optional, Union, Dict, Tuple

def plot_scatter(
    adata: sc.AnnData,
    xlims: Optional[tuple] = None,
    ylims: Optional[tuple] = None,
    color_by: Optional[str] = None,
    highlight_cell: Optional[str] = None,
    highlight_color: str = '#FC4B42', 
    genes: Optional[Union[str, List[str]]] = None,
    seed: int = 123,
    alpha: float = 0.8,
    cmap: Optional[str] = None,
    ptsize: float = 0.1,
    ticks: bool = False,
    dpi: int = 300,
    figsize: tuple = (5,5),
    legend_loc: str = "center left",
    legend_col: int = 2,
    bbox_to_anchor: tuple = (1.0, 0.5),
    save: Optional[str] = None,
) -> None:
    """
    Scatter plot of spatial omics data based on AnnData object.
    
    Parameters
    ----------
    adata : sc.AnnData
        An anndata object containing the data.
    xlims : tuple, optional
        The x-axis limits for cropping, by default None.
    ylims : tuple, optional
        The y-axis limits for cropping, by default None.
    color_by : str, optional
        The column to color by, by default None.
    highlight_cell : str, optional
        Cell type to highlight, by default None.
    highlight_color : str, optional
        Color for the highlighted cells, by default '#FC4B42'.
    genes : Union[str, List[str]], optional
        Gene to plot, by default None.
    seed : int, optional
        Random seed for reproducibility, by default 123.
    alpha : float, optional
        Transparency level of the points, by default 0.8.
    cmap : str, optional
        Colormap to use for plotting, by default None.
    ptsize : float, optional
        Size of the points, by default 0.1.
    ticks : bool, optional
        Whether to show ticks on the plot, by default False.
    dpi : int, optional
        Resolution for the plot, by default 300.
    figsize : tuple, optional
        Size of the figure, by default (5, 5).
    legend_loc : str, optional
        Location of the legend, by default "center left".
    legend_col : int, optional
        Number of columns in the legend, by default 2.
    bbox_to_anchor : tuple, optional
        Bounding box to anchor the legend, by default (1.0, 0.5).
    save : str, optional
        Filename to save the plot to, by default None.
    highlight_cells : List[str], optional
        List of cell types to highlight, by default None.
    
    Returns
    -------
    None
    """
    
    def plot(
        ax: plt.Axes,
        new_coord: pd.DataFrame, 
        colors, 
        alpha: float = 0.8, 
        ptsize: float = 0.01, 
        ticks: bool = False,
    ) -> None:

        if not ticks:
            plt.tick_params(
                left=False,
                right=False,
                labelleft=False,
                labelbottom=False,
                bottom=False,
            )
            ax.axis("off")
            
        ax.scatter(np.array(coords[:,0]), np.array(coords[:,1]), marker='o', linewidth=0,
                   alpha=alpha, color=np.array(colors), s=ptsize)

        ax.axis("scaled")
        ax.grid(False)
        ax.invert_yaxis()

    fig, ax = plt.subplots(
        figsize=figsize,
        dpi=dpi,
        layout="constrained",
    )

    if not genes:
        
        subset_idx, new_coord = crop(adata, xlims, ylims)
        map_dict = get_color_map(adata, seed = seed, color_by = color_by, cmap = cmap, genes = genes, subset_idx = subset_idx)
        coords = adata[
            new_coord.index,
        ].obsm["spatial"]
        colors = new_coord[color_by].map(map_dict)

        if highlight_cell:
            if isinstance(highlight_cell, str):
                highlight_cell = [highlight_cell]
            mask = new_coord[color_by].isin(highlight_cell)
            colors = pd.Series('#BFBFBF', index=new_coord.index)
            colors[mask] = highlight_color

        plot(ax, coords, colors, alpha, ptsize, ticks)
        
        if not highlight_cell:
            markers = [plt.Line2D([0, 0], [0, 0], color=color, marker="o", linestyle="") for color in map_dict.values()]
            labels = map_dict.keys()
        else:
            markers = [plt.Line2D([0, 0], [0, 0], color=highlight_color, marker="o", linestyle="")]
            labels = highlight_cell
        
        fig.legend(
                    markers,
                    labels,
                    numpoints=1,
                    loc=legend_loc,
                    bbox_to_anchor=bbox_to_anchor,
                    frameon=False,
                    ncol=legend_col,
                )
        
        if save:
            plt.savefig(save, transparent=True, bbox_inches="tight")

        plt.show()

    else:
        if isinstance(genes, str):
            genes = [genes]

        counts = sc.get.obs_df(adata, keys=genes[0]).to_list()
        subset_idx, new_coord = crop(adata, xlims, ylims)

        cmap = get_color_map(adata, seed = seed, color_by = color_by, cmap = cmap, genes = genes, subset_idx = subset_idx)

        c_max = np.quantile(counts, 0.99)
        c_min = min(counts)
        bar_colors = [cmap(c / c_max) for c in counts]
        bar_colors = np.clip(bar_colors, 0, 1)
        all_colors = [to_hex(i) for i in bar_colors]
        all_colors_list = [all_colors[i] for i in subset_idx]

        coords = adata[
            new_coord.index,
        ].obsm["spatial"]

        plot(ax, coords, all_colors_list, alpha, ptsize, ticks)
        
        norm = matplotlib.colors.Normalize(vmin=c_min, vmax=c_max)
        cbar_ax = fig.add_axes([1.05, 0.2, 0.02, 0.6])
        fig.colorbar(
            matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax
        )
        if save:
            plt.savefig(save)

        plt.show()

def plot_polygon(adata: Any, 
                 plot_type: str = "cell", 
                 xlims: Optional[Tuple[float, float]] = None, 
                 ylims: Optional[Tuple[float, float]] = None, 
                 color_by: Optional[str] = None, 
                 genes: Optional[Union[str, List[str]]] = None, 
                 seed: int = 123, 
                 alpha: float = 0.8,
                 cmap: Optional[Union[str, Dict[str, str]]] = None, 
                 ptsize: float = 0.01, 
                 ticks: bool = False, 
                 dpi: int = 300,
                 figsize: tuple = (5, 5), 
                 edgecolor: str = "#808080", 
                 linewidth: float = 0.2,
                 legend_col: int = 2,
                 save: Optional[str] = None) -> None:
    """
    Plot polygons based on the provided adata object.

    Parameters
    ----------
    adata : Any
        An AnnData object containing the data.
    plot_type : str, optional
        The type of plot ('cell', 'gene', or 'transcript'), by default "cell".
    xlims : tuple, optional
        The x-axis limits for cropping, by default None.
    ylims : tuple, optional
        The y-axis limits for cropping, by default None.
    color_by : str, optional
        The column from adata.obs to color by, by default None.
    genes : Union[str, List[str]], optional
        Genes to be considered for 'gene' or 'transcript' plots, by default None.
    seed : int, optional
        Seed for random color generation, by default 123.
    alpha : float, optional
        Opacity of the polygons, by default 0.8.
    cmap : Union[str, Dict[str, str]], optional
        Color map or a dictionary to map values to colors, by default None.
    ptsize : float, optional
        Size of points in the scatter plot, by default 0.01.
    ticks : bool, optional
        Whether to display axis ticks, by default False.
    dpi : int, optional
        Resolution of the plot, by default 300.
    figsize : tuple, optional
        Size of the figure, by default (5, 5).
    edgecolor : str, optional
        Color of the polygon edge, by default "#808080".
    linewidth : float, optional
        Width of the polygon edge line, by default 0.2.
    legend_col : int, optional
        Number of columns in the legend, by default 2.

    Returns
    -------
    None
        Displays the desired plot based on the given parameters.
    """
    if 'poly' not in adata.uns:
        print('Please create polygons first! To create polygons, please use the create_polygons function.')
        return

    poly_dict = {idx + 1: adata.uns['poly'][key] for idx, key in enumerate(adata.uns['poly'].keys())}
    new_lib = list(poly_dict.values())
    subset_idx, new_coord = crop(adata, xlims, ylims)
    print("Total number of polygons: ", len(new_coord))

    def add_scalebar(ax):
        scalebar = AnchoredSizeBar(ax.transData,
                                   50, ' ', 'lower left', 
                                   pad=0.5,
                                   sep=5,
                                   color='black',
                                   frameon=False,
                                   size_vertical=2,
                                   fontproperties=fm.FontProperties(size=12))
        ax.add_artist(scalebar)

    def plot_cells(ax, new_lib, subset_idx, cols, edgecolor, alpha, linewidth, ticks):
        if not ticks:
            ax.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
        for j, i in enumerate(subset_idx):
            ax.add_patch(PolygonPatch(new_lib[i], fc=cols[j], ec=edgecolor, alpha=alpha, zorder=0.5, linewidth=linewidth))
        ax.axis('scaled')
        ax.grid(False)
        ax.invert_yaxis()
        add_scalebar(ax)

    if plot_type == "cell":
        if color_by is None:
            cols = ['#6699cc'] * len(subset_idx)
        else:
            map_dict = get_color_map(adata, seed = seed, color_by = color_by, cmap = cmap, genes = genes, subset_idx = subset_idx)
            cols = new_coord[color_by].map(map_dict).tolist()

        fig, ax = plt.subplots(dpi=dpi, figsize=figsize)
        plot_cells(ax, new_lib, subset_idx, cols, edgecolor, alpha, linewidth, ticks)

        if color_by is not None:
            markers = [plt.Line2D([0, 0], [0, 0], color=color, marker='o', linestyle='') for color in map_dict.values()]
            ax.legend(markers, map_dict.keys(), numpoints=1, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, ncol=legend_col)
        if save:
            plt.savefig(save)
        plt.show()

    elif plot_type == "gene":
        if isinstance(genes, str):
            genes = [genes]

        def create_gene_plot(ax, gene, counts, cmap, subset_idx, alpha, edgecolor, ticks):
            c_max = np.quantile(counts, 0.99)
            bar_colors = [cmap(c / c_max) for c in counts]
            bar_colors = np.clip(bar_colors, 0, 1)
            all_colors = [to_hex(i) for i in bar_colors]
            cols = [all_colors[i] for i in subset_idx]
            if not ticks:
                ax.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)

            ax.set_title(gene)
            for j, i in enumerate(subset_idx):
                ax.add_patch(PolygonPatch(new_lib[i], fc=cols[j], ec=edgecolor, alpha=alpha, zorder=0.1, linewidth=linewidth))
            ax.axis('scaled')
            ax.grid(False)
            ax.invert_yaxis()
            add_scalebar(ax)
            norm = matplotlib.colors.Normalize(vmin=min(counts), vmax=c_max)
            return norm

        n_rows = math.ceil(len(genes) / 3)
        n_cols = min(len(genes), 3)
        if cmap is None:
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ['#f2f2f2', '#ffdbdb', '#fc0303'])

        fig, axs = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=figsize, dpi=dpi)

        if len(genes) == 1:
            counts = sc.get.obs_df(adata, keys=genes[0]).to_list()
            norm = create_gene_plot(axs, genes[0], counts, cmap, subset_idx, alpha, edgecolor, ticks)
            divider = make_axes_locatable(axs)
            cax = divider.append_axes("right", size="3%", pad=0.05)
            fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs, cax=cax)
        else:
            axs = axs.flatten()
            for gene, ax in zip(genes, axs):
                counts = sc.get.obs_df(adata, keys=gene).to_list()
                norm = create_gene_plot(ax, gene, counts, cmap, subset_idx, alpha, edgecolor, ticks)
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="3%", pad=0.05)
                plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, cax=cax)

            for idx in range(len(genes), n_rows * n_cols):
                fig.delaxes(axs[idx])
        
        plt.tight_layout()
        if save:
            plt.savefig(save)
        plt.show()

    elif plot_type == "transcript":
        mol_data = adata.uns['transcript']
        cell_id_keep = new_coord.cell_id.tolist()
        mol_data2 = mol_data[mol_data.cell_id.isin(cell_id_keep) & mol_data.feature_name.isin(genes)]
        if cmap is None:
            random.seed(seed)
            tx_col = ["#" + ''.join([random.choice('0123456789ABCDEF') for _ in range(6)]) for _ in genes]
            map_dict = dict(zip(genes, tx_col))
        else:
            map_dict = cmap

        all_tx = mol_data2['feature_name'].to_list()
        all_tx = pd.DataFrame(all_tx, columns=['colors']).replace({'colors': map_dict})

        fig, ax = plt.subplots(dpi=dpi, figsize=figsize)
        plot_cells(ax, new_lib, subset_idx, ['#fcfcfa'] * len(subset_idx), edgecolor, alpha, linewidth, ticks)
        ax.scatter(mol_data2.x_location, mol_data2.y_location, c=all_tx.colors, s=ptsize)
        markers = [plt.Line2D([0, 0], [0, 0], color=color, marker='o', linestyle='') for color in map_dict.values()]
        ax.legend(markers, map_dict.keys(), numpoints=1, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, ncol=legend_col)
        if save:
            plt.savefig(save)
        plt.show()

def plot_scatter_img(
    adata: sc.AnnData,
    xlims: Optional[tuple] = None,
    ylims: Optional[tuple] = None,
    color_by: Optional[str] = None,
    genes: Optional[Union[str, List[str]]] = None,
    seed: int = 123,
    alpha: float = 0.8,
    cmap: Optional[str] = None,
    ptsize: float = 0.1,
    ticks: bool = False,
    dpi: int = 300,
    legend_loc: str = "center left",
    legend_col: int = 2,
    image=None,
    image_key="spatial",
    img_alpha=1,
    bbox_to_anchor: tuple = (1.0, 0.5),
    save: Optional[str] = None,
    cutoff: float = 0.2, 
) -> None:
    """
    Scatter plot of spatial omics data based on AnnData object with optional registered background image.
    
    Parameters
    ----------
    adata : sc.AnnData
        An anndata object containing the data.
    xlims : tuple, optional
        The x-axis limits for cropping, by default None.
    ylims : tuple, optional
        The y-axis limits for cropping, by default None.
    color_by : str, optional
        The column to color by, by default None.
    genes : Union[str, List[str]], optional
        Gene(s) to plot, by default None.
    seed : int, optional
        Random seed for reproducibility, by default 123.
    alpha : float, optional
        Transparency level of the points, by default 0.8.
    cmap : str, optional
        Colormap to use for plotting, by default None.
    ptsize : float, optional
        Size of the points, by default 0.1.
    ticks : bool, optional
        Whether to show ticks on the plot, by default False.
    dpi : int, optional
        Dots per inch for the plot, by default 300.
    legend_loc : str, optional
        Location of the legend, by default "center left".
    legend_col : int, optional
        Number of columns in the legend, by default 2.
    image : optional
        Background image to plot, by default None.
    image_key : str, optional
        Key for accessing the image from `adata.uns`, by default "spatial".
    img_alpha : float, optional
        Transparency level of the image, by default 1.
    bbox_to_anchor : tuple, optional
        Bounding box to anchor the legend, by default (1.0, 0.5).

    Returns
    -------
    None
    """
    
    def plot(
        ax: plt.Axes,
        new_coord: pd.DataFrame,
        colors,
        alpha: float = 0.8,
        ptsize: float = 0.01,
        ticks: bool = False,
        xlims: Optional[tuple] = None,
        ylims: Optional[tuple] = None,
    ) -> None:

        if not ticks:
            plt.tick_params(
                left=False,
                right=False,
                labelleft=False,
                labelbottom=False,
                bottom=False,
            )
            ax.axis("off")

        ax.scatter(
            np.array(coords[:, 0]),
            np.array(coords[:, 1]),
            marker="o",
            linewidth=0,
            alpha=alpha,
            color=np.array(colors),
            s=ptsize,
        )

        ax.axis("scaled")
        ax.grid(False)
        ax.invert_yaxis()

        if image:
            ax.imshow(adata.uns[image_key]["image"], alpha=img_alpha)
            if xlims and ylims:
                ax.set_xlim(xlims[0], xlims[1])
                ax.set_ylim(ylims[1], ylims[0])

    fig, axs = plt.subplots(
        figsize=(5,5),
        dpi=dpi,
        layout="constrained",
    )

    if not genes:
        subset_idx, new_coord = crop(adata, xlims, ylims)
        map_dict = get_color_map(adata, seed = seed, color_by = color_by, cmap = cmap, genes = genes, subset_idx = subset_idx)
        coords = adata[
            new_coord.index,
        ].obsm["spatial"]
        colors = new_coord[color_by].map(map_dict)
        plot(axs, coords, colors, alpha, ptsize, ticks, xlims, ylims)

        markers = [
            plt.Line2D([0, 0], [0, 0], color=color, marker="o", linestyle="")
            for color in map_dict.values()
        ]
        fig.legend(
            markers,
            map_dict.keys(),
            numpoints=1,
            loc=legend_loc,
            bbox_to_anchor=bbox_to_anchor,
            frameon=False,
            ncol=legend_col,
        )
        if save:
            plt.savefig(save)
        plt.show()
    else:
        
        if isinstance(genes, str):
            genes = [genes]
            
        counts = sc.get.obs_df(adata, keys=genes[0])
        subset_idx, new_coord = crop(adata, xlims, ylims)

        filtered_idx = np.where(counts >= cutoff)[0]
        filtered_counts = [counts[i] for i in filtered_idx]
        filtered_subset_idx = np.intersect1d(subset_idx, filtered_idx)
        filtered_coords = adata.obs.iloc[filtered_subset_idx]

        cmap = get_color_map(adata, seed = seed, color_by = color_by, cmap = cmap, genes = genes, subset_idx = subset_idx)

        c_max = np.quantile(counts, 0.99)
        c_min = min(counts)
        bar_colors = [cmap(c / c_max) for c in counts]
        bar_colors = np.clip(bar_colors, 0, 1)
        all_colors = [to_hex(i) for i in bar_colors]
        all_colors_list = [all_colors[i] for i in filtered_subset_idx]

        coords = adata[
            filtered_coords.index,
        ].obsm["spatial"]

        plot(axs,coords,all_colors_list,alpha,ptsize, ticks, xlims, ylims)

        norm = matplotlib.colors.Normalize(vmin=c_min, vmax=c_max)
        cbar_ax = fig.add_axes([1.05, 0.2, 0.02, 0.6])
        cb = fig.colorbar(
            matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cbar_ax,
            orientation="vertical",
        )
        if save:
            plt.savefig(save)
        plt.show()


        
def plot_visiumHD(
    adata: sc.AnnData,
    genes: List[str] = [],
    color_by: Optional[str] = None,
    cmap: Optional[matplotlib.colors.Colormap] = None,
    alpha: float = 0.8,
    bin_size: int = 16,
    ticks: bool = False,
    legend_loc: str = "center left",
    bbox_to_anchor: Tuple[float, float] = (1.0, 0.5),
    scale_bar: bool = False,
    scale: int = 200,
    legend_col: int = 2,
    image: bool = True,
    library_id: str = "",
    crop: bool = False,
    xlims: Optional[Tuple[float, float]] = None,
    ylims: Optional[Tuple[float, float]] = None,
    dpi: int = 300,
    figsize: tuple = (5, 5),
    seed: int = 42,
    cutoff: int = 1,
    save: Optional[str] = None,
):
    """
    Plots spatial data from VisiumHD technology.

    Parameters
    ----------
    adata : sc.AnnData
            The annotated data matrix to plot.
    genes : List[str], optional
            List of genes to visualize.
    color_by : Optional[str], optional
            Column name to color cells by.
    cmap : Optional[matplotlib.colors.Colormap], optional
            Colormap for cell coloring.
    bin_size : int, optional
            Bin size for scaling spatial coordinates.
    ticks : bool, optional
            If True, display axis ticks.
    legend_loc : str, optional
            Location of the legend.
    bbox_to_anchor : Tuple[float, float], optional
            Bounding box to anchor the legend.
    scale_bar : bool, optional
            If True, add scalebar.
    legend_col : int, optional
            Number of columns in the legend.
    image : bool, optional
            If True, display underlying tissue image.
    library_id : str, optional
            Identifier for the library to plot.
    crop : bool, optional
            If True, crop the image to specified limits.
    xlims : Optional[Tuple[float, float]], optional
            X-axis limits for cropping.
    ylims : Optional[Tuple[float, float]], optional
            Y-axis limits for cropping.
    dpi : int, optional
            Dots per inch for the figure.
    figsize : tuple, optional
            Size of the figure.
    seed : int, optional
            Seed for random processes.
    cutoff : int, optional
            Minimum expression level to display.
    save : Optional[str], optional
            Path to save the plot image.

    Returns
    -------
    None
    """
    
    ### Helper functions
    def to_hex(rgb_tuple):
        return "#{:02x}{:02x}{:02x}".format(
            int(rgb_tuple[0] * 255), int(rgb_tuple[1] * 255), int(rgb_tuple[2] * 255)
        )
    
    def add_scalebar(ax, scale):
        scalebar = AnchoredSizeBar(
            ax.transData,
            scale,
            " ",
            "lower left",
            pad=0.5,
            sep=5,
            color="black",
            frameon=False,
            size_vertical=4,
            fontproperties=fm.FontProperties(size=12),
        )
        ax.add_artist(scalebar)
        
    def plot(ax, new_spatial, new_coord, px_width, colors, ticks, image, scale, alpha):
        for i in range(new_coord.shape[0]):
            new_x, new_y = new_coord[i, 0], new_coord[i, 1]
            rect = Rectangle(
                (new_x, new_y),
                px_width,
                px_width,
                edgecolor="none",
                facecolor='#0F4C5C',
                alpha = alpha,
            )
            ax.add_patch(rect)
            
        ax.axis("scaled")
        ax.grid(False)
        ax.invert_yaxis()
        ax.imshow(adata.uns["spatial"]["images"]["hires"])
        width, height, _ = adata.uns["spatial"]["images"]["hires"].shape
        #5995.30408607, 5397.05046309
        ax.set_xlim(1150, 6000)
        ax.set_ylim(5400, 0)
        #ax.set_xlim(new_spatial[:, 0].min(), new_spatial[:, 0].max())
        #ax.set_ylim(new_spatial[:, 1].max(), new_spatial[:, 1].min())
        
        if scale_bar:
            add_scalebar(ax, scale)
            
    scale_factor = adata.uns["spatial"]["scalefactors"]["tissue_hires_scalef"]
    spatial_coords = adata.obsm["spatial"] * scale_factor
    json_data = adata.uns["spatial"]["scalefactors"]
    px_width = (
        bin_size / json_data["microns_per_pixel"] * json_data["tissue_hires_scalef"]
    )
    scale = scale / json_data["microns_per_pixel"] * json_data["tissue_hires_scalef"]
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    if not ticks:
        plt.tick_params(
            left=False, right=False, labelleft=False, labelbottom=False, bottom=False
        )
        ax.axis("off")
        
    if image:
        if crop:
            subset_idx = np.where(
                (spatial_coords[:, 1] > xlims[0])
                & (spatial_coords[:, 1] < xlims[1])
                & (spatial_coords[:, 0] > ylims[0])
                & (spatial_coords[:, 0] < ylims[1])
            )[0]
        else:
            subset_idx = np.arange(spatial_coords.shape[0])
            
        new_spatial = spatial_coords[subset_idx]
        
        if not genes:
            map_dict = get_color_map(
                adata,
                seed=seed,
                color_by=color_by,
                cmap=cmap,
                genes=genes,
                subset_idx=subset_idx,
            )
            
            all_colors = adata.obs[color_by].map(map_dict).to_numpy()
            colors = all_colors[subset_idx]
            plot(ax, new_spatial, new_spatial, px_width, colors, ticks, image, scale, alpha)
            markers = [
                plt.Line2D([0, 0], [0, 0], color='#0F4C5C', marker="s", linestyle="")
                for color in map_dict.values()
            ]
            fig.legend(
                markers,
                map_dict.keys(),
                numpoints=1,
                loc=legend_loc,
                bbox_to_anchor=bbox_to_anchor,
                frameon=False,
                ncol=legend_col,
            )
            
        else:
            if isinstance(genes, str):
                genes = [genes]
                
            counts = sc.get.obs_df(adata, keys=genes[0]).to_numpy().flatten()
            filtered_idx = np.where(counts >= cutoff)[0]
            filtered_subset_idx = np.intersect1d(subset_idx, filtered_idx)
            new_coord = (
                adata[
                    filtered_subset_idx,
                ]
                .obsm["spatial"]
                .copy()
                * scale_factor
            )
            cmap = get_color_map(
                adata,
                seed=seed,
                color_by=color_by,
                cmap=cmap,
                genes=genes,
                subset_idx=subset_idx,
            )
            norm = matplotlib.colors.Normalize(
                vmin=min(counts), vmax=np.quantile(counts, 0.99)
            )
            mapper = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
            
            rgba_colors = mapper.to_rgba(counts)
            all_colors = (rgba_colors[:, :3] * 255).astype(int)
            hex_colors = np.apply_along_axis(
                lambda x: "#{:02x}{:02x}{:02x}".format(*x), 1, all_colors
            )
            colors = hex_colors[filtered_subset_idx]
            
            plot(ax, new_spatial, new_coord, px_width, colors, ticks, image, scale, alpha)
            cbar_ax = fig.add_axes([0.95, 0.25, 0.015, 0.5])
            fig.colorbar(mapper, cax=cbar_ax)
            
    else:
        pxl_rows = adata.obs["array_row"].values
        pxl_cols = adata.obs["array_col"].values
        max_row = int(np.max(pxl_rows)) + 1
        max_col = int(np.max(pxl_cols)) + 1
        expression_grid = np.full((max_row, max_col), np.nan)
        
        if genes:
            if isinstance(genes, str):
                genes = [genes]
            data = sc.get.obs_df(adata, keys=genes[0]).to_numpy().flatten()
        else:
            def process_data_labels(adata, color_by, cmap):
                data_series = adata.obs[color_by]
                is_string_data = data_series.dtype == 'object' or isinstance(data_series.iloc[0], str)
                data = data_series.astype(str).values if is_string_data else data_series.values
                unique_data_labels = np.unique(data)
                
                if cmap:
                    original_labels = list(cmap.keys())
                    missing_labels = set(unique_data_labels) - set(cmap.keys())
                    if missing_labels:
                        raise ValueError(f"The following labels are missing in cmap: {missing_labels}")
                else:
                    if is_string_data:
                        original_labels = sorted(unique_data_labels)  # Alphabetical order
                    else:
                        original_labels = np.sort(unique_data_labels)  # Numerical order
                        
                # Map labels to indices
                label_to_index = {label: idx for idx, label in enumerate(original_labels)}
                indexed_data = np.array([label_to_index[label] for label in data])
                
                return indexed_data, original_labels
            
            def create_colormap(original_labels, cmap):
                if cmap:
                    colors = [cmap[label] for label in original_labels]
                else:
                    num_categories = len(original_labels)
                    colors = (
                        sns.color_palette("tab20b", num_categories)
                        if num_categories <= 20
                        else sns.color_palette("tab20b", 20) + sns.color_palette("tab20c", num_categories - 20)
                    )
                cmap_listed = matplotlib.colors.ListedColormap(colors)
                norm = matplotlib.colors.BoundaryNorm(np.arange(len(original_labels) + 1) - 0.5, cmap_listed.N)
                return cmap_listed, norm, colors
            
            data, original_labels = process_data_labels(adata, color_by, cmap)
            
        for i in range(len(pxl_rows)):
            row = int(pxl_rows[i])
            col = int(pxl_cols[i])
            expression_grid[row : row + 1, col : col + 1] = data[i]
            
        if genes:
            norm = mcolors.Normalize(
                vmin=np.nanmin(data), vmax=np.nanquantile(data, 0.99)
            )
            cmap = get_color_map(
                adata, seed=seed, color_by=color_by, cmap=cmap, genes=genes
            )
            heatmap = ax.imshow(
                expression_grid,
                aspect="auto",
                origin="lower",
                cmap=cmap,
                norm=norm,
                interpolation="none",
            )
            cbar_ax = fig.add_axes([0.95, 0.25, 0.015, 0.5])
            fig.colorbar(heatmap, cax=cbar_ax)
            
        else:
            cmap_listed, norm, colors = create_colormap(original_labels, cmap)
            
            heatmap = ax.imshow(
                expression_grid,
                aspect="auto",
                origin="lower",
                cmap=cmap_listed,
                norm=norm,
                interpolation="none",
            )
            
            handles = [
                plt.Line2D(
                    [0],[0],
                    marker="s",
                    color="w",
                    label=str(original_labels[i]),
                    markersize=10,
                    markerfacecolor=colors[i],
                )
                for i in range(len(original_labels))
            ]
            ax.legend(
                handles=handles,
                loc=legend_loc,
                bbox_to_anchor=bbox_to_anchor,
                frameon=False,
                ncol=legend_col,
            )
            
        ax.axis("scaled")
        ax.grid(False)
        
        # Crop the axes if required
        if crop:
            print(xlims, ylims)
            ax.set_xlim(xlims[0], xlims[1])
            ax.set_ylim(ylims[0], ylims[1])
            
    if save:
        plt.savefig(save, dpi = 300, transparent=True, bbox_inches='tight')
        
    plt.show()    
def plot_integrate(adata_vis, 
                   adata_xe,
                    gene = "Nphs2", 
                    xlims = None, 
                    ylims = None, 
                    ptsize_vis = 400,
                    ptsize_xe = 10,
                    cmap_vis = mcolors.LinearSegmentedColormap.from_list("Greens", ["#e5e5e5", "#a1d99b", "#31a354", "#004a1b"]),
                    cmap_xe = mcolors.LinearSegmentedColormap.from_list("Reds", ['#f2f2f2', '#ffdbdb', '#fc0303']),
                    dpi = 200,
                    color_bar=False):
    
    ### Plot xenium expression 
    counts_xe = sc.get.obs_df(adata_xe, keys=gene, use_raw=True).to_list()
    subset_idx, new_coord = crop(adata_xe, xlims, ylims)
    new_spatial_xe = np.array(adata_xe.obs[["x_centroid","y_centroid"]].iloc[subset_idx])
    norm_xe = matplotlib.colors.Normalize(vmin=min(counts_xe), vmax=np.quantile(counts_xe, 0.99))
    mapper_xe = matplotlib.cm.ScalarMappable(norm=norm_xe, cmap=cmap_xe)
    all_colors_xe = [to_hex(mapper_xe.to_rgba(c)) for c in counts_xe]
    xe_colors = [all_colors_xe[i] for i in subset_idx]

    ### Plot Visium spot expression
    counts_vis = sc.get.obs_df(adata_vis, keys=gene, use_raw=True).to_list()
    subset_idx_vis, new_coord_vis = crop(adata_vis, xlims, ylims)
    new_spatial_vis = adata_vis.obsm["spatial"][subset_idx_vis]
    norm = matplotlib.colors.Normalize(vmin=min(counts_vis), vmax=np.quantile(counts_vis, 0.99))
    mapper_vis = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_vis)
    all_colors_vis = [to_hex(mapper_vis.to_rgba(c)) for c in counts_vis]
    vis_colors = [all_colors_vis[i] for i in subset_idx_vis]

    fig, ax = plt.subplots()
    ax.scatter(new_spatial_xe[:, 0], new_spatial_xe[:, 1], marker='o', linewidth=0, alpha=1, color=xe_colors, s=ptsize_xe)
    ax.scatter(new_spatial_vis[:, 0], new_spatial_vis[:, 1],  marker='o', alpha=0.5, color=all_colors_vis, s= ptsize_vis, edgecolors='black',linewidth=0.2)
    if xlims and ylims:
        ax.set_xlim(xlims[0], xlims[1])
        ax.set_ylim(ylims[0], ylims[1])
    ax.axis('scaled')
    ax.invert_yaxis()
    ax.axis("off")

    if color_bar:
        cbar1 = fig.add_axes([0.90, 0.2, 0.01, 0.6])  
        fig.colorbar(mapper_xe, ax=ax, orientation='vertical', fraction=0.02, pad=0.02, cax=cbar1)
        cbar2 = fig.add_axes([1.0, 0.2, 0.01, 0.6])  
        fig.colorbar(mapper_vis, ax=ax, orientation='vertical', fraction=0.02, pad=0.02, cax=cbar2)
        
    plt.tight_layout()
    plt.show()


def plot_compare(adata_xe, Xenium_agg_gene_df, adata_vis, gene, cmap, save = None):

    counts_xe = sc.get.obs_df(adata_xe, keys = gene, use_raw=True).to_numpy().flatten()
    norm_xe = matplotlib.colors.Normalize(vmin=min(counts_xe), vmax=np.quantile(counts_xe, 0.99))
    mapper_xe = matplotlib.cm.ScalarMappable(norm=norm_xe, cmap=cmap)
    all_colors_xe = [to_hex(mapper_xe.to_rgba(c)) for c in counts_xe]
    
    counts = Xenium_agg_gene_df[gene].to_numpy()
    norm = matplotlib.colors.Normalize(vmin=min(counts), vmax=np.quantile(counts, 0.99))
    mapper = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    all_colors = [to_hex(mapper.to_rgba(c)) for c in counts]
    
    ### spot expression of True Visium
    counts_vis = sc.get.obs_df(adata_vis, keys = gene, use_raw=True).to_numpy().flatten()
    norm_vis = matplotlib.colors.Normalize(vmin=min(counts_vis), vmax=np.quantile(counts_vis, 0.99))
    mapper_vis = matplotlib.cm.ScalarMappable(norm=norm_vis, cmap=cmap)
    all_colors_vis = [to_hex(mapper_vis.to_rgba(c)) for c in counts_vis]
    
    x = adata_vis.obsm["spatial"][:,0]
    y = adata_vis.obsm["spatial"][:,1]
    
    fig, axs = plt.subplots(1, 3, dpi=150, figsize=(18, 6))  

    # Xenium plot
    axs[0].scatter(np.array(adata_xe.obs.x_centroid), np.array(adata_xe.obs.y_centroid), marker='o', linewidth=0, alpha=1, color=all_colors_xe, s=1)
    axs[0].axis('scaled')
    axs[0].axis("off")
    axs[0].invert_yaxis()
    axs[0].set_title("Xenium (observed)")
    
    # Xenium Pseudobulk
    axs[1].scatter(x, y, marker='o', linewidth=0, alpha=1, color=all_colors, s=40)
    axs[1].axis('scaled')
    axs[1].axis("off")
    axs[1].invert_yaxis()
    axs[1].set_title("Xenium (pesudobulk)")
    
    # Visium plot 
    axs[2].scatter(x, y, marker='o', linewidth=0, alpha=1, color=all_colors_vis, s=40)
    axs[2].axis('scaled')
    axs[2].axis("off")
    axs[2].invert_yaxis()
    axs[2].set_title("Visium")
    
    # Colorbar for the first plot (Xenium)
    cbar_ax1 = fig.add_axes([0.35, 0.2, 0.01, 0.6])  
    fig.colorbar(mapper_xe, ax=axs[2], orientation='vertical', fraction=0.02, pad=0.02, cax=cbar_ax1)
    
    # Colorbar for the second plot (Xenium Pseudobulk)
    cbar_ax2 = fig.add_axes([0.65, 0.2, 0.01, 0.6])  
    fig.colorbar(mapper, ax=axs[2], orientation='vertical', fraction=0.02, pad=0.02, cax=cbar_ax2)
    
    # Colorbar for the third plot (Visium)
    cbar_ax3 = fig.add_axes([0.90, 0.2, 0.01, 0.6])  
    fig.colorbar(mapper_vis, ax=axs[2], orientation='vertical', fraction=0.02, pad=0.02, cax=cbar_ax3)

    
    if save:
        plt.tight_layout()
        fig.savefig(save, dpi=500, transparent=True,  bbox_inches="tight", pad_inches=0)
        
    plt.show()
    

def plot_fov(
    data: "AnnData", 
    n_fields_x: int, 
    n_fields_y: int, 
    x_col: str = "x_centroid", 
    y_col: str = "y_centroid",
    group_label: Optional[str] = None, 
    highlight_cell: Optional[str] = None, 
    highlight_color: str = '#FC4B42',
    fill: bool = True, 
    point_size: float = 0.8, 
    alpha: float = 0.3,
    font_size: Optional[float] = None, 
    resolution: int = 200, 
    plot_width: int = 5, 
    plot_height: int = 5
) -> None:
    """
    Plot field of view (FoV) with optional cell highlighting and field partitioning.

    Parameters
    ----------
    data : AnnData
        An anndata object containing the data.
    n_fields_x : int
        Number of fields in the x dimension to partition the plot.
    n_fields_y : int
        Number of fields in the y dimension to partition the plot.
    x_col : str, optional
        Column name indicating x coordinates in the data, by default "x_centroid".
    y_col : str, optional
        Column name indicating y coordinates in the data, by default "y_centroid".
    group_label : str, optional
        Column name to use for grouping and highlighting cells, by default None.
    highlight_cell : str, optional
        Cell group label to be highlighted, by default None.
    highlight_color : str, optional
        Color for highlighted cells, by default '#FC4B42'.
    fill : bool, optional
        Whether to fill the polygons representing fields, by default True.
    point_size : float, optional
        Size of the plotted cell points, by default 0.8.
    alpha : float, optional
        Transparency level for the polygons and cell points, by default 0.3.
    font_size : float, optional
        Font size for the field numbers. Calculated based on n_fields_x and n_fields_y if not provided, by default None.
    resolution : int, optional
        Resolution of the plot in dpi, by default 200.
    plot_width : int, optional
        Width of the plot, by default 5.
    plot_height : int, optional
        Height of the plot, by default 5.

    Returns
    -------
    None
    """
    
    coord = data.obs.copy()
    
    rectangles, centroids = split_field(coord, n_fields_x, n_fields_y, x_col=x_col, y_col=y_col)
    
    # Convert rectangles to DataFrame for easier manipulation
    rectangle_df = pd.concat([pd.DataFrame(rect, columns=['x', 'y']).assign(id=idx) for idx, rect in enumerate(rectangles)])
    
    # Create polygons for plotting
    polygons = [{'coordinates': [group[['x', 'y']].values.tolist()], 'type': 'Polygon'} for _, group in rectangle_df.groupby('id')]
    
    if fill:
        fill_color = "#D3D3D3"
    else:
        fill_color = 'None'
    
    if font_size is None:
        font_size = 120 / (n_fields_x + n_fields_y)
    
    fig, ax = plt.subplots(dpi=resolution, figsize=(plot_width, plot_height))
    plt.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    
    if isinstance(group_label, str) and highlight_cell:
        coord['highlight'] = coord[group_label].apply(lambda x: highlight_color if x == highlight_cell else '#BFBFBF')
        [ax.plot(group_data[x_col], group_data[y_col], 'o', markersize=point_size, markerfacecolor=color, markeredgecolor='none', alpha=alpha) for color, group_data in coord.groupby('highlight')]
    else:
        ax.plot(coord[x_col].tolist(), coord[y_col].tolist(), 'o', markersize=point_size, markerfacecolor='#bfbfbf', markeredgecolor='none', alpha=alpha)
    
    for polygon in polygons:
        ax.add_patch(PolygonPatch(polygon, fc=fill_color, ec='black', alpha=alpha, linewidth=0.5))
    
    for idx, (x, y) in enumerate(centroids):
        ax.text(x, y, str(idx), fontsize=font_size, ha='center', va='center', weight='bold')
    
    ax.axis('scaled')
    ax.grid(False)
    ax.invert_yaxis()
    ax.set_frame_on(False)
    plt.show()

    
def subset_fov(adata: "AnnData", 
               fov: List[int], 
               n_fields_x: int, 
               n_fields_y: int, 
               x_col: str = 'x_centroid', 
               y_col: str = 'y_centroid') -> Tuple[float, float, float, float]:
    """
    Extract the coordinates that define the boundary of the specified field of view (FoV) from an AnnData object.

    Parameters:
    - adata: An anndata object.
    - fov: List of field of view indices to be considered.
    - n_fields_x, n_fields_y: Number of fields in x and y dimensions used to partition the plot.
    - x_col, y_col: Column names indicating x and y coordinates in the data.

    Returns:
    (xmin, xmax, ymin, ymax): Boundary coordinates in the format.
    """

    df = adata.obs.copy()
    rectangles, _ = split_field(df, n_fields_x, n_fields_y, x_col=x_col, y_col=y_col)
    
    selected_rectangles = [rectangles[i] for i in fov]
    min_pt = np.min(np.min(selected_rectangles, axis=0), axis=0)
    max_pt = np.max(np.max(selected_rectangles, axis=0), axis=0)
    
    return min_pt[0], max_pt[0], min_pt[1], max_pt[1]

def plot_network(
    adata: "AnnData", 
    df: pd.DataFrame, 
    color_map: Optional[dict] = None, 
    groupby: str = "celltype", 
    highlight: Optional[List[str]] = None, 
    seed: int = 1
) -> None:
    """
    Plot a network graph based on cell type proximity enrichment frequency.

    Parameters
    ----------
    adata : AnnData
        An anndata object containing the data.
    df : DataFrame
        DataFrame containing interaction scores between cell types.
    color_map : dict, optional
        Dictionary mapping cell types to colors, by default None.
    groupby : str, optional
        Column name to group by, by default "celltype".
    highlight : list, optional
        List of cell types to highlight, by default None.
    seed : int, optional
        Random seed for layout reproducibility, by default 1.

    Returns
    -------
    None
    """
    
    celltype_counts = pd.DataFrame(adata.obs[groupby].value_counts()).reset_index()
    celltype_counts.columns = ['CellType', "Count"]
    
    df.replace(0, 0.001, inplace=True)  # Add small values to facilitate the network visualization
    df = df.reset_index()
    zscore_long = pd.melt(df, id_vars=["index"], value_vars=df.columns)
    zscore_long.columns = ['CellType1', 'CellType2', 'zscore']
    zscore_long = zscore_long[zscore_long["CellType1"] != zscore_long["CellType2"]]
    zscore_long = zscore_long.sort_values(by=['CellType1', 'CellType2']).reset_index(drop=True)
    
    if color_map is None:
        celltypes = adata.obs[groupby].cat.categories
        n = len(celltypes)
        cmap = plt.get_cmap("tab20", lut=n)
        ct_colors = [plt.cm.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
        color_map = dict(zip(celltypes, ct_colors))
        
    if highlight:
        for ct, color in color_map.items():
            if ct in highlight:
                color_map[ct] = color
            else:
                color_map[ct] = '#d3d3d3'
        
    edge_weights_scale_factor = zscore_long['zscore'].max() / 5
    node_size_scale_factor = celltype_counts["Count"].max() / 500
    
    import networkx as nx
    G = nx.Graph()
    
    for _, row in zscore_long.iterrows():
        G.add_edge(row['CellType1'], row['CellType2'], weight=row['zscore'])
        
    for _, row in celltype_counts.iterrows():
        G.nodes[row['CellType']]['count'] = row['Count']
        G.nodes[row['CellType']]['color'] = color_map.get(row['CellType'])
        
    node_sizes = [G.nodes[node]['count'] / node_size_scale_factor for node in G.nodes]
    node_colors = [G.nodes[node]['color'] for node in G.nodes]
    edge_weights = [G[u][v]['weight'] for u, v in G.edges]
    edge_weights = np.array(edge_weights) / edge_weights_scale_factor
    
    pos = nx.spring_layout(G, seed=seed, iterations=1000)
    label_pos = {k: [v[0], v[1] - 0.01] for k, v in pos.items()}  
    
    fig, ax = plt.subplots(figsize=(4, 4), dpi=150)
    nx.draw_networkx_edges(G, pos, width=edge_weights, edge_color='#a7a7a7')
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, alpha=0.8)
    for node, coordinates in label_pos.items():
        plt.text(coordinates[0], coordinates[1], s=node, fontsize=8, ha='center')
    
    # Create discrete legend for edge width
    edge_bins = [1, 10, 20, 40]
    edge_labels = ['<1', '10', '20', '40']
    edge_size = edge_bins / edge_weights_scale_factor
    edge_handles = [Line2D([0], [0], color='#a7a7a7', lw=w, label=label) for w, label in zip(edge_size, edge_labels)]
    
    # Create discrete legend for node size
    size_bins = [500, 2500, 5000, 10000]
    size_labels = ['500', '2500', '5000', '10000']
    node_size = size_bins / node_size_scale_factor
    size_handles = [Line2D([0], [0], marker='o', color='#808080', markersize=np.sqrt(size), label=label, linestyle='None') for size, label in zip(node_size, size_labels)]
    
    # Add legends to the plot
    plt.legend(handles=edge_handles + size_handles, bbox_to_anchor=(1, 0.8), labelspacing=0.8, frameon=False, fontsize=8)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.show()
    