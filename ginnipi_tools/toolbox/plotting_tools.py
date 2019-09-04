#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 2018

Useful functions for plotting using nilearn plotting.
Some functions from plotting_tools by Ami Tsuchida (atsuch@gmail.com)

@author: tsuchida
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')

import bokeh.io
from bokeh.models import (Label, Range1d, ColumnDataSource,
                          CategoricalColorMapper, GlyphRenderer, HoverTool,
                          TapTool, OpenURL)
from bokeh.models.glyphs import HBar, VBar, Circle
from bokeh.transform import jitter
from bokeh.plotting import figure, output_file, show, save
from bokeh.layouts import gridplot, column, row
from bokeh.embed import file_html, autoload_static, components

from scipy.odr import Model, Data, ODR
import scipy.stats as stats


def get_lims(arr, delta=0.01):
    max_val = np.nanmax(arr)
    min_val = np.nanmin(arr)
    margin = np.absolute((max_val - min_val))* delta
    upper = max_val + margin
    lower = min_val - margin
    return (lower, upper)

def orthoregress(x, y):
    """Perform an Orthogonal Distance Regression on the given data,
    using the same interface as the standard scipy.stats.linregress function.
    Adapted from https://gist.github.com/robintw/d94eb527c44966fbc8b9#file-orthoregress-py
    
    Arguments:
    x: x data
    y: y data

    Returns:
    [slope, intercept, residual]

    Uses standard ordinary least squares to estimate the starting parameters
    then uses the scipy.odr interface to the ODRPACK Fortran code to do the
    orthogonal distance calculations.
    """
      
    def f(p, x):
        """Basic linear regression 'model' for use with ODR"""
        return (p[0] * x) + p[1]
    
    linreg = stats.linregress(x, y)
    mod = Model(f)
    dat = Data(x, y)
    od = ODR(dat, mod, beta0=linreg[0:2])
    out = od.run()

    return list(out.beta) + [out.res_var]


def stack_df(summary_df, stack_cols, stacked_col_name, val_name):
    """
    Stack a dataframe (summary_df)
    """
    # index columns that should not be stacked
    indices = [col for col in summary_df.columns if col not in stack_cols]
    if len(indices) == 0:
        stacked = summary_df.stack().reset_index(level=1)
        stacked.columns = [stacked_col_name, val_name]
    else:
        indexed = summary_df.set_index(indices)
        stacked = indexed.stack().reset_index()
        stacked.columns = indices + [stacked_col_name, val_name]

    return stacked


def get_bokehpalette(cmap, N):
    """
    Create bokeh palette of size N

    cmap: any matplotlib colormap
    """
    colormap =cm.get_cmap(cmap, lut=N)
    bokehpalette = [mpl.colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]

    return bokehpalette


def str_format_val(val, num_digits=3, signed=False):
    """
    Function to return float or scientific format of the val.
    
    Returns float format if absolute(val) >= 0.1
    Returns scientifc format if absolute(val) < 0.1
    """
    
    if np.absolute(val) >= 0.1:
        return '{1:+0.{0:d}f}'.format(num_digits, val) if signed else '{1:0.{0:d}f}'.format(num_digits, val)
    
    else:
        return '{1:+0.{0:d}e}'.format(num_digits, val) if signed else '{1:0.{0:d}e}'.format(num_digits, val)
    
        
    
def plot_dist_box_by_cols(summary_df, measure_name,
                          cols_to_plot, col_groupname, colnames=None,
                          title=None, horizontal=True,
                          cmap="Spectral", palette=None,
                          urlbase=None, part="anatomical", bgcolor="#D9EEF1",
                          plot_range=None, plot_size=None, 
                          linked_plot=None, out_html=None, nolegend=False):
    '''
    Build boxplot
    '''

    # We only need cols_to_plot
    # Note that it uses idx values and name to label scatter points
    df_toplot = summary_df[cols_to_plot]

    # default title
    if title is None:
        title = '%s by %s' % (measure_name, col_groupname)

    # Stack them to prepare for plotting
    stacked = stack_df(df_toplot, cols_to_plot, col_groupname, measure_name)

    # find the quartiles and IQR for each category
    groups = stacked.groupby(col_groupname)
    mean = groups.mean()
    count = groups.count()
    q1 = groups.quantile(q=0.25)
    q2 = groups.quantile(q=0.5)
    q3 = groups.quantile(q=0.75)
    iqr = q3 - q1
    upper = q3 + 1.5*iqr
    lower = q1 - 1.5*iqr

    # find the outliers for each category
    def outliers(group):
        cat = group.name
        return group[(group[measure_name] > upper.loc[cat][measure_name]) | (group[measure_name] < lower.loc[cat][measure_name])][measure_name]
    out = groups.apply(outliers).dropna()

    # if no outliers, shrink lengths of stems to be no longer than the minimums or maximums
    qmin = groups.quantile(q=0.00)
    qmax = groups.quantile(q=1.00)
    upper.loc[:,measure_name] = [min([x,y]) for (x,y) in zip(list(qmax.loc[:,measure_name]),upper.loc[:,measure_name])]
    lower.loc[:,measure_name] = [max([x,y]) for (x,y) in zip(list(qmin.loc[:,measure_name]),lower.loc[:,measure_name])]

    # Put them in a single df to create a datasource for box plots
    box_source = pd.concat([mean, count, q1, q2, q3, iqr, upper, lower], axis=1)
    box_source.columns = ['mean', 'count', '25p', '50p', '75p', 'IQR', 'upper', 'lower']
    box_source = ColumnDataSource(box_source)
    
    # if there are outliers, put them in a single df to create a datasource for
    # outlier scatter plots
    if not out.empty:
        out_dat = {}
        for cat in cols_to_plot:
            if not out.loc[cat].empty:
                num_out = len(out.loc[cat])
                out_dat[col_groupname] = out_dat.get(col_groupname, []) + [cat]*num_out
                out_dat[out.loc[cat].index.name] = out_dat.get(out.loc[cat].index.name, []) + (out.loc[cat].index.tolist())
                out_dat[measure_name] = out_dat.get(measure_name, []) + (out.loc[cat].values.tolist())

        out_source = ColumnDataSource(pd.DataFrame(out_dat))

    # Categorical color mapper
    # If palette is provided, it overrides cmap
    if palette is not None:
        my_palette = palette
    else:
        my_palette = get_bokehpalette(cmap, len(cols_to_plot))

    color_mapper = CategoricalColorMapper(
        factors=cols_to_plot, palette=my_palette)

    # properties for color etc.
    style_d = dict(fill_color={'field': col_groupname, 'transform': color_mapper},
                   fill_alpha=0.5, line_color='slategray')

    hover_style_d = dict(fill_color={'field': col_groupname, 'transform': color_mapper},
                        fill_alpha=0.9, line_color='black')
    
    # Set the plot range on the measure dimension if plot_range is not given
    if plot_range is None:
        plot_range = get_lims(stacked[measure_name].values)
    measure_range = Range1d(plot_range[0], plot_range[1])
    very_thin_range = 0.001*(plot_range[1] - plot_range[0])
     
    # plot top to bottom if horizontal
    if horizontal:
        y_range = cols_to_plot[::-1]
        if plot_size is None:
            plot_size = (800, 80+50*len(y_range))
        p = figure(plot_width=plot_size[0],
                   plot_height=plot_size[1],
                   x_range=measure_range, 
                   y_range=y_range, 
                   title=title,
                   background_fill_color=bgcolor)
        
        # components for vertically arranged horizontal boxplots
        seg1_d = dict(x0='25p', y0=col_groupname, x1='lower', y1=col_groupname)
        seg2_d = dict(x0='75p', y0=col_groupname, x1='upper', y1=col_groupname)
        box = HBar(y=col_groupname, right='75p', left='25p', height=0.5,
                   line_width=2, **style_d)
        box_h = HBar(y=col_groupname, right='75p', left='25p', height=0.5,
                     line_width=2, **hover_style_d)
        whisk1_d = dict(x='lower', y=col_groupname, width=very_thin_range, height=0.2)
        whisk2_d = dict(x='upper', y=col_groupname, width=very_thin_range, height=0.2)
        whisk3_d = dict(x='mean', y=col_groupname, width=very_thin_range, height=0.5)
        whisk4_d = dict(x='50p', y=col_groupname, width=very_thin_range, height=0.5)
        scatter = Circle(x=measure_name, size=8,
                         y=jitter(col_groupname, width=0.6, range=p.y_range),
                         **style_d)
        scatter_h = Circle(x=measure_name, size=8,
                         y=jitter(col_groupname, width=0.6, range=p.y_range),
                         **hover_style_d)
        p.xaxis.axis_label = "%s" % measure_name

    else:
        x_range = cols_to_plot
        if plot_size is None:
            plot_size = (80+50*len(x_range), 800)
        p = figure(plot_width=plot_size[0], 
                   plot_height=plot_size[1],
                   x_range=x_range, 
                   y_range=measure_range, 
                   title=title,
                   background_fill_color=bgcolor)

        # components for horizontally arranged vertical boxplots
        seg1_d = dict(x0=col_groupname, y0='25p', x1=col_groupname, y1='lower')
        seg2_d = dict(x0=col_groupname, y0='75p', x1=col_groupname, y1='upper')
        box = VBar(x=col_groupname, top='75p', bottom='25p', width=0.5,
                   line_width=2, **style_d)
        box_h = VBar(x=col_groupname, top='75p', bottom='25p', width=0.5,
                     line_width=2, **hover_style_d)
        whisk1_d = dict(x=col_groupname, y='lower', width=0.2, height=very_thin_range)
        whisk2_d = dict(x=col_groupname, y='upper', width=0.2, height=very_thin_range)
        whisk3_d = dict(x=col_groupname, y='mean', width=0.5, height=very_thin_range)
        whisk4_d = dict(x=col_groupname, y='50p', width=0.5, height=very_thin_range)
        scatter = Circle(x=jitter(col_groupname, width=0.6, range=p.y_range),
                         y=measure_name, size=8,
                         **style_d)
        scatter_h = Circle(x=jitter(col_groupname, width=0.6, range=p.y_range),
                         y=measure_name, size=8,
                         **hover_style_d)
        p.yaxis.axis_label = "%s" % measure_name

    # Now plot the box by each component
    p.segment(line_color='slategray', line_width=4, source=box_source, **seg1_d)
    p.segment(line_color='slategray', line_width=4, source=box_source, **seg2_d)
    p.rect(line_color='slategray', line_width=4, source=box_source, **whisk1_d)
    p.rect(line_color='slategray', line_width=4, source=box_source, **whisk2_d)
    p.rect(line_color='slategray', source=box_source, **whisk3_d)
    p.rect(line_color='slategray', line_width=4, source=box_source, **whisk4_d)
    box_glyph = GlyphRenderer(data_source=box_source, glyph=box, hover_glyph=box_h)

    # Make and add the hover tool
    box_tooltips = [(col_groupname, '@%s' % col_groupname),
                    ('count', '@count'),
                    ('mean', '@mean'),
                    ('median', '@50p'),
                    ('IQR', '@IQR'),
                    ('upper', '@upper'),
                    ('lower', '@lower')]
    box_hover = HoverTool(renderers=[box_glyph], tooltips=box_tooltips)
    p.add_tools(box_hover)
    p.renderers.extend([box_glyph])

    # Plot outliers if any
    if not out.empty:
        out_scatter = GlyphRenderer(data_source=out_source, glyph=scatter, hover_glyph=scatter_h)
        scatter_tooltips = [(df_toplot.index.name, '@%s' % df_toplot.index.name),
                            (col_groupname, '@%s' % col_groupname),
                            ('%s value' % measure_name, '@%s' % measure_name)]
        scat_hover = HoverTool(renderers=[out_scatter], tooltips=scatter_tooltips)
        
        
        url = urlbase + '/QCindividual/QCreport_' + '@%s'% df_toplot.index.name + '/indexed_{}.html'.format(part)
        tap_hover = TapTool(renderers=[out_scatter],callback=OpenURL(url=url))
        
        p.add_tools(scat_hover,tap_hover)
        p.renderers.extend([out_scatter])

    p.min_border = 10 #5
    #p.min_border_bottom = 100
    if nolegend:
        p.yaxis.visible=False
        p.min_border_left = 60
        p.min_border_right = 60
    
    if out_html is not None:
        bokeh.io.save(p, filename=out_html, title=title)
    return p


def build_histogram(df, measure_name, 
                    cols_to_plot, 
                    col_groupname,
                    palette=None, cmap="Spectral", bgcolor="#D9EEF1",
                    num_bins=100, bin_range=None, plot_size=(750, 450),
                    title=None,linked_plot=None, out_html=None):
    '''
    Build histogram
    '''
    # We only need cols_to_plot
    df = df[cols_to_plot]
    
    # default title
    if title is None:
        title = '%s by %s' % (measure_name, col_groupname)
        
    #p = figure(tools=tools, x_range=linked_plot.x_range)
    if linked_plot:
        p = figure(plot_width=plot_size[0],
                   plot_height=plot_size[1],
                   title=title,
                   x_range=linked_plot.x_range,
                   background_fill_color=bgcolor)
    else:
        p = figure(plot_width=plot_size[0],
                   plot_height=plot_size[1],
                   title=title,
                   background_fill_color=bgcolor)
    
    # Get histogram edges for the whole dataset if range is not provided
    if bin_range is None:
        hist, edges = np.histogram(df.dropna(), bins=num_bins)
    else:
        edges = np.linspace(bin_range[0], bin_range[1], num=num_bins, endpoint=False)
        np.append(edges, bin_range[1])
        if not linked_plot:
            p.x_range = Range1d(bin_range[0], bin_range[1])
    
    # Compute histogram values for each column in cols_to_plot and keep them in a df
    hist_df = pd.DataFrame({'bin_l': edges[:-1], 'bin_u': edges[1:]})
    for col in cols_to_plot:
        hist_df[col] = np.histogram(df[col].dropna(), bins=edges)[0]    
        
    # Stack them to create a source
    stacked = stack_df(hist_df, cols_to_plot, col_groupname, 'count')
    source = ColumnDataSource(stacked)
    
    # Categorical color mapper
    # If palette is provided, it overrides cmap
    if palette is not None:
        my_palette = palette
    else:
        my_palette = get_bokehpalette(cmap, len(cols_to_plot))

    color_mapper = CategoricalColorMapper(factors=cols_to_plot, palette=my_palette)
    
    # Tooltips
    tooltips = [('%s bin' % measure_name, '%s - %s' %('@bin_l', '@bin_u')),
                ('%s' % col_groupname, '@%s' % col_groupname),
                ('count', '@count')]
    hover = HoverTool(tooltips=tooltips)
    p.add_tools(hover)

    # plot histograms
    p.quad(top='count', bottom=0, left='bin_l', right='bin_u',
           fill_color={'field': col_groupname, 'transform': color_mapper},
           line_color='slategray', alpha=0.6,
           hover_fill_color={'field': col_groupname, 'transform': color_mapper},
           hover_line_color='black', hover_alpha=1.0,
           legend=col_groupname, source=source)

    # Save as html
    if out_html is not None:
        bokeh.io.save(p, filename=out_html, title=title, template="file.html")
    
    return p

def plot_hist_box(summary_df, 
                  measure_name,
                  col_groupname, cols_to_plot,
                  title=None, bgcolor=None,
                  plot_range=None, hist_bins=100,
                  plot_width=None,
                  palette=None, out_html=None,
                  nolegend = False,
                  urlbase='file:///beegfs_data/repos',
                  part="anatomical"):
    '''
    Stack histogram + boxplot in one figure
    + plot as a HTML if out_html != None
    '''
    
    # Set common range
    if plot_range is None:
        plot_range = get_lims(summary_df[cols_to_plot].values.flatten())
        
    # Default title
    if title is None:
        title = '%s distributions by %s' % (measure_name, col_groupname)
    
    # set proper plot size if plot_width is not None
    if plot_width is not None:
        hist_size = (plot_width, int(round(plot_width*(5.0/8.0), 0)))
        distbox_size = (plot_width, 80+50*len(cols_to_plot))
    else:
        hist_size = (800, 500)
        distbox_size = None
        
    bottom = plot_dist_box_by_cols(summary_df,
                                measure_name = measure_name,
                                col_groupname = col_groupname,
                                cols_to_plot = cols_to_plot,
                                plot_range = plot_range,
                                title = '',
                                palette = palette,
                                nolegend = nolegend,
                                urlbase = urlbase,
                                plot_size = distbox_size,
                                bgcolor = bgcolor,
                                part = part)
    
    top = build_histogram(summary_df,
                        measure_name = measure_name,
                        col_groupname = col_groupname,
                        cols_to_plot = cols_to_plot,
                        bin_range = plot_range,
                        num_bins = hist_bins,
                        title = title,
                        palette = palette,
                        linked_plot = bottom,
                        plot_size=hist_size,
                        bgcolor=bgcolor)
    
    plots = column(top, bottom)
    
    # Save as html
    if out_html is not None:
        save(plots, filename=out_html, title=title)

    return plots


def pairplots_by_region(summary_df, measure_name,
                        col1, col2,
                        title=None, color='navy', bgcolor="#D9EEF1", plot_size=(500, 500),
                        line_fit='ortho', out_html=None, urlbase=None, part = "anatomical"):

    # We only need cols_to_plot
    # Note that it uses idx values and name to label scatter points
    df_toplot = summary_df[[col1, col2]]

    # default title
    if title is None:
        title = 'Comparison of %s values between %s and %s' \
                % (measure_name, col1, col2)

    # get lims from range of values in col1 and col2
    lims = get_lims(df_toplot.values)
    p = bokeh.plotting.figure(plot_width=plot_size[0], plot_height=plot_size[1],
                              x_range=lims, y_range=lims,
                              title=title, background_fill_color=bgcolor)

    # Add a diagnal line in the background
    p.line(lims, lims, line_width=2, line_color='slategray',
           line_alpha=0.5, line_dash='dashed')
    source = bokeh.models.ColumnDataSource(df_toplot)
    points=p.circle(x=col1, y=col2, source=source,
                    size=10, fill_color=color,
                    hover_fill_color="firebrick",
                    fill_alpha=0.2, hover_alpha=0.5,
                    line_color='slategray', hover_line_color="firebrick")
    # Make and add the hover tool
    tooltips = [(df_toplot.index.name, '@%s' % df_toplot.index.name),
                ('%s' % col1, '@%s' % col1),
                ('%s' % col2, '@%s' % col2)]
    hover = bokeh.models.HoverTool(tooltips=tooltips, renderers=[points])
    p.add_tools(hover)
    p.renderers.extend([points])
    
    url = urlbase + '/QCindividual/QCreport_' + '@%s'% df_toplot.index.name + '/indexed_{}.html'.format(part)
    tap_hover = TapTool(renderers=[points],callback=OpenURL(url=url))
    p.add_tools(tap_hover)
    
    p.xaxis.axis_label = col1
    p.yaxis.axis_label = col2

    # Add line fit if line_fit is not None
    if line_fit is not None:
        # Drop na from data
        nonan_dat = df_toplot.dropna()
        x, y = nonan_dat[col1].values, nonan_dat[col2].values
        # Pearson correlation
        corr, corr_p = stats.pearsonr(x, y)
        delta = (lims[1] - lims[0])*0.05
        x_0 = lims[0] + delta
        x_1 = lims[1] - delta
        
        if line_fit == 'lin':  
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            txt2 = 'Rsquare=%0.3f' % (r_value**2)
        elif line_fit == 'ortho':
            slope, intercept, res = orthoregress(x, y) 
            txt2 = 'Residual Var=%0.3f' % (res)

        y_0 = slope * x_0 + intercept
        y_1 = slope * x_1 + intercept
        
        # Add line
        p.line([x_0, x_1], [y_0, y_1], line_width=1, line_color='firebrick')
        # Add txt
        if intercept >= 0:
            equation_txt = 'Y=%0.3f*X+%0.3e' % (slope, intercept) 
        else:
            equation_txt = 'Y=%0.3f*X-%0.3e' % (slope, np.absolute(intercept))
        corr_txt = 'Pearson r=%0.3f' % (corr)
        
        for i, txt in enumerate([equation_txt, txt2, corr_txt]):
            txt_label = Label(x=x_0, y=(x_1 - (i + 1)*delta), text=txt)
            p.add_layout(txt_label)

    if out_html is not None:
        title = out_html.split('.')[0]
        bokeh.io.save(p, filename=out_html, title=title)
    return p


def plot_surf_stat(lh_surf, lh_stat_map, lh_bg_map,
                   rh_surf, rh_stat_map, rh_bg_map, out_fname,
                   cmap='coolwarm', symmetric_cbar='auto',
                   upper_lim=None, threshold=None):
    '''Use Nilearn to plot statistical surface map for L lat, med, R lat, med view.
    '''
    import os.path as op
    import numpy as np
    import matplotlib.pyplot as plt
    import nibabel.freesurfer.io as fsio
    import nibabel.freesurfer.mghformat as fsmgh
    from nilearn import plotting
    
    # Get max and min value across stat maps from the two hemi
    lh_stat_dat = fsmgh.load(lh_stat_map).get_data()
    rh_stat_dat = fsmgh.load(rh_stat_map).get_data()
    flat_dat = np.hstack((lh_stat_dat.flatten(), rh_stat_dat.flatten()))
    max_val = np.maximum(np.abs(np.nanmax(flat_dat)), np.abs(np.nanmin(flat_dat)))
    if upper_lim is not None:
        vmax = upper_lim if max_val > upper_lim else max_val
    else:
        vmax = max_val
    
    fig, axs = plt.subplots(2, 2, figsize=(8,6), subplot_kw={'projection': '3d'})
    
    # Get threshold if txt is specified
    if isinstance(threshold, str):
        if op.exists(threshold):
            thresh_arr = np.loadtxt(threshold)
            thresh = thresh_arr if thresh_arr.shape == () and thresh_arr != np.inf else None
        else:
            thresh = None
    
    elif isinstance(threshold, int) or isinstance(threshold, float):
        thresh = threshold
        
    else:
        thresh = None
        
    # Add threshold to title if not None
    if thresh is not None:
        if float(thresh) >= 1e-2 and float(thresh) < 1e2:
            title_txt = 'Threshold = {:.2f}'.format(float(thresh))
        else:
            title_txt = 'Threshold = {:.2e}'.format(float(thresh))
    else:
        title_txt = ''
         
    for i, ax in enumerate(fig.axes):
        if i <= 1:
            hemi='left'
            surf = lh_surf 
            stat_map = lh_stat_dat
            bg = lh_bg_map
            
        else:
            hemi='right'
            surf = rh_surf 
            stat_map = rh_stat_dat
            bg = rh_bg_map
            
        view = 'lateral' if i % 2 == 0 else 'medial'
        colorbar = True if i == 3 else False  
        title = title_txt if i == 0 else ''
        
        plotting.plot_surf_stat_map(surf,
                                    stat_map,
                                    hemi=hemi,
                                    bg_map=bg,
                                    view=view,
                                    vmax=vmax,
                                    threshold=thresh,
                                    title=title,
                                    cmap=cmap,
                                    symmetric_cbar=symmetric_cbar,
                                    colorbar=colorbar,
                                    axes=ax,
                                    figure=fig)
    
    fig.savefig(out_fname, dpi=200, bbox_inches='tight')
    
    return op.abspath(out_fname)
    

def plot_surf_map(lh_surf, lh_surf_map, lh_bg_map,
                  rh_surf, rh_surf_map, rh_bg_map, out_fname,
                  cmap='jet', vmin=None, vmax=None):
    '''Use Nilearn to plot non-statistical surface map for L lat, med, R lat, med view.
    '''
    import os.path as op
    import numpy as np
    import matplotlib.pyplot as plt
    import nibabel.freesurfer.io as fsio
    import nibabel.freesurfer.mghformat as fsmgh
    from nilearn import plotting
    
    # Get max and min value across stat maps from the two hemi
    lh_surf_dat = fsmgh.load(lh_surf_map).get_data()
    rh_surf_dat = fsmgh.load(rh_surf_map).get_data()
    if vmin is None or vmax is None:
        flat_dat = np.hstack((lh_surf_dat.flatten(), rh_surf_dat.flatten()))
    if vmin is None:
        vmin = np.nanmin(flat_dat)
    if vmax is None:
        vmax = np.nanmax(flat_dat)
    
    fig, axs = plt.subplots(2, 2, figsize=(8,6), subplot_kw={'projection': '3d'})
    
    for i, ax in enumerate(fig.axes):
        if i <= 1:
            hemi='left'
            surf = lh_surf 
            surf_map = lh_surf_dat
            bg = lh_bg_map
            
        else:
            hemi='right'
            surf = rh_surf 
            surf_map = rh_surf_dat
            bg = rh_bg_map
            
        view = 'lateral' if i % 2 == 0 else 'medial'
        colorbar = True if i == 3 else False   
        
        plotting.plot_surf(surf,
                           surf_map,
                           hemi=hemi,
                           bg_map=bg,
                           view=view,
                           vmin=vmin,
                           vmax=vmax,
                           cmap=cmap,
                           colorbar=colorbar,
                           axes=ax,
                           figure=fig)
    
    fig.savefig(out_fname, dpi=200, bbox_inches='tight')
    
    return op.abspath(out_fname)

def create_html_view(stat_map_img, bg_img, out_fname='viewer.html', title=None,
                     opacity=1, symmetric_cmap=True, cmap='coolwarm', 
                     threshold=None, vmin=None, vmax=None):
    '''Use Nilearn to create html of volumetric map.
    '''
    import os.path as op
    import numpy as np
    from nilearn import plotting
    import nibabel as nib
    
    # Check if stat_map_img is all 0, in which case the plotting will fail...
    stat_map_dat = nib.load(stat_map_img).dataobj
    if np.min(stat_map_dat) == 0 and np.max(stat_map_dat) == 0:
        vmin, vmax = (0, 0)
    
    html_view = plotting.view_img(stat_map_img, bg_img=bg_img, title=title,
                                  threshold=threshold, symmetric_cmap=symmetric_cmap,
                                  black_bg=False, cmap=cmap, vmin=vmin, vmax=vmax,
                                  opacity=opacity)
     
    html_view.save_as_html(out_fname)
     
    return op.abspath(out_fname)
     
     