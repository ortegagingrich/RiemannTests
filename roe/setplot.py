
#import clawpack.geoclaw.shallow_1d.plot as geoplot

#let's try something
import sys
sys.path.insert(0,'/home/jog/Software_Development/uw/Riemann/roe/src/python')
import plot as geoplot
#from src.python import plot as geoplot

import setrun
rundata=setrun.setrun()


def setplot(plotdata):
    
    xlims=[-50.0,100.0]
    ylims=[-1,1.5]

    plotdata.clearfigures()

    plotfigure = plotdata.new_plotfigure(name='surface height', figno=0)
    
    #first plot the results from the full solver
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd='subplot(3,1,1)'
    plotaxes.xlimits = xlims
    plotaxes.ylimits = ylims
    plotaxes.title = 'full solver'

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = geoplot.surface_full
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = 'b'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'k'
    
    #second, plot the results from the roe solver
    plotaxes=plotfigure.new_plotaxes()
    plotaxes.axescmd='subplot(3,1,2)'
    plotaxes.xlimits=xlims
    plotaxes.ylimits=ylims
    plotaxes.title='roe solver'
    
    plotitem=plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var=geoplot.surface_roe
    plotitem.plot_var2=geoplot.topo
    plotitem.color='m'
    
    plotitem=plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var=geoplot.topo
    plotitem.color='k'
    
    #plot for the error
    plotaxes=plotfigure.new_plotaxes()
    plotaxes.axescmd='subplot(3,1,3)'
    plotaxes.xlimits=xlims
    plotaxes.ylimits='auto'
    plotaxes.title='roe error'
    
    plotitem=plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var=geoplot.roe_error
    plotitem.color='r'
    
    
    
    #moving on

    plotdata.printfigs = True          # Whether to output figures
    plotdata.print_format = 'png'      # What type of output format
    plotdata.print_framenos = 'all'      # Which frames to output
    plotdata.print_fignos = 'all'      # Which figures to print
    plotdata.html = True               # Whether to create HTML files
    plotdata.latex = False             # Whether to make LaTeX output

    return plotdata

