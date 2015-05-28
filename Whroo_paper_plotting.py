# -*- coding: utf-8 -*-
"""
Created on Mon May 25 14:27:04 2015

@author: imchugh
"""

import plots_generic as pg

def daily_NEE_years_comparison(df):
    pg.line_plot_test(df.index, df, title='$Whroo\/NEE$', hor_line=[0],
                      var_labels=['','','','','2011','2012','2013','2014'],
                      line_width=[0.5,0.5,0.5,0.5,3,3,3,3],
                      colors=['r','b','g','m'],
                      xlab='$DOY$',
                      ylab='$NEE\/(gC\/m^{-2}d^{-1})$',
                      xlim=[0,366],
                      ylim=[-6,4])
                      
def cuml_NEE_years_comparison(df):
    reload(pg)
    pg.line_plot_test(df.index, df, title='$Whroo\/cumulative\/NEE$', hor_line=[0],
                      var_labels=['2012','2013','2014'],
                      line_width=[3,3,3],
                      xlab='$DOY$',
                      ylab='$NEE\/(gC\/m^{-2})$',
                      xlim=[0,366])
    
                      
