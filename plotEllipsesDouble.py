from matplotlib import pyplot as plt
from matplotlib import patches as ptc
from matplotlib.ticker import MaxNLocator
import sys
from math import *

filename = sys.argv[1]
filename2 = sys.argv[2]
param_filename = sys.argv[3]
normalized_As = sys.argv[4]
with open(filename) as f:
    content = [float(x) for x in f if x != '\n']
with open(filename2) as f:
    content2 = [float(x) for x in f if x != '\n']

with open(param_filename) as f:
    params = [line.rstrip('\n') for line in f]

num_params = int(content[0])
fig, ax = plt.subplots()

if normalized_As == '0':
    latex_params = {'ombh2':'\Omega_b h^2', 'omch2':'\Omega_{CDM} h^2', 'omk':'\Omega_k',\
                'hubble':'H_0', 'fesc':'f_{esc}', 'fstar':'f_*', 'w_DE':'w',\
                'omnuh2':'\Omega_\\nu h^2', 'T_CMB': 't_{CMB}', 'sigma8': '\sigma_8',\
                'n_s':'n_s', 'A_s':'A_s','100*theta_s':'100 \\theta_s','nion':'N_{ion}',\
                'omega_lambda':'\Omega_\Lambda', 'alpha':'\\alpha', 'beta':'\\beta',\
                'gamma':'\gamma', 'RLy':'R_{Ly}'}
else:
    latex_params = {'ombh2':'\Omega_b h^2', 'omch2':'\Omega_{CDM} h^2', 'omk':'\Omega_k',\
                'hubble':'H_0', 'fesc':'f_{esc}', 'fstar':'f_*', 'w_DE':'w',\
                'omnuh2':'\Omega_\\nu h^2', 'T_CMB': 't_{CMB}', 'sigma8': '\sigma_8',\
                'n_s':'n_s', 'A_s':'A_s * 10^9','100*theta_s':'100 \\theta_s','nion':'N_{ion}',\
                'omega_lambda':'\Omega_\Lambda', 'alpha':'\\alpha', 'beta':'\\beta',\
                'gamma':'\gamma', 'RLy':'R_{Ly}'}

ellipse_number = 0
for i in range(0,num_params - 1):
    for j in range(i+1, num_params):

        #a^2 has to correspond to the height because we're plotting the x parameter on the y axis
        w = 2*sqrt(content[ellipse_number*7 + 1])
        h = 2*sqrt(content[ellipse_number*7 + 2])
        w2 = 2*sqrt(content2[ellipse_number*7 + 1])
        h2 = 2*sqrt(content2[ellipse_number*7 + 2])

        #Here we need to convert radiants into degrees
        theta = 180*content[ellipse_number*7 + 3]/pi
        x = content[ellipse_number*7 + 4]
        y = content[ellipse_number*7 + 5]
        sig_x = content[ellipse_number*7 + 6]
        sig_y = content[ellipse_number*7 + 7]
        theta2 = 180*content2[ellipse_number*7 + 3]/pi
        x2 = content2[ellipse_number*7 + 4]
        y2 = content2[ellipse_number*7 + 5]
        sig_x2 = content2[ellipse_number*7 + 6]
        sig_y2 = content2[ellipse_number*7 + 7]


        frame_index = (num_params-1)*i + j
        ax1 = plt.subplot(num_params-1,num_params-1, frame_index)
        #ax1.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
        #ax1.ticklabel_format(axis='x',style='sci',scilimits=(-2,2))
        #ax1.yaxis.get_major_formatter().set_powerlimits((0,1))
        plt.subplots_adjust(hspace = 0., wspace = 0.)
        if (j-i) > 1:
            ax1.xaxis.set_ticklabels([])
            ax1.yaxis.set_ticklabels([])
        if (j-i) == 1:
    
            if params[i] in latex_params:
                plt.ylabel(r'$'+latex_params[params[i]]+'$', rotation='horizontal', fontsize = 20)
            else:
                plt.ylabel(params[i],rotation='horizontal', fontsize = 20) 
            ax1.yaxis.labelpad = 40

            #labels[-1].label1.set_visible(False)
            plt.gca().yaxis.set_major_locator(MaxNLocator(nbins = 5,prune='upper'))
            plt.gca().xaxis.set_major_locator(MaxNLocator(nbins = 5))

            if (i == num_params - 2):
                if params[i+1] in latex_params:
                    plt.xlabel(r'$'+latex_params[params[i+1]]+'$', fontsize = 20)
                else:
                    plt.xlabel(params[i+1], fontsize = 20)
                ax1.xaxis.labelpad = 10
        alpha = 1.52
        contour_1sig = ptc.Arc(xy = (x,y), width=alpha * w, height=alpha * h, angle=theta)
        ellipse_1sig = ptc.Ellipse(xy = (x,y), width=alpha * w, height=alpha * h, angle=theta)
        contour_1sig2 = ptc.Arc(xy = (x2,y2), width=alpha * w2, height=alpha * h2, angle=theta2)
        ellipse_1sig2 = ptc.Ellipse(xy = (x2,y2), width=alpha * w2, height=alpha * h2, angle=theta2,\
                alpha = 1, fc = 'r')


        alpha = 2.48
        contour_2sig = ptc.Arc(xy = (x,y), width=alpha * w, height=alpha * h, angle=theta)
        ellipse_2sig = ptc.Ellipse(xy = (x,y), width=alpha * w, height=alpha * h, angle=theta, alpha = 0.5)
        contour_2sig2 = ptc.Arc(xy = (x2,y2), width=alpha * w2, height=alpha * h2, angle=theta2)
        ellipse_2sig2 = ptc.Ellipse(xy = (x2,y2), width=alpha * w2, height=alpha * h2, angle=theta2,\
                alpha = 0.7, fc = 'r')

        ax1.add_artist(contour_1sig)
        ax1.add_artist(contour_2sig)
        ax1.add_artist(ellipse_1sig)
        ax1.add_artist(ellipse_2sig)
        
        ax1.add_artist(contour_1sig2)
        ax1.add_artist(contour_2sig2)
        ax1.add_artist(ellipse_1sig2)
        ax1.add_artist(ellipse_2sig2)
        
        ax1.plot([x],[y],marker='x',color='black', markersize=8, mew=2)

        
        #alpha = 1.
        #contour = ptc.Arc(xy = (x,y), width=alpha * w, height=alpha * h, angle=theta)
        #ax1.add_artist(contour)
        #r = ptc.Rectangle((x-alpha*sig_x,y-alpha*sig_y),2*alpha*sig_x,2*alpha*sig_y,fill=False)
        #ax1.add_artist(r)

        alpha = 3
        ax1.set_xlim([ x - alpha*sig_x,  x + alpha*sig_x])
        ax1.set_ylim([ y - alpha*sig_y,  y + alpha*sig_y])

        ellipse_number += 1

legend_elements = [ptc.Patch(facecolor = 'blue', label ='Power spectrum'),
                    ptc.Patch(facecolor = 'red', label ='Power spectrum + Bispectrum')]
plt.legend(handles=legend_elements, loc=(-3.5,0.5))

plt.show()
