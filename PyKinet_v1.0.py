#-------------------------------------------------------------------------------
# Name:        PyKinet
# Purpose:     This script is used to fit the model by Ballauf et al.
#              (https://doi.org/10.1021/jp101125j) to the kinetics data
#              for 4-nitrophenol reduction with sodium borohydride catalyzed
#              by metal nanoparticles
#
# Author:      Andrey Romanyuk
#-------------------------------------------------------------------------------


"""
The script works in several steps:
1) Necessary libraries are imported, the default directories are set up
2) Default values for the parameters are set up
3) Looks for data files in the subdirectory /raw, the names of the files should end with _SS.txt to get recognized and processed
4) Reads the data files line by line and extracts the data according to the pattern of data files:
    a) Block 1: parameter name and its value separated by a tab or a space
    b) Separator line: 20 dash symbols
    c) Block 2: the data starting with the line of headers
5) Groups the data sets by samples (all replicats of the same sample belong to one group)
6) The system of ODEs and the objective function (the one that gets minimized by the solver) are described
7) For each data set fitting is set up and the optimal values of the parameters are found
8) Both input and fitted data are plotted
    A) For individual samples
    B) For grouped samples
9) Both input and fitted data are exported to a txt file
    A) For individual samples
    B) For grouped samples
10) All rate constants are exported to a txt file

"""

# STEP_1
import os, re
import numpy as np
from scipy import integrate
import pylab as plt
from lmfit import Parameters, minimize, report_fit

# the directory of the script is the working directory
script_folder = os.path.dirname(os.path.abspath(__file__))
# the data files are expected in the subdirectory /raw
if os.path.isdir(script_folder + '/raw'):
    print('The directory for files with raw data was found')
else:
    print('No directory for files with raw data was found')

# the directory where to look for files
list_f = os.listdir(script_folder + '/raw')


# STEP_2
"""
Definition of the default values for the parameters
These values are used if not imported from the files containing data
"""
Knp_init, Kbh_init, k_init, n_init, m_init = 6922, 69, 1.6e-4, 0.6, 1.0
Cbh, S = 0.03, 1.0  # these default values will be replaced for the ones read from the data files

# STEP_3
# dictionary to store the input data
samples_dict = {'Filenames':[], 'Samples':[], 'Groups':[]}

"""
Files ending with "SS_raw.txt" will be recognized as data files and then processed (SS stands for steady state)
Also, the file name is expected to start with sX where X is one or several digits.
Therefore, the file name should have the following format: "sXX(any characters)_SS_raw.txt"
"""
for file in list_f:
    if file.endswith("_SS_raw.txt") or file.endswith("_SS_raw.TXT"):
        filename_w_ext = os.path.basename(file)
        samples_dict['Filenames'].append(filename_w_ext)

print('The following data files will be processed using the model for the stationary kinetics:')
for f_name in samples_dict['Filenames']:
    print(f_name)

# STEP_4
# reading all the data files and extracting the data
for file in samples_dict['Filenames']:
    sample = {}  # dictionary to store data of a single curve
    print('Processing file...\t', file)
    ts = []  # a list of time points in seconds
    c = []  # a list of experimental C(NP) values
    param_list = [Cbh, S] # a list of default values
    param_upd = [0, 0]  # flags: 0 -- updated, 1 -- not updated
    param_names = ['Cbh', 'S']  # names of the parameters that will be searched for in the file
    # the full list of parameters can be read from the file:
    # param_names = ['Knp_init', 'Kbh_init', 'n_init', 'm_init', 'k_init', 'Cbh', 'S']

    # opening files with UTF-8 enconding
    with open(script_folder + '/raw/' + file, 'r', encoding="utf8") as input_file:
        lines = input_file.readlines()  # read all lines and make them a list of lines

    # looking for the line separating the block of parameters from the block of data points
    for j in range(len(lines)):
        if re.fullmatch('-{20}\s+', lines[j]):  # the line should match the pattern precisely
            file_divider = j
            break
    # extracting parameters from the block before the separating line
    for j in range(0, file_divider, 1):
        for i in range(len(param_list)):
            if re.match(param_names[i] + '\s+', lines[j]) and param_upd[i] == 0:
                param_list[i] = float(re.split('\s+', lines[j])[1])
                param_upd[i] = 1
    # use the values from the file to update the parameters
    Cbh, S = [param_list[i] for i in range(len(param_list))]

    # Extraction of data after the divider and the headers line.
    # Each line with data should have the values for time point and concentration devided by a tab or a space.
    for l in lines[file_divider+2:]:
        if re.match('\d+(\.\d+)?\s\d+[\.\,]\d+E?-?\d+\s?', l):
            if l.find(',') != -1:  # if commas are used as a decimal point, they are replaced for dots
                l = l.replace(',', '.')
            split_line = re.split("\s", l)
            ts.append(float(split_line[0]))
            c.append(float(split_line[1]))
    # the data for each sample is stored as a separate dictionary
    sample = {'Title': file, 't': ts, 'c': c, 'X0': c[0], 'Knp': Knp_init, 'Kbh': Kbh_init, 'n': n_init, 'm': m_init,\
              'Cbh': Cbh, 'S': S}
    samples_dict['Samples'].append(sample)
    print('Parameters used (Knp, Kbh, n, m, Cbh, S, X0): ',\
          sample['Knp'], sample['Kbh'], sample['n'], sample['m'], sample['Cbh'], sample['S'], sample['X0'])

# STEP_5
# grouping the samples
for s in samples_dict['Samples']:
    # the file name is expected to have a name of the format sX where X is one or several digits
    if re.match('s\d\S+', s['Title']):
        g_name = re.split('_', s['Title'])[0]
        if g_name not in [g['Group name'] for g in samples_dict['Groups']]:
            samples_dict['Groups'].append({'Group name': g_name, 'Files': []})
            samples_dict['Groups'][-1]['Files'].append(s)
        else:
            samples_dict['Groups'][-1]['Files'].append(s)

# calculating the mean and SD values within each group
for group in samples_dict['Groups']:
    group['t'] = sorted(list(set(np.concatenate([file['t'] for file in group['Files']]))))
    for par in ['Knp', 'Kbh', 'n', 'm', 'Cbh', 'S']:
        group[par] = list(set([file[par] for file in group['Files']]))
        assert(len(group[par]) == 1), f'{par} values are not consistent within a group of replicates'
        group[par] = group[par][0]
    group['Data'] = {}  # will contain data  from replicates for each time point
    group['k_all'] = []  # list of rate constants obtained for replicates within one group
    group['Mean+SD'] = {'t': [], 'Mean': [], 'SD': []}
    for tt in group['t']:
        group['Data'][tt] = []  # Conc go into the 0th list, SDs go into the 1st one
        for file in group['Files']:
            if tt in file['t']:  # if there is a conc value at this time point
                group['Data'][tt].append(file['c'][list(file['t']).index(tt)])  # add this conc value to the list
        # if there are more than 1 conc value at this time point,
        # then we can average them out  and calculate SD
        # single values are ignored
        if len(group['Data'][tt]) > 1:
            group['Mean+SD']['t'].append(tt)
            group['Mean+SD']['Mean'].append(np.mean(group['Data'][tt]))
            group['Mean+SD']['SD'].append(np.std(group['Data'][tt]))
g_num = len(samples_dict['Groups'])  # the number of groups


# STEP_6
# the system of ODEs
def dX_dt(X, t, Knp, Kbh, k, n, m, Cbh, S):
    Cnp = X
    dxdt = -(k * S * ((Kbh * Cbh) ** m) * ((Knp * Cnp) ** n)) / (1 + ((Kbh * Cbh) ** m) + ((Knp * Cnp) ** n)) ** 2
    return dxdt

# the objective function
def ode_fit(params, t, data):
    solutions = integrate.odeint(dX_dt, y0=F['X0'], t=F['t'], args=(
    params['Knp'], params['Kbh'], params['k'], params['n'], params['m'], params['Cbh'], params['S']))[:, 0]
    resids = solutions - data  # return a 1D-array of residuals
    return resids

# STEP_8A
# plotting the results
plt.figure('Data fitting: Individual datasets')  # figure #1
plt.clf()
plt.suptitle('$K_{NP} = ' + f'{Knp_init:.2E}' + '\ [L/mol],\ K_{BH} = ' + f'{Kbh_init:.2E}' + '\ [L/mol]$', usetex = True, fontsize='18') # title for the plot
plt.ylabel('Cnp')  # title for the Y-axis
marker_style = dict(linestyle=' ', marker='o', markersize=6, markeredgewidth=1.5, fillstyle='none')
col_list = ['b', 'g', 'k', 'm', 'c']  # list of colors used to plot replicates

for g, group in zip(range(0, g_num), samples_dict['Groups']):
    print('Group\t' + group['Group name'] + '\n')
    plt.subplot(1, g_num + 1, g + 1)  # plot different concentrations in a separate subplot

    for F, clr in zip(group['Files'], col_list):
        # STEP_7
        # setting up the fitting
        fit_params = Parameters()
        # each parameter can take only positive values
        min_bnd = 0.0  # the lower bound for parameters
        max_bnd = np.inf  # the upper bound for parameters
        fit_params.add('Knp', value=F['Knp'], vary=False, min=min_bnd, max=max_bnd)
        fit_params.add('Kbh', value=F['Kbh'], vary=False, min=min_bnd, max=max_bnd)
        fit_params.add('k',  min=min_bnd, max=max_bnd, vary=True)
        fit_params.add('n', value=F['n'], min=min_bnd, max=max_bnd, vary=False)
        fit_params.add('m', value=F['m'], min=min_bnd, max=max_bnd, vary=False)
        fit_params.add('Cbh', value=F['Cbh'], vary=False)
        fit_params.add('S', value=F['S'], vary=False)

        # minimize the residuals
        result = minimize(ode_fit, params=fit_params, args=(F['t'], F['c']), method='least_squares', verbose=1)
        # show the report
        result.params.pretty_print()  # detailed fitting report
        report_fit(result)

        lmfit_params = [param.value for name, param in result.params.items()]  # a list of values of the parameters
        F['k'] = result.params['k'].value
        group['k_all'].append(F['k'])

        legend_params = f'$k_{(g+1)}' + f" = {F['k']:.2E}" + '\ [mol/m^2*s]$'
        F['Model'] = integrate.odeint(dX_dt, y0=F['X0'], t=F['t'], args=tuple(lmfit_params))[:, 0]

        # STEP_8A
        plt.plot(F['t'], np.array(F['c']), color=clr, **marker_style)
        plt.plot(F['t'], np.array(F['Model']), '--',  color=clr, label=legend_params,linewidth=2, markersize=5)

        # STEP_9A
        if not os.path.isdir(script_folder + '/model'):  # check if the subdirectory /model exists
            os.mkdir(script_folder + '/model')  # create the subdirectory /model for the output data
            print(r'The directory /model for output files was created')

        with open(script_folder + '/model/' + F['Title'][:-7] + 'model.txt', 'w', encoding="utf8") as output_file:
            output_line = 'Time [s]' + '\t' + 'Conc [M] (experiment)'+ '\t' + 'Conc [M] (Model)' + '\n'
            for e1, e2, e3 in zip(F['t'], F['c'], F['Model']):
                output_line += str(e1) + '\t' + str(e2) + '\t' + str(e3) + '\n'
            output_file.write(output_line)

    # STEP_8A
    group['k_mean'] = np.mean(group['k_all'])
    group['k_SD'] = np.std(group['k_all'])
    g_params = [group['Knp'], group['Kbh'], group['k_mean'], group['n'], group['m'], group['Cbh'], group['S']]
    group['Model'] = integrate.odeint(dX_dt, y0=group['Mean+SD']['Mean'][0], t=group['Mean+SD']['t'],
                                  args=tuple(g_params))[:, 0]
    label = f"$Mean\ k={group['k_mean']:.2e} \pm {group['k_SD']:.2e}$"
    plt.plot(group['Mean+SD']['t'], group['Model'], 'r--', label=label, linewidth=2, markersize=5)

    # settings for the figure #1
    plt.ylim((1e-6, 5.0e-5))
    plt.xlim((-5.0, 150.0))
    plt.xticks(np.arange(0.0, 151.0, 50))
    if g == 0:
        plt.ylabel('$C_{NP},\ [mol/L]$', fontsize=18)
    else:
        plt.tick_params(labelleft=False)
    # plt.yscale('log')  # optional log scale
    plt.tick_params(labelsize=14)
    if g == g_num//2:
        plt.xlabel('$Time\ [s]$', fontsize=18)  # title for the X-axis
    plt.legend(loc='upper right', prop={'size': 12})
# plt.show()  # show the graph

# STEP_8B
plt.figure('Data fitting: Grouped datasets')  # figure #2
plt.clf()
plt.suptitle('$K_{NP} = ' + f'{Knp_init:.2E}' + '\ [L/mol],\ K_{BH} = ' + f'{Kbh_init:.2E}' + '\ [L/mol]$', usetex=True,
             fontsize='18')  # title for the plot
plt.ylabel('Cnp')  # title for the Y-axis
for g, group in zip(range(0, g_num), samples_dict['Groups']):
    plt.subplot(1, g_num + 1, g + 1)  # plot different concentrations in a separate subplot
    group_legend = ''
    g_params = [group['Knp'], group['Kbh'], group['k_mean'], group['n'], group['m'], group['Cbh'], group['S']]
    for i, par in zip(g_params, ['Knp', 'Kbh', 'k', 'n', 'm', 'Cbh', 'S']):
        group_legend += par + ':  ' + f'{i:.2E}\n'
    plt.errorbar(group['Mean+SD']['t'], group['Mean+SD']['Mean'], yerr=group['Mean+SD']['SD'],
                 label=group['Group name'] + ' Mean', c='c', marker='o')
    plt.plot(group['Mean+SD']['t'], group['Model'], 'r--', label='Fit\n' + group_legend, linewidth=2, markersize=5)
    # settings for the figure 2
    plt.ylim((1e-6, 5.0e-5))
    plt.xlim((-5.0, 150.0))
    plt.xticks(np.arange(0.0, 151.0, 50))
    if g == 0:
        plt.ylabel('$C_{NP},\ [mol/L]$', fontsize=18)
    else:
        plt.tick_params(labelleft=False)
    # plt.yscale('log')  # optional log scale
    plt.tick_params(labelsize=14)
    if g == g_num//2:
        plt.xlabel('$Time\ [s]$', fontsize=18)  # title for the X-axis
    plt.legend(loc='upper right', prop={'size': 12})

    # STEP_9B
    with open(script_folder + '/model/' + group['Group name'] + '_model.txt', 'w', encoding="utf8") as data_output:
        data_lines = 'Time [s]' + '\t' + 'Conc [M] (mean)' + '\t' + 'SD' + '\t' + 'Conc [M] (model)' + '\n'
        for e1, e2, e3, e4 in zip(group['Mean+SD']['t'], group['Mean+SD']['Mean'], group['Mean+SD']['SD'], group['Model']):
            data_lines += str(e1) + '\t' + str(e2) + '\t' + str(e3) + '\t' + str(e4) + '\n'
        data_output.write(data_lines)
plt.show()  # show the graph

# STEP_10
# exporting the rate constants into a .txt file
rate_constants = [gr['k_all'] for gr in samples_dict['Groups']]
with open(script_folder + '/model/' + 'Individual rate constants.txt', 'w', encoding="utf8") as rate_output:
    rate_lines = 'Sample\tReplicates (k)\n'
    for g, group in zip(rate_constants, samples_dict['Groups']):
        rate_lines += group['Group name'] + '\t'
        for r in g:
            rate_lines += str(r) + '\t'
        rate_lines += '\n'
    rate_output.write(rate_lines)
