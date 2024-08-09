from __future__ import print_function
import matplotlib

matplotlib.use('agg')
import csv
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from pyteomics import mzml, mgf
import seaborn as sns
import pkg_resources

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
import os
import argparse
import logging
import re
import matplotlib.gridspec as gridspec
from sklearn.metrics import mean_squared_error
from scipy import optimize
from scipy.stats import norm
from collections import Counter
from scipy.signal import find_peaks

COLORS = ["#b84c7d", "#4d8ac9", "#4bc490", "#7f63b8", "b", "g", '#edd630'] + ['k'] * 50


class PyplotContext:
    def __init__(
        self, separate_figures: bool,
        wide_figure: bool = False, output_path=None,
        grid_position=None, grid_size=None
    ):
        self.separate_figures = separate_figures
        self.output_path = output_path
        self.wide_figure = wide_figure
        self.grid_position = grid_position
        self.grid_size = grid_size

    def __enter__(self):
        if self.separate_figures:
            figsize = ((15. if self.wide_figure else 7.5), 7.)
            plt.figure(figsize=figsize)
        else:
            colspan = (2 if self.wide_figure else 1)
            plt.subplot2grid(self.grid_size, self.grid_position, colspan=colspan)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.separate_figures:
            plt.savefig(self.output_path)
            plt.close()       


class PyplotContextConstructor:
    def __init__(self, separate_figures: bool = False, out_dir=None, grid_size=None) -> None:
        self.separate_figures = separate_figures
        self.out_dir = out_dir
        self.grid_size = grid_size
        if separate_figures and not os.path.exists(out_dir):
            os.makedirs(out_dir)

    def __call__(self, name, grid_position=None, grid_size=None, wide_figure=False):
        return PyplotContext(
            separate_figures=self.separate_figures,
            output_path=os.path.join(self.out_dir, f'{name}.png'),
            grid_position=grid_position,
            grid_size=self.grid_size,
            wide_figure=wide_figure
        )


class GridOrSaveContext:
    def __init__(
        self, separate_figures: bool = False, out_path=None,
        grid_position=0, gs0=None, f: Figure = None):
        self.separate_figures = separate_figures
        self.out_path = out_path
        self.grid_position = grid_position
        self.gs0 = gs0
        self.f = f
        self.sps = self.gs0[self.grid_position]

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        if self.separate_figures:
            self.f.savefig(self.out_path)
            plt.close(self.f)


class GridOrSaveContextConstructor:
    def __init__(self, separate_figures: bool = False, out_dir=None, grid_size=None, figsize=None) -> None:
        self.separate_figures = separate_figures
        self.out_dir = out_dir
        self.grid_size = grid_size
        self.figsize = figsize
        if separate_figures and not os.path.exists(out_dir):
            os.makedirs(out_dir)
        if not self.separate_figures:
            self.f = plt.figure(figsize=figsize)
            self.gs0 = gridspec.GridSpec(*self.grid_size, figure=self.f)

    def __call__(self, name, grid_position=0):
        f = (plt.figure(figsize=self.figsize) if self.separate_figures else self.f)
        return GridOrSaveContext(
            separate_figures=self.separate_figures,
            out_path=os.path.join(self.out_dir, f'{name}.png'),
            grid_position=(0 if self.separate_figures else grid_position),
            gs0=(gridspec.GridSpec(1, 1, f) if self.separate_figures else self.gs0),
            f=f
        )


def read_data(name_mzml, name_calibration, extension, name_psm, delim, colname):
    name = os.path.split(name_mzml)[1].split('.')[0]
    logging.info('Reading data for %s... ' % name)
    injtime_ms1 = []
    injtime_ms2 = []
    starttime_ms1 = []
    starttime_ms2 = []
    indexms1 = []
    charge_ms2 = []
    mz_ms2 = []
    angle_x = []
    angle_y = []
    prec_int = []
    prec_isolated_mz = []
    x = []
    y = []
    errors = set()
    a = 0
    for i, sc in enumerate(mzml.read(name_mzml)):
        if 'injection time' not in errors:
            try:
                injtime = sc['scanList']['scan'][0]['ion injection time']
            except KeyError:
                errors.add('injection time')
                injtime_ms1 = None
                injtime_ms2 = None
        if 'start time' not in errors:
            try:
                starttime = sc['scanList']['scan'][0]['scan start time']
            except KeyError:
                errors.add('start time')
                starttime_ms1 = None
                starttime_ms2 = None
        if sc['ms level'] == 1:
            try:
                injtime_ms1.append(injtime)
            except AttributeError:
                pass
            try:
                starttime_ms1.append(starttime)
            except AttributeError:
                pass
            indexms1.append(i)
        if sc['ms level'] == 2:
            try:
                injtime_ms2.append(injtime)
            except AttributeError:
                pass
            try:
                starttime_ms2.append(starttime)
            except AttributeError:
                pass
            try:
                charge = int(sc['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
            except (KeyError, ValueError):
                charge = 0
                errors.add('charge')
            mz = sc['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
            try:
                isolated = sc["precursorList"]['precursor'][0]['isolationWindow']['isolation window target m/z']
            except (KeyError, ValueError):
                errors.add('precursor isolated m/z')
                isolated = None
            try:
                prec_intensity = sc['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0][
                    'peak intensity']
            except (KeyError, ValueError):
                errors.add('precursor intensity')
                prec_intensity = None
                a += 1
            if sc['intensity array'].size:
                intensity = sc['intensity array'].sum(dtype=float) / sc['intensity array'].size
            else:
                intensity = 0
            prec_isolated_mz.append(isolated)
            prec_int.append(prec_intensity)
            angle_y.append(intensity)
            angle_x.append(len(sc['intensity array']))
            mz_ms2.append(mz)
            charge_ms2.append(charge)

    def convert_times(arr):
        if arr:
            uinfo = arr[0].unit_info
            arr = np.array(arr)
            if uinfo == 'second':
                arr /= 60.
            elif uinfo != 'minute':
                logging.warning('Unexpected unit: %s', uinfo)
            return arr

    logging.debug('starttime_ms1: %s', starttime_ms1)
    starttime_ms1 = convert_times(starttime_ms1)
    logging.debug('starttime_ms1: %s', starttime_ms1)
    starttime_ms2 = convert_times(starttime_ms2)

    coef = 10 ** (len(str(int(np.mean(angle_y)))) - 1)

    # read_for_angle_calibration
    if name_psm is not None:
        logging.info('Reading reference files for angle score calculation...')
        d = {'mgf': mgf.IndexedMGF, 'mzml': mzml.MzML}
        reader = d[extension.lower()]
        f = reader(name_calibration)
        pattern = r'scan=\d+'
        replace = 'scan={}'
        if extension.lower() == 'mzml':
            id_format = re.sub(pattern, replace, f[0]['id'])
        if extension.lower() == 'mgf':
            id_format = re.sub(pattern, replace, f[0]['params']['title'])
        with open(name_psm) as fIn:
            reader = csv.DictReader(fIn, delimiter=delim)
            next(reader)
            for row in reader:
                psm = row[colname].split(' RTINSECONDS')[0]
                if psm in f:
                    match = True
                else:
                    if extension.lower() == 'mgf':
                        logging.info('Cannot match spectrum names, please check RefFile name for spaces')
                        f_split = dict()
                        for k in f.index.keys():
                            f_split[k.split(' ')[0]] = f[k]
                        f = f_split
                    if psm in f:
                        match = True
                    else:
                        match = False
                break
            for row in reader:
                psm = row[colname].split(' RTINSECONDS')[0]
                if match:
                    sc = f[psm]
                else:
                    sc_number = re.search(r'(\d+)\.\1', psm).group(1).lstrip('0')
                    sc = f[id_format.format(sc_number)]
                x1 = len(sc['m/z array'])
                x.append(x1)
                y1 = sc['intensity array'].sum(dtype=float) / sc['intensity array'].size
                y.append(y1)
        coef = 10 ** (len(str(int(np.mean(y)))) - 1)
        y = np.array(y) / coef
    prec_without_intense = round (100 * prec_int.count(None) / len(prec_int), 2)
    logging.info('Reading is complete.')
    for error in errors:
        logging.warning('There was an error extracting %s information. Some figures will not be produced.', error)
        logging.info('%s %% of precursor ions have no intensity information', prec_without_intense)
    return (injtime_ms1, injtime_ms2, starttime_ms1, starttime_ms2, indexms1, charge_ms2,
            mz_ms2, angle_x, angle_y, prec_int, prec_isolated_mz, x, y, coef)


def ms1_ms2(starttime_ms1, starttime_ms2):
    if starttime_ms1 is None or starttime_ms2 is None:
        plt.text(0.5, 0.5, 'Start time information missing', ha='center')
        return None, None
    width = 0.5
    ms1 = starttime_ms1.size
    ms2 = starttime_ms2.size
    plt.bar([0], ms1, width, alpha=1, color=COLORS[0])
    plt.text(0, ms1 / 2, ms1, ha='center', fontsize=20)
    plt.bar([1], ms2, width, alpha=1, color=COLORS[1])
    plt.text(1, ms2 / 2, ms2, ha='center', fontsize=20)
    plt.xticks(np.arange(2))
    xtick_marks = ['MS1', 'MS2']
    _, xtick_names = plt.xticks(np.arange(2), xtick_marks)
    plt.setp(xtick_names, rotation=0, fontsize=15)
    plt.title('MS1/MS2')
    return ms1, ms2


def aqtime(starttime_ms1, starttime_ms2):
    if starttime_ms1 is None or starttime_ms2 is None:
        plt.text(0.5, 0.5, 'Start time information missing', ha='center')
        return

    starttime_all = sorted([(x, 1) for x in starttime_ms1] + [(x, 2) for x in starttime_ms2])
    aqtime_ms1 = []
    aqtime_ms2 = []
    for i, k in enumerate(starttime_all[:-1]):
        aq = starttime_all[i + 1][0] - k[0]
        if k[1] == 1:
            aqtime_ms1.append(aq)
        if k[1] == 2:
            aqtime_ms2.append(aq)
    plt.hist(np.array(aqtime_ms1) * 60, histtype='step', lw=2, density=True,
               label='MS1, sum=%.2f min' % sum(aqtime_ms1), color=COLORS[0])
    plt.hist(np.array(aqtime_ms2) * 60, histtype='step', lw=2, density=True,
               label='MS2, sum=%.2f min' % sum(aqtime_ms2), color=COLORS[1])
    plt.legend(loc=1, fontsize=12)
    plt.xlabel("AT, sec", fontsize=15)
    plt.ylabel("# of spectra, normalized", fontsize=15)
    plt.title('Acquisition time')


def it_ms1(starttime_ms1, injtime_ms1, start, finish, mult):
    if starttime_ms1 is None:
        plt.text(0.5, 0.5, 'Start time information missing', ha='center')
        return None, None
    if injtime_ms1 is None:
        plt.text(0.5, 0.5, 'Injection time information missing', ha='center')
        return None, None
    plt.scatter(starttime_ms1, injtime_ms1, s=10, alpha=0.7, color=COLORS[0])
    plt.xlabel("Start time, min", fontsize=15)
    plt.ylabel("Injection time, ms", fontsize=15)
    plt.title('MS1')
    if mult:
        ind = np.logical_and(np.array(starttime_ms1) > start, np.array(starttime_ms1) < finish)
        mean = np.mean(np.array(injtime_ms1)[ind])
        # split = np.split(np.array(injtime_ms1)[ind][100:-(ind.sum()%100)], 100)
        # stack = np.vstack(split)
        # avg = np.sqrt((((stack-mean)**2).sum(axis=1)))
        # sliding_std = avg.std()
        perc_95 = np.percentile(np.array(injtime_ms1)[ind], 95)
        return mean, perc_95
    else:
        return None, None


def monoisotopic_error(charge_ms2, mz_ms2, prec_isolated_mz, mult):
    if set(prec_isolated_mz) == {None}:
        plt.text(0.5, 0.5, 'Prec. isolation m/z information missing', ha='center')
        return None, None, None

    mask = [i != 0 for i in charge_ms2]
    zeros = charge_ms2.count(0)
    isolated = np.array(prec_isolated_mz)[mask]
    mono = np.array(mz_ms2)[mask]
    ch = np.array(charge_ms2)[mask]
    diff = (isolated - mono) * ch

    plt.hist(np.round(diff, 0), bins=np.arange(-0.25, 5, 1), width=0.5, align='mid', color=COLORS[2])
    a, b = np.histogram(diff, bins=np.arange(-0.25, 5, 1))
    mod = 100 * sum(a[1:3]) / a[0]
    t = str(int(mod)) + '% $1^{st}$ and $2^{nd}$ isotopes \n' + str(zeros) + ' prec. with zero charges'
    ax = plt.gca()
    plt.text(0.7, 0.9, t, ha='center', va='center', transform=ax.transAxes)
    plt.title('Monoisotopic error')
    plt.ylabel('# scans')
    plt.xlabel('Isolated - Mono, Da')
    if mult:
        return mod


def inten_prec(starttime_ms2, start, finish, prec_int, mult):
    logging.debug('inten_prec received arguments: %s, %s, %s, %s, %s', starttime_ms2, start, finish, prec_int, mult)
    if starttime_ms2 is None:
        plt.text(0.5, 0.5, 'Start time information missing', ha='center')
        return None, None, None
    if set(prec_int) == {None}:
        plt.text(0.5, 0.5, 'Prec_intensity information missing', ha='center')
        return None, None, None
    ind = np.logical_and(starttime_ms2 > start, starttime_ms2 < finish)
    prec = np.array(prec_int, dtype=float)[ind]
    non0 = (prec > 0)
    logging.debug('%d of %d precursors have zero intensity, excluding them from analysis.',
        prec.size - non0.sum(), prec.size)

    prec = np.log10(prec[(~np.isnan(prec)) & non0])
    logging.debug('prec: %s', prec)
    logging.debug('prec min: %s, prec max: %s', prec.min(), prec.max())

    if mult:
        def gaussian(x, a, x0, sigma):
            return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

        counts, bins, bars = plt.hist(prec, bins=np.linspace(0, max(prec), 100), density=True, color=COLORS[0],
                                        alpha=0.5, lw=1, edgecolor='k')
        popt, pcov = optimize.curve_fit(gaussian, bins[:-1], counts, bounds=(0, np.nanmax(prec)))
        x = np.linspace(0, max(prec), 100)
        y = norm.pdf(x, popt[1], popt[2])
        plt.plot(x, gaussian(x, *popt), color='black')
        plt.xlim(np.nanpercentile(prec, 0.1) - 1, np.nanpercentile(prec, 99.9) + 1)
        plt.ylabel("# of spectra, normalized", fontsize=15)
        plt.xlabel("Log10(intensity)", fontsize=15)
        plt.title('Intensity of precursor ions')
        rsme = mean_squared_error(y[:-1], counts)
        return np.mean(prec), rsme, popt[2]
    else:
        a = np.nanpercentile(prec, 0.1)
        b = np.nanpercentile(prec, 99.9)
        logging.debug('a = %s, b = %s', a, b)
        plt.hist(prec, bins=np.linspace(0, max(prec), 100), color=COLORS[0], alpha=0.5, lw=1, edgecolor='k')
        plt.xlim(a, b)
        plt.ylabel("# of spectra", fontsize=15)
        plt.xlabel("Log10(intensity)", fontsize=15)
        plt.title('Intensity of precursor ions')
        return None, None, None


def it_ms2(starttime_ms2, start, finish, injtime_ms2):
    if starttime_ms2 is None:
        plt.text(0.5, 0.5, 'Start time information missing', ha='center')
        return
    if injtime_ms2 is None:
        plt.text(0.5, 0.5, 'Injection time information missing', ha='center')
        return
    ind = np.logical_and(starttime_ms2 > start, starttime_ms2 < finish)
    plt.hist(np.array(injtime_ms2)[ind], bins=np.linspace(0, max(injtime_ms2), 100),
               color=COLORS[1], alpha=0.5, lw=1, edgecolor='k')
    plt.xlabel("Injection time, ms", fontsize=15)
    plt.ylabel("# of spectra", fontsize=15)
    plt.title('MS2 injection time')


def realtop(starttime_ms1, indexms1, start, finish, mult):
    if starttime_ms1 is None:
        plt.text(0.5, 0.5, 'Start time information missing', ha='center')
        return
    plt.scatter(starttime_ms1[:-1], np.ediff1d(indexms1) - 1, alpha=0.7, s=15,
                  label='TopN', color=COLORS[1])
    fit = lowess(np.ediff1d(indexms1) - 1, starttime_ms1[:-1], frac=0.05, it=0)
    plt.plot(fit[:, 0], fit[:, 1], "r-", label='Average TopN')
    plt.ylim(0, max(np.ediff1d(indexms1)) * 1.1)
    plt.xlim(starttime_ms1.min() * .95, starttime_ms1.max() * 1.07)
    plt.legend(loc=1, markerscale=2, fontsize=15)
    plt.xlabel('RT, min', fontsize=15)
    plt.ylabel('TopN', fontsize=15)
    plt.title('Real TopN', fontsize=15)

    if mult:
        indexms1 = np.array(indexms1)
        ind = np.logical_and(np.array(starttime_ms1) > start, np.array(starttime_ms1) < finish)
        fit1 = lowess(np.ediff1d(indexms1[ind]) - 1, np.array(starttime_ms1[ind])[:-1], frac=0.05, it=0)
        # median_top_fit = np.median(fit1[:, 1])
        return fit1[:, 1]


def charge(maxcharge, charge_ms2, starttime_ms2, mz_ms2):
    ch_scans = []
    if starttime_ms2 is None:
        plt.text(0.5, 0.5, 'Start time information missing', ha='center')
        return
    for i in range(0, maxcharge + 1):
        charge_ms2 = np.array(charge_ms2)
        starttime_ms2 = np.array(starttime_ms2)
        mz_ms2 = np.array(mz_ms2)
        msk = [k == i for k in charge_ms2]
        plt.scatter(np.array(starttime_ms2)[msk], np.array(mz_ms2)[msk],
                      s=20, alpha=0.5, color=COLORS[i - 1],
                      label='charge_%s+(%i scans)' % (i, len(np.array(mz_ms2)[msk])))
        ch_scans.append(len(np.array(mz_ms2)[msk]))
        plt.title('Precursor ions', fontsize=15)
        plt.legend(markerscale=1.5, fontsize=15, facecolor='white')
        plt.xlabel('RT, min', fontsize=15)
        plt.ylabel('m/z', fontsize=15)
    return ch_scans


def angle_calculation(x, y, angle_x, angle_y, coef):
    x_mod = []
    y_mod = []
    per_1_x = np.percentile(x, 2)
    per_99_x = np.percentile(x, 98)
    per_1_y = np.percentile(y, 2)
    per_99_y = np.percentile(y, 98)
    for i in zip(x, y):
        if per_1_x < i[0] < per_99_x and per_1_y < i[1] < per_99_y:
            x_mod.append(i[0])
            y_mod.append(i[1])
    for i in np.arange(0, 1, 0.0001):
        dif = []
        for k in zip(x_mod, y_mod):
            dif.append(i * (k[0] - per_1_x) + per_1_y - k[1])
        under_line = sum(x > 0 for x in dif)
        if float(under_line) / len(x_mod) > 0.01:
            under = i
            break
        else:
            continue
    for j in np.arange(1, 0, -0.001):
        dif = []
        for k in zip(x_mod, y_mod):
            dif.append(j * (k[0] - per_1_x) + per_1_y - k[1])
        above_line = sum(x < 0 for x in dif)
        if float(above_line) / len(x_mod) > 0.01:
            above = j
            break
        else:
            continue
    a = 0
    angle_y_mod = np.array(angle_y) / coef
    for x1, y1 in zip(angle_x, angle_y_mod):
        if y1 - (x1 - per_1_x) * above < per_1_y and y1 - (x1 - per_1_x) * under > per_1_y:
            a += 1
    per = 100 * float(a) / len(angle_x)
    return under, above, per_1_x, per_1_y, per, angle_y_mod, coef


def angle(under, above, per_1_x, per_1_y, per, angle_x, angle_y_mod, coef):
    a = np.arange(per_1_x, max(angle_x), 100.)
    b = under * (a - per_1_x) + per_1_y
    plt.plot(a, b, color='black')
    c = above * (a - per_1_x) + per_1_y
    plt.plot(a, c, color='black')
    plt.scatter(angle_x, angle_y_mod, alpha=0.5, color=COLORS[1], label='score=%.2f' % per)
    plt.xlim(0, np.percentile(angle_x, 99.5))
    plt.ylim(0, np.percentile(angle_y_mod, 99.5))
    plt.xlabel('# peaks', fontsize=15)
    plt.legend(fontsize=15, loc=1)
    plt.ylabel('avg intensity, 10^%i' % np.log10(coef), fontsize=15)


def inten_number_peaks_ms1(angle_x, angle_y):
    coef = 10 ** (len(str(int(np.mean(angle_y)))) - 1)
    angle_y_mod = np.array(angle_y) / coef
    plt.scatter(angle_x, angle_y_mod, alpha=0.5, color=COLORS[1])
    plt.xlim(0, np.percentile(angle_x, 99.5))
    plt.ylim(0, np.percentile(angle_y_mod, 99.5))
    plt.xlabel('# peaks', fontsize=15)
    plt.ylabel('avg intensity, 10^%i' % np.log10(coef), fontsize=15)
    plt.title('MS/MS', fontsize=15)
    return np.nanmedian(angle_x), np.nanmedian(angle_y)


# Functions_for_multi_files

def name_red(names):
    stop = 0
    for i in range(len(names[0]) + 1):
        if len(set([x[:i] for x in names])) == 1:
            continue
        else:
            stop = i
            break
    names_red = [x[stop - 1:] for x in names]
    return names_red


def start_finish(indexms1, starttime_ms1):
    fit = lowess(np.ediff1d(indexms1) - 1, np.array(starttime_ms1)[:-1], frac=0.1, it=0)
    deriv = np.gradient(fit[:, 1], fit[:, 0])
    great_peaks_pos, _ = find_peaks(deriv, prominence=0.5)
    great_peaks_neg, _ = find_peaks(-deriv, prominence=0.5)
    if deriv.argmax() in great_peaks_pos:
        start = fit[:, 0][deriv.argmax()] + max(starttime_ms1) / 50
    else:
        start = min(starttime_ms1)
        logging.info('Cannot find start time, use %.2f', start)

    if deriv.argmin() in great_peaks_neg:
        finish = fit[:, 0][deriv.argmin()] - max(starttime_ms1) / 50
    else:
        finish = max(starttime_ms1)
        logging.info('Cannot find stop time, use %.2f', finish)

    return start, finish


def graph_with_break(names, ms1, ms2, title1, title2, ytitle, sps, f: Figure):
    '''gs0[0]'''
    gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=sps)
    n_files = len(names)
    ms1 = np.true_divide(ms1, 1000)
    ms2 = np.true_divide(ms2, 1000)
    ax2 = f.add_subplot(gs00[:-1, -1])
    ax1 = f.add_subplot(gs00[-1, -1])
    #     f, (ax2, ax1) = plt.subplots(2, 1, sharex=True)
    a = np.arange(0, n_files)
    ax1.plot(a, ms1, 'o:', color=COLORS[0], label='MS1')
    ax2.plot(a, ms2, 'o:', color=COLORS[1], label='MS2')
    ax1.set_title(title1 + ',10^3')
    ax2.set_title(title2 + ',10^3')
    ax1.set_ylim(min(ms1) - max(ms1) / 10, max(ms1) + max(ms1) / 10)
    ax2.set_ylim(min(ms2) - max(ms2) / 10, max(ms2) + max(ms2) / 10)
    ax2.set_xticks(range(len(names)))
    ax2.spines['bottom'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax2.xaxis.tick_top()
    ax2.tick_params(labeltop=False)
    ax1.xaxis.tick_bottom()
    plt.xticks(a, names, rotation=60)
    d = .015
    kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
    ax2.plot((-d, +d), (-d, +d), **kwargs)
    ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    kwargs.update(transform=ax1.transAxes)
    ax1.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)


def it_ms1_mult(names, mean_it_ms1, perc_95, sps, f: Figure):
    n_files = len(names)
    a = np.arange(0, n_files)
    g00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=sps)
    ax = f.add_subplot(g00[:, :])
    ax.errorbar(a, mean_it_ms1, perc_95, linestyle='None', capsize=10, fmt='-o', color=COLORS[0], ecolor='black',
                   elinewidth=1)
    ax.set_xticks(np.arange(n_files), names, rotation=60, fontsize=15)
    ax.set_ylim(min(np.array(mean_it_ms1) - np.array(perc_95)) - 1, max(np.array(mean_it_ms1) + np.array(perc_95)) + 1)
    #ax.legend(loc=1, fontsize=15)
    ax.set_ylabel('Injection time, ms', fontsize=15)
    ax.set_title('Mean IT MS1, 95th percentile')


def prec_int_mult(names, input_data, sps, f: Figure):
    '''gs0[2]'''
    input_names = ['mean', 'rsme', 'std']
    a = np.arange(0, len(names))
    gs11 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=sps)
    for i in range(3):
        ax = f.add_subplot(gs11[i, :])
        ax.plot(a, input_data[i], 'o:', color=COLORS[i], label=input_names[i])
        ax.set_ylim(min(input_data[i]) - max(input_data[i]) / 10, max(input_data[i]) + max(input_data[i]) / 10)
        ax.legend()
        ax.set_xticks(range(len(names)))
        ax.tick_params(axis='x', colors='w')
        if i == 0:
            ax.set_title('Precursor intensity')
    plt.xticks(np.arange(len(names)), names, rotation=60, color='black')


def charge_mult(names, ch_scans, sps, f: Figure):
    pad = len(max(ch_scans, key=len))
    charges = np.array([i + [0] * (pad - len(i)) for i in ch_scans])
    msk = ~np.all(charges == 0, axis=0)
    gs00 = gridspec.GridSpecFromSubplotSpec(msk.sum(), 1, subplot_spec=sps)
    for i, j in zip(np.arange(msk.sum()), np.arange(pad)[msk]):
        a = np.arange(0, len(names))
        ch_state = charges[:, j]
        ax = f.add_subplot(gs00[i, :])
        ax.plot(a, ch_state, 'o:', color=COLORS[i], label='ch_state=%i+' % j)
        ax.legend()
        ax.set_xticks(range(len(names)))
        ax.tick_params(axis='x', colors='w')
        if i == 0:
            ax.set_title('#scans with such prec.')
    plt.xticks(np.arange(len(names)), names, rotation=60, color='black')


def real_top_mult(names, fit, sps, f: Figure):
    n_files = len(names)
    gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=sps)
    ax = f.add_subplot(gs00[:, :])
    sns.boxplot(data = fit, orient='v', color=COLORS[2], whis=1, notch=False, ax=ax)
    ax.set_xticks(np.arange(n_files), names, rotation=60)
    ax.set_ylabel('smooth topN')


def simple_graph_mult(names, values, title, sps, f: Figure):
    a = range(len(names))
    gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=sps)
    ax = f.add_subplot(gs00[:, :])
    ax.plot(a, values, 'o:', color=COLORS[3])
    ax.set_title(title)
    r = (max(values) - min(values)) / 10
    ax.set_ylim(min(values) - r, max(values) + r)
    ax.set_xticks(a, names, rotation=60)
    #plt.legend()


def process_file(name_mzml, args):
    if len(args.input) > 1:
        mult = True
    else:
        mult = False
    if args.output is None:
        output = os.path.split(name_mzml)[0]
    else:
        output = args.output

    if (args.refPSM is None) + (args.refFile is None) == 1:
        logging.error('Not enough files for angle score calculation provided, skipping angle score.')


    name = os.path.split(name_mzml)[1].split('.')[0]
    name_psm = args.refPSM
    name_calibr = args.refFile
    extension = os.path.split(name_calibr)[1].split('.')[1] if name_calibr else None
    delim = str(args.d)
    colname = str(args.cn)
    plt_manager = PyplotContextConstructor(args.sf, os.path.join(output, name + '_viQC'), grid_size=(6,2))

    if (args.refPSM is not None) and (os.path.split(name_calibr)[1].split('.')[0] not in name_psm):
        logging.warning('File names for angle score calibration don\'t match!')

    injtime_ms1, injtime_ms2, starttime_ms1, starttime_ms2, indexms1, charge_ms2, \
    mz_ms2, angle_x, angle_y, prec_int, prec_isolated_mz, x, y, coef = read_data(
        name_mzml, name_calibr, extension, name_psm, delim, colname)

    logging.debug('starttime_ms1: %s', starttime_ms1)

    if args.stop is None:
        start, finish = start_finish(indexms1, starttime_ms1) if starttime_ms1 is not None else (None, None)
    else:
        finish = args.stop
        start = args.start
    logging.debug('start = %s, finish = %s', start, finish)
    if args.charge is None:
        maxcharge = max(charge_ms2)
        logging.info('Maximum charge in file: %s', maxcharge)
    else:
        maxcharge = args.charge

    # for pretty pictures
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams['ytick.labelsize'] = 15
    plt.rcParams['xtick.labelsize'] = 15
    plt.rcParams['axes.titlesize'] = 15
    plt.rcParams['axes.titlesize'] = 15
    plt.rcParams['font.size'] = 15

    # build pictures
    if not args.sf:
        fig = plt.figure(figsize=(15, 40))
    with plt_manager("ms1_ms2", (0, 0), False):
        ms1_f, ms2_f = ms1_ms2(starttime_ms1, starttime_ms2)
    with plt_manager("aqtime", (0, 1), False):
        aqtime(starttime_ms1, starttime_ms2)

    with plt_manager("inten_prec", (1, 0), False):
        mean_prec_f, rsme_prec_f, std_prec_f = inten_prec(starttime_ms2, start, finish, prec_int, mult)
    with plt_manager("monoisotopic_error", (1, 1), False):
        MI_error_percent = monoisotopic_error(charge_ms2, mz_ms2, prec_isolated_mz, mult)

    with plt_manager("it_ms1", (2, 0), True):
        mean_it_f, perc_95_it_f = it_ms1(starttime_ms1, injtime_ms1, start, finish, mult)

    with plt_manager("it_ms2", (3, 0), True):
        it_ms2(starttime_ms2, start, finish, injtime_ms2)

    with plt_manager("charge", (4, 0), True):
        ch_state_numbers = charge(maxcharge, charge_ms2, starttime_ms2, mz_ms2)

    with plt_manager("realtop", (5, 0), False):
        fit_realtop = realtop(starttime_ms1, indexms1, start, finish, mult)

    if name_psm is not None:
        under, above, per_1_x, per_1_y, per, angle_y_mod, coef = angle_calculation(x, y, angle_x, angle_y, coef)
        with plt_manager("angle", (5, 1), False):
            angle(under, above, per_1_x, per_1_y, per, angle_x, angle_y_mod, coef)
        median_peaks_ms2, median_intens_ms2 = None, None
    else:
        with plt_manager("inten_number_peaks_ms1", (5, 1), False):
            median_peaks_ms2, median_intens_ms2 = inten_number_peaks_ms1(angle_x, angle_y)
        per = None

    outname = os.path.join(output, name + '_viQC.%s' % args.pic)
    if not args.sf:
        plt.savefig(outname)
        # outname = os.path.join(output, name + '_viQC.svg')
        # plt.savefig(outname)
        plt.close(fig)
    logging.info("Calculating metrics for %s", name)
    return ms1_f, ms2_f, mean_it_f, perc_95_it_f, mean_prec_f, rsme_prec_f, std_prec_f, ch_state_numbers, fit_realtop, median_peaks_ms2, median_intens_ms2, per, MI_error_percent


def mult_process(name_files, args):
    results_mult = [[] for i in range(13)]
    for name_mzml in name_files:
        for value, value_list in zip(process_file(name_mzml, args), results_mult):
            value_list.append(value)
            
    if args.output is None:
        output = os.path.split(name_mzml)[0]
    else:
        output = args.output
    outname = os.path.join(output, 'Common_%i_files_viQC.%s' % (len(name_files), args.pic))

    names_full = [os.path.split(x)[1].split('.')[0] for x in name_files]
    names = name_red(names_full)
    logging.info("Combine all together...")
    #  build common pictures
    grid_manager = GridOrSaveContextConstructor(
        args.sf, os.path.splitext(outname)[0], grid_size=(4, 2),
        figsize=((15, 10) if args.sf else(30, 40))
    )
    with grid_manager("N of MS scans", 0) as gm:
        graph_with_break(names, results_mult[0], results_mult[1], '# of MS1 scans', '# of MS/MS scans', '#scans, 10^3',
                        gm.sps, gm.f)  # MS1/MS2 graph
    with grid_manager("Injection time", 1) as gm:
        it_ms1_mult(names, results_mult[2], results_mult[3], gm.sps, gm.f)
    with grid_manager("Precursor intensity", 2) as gm:
        prec_int_mult(names, [results_mult[4], results_mult[5], results_mult[6]], gm.sps, gm.f)
    with grid_manager("N scans with such prec.", 3) as gm:
        charge_mult(names, results_mult[7], gm.sps, gm.f)
    with grid_manager("smooth topN", 4) as gm:
        real_top_mult(names, results_mult[8], gm.sps, gm.f)

    if set(results_mult[11]) == {None}:
        with grid_manager("smooth topN", 5) as gm:
            simple_graph_mult(names, results_mult[9], 'median #peaks in MS/MS', gm.sps, gm.f)
        with grid_manager("smooth topN", 7) as gm:
            simple_graph_mult(names, results_mult[10], 'median intensity in MS/MS', gm.sps, gm.f)
    else:
        with grid_manager("smooth topN", 5) as gm:
            simple_graph_mult(names, results_mult[11], 'Angle score', gm.sps, gm.f)

    with grid_manager("smooth topN", 6) as gm:
        simple_graph_mult(names, results_mult[12], 'Monoisotopic error, %', gm.sps, gm.f)

    if not args.sf:
        grid_manager.f.savefig(outname)

    logging.info('Output common figure saved to %s', outname)


class StatsHandler(logging.Handler):
    # inspired by https://stackoverflow.com/a/31142078/1258041

    def __init__(self, *args, **kwargs):
        super(StatsHandler, self).__init__(*args, **kwargs)
        self.level2count = Counter()

    def emit(self, record):
        l = record.levelno
        self.level2count[l] += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs='+',
                        help='mzML file(s) with path(s), if more than 1 file specified the common picture is built')
    parser.add_argument('-o', '--output',
                        help='path to save result, by default save in the folder of the input file')
    parser.add_argument('-refPSM',
                        help='CSV file with PSM identifications for angle score calculation')
    parser.add_argument('-refFile',
                        help='MGF or mzML file for angle score calculation')
    parser.add_argument('-d',
                        help='delimiter in CSV file with PSM identifications for angle score calculation; '
                             'tab by default', default='\t')
    parser.add_argument('-cn',
                        help='column name with spectrum titles in CSV file with PSM identifications '
                             'for angle score calculation; "spectrum" by default',
                        default='spectrum')
    parser.add_argument('-start', type=float,
                        help='delay time before sample actually comes to mass spec; '
                             'if cannot be identified 0 is used',
                        default=0)
    parser.add_argument('-stop', type=float,
                        help='time of wash start; if cannot find be identified  '
                             'maximum analysis time is used')
    parser.add_argument('-charge', type=int,
                        help='max charge of precursor ions. By default, all charges are considered')
    parser.add_argument('-sf', '--separate-figures', action='store_true',
                        help='save figures as separate files', dest="sf")
    parser.add_argument('-pic',
                        help='the output figure type (png or svg for vector graphic). Default png',
                        default='png')
    parser.add_argument('-log', '--logname', help='log file name. By default, log to stdout')
    parser.add_argument('--debug', action='store_true')
    try:
        parser.add_argument('-V', '--version', action='version',
            version='%s' % (pkg_resources.require("viQC")[0], ))
    except:
        print('Version information is not available, please clone the whole directory from https://github.com/lisavetasol/viQC')
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)7s: %(asctime)s %(message)s',
                        datefmt='[%H:%M:%S]', level=(logging.INFO, logging.DEBUG)[args.debug], filename=args.logname)
    stats = StatsHandler()
    logging.getLogger().addHandler(stats)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)

    exist = True
    for infile in args.input:
        if not os.path.exists(infile):
            logging.error("Could not find the input file %s", infile)
            exist = False

    if not exist:
        return

    if len(args.input) == 1:
        infile = args.input[0]
        process_file(infile, args)
    else:
        mult_process(args.input, args)

    logging.info('There were %s error(s), %s warning(s).', stats.level2count[logging.ERROR],
                 stats.level2count[logging.WARNING])
    logging.info('Enjoy your QC!')


if __name__ == '__main__':
    main()
