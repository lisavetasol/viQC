from __future__ import print_function
import matplotlib
matplotlib.use('agg')
import csv
import pylab
from pyteomics import mzml, mgf
try:
    import seaborn as sns
except ImportError:
    pass
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
import os
import argparse
import logging
from collections import Counter

COLORS = ["#b84c7d", "#4d8ac9", "#4bc490", "#7f63b8", "b", "g", '#edd630'] + ['k'] * 50

def read_data(name_mzml, name_calibration, extension, name_psm, delim, colname):
    logging.info('Reading data...')
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
    x = []
    y = []
    errors = set()

    for sc in mzml.read(name_mzml):
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
            indexms1.append(int(sc['id'].split(' ')[2].split('=')[1]))
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
                prec_intensity = sc['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['peak intensity']
            except (KeyError, ValueError):
                errors.add('peak intensity')
            if sc['intensity array'].size:
                intensity = sc['intensity array'].sum(dtype=float) / sc['intensity array'].size
            else:
                intensity = 0
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

    starttime_ms1 = convert_times(starttime_ms1)
    starttime_ms2 = convert_times(starttime_ms2)

    coef = 10**(len(str(int(np.mean(angle_y))))-1)

    #read_for_angle_calibration
    if name_psm is not None:
        logging.info('Reading reference files for angle score calculation...')
        d = {'mgf': mgf.IndexedMGF, 'mzml': mzml.MzML}
        reader = d[extension.lower()]
        f = reader(name_calibration)
        with open(name_psm) as fIn:
            reader = csv.DictReader(fIn, delimiter=delim)
            next(reader)
            for row in reader:
                psm = row[colname].split(' RTINSECONDS')[0]
                sc = f[psm]
                x1 = len(sc['m/z array'])
                x.append(x1)
                y1 = sc['intensity array'].sum(dtype=float) / sc['intensity array'].size
                y.append(y1)
        coef = 10**(len(str(int(np.mean(y))))-1)
        y = np.array(y) / coef
    logging.info('Reading is complete.')
    for error in errors:
        logging.warning('There was an error extracting %s information. Some figures will not be produced.', error)
    return (injtime_ms1, injtime_ms2, starttime_ms1, starttime_ms2, indexms1, charge_ms2,
        mz_ms2, angle_x, angle_y, prec_int, x, y, coef)

def ms1_ms2(starttime_ms1, starttime_ms2):
    if starttime_ms1 is None or starttime_ms2 is None:
        pylab.text(0.5, 0.5, 'Start time information missing', ha='center')
        return
    width = 0.5
    ms1 = starttime_ms1.size
    ms2 = starttime_ms2.size
    pylab.bar([0], ms1, width, alpha=1, color=COLORS[0])
    pylab.text(0, ms1/2, ms1, ha='center', fontsize=20)
    pylab.bar([1], ms2, width, alpha=1, color=COLORS[1])
    pylab.text(1, ms2/2, ms2, ha='center', fontsize=20)
    pylab.xticks(np.arange(2))
    xtick_marks = ['MS1', 'MS2']
    _, xtick_names = pylab.xticks(np.arange(2), xtick_marks)
    pylab.setp(xtick_names, rotation=0, fontsize=15)
    pylab.title('MS1/MS2')

def aqtime(starttime_ms1, starttime_ms2):
    if starttime_ms1 is None or starttime_ms2 is None:
        pylab.text(0.5, 0.5, 'Start time information missing', ha='center')
        return

    starttime_all = sorted([(x, 1) for x in starttime_ms1] + [(x, 2) for x in starttime_ms2])
    aqtime_ms1 = []
    aqtime_ms2 = []
    for i, k in enumerate(starttime_all[:-1]):
        aq = starttime_all[i+1][0] - k[0]
        if k[1] == 1:
            aqtime_ms1.append(aq)
        if k[1] == 2:
            aqtime_ms2.append(aq)
    pylab.hist(np.array(aqtime_ms1)*60, histtype='step', lw=2, density=True, label='MS1, sum=%.2f min'%sum(aqtime_ms1), color=COLORS[0])
    pylab.hist(np.array(aqtime_ms2)*60, histtype='step', lw=2, density=True, label='MS2, sum=%.2f min'%sum(aqtime_ms2), color=COLORS[1])
    pylab.legend(loc=1, fontsize=12)
    pylab.xlabel("AT, sec", fontsize=15)
    pylab.ylabel("# of spectra, normalized", fontsize=15)
    pylab.title('Acquisition time')

def it_ms1(starttime_ms1, injtime_ms1):
    if starttime_ms1 is None:
        pylab.text(0.5, 0.5, 'Start time information missing', ha='center')
        return
    if injtime_ms1 is None:
        pylab.text(0.5, 0.5, 'Injection time information missing', ha='center')
        return
    pylab.scatter(starttime_ms1, injtime_ms1, s=10, alpha=0.7, color=COLORS[0])
    pylab.legend(markerscale=2)
    pylab.xlabel("Start time, min", fontsize=15)
    pylab.ylabel("Injection time, ms", fontsize=15)
    pylab.title('MS1')

def inten_prec(starttime_ms2, start, finish, prec_int):
    if starttime_ms2 is None:
        pylab.text(0.5, 0.5, 'Start time information missing', ha='center')
        return
    ind = np.logical_and(starttime_ms2 > start, starttime_ms2 < finish)
    prec = np.log10(np.array(prec_int))[ind]
    b = np.percentile(prec, 99.9)
    pylab.hist(prec, bins=np.linspace(0, max(prec), 100), color=COLORS[0], alpha=0.5, lw=1, edgecolor='k')
    pylab.xlabel("Log10(intensity)", fontsize=15)
    pylab.xlim(np.percentile(prec, 0.1), b)
    pylab.ylabel("# of spectra", fontsize=15)
    pylab.title('Intensity of precursor ions')

def it_ms2(starttime_ms2, start, finish, injtime_ms2):
    if starttime_ms2 is None:
        pylab.text(0.5, 0.5, 'Start time information missing', ha='center')
        return
    if injtime_ms2 is None:
        pylab.text(0.5, 0.5, 'Injection time information missing', ha='center')
        return
    ind = np.logical_and(starttime_ms2 > start, starttime_ms2 < finish)
    pylab.hist(np.array(injtime_ms2)[ind], bins=np.linspace(0, max(injtime_ms2), 100),
        color=COLORS[1], alpha=0.5, lw=1, edgecolor='k')
    pylab.xlabel("Injection time, ms", fontsize=15)
    pylab.ylabel("# of spectra", fontsize=15)
    pylab.title('MS2 injection time')

def realtop(starttime_ms1, indexms1):
    if starttime_ms1 is None:
        pylab.text(0.5, 0.5, 'Start time information missing', ha='center')
        return
    pylab.scatter(starttime_ms1[:-1], np.ediff1d(indexms1)-1, alpha=0.7, s=15,
        label='TopN', color=COLORS[1])
    fit = lowess(np.ediff1d(indexms1)-1, starttime_ms1[:-1], frac=0.05, it=0)
    pylab.plot(fit[:,0],fit[:,1], "r-" ,label='Average TopN')
    pylab.ylim(0, max(np.ediff1d(indexms1)) * 1.1)
    pylab.xlim(starttime_ms1.min() * .95, starttime_ms1.max() * 1.07)
    pylab.legend(loc=1,markerscale=2,fontsize=15)
    pylab.xlabel('RT, min',fontsize=15)
    pylab.ylabel('TopN',fontsize=15)
    pylab.title('Real TopN',fontsize=15)

def charge(maxcharge, charge_ms2, starttime_ms2, mz_ms2):
    if starttime_ms2 is None:
        pylab.text(0.5, 0.5, 'Start time information missing', ha='center')
        return
    for i in range(0, maxcharge+1):
        charge_ms2 = np.array(charge_ms2)
        starttime_ms2 = np.array(starttime_ms2)
        mz_ms2 = np.array(mz_ms2)
        msk = [k==i for k in charge_ms2]
        pylab.scatter(np.array(starttime_ms2)[msk], np.array(mz_ms2)[msk],
            s=20, alpha=0.5, color=COLORS[i-1], label='charge_%s+(%i scans)'%(i, len(np.array(mz_ms2)[msk])))
        pylab.title('Precursor ions', fontsize=15)
        pylab.legend(markerscale=1.5, fontsize=15)
        pylab.xlabel('RT, min', fontsize=15)
        pylab.ylabel('m/z', fontsize=15)

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
    for i in np.arange(0.001, 1, 0.001):
        dif = []
        for k in zip(x_mod, y_mod):
            dif.append(i * (k[0]-per_1_x) + per_1_y - k[1])
        under_line = sum(x > 0 for x in dif)
        if float(under_line) / len(x_mod) > 0.01:
            under = i
            break
        else:
            continue
    for j in np.arange(1, 0.001, -0.01):
        dif = []
        for k in zip(x_mod, y_mod):
            dif.append(j * (k[0]-per_1_x) + per_1_y - k[1])
        above_line = sum(x < 0 for x in dif)
        if float(above_line) / len(x_mod) > 0.01:
            above = j
            break
        else:
            continue
    a = 0
    angle_y_mod = np.array(angle_y) / coef
    for x1, y1 in zip(angle_x, angle_y_mod):
        if y1 - (x1-per_1_x) * above < per_1_y and y1 - (x1-per_1_x) * under > per_1_y:
            a += 1
    per = 100 * float(a) / len(angle_x)
    return under, above, per_1_x, per_1_y, per, angle_y_mod, coef

def angle(under, above, per_1_x, per_1_y, per, angle_x, angle_y_mod, coef):
    a = np.arange(per_1_x, 1000, 100.)
    b = under * (a-per_1_x) + per_1_y
    pylab.plot(a, b, color='black')
    c = above * (a-per_1_x) + per_1_y
    pylab.plot(a, c, color='black')
    pylab.scatter(angle_x, angle_y_mod, alpha=0.5, color=COLORS[1], label='score=%.2f' % per)
    pylab.xlim(0, np.percentile(angle_x, 99.5))
    pylab.ylim(0, np.percentile(angle_y_mod, 99.5))
    pylab.xlabel('# peaks', fontsize=15)
    pylab.legend(fontsize=15, loc=1)
    pylab.ylabel('avg intensity, 10^%i'%np.log10(coef), fontsize=15)

def inten_number_peaks_ms1(angle_x, angle_y):
    coef = 10**(len(str(int(np.mean(angle_y))))-1)
    angle_y_mod = np.array(angle_y) / coef
    pylab.scatter(angle_x, angle_y_mod, alpha=0.5, color=COLORS[1])
    pylab.xlim(0, np.percentile(angle_x, 99.5))
    pylab.ylim(0, np.percentile(angle_y_mod, 99.5))
    pylab.xlabel('# peaks', fontsize=15)
    pylab.ylabel('avg intensity, 10^%i' % np.log10(coef), fontsize=15)
    pylab.title('MS/MS', fontsize=15)

def process_file(name_mzml, args):
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
    start = args.start
    delim = str(args.d)
    colname = str(args.cn)

    if (args.refPSM is not None) and (os.path.split(name_calibr)[1].split('.')[0] not in name_psm):
        logging.warning('File names for angle score calibration don\'t match!')

    injtime_ms1, injtime_ms2, starttime_ms1, starttime_ms2, indexms1, charge_ms2, \
        mz_ms2, angle_x, angle_y, prec_int, x, y, coef = read_data(
            name_mzml, name_calibr, extension, name_psm, delim, colname)

    if args.stop is None:
        finish = starttime_ms1.max() if starttime_ms1 is not None else None
    else:
        finish = args.stop

    if args.charge is None:
        maxcharge = max(charge_ms2)
        logging.info('Maximum charge in file: %s', maxcharge)
    else:
        maxcharge = args.charge

    # for pretty pictures
    pylab.style.use('seaborn-whitegrid')
    pylab.rcParams['ytick.labelsize'] = 15
    pylab.rcParams['xtick.labelsize'] = 15
    pylab.rcParams['axes.titlesize']  = 15
    pylab.rcParams['axes.titlesize']  = 15
    pylab.rcParams['font.size']  = 15

    pylab.figure(figsize=(15, 40))
    pylab.subplot2grid((6, 2), (0, 0))
    ms1_ms2(starttime_ms1, starttime_ms2)
    pylab.subplot2grid((6, 2), (0, 1))
    aqtime(starttime_ms1, starttime_ms2)
    pylab.subplot2grid((6, 2), (1, 0), colspan=2)
    it_ms1(starttime_ms1, injtime_ms1)
    pylab.subplot2grid((6, 2), (2, 0), colspan=2)
    inten_prec(starttime_ms2, start, finish, prec_int)
    pylab.subplot2grid((6, 2), (3, 0), colspan=2)
    it_ms2(starttime_ms2, start, finish, injtime_ms2)
    pylab.subplot2grid((6, 2), (4, 0), colspan=2)
    charge(maxcharge, charge_ms2, starttime_ms2, mz_ms2)
    pylab.subplot2grid((6, 2), (5, 0))
    realtop(starttime_ms1, indexms1)
    if name_psm is not None:
        under, above, per_1_x, per_1_y, per, angle_y_mod, coef = angle_calculation(x, y, angle_x, angle_y, coef)
        pylab.subplot2grid((6, 2), (5, 1))
        angle(under, above, per_1_x, per_1_y, per, angle_x, angle_y_mod, coef)
    else:
        pylab.subplot2grid((6, 2), (5, 1))
        inten_number_peaks_ms1(angle_x, angle_y)
    logging.info("Saving results...")

    outname = os.path.join(output, name + '_viQC.png')
    pylab.savefig(outname)
    logging.info('Output figure saved to %s', outname)


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
    parser.add_argument('input', nargs='+', help='mzML file with path(s)')
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
        'used for precursor intensity and injection time (MS/MS) calculation; 0 by default',
        default=0)
    parser.add_argument('-stop', type=float,
        help='time of wash start; used for precursor intensity and '
        'injection time (MS/MS) calculation. By default, maximum analysis time')
    parser.add_argument('-charge', type=int,
        help='max charge of precursor ions. By default, all charges are considered')
    parser.add_argument('-log', '--logname',
        help='log file name. By default, log to stdout')
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)7s: %(asctime)s %(message)s',
            datefmt='[%H:%M:%S]', level=logging.INFO, filename=args.logname)
    stats = StatsHandler()
    logging.getLogger().addHandler(stats)

    for infile in args.input:
        if os.path.exists(infile):
            process_file(infile, args)
        else:
            logging.error("Could not find the input file %s", infile)
    logging.info('There were %s error(s), %s warning(s).', stats.level2count[logging.ERROR], stats.level2count[logging.WARNING])
    logging.info('Enjoy your QC!')


if __name__ == '__main__':
    main()
