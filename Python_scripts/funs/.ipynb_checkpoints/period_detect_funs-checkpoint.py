import numpy as np
from scipy import signal
from scipy.stats import linregress

import time
from timeit import default_timer as timer
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# moving average
def moving_mean_fun(ts, lag = 3):
    if (lag % 2) == 0: # even
        q = int(lag/2)
    elif (lag % 2) == 1: # odd
        q = int((lag-1)/2)
    
    len_ts = len(ts)
    
    mt = np.full(len_ts, np.nan)
    if (lag % 2) == 1: # odd
        sum = ts[0]*(q+1)
        for i in range(q):
            sum = sum + ts[i]
        for i in range(q+1):
            sum = sum + ts[i+q] - ts[0]
            mt[i] = sum / lag
        for i in range(q+1, len_ts-q):
            sum = sum + ts[i+q] - ts[i-q-1]
            mt[i] = sum / lag
        for i in range(len_ts-q, len_ts):
            sum = sum + ts[len_ts-1] - ts[i-q-1]
            mt[i] = sum / lag
    elif (lag % 2) == 0: # even
        sum = ts[0]*(q+0.5)
        for i in range(q-1):
            sum = sum + ts[i]
        sum = sum + ts[q-1]/2
        for i in range(q+1):
            sum = sum + ts[i + q] / 2 + ts[i + q -1] / 2 - ts[0]
            mt[i] = sum / lag
        for i in range(q+1, len_ts-q):
            sum = sum + ts[i + q] / 2 + ts[i + q -1] / 2 - ts[i - q - 1] / 2 - ts[i - q] / 2
            mt[i] = sum / lag
        for i in range(len_ts-q, len_ts):
            sum = sum + ts[len_ts - 1] - ts[i - q - 1] / 2 - ts[i - q] / 2
            mt[i] = sum / lag
        
    
    return mt

# loess
def weight_extend_fun2(ts, k, lag, current_value_flag = False):
    num = len(ts)
    
    idx_set = np.arange(num)
    if (lag % 2) == 0: # even
        q = int(lag/2)
        new_lag = lag + 1
    elif (lag % 2) == 1: # odd
        q = int((lag-1)/2)
        new_lag = lag
    
    if k >= num+1 or k < -1:
        print('exceed error')
        return
    if new_lag >= num:
        window_idx_set = idx_set
        window_val_set = ts
    else:
        if k>=q and k<num-q:
            window_idx_set = idx_set[k-q:k+q+1]
            window_val_set = ts[window_idx_set]
        elif k<q:
            window_idx_set = idx_set[0:new_lag]
            window_val_set = ts[window_idx_set]
        elif k>=num-q:
            window_idx_set = idx_set[-new_lag:]
            window_val_set = ts[window_idx_set]
    
    if current_value_flag == False:
        if k >= 0 and k < num:
            window_idx_set = window_idx_set[window_idx_set != idx_set[k]]
            window_val_set = ts[window_idx_set]
    
    # cal weight
    dis_set = np.abs(window_idx_set-k)
    max_dis = np.max(dis_set)
    if new_lag > num:
        max_dis = max_dis + (new_lag - num)//2
    delta = 1.001*max_dis
    d_norm_set = dis_set/delta
    weight_set = (1-d_norm_set**3)**3
    weight_set = weight_set/np.sum(weight_set)
            
    return window_idx_set, window_val_set, weight_set

def stats_lowess_fit_fun(x, y, xval, weights):
    num = len(x)
    sum_weighted_x = (np.array(weights) * np.array(x)).sum()
    weighted_sqdev_x = (np.array(weights) * ((np.array(x) - sum_weighted_x) ** 2)).sum()
    # print(weights)
    # print(x)
    # print(sum_weighted_x)
    # print(weighted_sqdev_x)
    p = np.array(weights) * (1.0 + (xval - sum_weighted_x) * (np.array(x) - sum_weighted_x) / weighted_sqdev_x)
    y_fit = (np.array(y) * p).sum()
    return y_fit

def resid_weight_fun(y, y_fit):
    abs_resid = np.abs(y-y_fit)
    median_abs_resid = np.median(abs_resid)
    if median_abs_resid == 0:
        abs_resid[abs_resid > 0] = 1
    else:
        scale = 6.0 * median_abs_resid
        abs_resid = abs_resid/scale
        abs_resid[abs_resid > 1] = 1
    resid_weight = (1-abs_resid**2)**2
    return resid_weight

def weight_with_resid_fun(window_idx_set, weight_set, resid_weight):
    resid_weight_set = resid_weight[window_idx_set]
    new_weight_set = resid_weight_set * weight_set
    tmp = new_weight_set.copy()
    tmp[tmp>1e-12]=1
    non_zero_count = np.sum(tmp)
    if non_zero_count < 2:
        return weight_set
    else:
        sum_weight = np.sum(new_weight_set)
        new_weight_set = new_weight_set/sum_weight
        return new_weight_set
    
def loess_stats_with_resid_weight_fun(ts, lag, current_value_flag = False):
    len_ts = len(ts)
    def get_all_tuple(ts, k, lag, current_value_flag):
        # window_idx_set, window_val_set, weight_set = weight_fun(ts, k, lag, current_value_flag)
        window_idx_set, window_val_set, weight_set = weight_extend_fun2(ts, k, lag, current_value_flag)
        mtv = stats_lowess_fit_fun(window_idx_set, window_val_set, k, weight_set)
        return window_idx_set, window_val_set, weight_set, mtv
    all_tuples = [get_all_tuple(ts, k, lag, current_value_flag) for k in range(len_ts)]
    mt = np.array([tp[3] for tp in all_tuples])
    
    # consider resid weight
    resid_weight = resid_weight_fun(ts, mt)
    new_weight_set = [weight_with_resid_fun(tp[0], tp[2], resid_weight) for tp in all_tuples]
    final_mt = np.array([stats_lowess_fit_fun(tp[0], tp[1], k, new_weight_set[k]) for k, tp in enumerate(all_tuples)])
    
    return final_mt

# period detection    
def getPeriodHints(ts, num_permutations = 100, confidence_interval = 0.99, random_seed = None):
    
    num = len(ts)
    np.random.seed(random_seed)
    maxPowers = np.array([])
    fft_ks = np.array([], dtype = int)
    
    # get threshold
    for i in range(num_permutations):
        ts_p = ts.copy()
        np.random.shuffle(ts_p)
        f, Pxx_den = signal.periodogram(ts_p)
        # maxPowers = np.append(maxPowers, np.max(Pxx_den))
        maxPowers = np.append(maxPowers, np.max(Pxx_den[2:]))
    
    maxPowers_sort = np.sort(maxPowers)
    P_threshold = maxPowers_sort[int(len(maxPowers_sort)*confidence_interval)-1]
    
    # get period hints
    f, Pxx_den = signal.periodogram(ts)
    for i in range(1, len(f)):
        if Pxx_den[i] > P_threshold:
            if 2 <= 1/f[i] <= np.ceil(num/2):
                fft_ks = np.append(fft_ks, i)
    return fft_ks

def getSearchRange(k, num):
    low = 0.5*(num/(k+1) + num/k)-1
    high = 0.5*(num/(k-1) + num/k)+1
    low_int = int(np.ceil(low))
    high_int = int(high)
    low_int = np.max([low_int, 2])
    high_int = np.min([high_int, int(num/2)])
    
#     while high_int - low_int + 1 < 4:
#         if low_int > 2:
#             low_int = low_int - 1
#         if high_int < int(num/2) - 1:
#             high_int = high_int + 1
    # the range should be larger than 6
    min_range = np.min([int(num/2) - 2, 5])
    while high_int - low_int < min_range:
        if low_int > 2:
            low_int = low_int - 1
        if high_int < int(num/2) - 1:
            high_int = high_int + 1
    
    return np.arange(low_int, high_int+1)

def validation_hint(k, num, acfs):
    period = None
    # acfs begins with period = 2
    x_range = getSearchRange(k, num)
    y_range = acfs[x_range-2]
    
    range_len = len(x_range)
    min_err = float("inf")
    min_slope1 = 0
    min_slope2 = 0
    
    for i in range(1, range_len-1):
        seg1_x = x_range[:i+1]
        seg1_y = y_range[:i+1]
        seg2_x = x_range[i:]
        seg2_y = y_range[i:]

        slope1, _, _, _, stderr1 = linregress(seg1_x, seg1_y)
        slope2, _, _, _, stderr2 = linregress(seg2_x, seg2_y)
        # 20230609 include and len(seg1_x) > 2 and len(seg2_x) > 2
        if stderr1 + stderr2 < min_err and len(seg1_x) > 2 and len(seg2_x) > 2:
            min_err = stderr1 + stderr2
            min_slope1 = slope1
            min_slope2 = slope2
    
    # angle1 = np.arctan(min_slope1) / (np.pi / 2)
    # angle2 = np.arctan(min_slope2) / (np.pi / 2)
    # valid = min_slope1 > min_slope2 and not np.isclose(np.abs(angle1 - angle2), 0, atol=0.01)
    # valid = min_slope1 > 0 and min_slope2 < 0
    valid = min_slope1 > min_slope2
    if valid:
        period = np.argmax(y_range) + x_range[0]
    
    return valid, period

def validation_all_hints(fft_ks, num, acfs):
    periods = np.array([],dtype = int)
    # print(fft_ks)
    for k in fft_ks:
        valid, period = validation_hint(k, num, acfs)
        if valid:
            # print(period)
            periods = np.append(periods, period)
            
    return periods

def validation_hint_with_acfs(k, x_range, y_range):
    # x_range: period hints
    # y_range: acfs of period hints
    period = None 
    acf_max = np.max(y_range)
    range_len = len(x_range)
    min_err = float("inf")
    min_slope1 = 0
    min_slope2 = 0
    
    for i in range(1, range_len-1):
        seg1_x = x_range[:i+1]
        seg1_y = y_range[:i+1]
        seg2_x = x_range[i:]
        seg2_y = y_range[i:]

        slope1, _, _, _, stderr1 = linregress(seg1_x, seg1_y)
        slope2, _, _, _, stderr2 = linregress(seg2_x, seg2_y)
        # 20230609 include and len(seg1_x) > 2 and len(seg2_x) > 2
        if stderr1 + stderr2 < min_err and len(seg1_x) > 2 and len(seg2_x) > 2:
            min_err = stderr1 + stderr2
            min_slope1 = slope1
            min_slope2 = slope2
    
    # angle1 = np.arctan(min_slope1) / (np.pi / 2)
    # angle2 = np.arctan(min_slope2) / (np.pi / 2)
    # valid = min_slope1 > min_slope2 and not np.isclose(np.abs(angle1 - angle2), 0, atol=0.01)
    # valid = min_slope1 > 0 and min_slope2 < 0
    valid = min_slope1 > min_slope2
    if valid:
        period = np.argmax(y_range) + x_range[0]
    
    return valid, period, acf_max

# autocorrelation
def raw_acf_fun(ts, period):
    num = len(ts)
    num_acf = num - period
    ts_mean = np.mean(ts)
    c0_no = np.mean((ts - ts_mean*np.ones(num))**2)
    ch_set_no = np.array([(ts[i]-ts_mean)*(ts[i+period]-ts_mean) for i in range(num_acf)])
    ch_no = np.sum(ch_set_no) / num
    acf_no = ch_no / c0_no
    
    return {'period':period, 'acf_no':acf_no}

def parallel_raw_acf_fun(ts, period_set, ncpus = cpu_count()):
    print(f'starting computations on {ncpus} cores')

    with Pool(ncpus) as p, tqdm(total=len(period_set), miniters=len(period_set)/50) as pbar:
        res = [p.apply_async(raw_acf_fun, args=(ts, period), callback=lambda _: pbar.update(1)) for period in period_set]
        results = [p.get() for p in res]
     
    return results

def get_raw_acfs(ts):
    num = len(ts)
    period_set = np.arange(int(num/2)-1)+2
    num_acf = len(period_set)
    
    raw_results = parallel_raw_acf_fun(ts, period_set)
    raw_acf_set_no = np.full(num_acf, np.nan)
    for i in range(num_acf):
        period = raw_results[i]['period']
        raw_acf_set_no[period - 2] = raw_results[i]['acf_no']
    return period_set, raw_acf_set_no

def weight_fun(detrend, period, threshold_low, threshold_range):
    num = len(detrend)
    num_acf = num - period
    
    # acf weights
    ts_delta = np.array([detrend[i+period] - detrend[i] for i in range(num_acf)])
    abs_ts_delta = np.abs(ts_delta)
    
    mean = np.mean(abs_ts_delta)
    sigma = np.std(abs_ts_delta, ddof = 1)
    
    delta_set = (abs_ts_delta - mean - threshold_low*sigma) / (threshold_range*sigma)
    
    delta_set[delta_set < 0] = 0
    delta_set[delta_set > 1] = 1
    
    weights = (1-delta_set**2)**2
    
    weights = weights / np.mean(weights)
    
    # c0 weights
    c0_weights = np.full(num, np.nan)
    for i in range(period, num_acf):
        c0_weights[i] = (weights[i - period] + weights[i])/2
    
    end_sum = np.sum(weights[:period])
    end_sum = end_sum + np.sum(weights[-period:])
    # print(end_sum)
    
    if end_sum > 2*period:
        # print('end_sum > 2*period')
        alpha = period / end_sum
        c0_weights[:period] = (0.5 + alpha) * weights[:period]
        c0_weights[num - period:] = (0.5 + alpha) * weights[num_acf - period:]
    else:
        # print('end_sum <= 2*period')
        c0_weights[:period] = weights[:period]
        c0_weights[num - period:] = weights[num_acf - period:]
        c0_weights = c0_weights / np.mean(c0_weights)
    
    return weights, c0_weights

def weighted_acf_fun(detrend, period, threshold_low = 3, threshold_range = 3):
    # acf weights and c0 weights
    weights, c0_weights = weight_fun(detrend, period, threshold_low, threshold_range)
    
    # weighted acf
    detrend_mean = np.mean(detrend)
    # c0 = np.mean((detrend - detrend_mean*np.ones(num))**2)
    c0 = np.mean(c0_weights*(detrend - detrend_mean*np.ones(num))**2)
    ch_set = np.array([weights[i]*(detrend[i]-detrend_mean)*(detrend[i+period]-detrend_mean) for i in range(num_acf)])
    # c0 = np.mean((detrend)**2)
    # ch_set = np.array([weights[i]*(detrend[i])*(detrend[i+period]) for i in range(num_acf)])
    ch = np.sum(ch_set) / num
    acf = ch / c0
    
    return {'period':period, 'acf':acf}

def mean_acf_fun(ts, period, threshold_low = 3, threshold_range = 3):
    num = len(ts)
    num_acf = num - period
    
    # moving average to detrend
    ts_smooth = moving_mean_fun(ts, period)
    detrend = ts - ts_smooth
    
    # acf weights and c0 weights
    weights, c0_weights = weight_fun(detrend, period, threshold_low, threshold_range)
    
    # weighted acf
    detrend_mean = np.mean(detrend)
    # c0 = np.mean((detrend - detrend_mean*np.ones(num))**2)
    c0 = np.mean(c0_weights*(detrend - detrend_mean*np.ones(num))**2)
    ch_set = np.array([weights[i]*(detrend[i]-detrend_mean)*(detrend[i+period]-detrend_mean) for i in range(num_acf)])
    # c0 = np.mean((detrend)**2)
    # ch_set = np.array([weights[i]*(detrend[i])*(detrend[i+period]) for i in range(num_acf)])
    ch = np.sum(ch_set) / num
    acf = ch / c0
    
    # acf without weight
    c0_no = np.mean((detrend - detrend_mean*np.ones(num))**2)
    ch_set_no = np.array([(detrend[i]-detrend_mean)*(detrend[i+period]-detrend_mean) for i in range(num_acf)])
    ch_no = np.sum(ch_set_no) / num
    acf_no = ch_no / c0_no
        
    return {'period':period, 'acf':acf, 'acf_no':acf_no}

def loess_acf_fun(ts, period, threshold_low = 3, threshold_range = 3, current_value_flag = True):
    num = len(ts)
    num_acf = num - period
    # loess smooth to detrend
    loess_lag = 2 * period
    ts_smooth = loess_stats_with_resid_weight_fun(ts, loess_lag, current_value_flag)
    detrend = ts - ts_smooth
    
    # acf weights and c0 weights
    weights, c0_weights = weight_fun(detrend, period, threshold_low, threshold_range)
    
    # weighted acf
    detrend_mean = np.mean(detrend)
    # c0 = np.mean((detrend - detrend_mean*np.ones(num))**2)
    c0 = np.mean(c0_weights*(detrend - detrend_mean*np.ones(num))**2)
    ch_set = np.array([weights[i]*(detrend[i]-detrend_mean)*(detrend[i+period]-detrend_mean) for i in range(num_acf)])
    # c0 = np.mean((detrend)**2)
    # ch_set = np.array([weights[i]*(detrend[i])*(detrend[i+period]) for i in range(num_acf)])
    ch = np.sum(ch_set) / num
    acf = ch / c0
    
    # acf without weight
    c0_no = np.mean((detrend - detrend_mean*np.ones(num))**2)
    ch_set_no = np.array([(detrend[i]-detrend_mean)*(detrend[i+period]-detrend_mean) for i in range(num_acf)])
    ch_no = np.sum(ch_set_no) / num
    acf_no = ch_no / c0_no
    
    return {'period':period, 'acf':acf, 'acf_no':acf_no}
        
def parallel_loess_acf_fun(ts, period_set, threshold_low = 3, threshold_range = 3, current_value_flag = True, ncpus = cpu_count()):
    print(f'starting computations on {ncpus} cores')

    with Pool(ncpus) as p, tqdm(total=len(period_set), miniters=len(period_set)/50) as pbar:
        res = [p.apply_async(loess_acf_fun, args=(ts, period, threshold_low, threshold_range, current_value_flag), callback=lambda _: pbar.update(1)) for period in period_set]
        results = [p.get() for p in res]
     
    return results

def parallel_mean_acf_fun(ts, period_set, threshold_low = 3, threshold_range = 3, ncpus = cpu_count()):
    print(f'starting computations on {ncpus} cores')

    with Pool(ncpus) as p, tqdm(total=len(period_set), miniters=len(period_set)/50) as pbar:
        res = [p.apply_async(mean_acf_fun, args=(ts, period, threshold_low, threshold_range), callback=lambda _: pbar.update(1)) for period in period_set]
        results = [p.get() for p in res]
     
    return results  

def get_loess_acfs(ts, threshold_low = 3, threshold_range = 3, current_value_flag = True):
    num = len(ts)
    period_set = np.arange(int(num/2)-1)+2
    num_acf = len(period_set)
    
    # loess acf
    loess_results = parallel_loess_acf_fun(ts, period_set, threshold_low, threshold_range, current_value_flag)
    loess_acf_set = np.full(num_acf, np.nan)
    loess_acf_set_no = np.full(num_acf, np.nan)
    for i in range(num_acf):
        period = loess_results[i]['period']
        loess_acf_set[period - 2] = loess_results[i]['acf']
        loess_acf_set_no[period - 2] = loess_results[i]['acf_no']
    return period_set, loess_acf_set, loess_acf_set_no

def get_mean_acfs(ts, threshold_low = 3, threshold_range = 3):
    num = len(ts)
    period_set = np.arange(int(num/2)-1)+2
    num_acf = len(period_set)
    
    # mean acf
    mean_results = parallel_mean_acf_fun(ts, period_set, threshold_low, threshold_range)
    mean_acf_set = np.full(num_acf, np.nan)
    mean_acf_set_no = np.full(num_acf, np.nan)
    for i in range(num_acf):
        period = mean_results[i]['period']
        mean_acf_set[period - 2] = mean_results[i]['acf']
        mean_acf_set_no[period - 2] = mean_results[i]['acf_no']
    return period_set, mean_acf_set, mean_acf_set_no
