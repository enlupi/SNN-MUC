import random as rnd
import numpy as np
import math
import tqdm
from scipy.stats import binom
import Params as ps 
import EventData as ed


def correct_hit(psi, dx):
    ## corrections
    slope = 4
    offset = 6#4
    linear_factor = 15

    # region near the wire
    x0 = (offset+linear_factor*math.pow(math.tan(psi),1.5)) / slope

    if dx<=x0:
        # linear relation
        return offset+linear_factor*math.pow(math.tan(psi),1.5) - slope*dx
    else:
        return  0#-5*math.tan(psi)*np.sin(math.tan(psi)*4*math.pi/15*(dx-x0))
    


# SIMULATE MUON

def generate_muon(bx0, useIdealEF=False):

    tdc0 = np.random.randint(0,30)
    angle = rnd.uniform(math.pi*1./4., math.pi*3./4.)

    m = math.tan(angle)
    xmin = -(1+(ps.NWIRES-1)//2)*ps.XCELL
    xmax = (0.5+ps.NWIRES//2)*ps.XCELL
    entry_point = rnd.uniform(xmin, xmax)
    q = - m*entry_point

    psi = math.pi/2-angle

    muon_hits = []

    wire_pattern = ''
    side_pattern = ''

    for l in range(ps.NLAYERS):
        x = (ps.pos_shift_z[l] - q)/m

        # apply smearing
        x += 0.25*np.random.randn()

        # find the wire
        nw = math.floor((x-xmin-ps.is_shifted_right[l]*0.5*ps.XCELL)/ps.XCELL)+1

        # remove hits outside the macrocell
        if nw < 1 or nw > ps.NWIRES: continue

        # find the distance from the wire
        dx = x - xmin - (ps.is_shifted_right[l]+1)*0.5*ps.XCELL - (nw-1)*ps.XCELL

        # remove hits in the I
        if abs(dx) >= (21 - 1.4/2):
            continue

        ## apply EF dishomogeneity
        # drift time
        t = abs(dx) / ps.VDRIFT

        if not useIdealEF:
            t -= 19.79*math.pow(math.tan(psi),2)
            t += correct_hit(abs(psi), abs(dx))

        # compute BX and TDC
        bx = t // ps.DURATION['bx']
        dt = t % ps.DURATION['bx']


        tdc_meas = (int(np.floor(dt/ps.DURATION['tdc'])) + tdc0) % 30
        bx_counter = int(bx0 + bx) + (int(np.floor(dt/ps.DURATION['tdc'])) + tdc0)//30

        wire_pattern += f"{l+1}{ps.WIRE_MAP[nw]}"
        side_pattern += f"{ps.SIDE_MAP[1 if dx>0 else -1]}"

        muon_hits.append({
            'layer': l + 1,
            'wire_num': nw,
            'bx': bx_counter,
            'tdc': tdc_meas,
            'label': 1 if dx>0 else -1,
            't0': bx0+tdc0/30,
            'psi': psi,
            'x0': entry_point,
            'signal': True
        })

    # get expected eq label
    """ wire_pos, side = wire_pattern, side_pattern
    wire_nums = [int(wire_pos[i])-1 for i in range(0,len(wire_pos),2)]

    if (side[:3] in ['LLL', 'RRR']) and (len(side)==4):
        wire_pos = wire_pos[2:]
        side = side[1:]
    else:
        wire_pos = wire_pos[:6]
        side = side[:3]

    selected_pattern = f'{wire_pos}-{side}'
    eq_label = find_label_from_pattern(selected_pattern)"""

    return muon_hits, f'{wire_pattern}-{side_pattern}' #,eq_label



# GENERATE PURE SIGNAL

def get_event(bx0):
    valid_event_flag = False
    num_muon_hits = 0
    while not valid_event_flag:
        muon_hits, gen_pattern = generate_muon(bx0)
        if len(muon_hits) >= (ps.NLAYERS-1):
            lat = gen_pattern.split('-')[1]
            if lat.find('LLL') == -1 and lat.find('RRR') == -1:
                valid_event_flag = True
                num_muon_hits = len(muon_hits)

    # add noise
    num_hits = len(muon_hits)

    return muon_hits, gen_pattern, num_muon_hits, num_hits


def generate_clean_evts(num_events, bx0=500):
    events_arr_no_noise = np.zeros(num_events, dtype=ed.gen_event_dtype)
    # initialize array
    events_arr_no_noise['mc']['bx'] = -1
    events_arr_no_noise['mc']['tdc'] = -1

    muon_list = []
    max_n_hit = 0
    for ev_id in tqdm.tqdm(range(num_events)):
        muon_hits, pattern, num_muon_hits, num_hits = get_event(bx0)

        events_arr_no_noise[ev_id]['id'] = ev_id
        ed.hits_to_numpy(events_arr_no_noise[ev_id], muon_hits)

        muon_list.append(muon_hits)

        if len(muon_hits) > max_n_hit:
            max_n_hit = len(muon_hits)

    return events_arr_no_noise, muon_list, max_n_hit


# GENERATE NOISY EVENTS

def noise_distributions(bx0=500, bx_oot=ps.bx_oot):
    return np.concatenate([
                np.random.triangular(bx0-bx_oot, bx0, bx0,10000),
                np.random.uniform(bx0,bx0+16,30000),
                np.random.triangular(bx0+16,bx0+16,bx0+16+bx_oot,10000)
            ])

def get_event_noise(bx0, noise_frac=0, bkg_frac=0.2):
    # simulate true event
    if np.random.rand() <= (1 - bkg_frac):
        valid_event_flag = False
        num_muon_hits = 0
        while not valid_event_flag:
            muon_hits, gen_pattern = generate_muon(bx0)
            if len(muon_hits) >= (ps.NLAYERS-1):
                lat = gen_pattern.split('-')[1]
                if lat.find('LLL') == -1 and lat.find('RRR') == -1:
                    valid_event_flag = True
                    num_muon_hits = len(muon_hits)

        t0 = muon_hits[0]['t0']
        angle = muon_hits[0]['psi']
        x0 = muon_hits[0]['x0']

        # simulate cell inefficiency
        dead_cells = binom.rvs(num_muon_hits, ps.cell_ineff) #if (num_muon_hits==NLAYERS) and (np.random.rand()<=0.2):
        for i in range(dead_cells):
            # remove one hit
            _ = muon_hits.pop(np.random.randint(len(muon_hits)))

        # add noise
        if np.random.rand() < noise_frac:
            # number of noise hits
            n_noise = np.random.choice([1,2,3,4], p=[0.45,0.3,0.2,0.05])
            for _ in range(n_noise):
                layer,wire_num = np.random.randint(1,ps.NLAYERS+1), np.random.randint(1,ps.NWIRES+1)
                #bx = bx0+np.random.randint(-10,20)
                bx = round(np.random.choice(noise_distributions(bx0), 1)[0])
                tdc = np.random.randint(0,31)
                label = 0
                # check if it can be a real hit..
                """if bx>bx0:
                    #tdrift = (bx-bx0+tdc/30)*25
                    tdrift = (bx+tdc/30-t0)*25
                    if tdrift<TDRIFT-10:
                        dx = tdrift*VDRIFT
                        wire_pos = pos_shift_x[layer-1] + (wire_num-1)*XCELL
                        x_l = wire_pos - dx
                        x_r = wire_pos + dx

                        x_th = math.tan(angle)*pos_shift_z[layer-1]+x0

                        res_l, res_r = abs(x_l-x_th), abs(x_r-x_th)

                        if (res_l <= res_r) and (res_l<=1):
                            label = -1
                        elif (res_l >= res_r) and (res_r<=1):
                            label = +1
                        else:
                            label = 0"""

                muon_hits.append({
                        'layer': layer,
                        'wire_num': wire_num,
                        'bx': bx,
                        'tdc': tdc,
                        'label': label,
                        't0': t0,
                        'psi': angle,
                        'x0': x0,
                        'signal': False
                    })

    # simulate noise
    else:
        muon_hits = []
        num_muon_hits = 0
        n_noise = np.random.choice([1,2,3,4], p=[0.45,0.3,0.2,0.05])
        for _ in range(n_noise):
            layer,wire_num = np.random.randint(1,ps.NLAYERS+1), np.random.randint(1,ps.NWIRES+1)
            t0 = bx0
            angle = -9
            x0 = -9
            bx = round(np.random.choice(noise_distributions(bx0), 1)[0])
            tdc = np.random.randint(0,31)
            label = 0
            gen_pattern = ''
            muon_hits.append({
                'layer': layer,
                'wire_num': wire_num,
                'bx': bx,
                'tdc': tdc,
                'label': label,
                't0': t0,
                'psi': angle,
                'x0': x0,
                'signal': False
            })

    # shuffle list
    rnd.shuffle(muon_hits)

    num_hits = len(muon_hits)

    return muon_hits, gen_pattern, num_muon_hits, num_hits


def generate_noisy_evts(num_events, bx0=500, noise_frac=0.1, bkg_frac=0.5):
    events_arr = np.zeros(num_events, dtype=ed.gen_event_dtype)
    # initialize array
    events_arr['mc']['bx'] = -1
    events_arr['mc']['tdc'] = -1

    muon_list = []
    max_n_hit = 0
    for ev_id in tqdm.tqdm(range(num_events)):
        muon_hits, pattern, num_muon_hits, num_hits = get_event_noise(bx0, noise_frac=noise_frac, bkg_frac=bkg_frac)

        events_arr[ev_id]['id'] = ev_id
        ed.hits_to_numpy(events_arr[ev_id], muon_hits)
        muon_list.append(muon_hits)
        
        if len(muon_hits) > max_n_hit:
            max_n_hit = len(muon_hits)

    return events_arr, muon_list, max_n_hit



# EVENT TIMING

def get_timebox(events_arr, nhits=[ps.NLAYERS-1,ps.NLAYERS]):
    curr_events = events_arr[np.isin(events_arr['n_true_hits'], nhits)]
    idxs = np.where((curr_events['mc']['bx']!=0))
    a = (curr_events['mc']['bx']+curr_events['mc']['tdc']/30)
    b = curr_events['t0']

    return (a-b[:,np.newaxis,np.newaxis])[idxs]*25