import numpy as np
from Params import NLAYERS
from Params import NWIRES

hit_dtype = np.dtype([
    ('bx', np.int16),
    ('tdc', np.int16),
    ('label', np.int8)
])

gen_event_dtype = np.dtype([
    ('id', np.int16),
    ('mc', hit_dtype, (NLAYERS,NWIRES)),
    ('t0', np.float32),
    ('angle', np.float16),
    ('x0', np.float16),
    ('n_true_hits', np.int8),
    ('n_hits', np.int8),
    ('signal', np.bool_)
])



# fill event structure
def hits_to_numpy(event, muon_hits, signal_type):
    for hit in muon_hits:
        layer, wire = hit['layer']-1, hit['wire_num']-1

        curr_mc = event['mc'][::-1]
        curr_mc['bx'][layer, wire] = hit['bx']
        curr_mc['tdc'][layer, wire] = hit['tdc']
        curr_mc['label'][layer, wire] = hit['label']

    event['t0'] = muon_hits[0]['t0']
    event['angle'] = muon_hits[0]['psi']
    event['x0'] = muon_hits[0]['x0']

    event['n_hits'] = np.count_nonzero(event['mc']['bx']!=-1)
    event['n_true_hits'] = np.count_nonzero(event['mc']['label']!=0)
    event['signal'] = signal_type


# get hits list from event
def numpy_to_hits(event):
    muon_hits = []

    # find hits in macrocell
    curr_mc = event['mc'][::-1]

    hits_idx = np.where(curr_mc['bx']!=-1)
    for l, w in zip(*hits_idx):
        muon_hits.append({
            'layer': l+1,
            'wire_num': w+1,
            'bx': curr_mc['bx'][l, w],
            'tdc': curr_mc['tdc'][l, w],
            'label': curr_mc['label'][l,w],
            't0': event['t0'],
            'psi': event['angle'],
            'x0': event['x0']
        })

    return muon_hits