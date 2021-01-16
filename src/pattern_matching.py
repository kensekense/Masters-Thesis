import numpy as np

#PREFIX-SPAN variant
def create_init_prefix(trace, min_supp):

    assert type(trace) == list
    assert type(min_supp) == int

    labels = list(set(trace)) #create the list of symbols

    #create initial prefix

    #chunk trace for sequences
    sequences = []
    for i in range(0, len(trace), 30):
        sequences.append((i,i+30)) #REFERENCE: sequences[i] = (start, end)

    #create bucket for symbols
    bucket = {}
    for label in labels:
        bucket["{0}".format(label)] = 0

    #check for support
    for prefix in list(bucket.keys()):

        #check all sequences
        for chunk in sequences:
            if prefix in trace[chunk[0]:chunk[1]]: #in the sequence
                bucket["{0}".format(prefix)] += 1

    #drop lower than min_supp
    for key in list(bucket.keys()):
        if bucket[key] < min_supp:
            bucket.pop(key)

    #create initial projections
    projections = {}

    #initialize the list
    for key in list(bucket.keys()):
        projections[key] = []

    #populate by searching each sequence
    for item in list(projections.keys()):
        for seq_id in sequences: #should be in the form (start, end)
            #look to replace the start index
            sequence = trace[seq_id[0]:seq_id[1]]
            for i in range(0, len(sequence), len(item)):
                if item == sequence[i]:
                    projections[item].append((seq_id[0]+i, seq_id[1]))
                    break

    #drop entry if projections <= 1
    for key in list(projections.keys()):
        if len(projections[key]) <= 1:
            projections.pop(key)

    return bucket, projections

def generate_projections(bucket, trace):

    assert type(bucket) == dict #should be in the form {'use.social.convention': 23}
    assert type(trace) == list #should be in the form ['misc', ..., 'use.social.convention']

    projections = {}
    sequences = []
    for i in range(0, len(trace), 30):
        sequences.append((i,i+30)) #REFERENCE: sequences[i] = (start, end)

    #generate a projection for each bucket item

    #initialize the list
    for key in list(bucket.keys()):
        projections[key] = []

    #populate by searching each sequence
    for item in list(projections.keys()):
        pref = item.split(", ")
        for seq_id in sequences: #should be in the form (start, end)
            #look to replace the start index
            sequence = trace[seq_id[0]:seq_id[1]]
            for i in range(0, len(sequence), len(pref)):
                if pref == sequence[i:i+len(pref)]:
                    projections["{0}".format((", ").join(pref))].append((seq_id[0]+i+len(pref)-1, seq_id[1]))
                    break

    #drop if not enough
    for k in list(projections.keys()):
        if len(projections[k]) <= 1:
            projections.pop(k)

    return projections

def compute_and_drop(projections, trace, min_supp):

    assert type(projections) == dict #should be in the form ("alpha": [(start, end), ... (startN, endN)])
    assert type(trace) == list

    labels = list(set(trace))
    proj_list = []

    #for every projection, count the i+1 patterns
    for proj in list(projections.keys()):
        #append from the labels to initialize
        iplus1_count = {}
        for L in labels:
            prefix = proj.split(", ") + [L]
            iplus1_count["{0}".format((", ").join(prefix))] = 0

            #count the labels
            sequences = projections[proj]
            for seq_id in sequences: #should be in the form (start, end)
                sequence = trace[seq_id[0]:seq_id[1]]
                for i in range(0, len(sequence), len(prefix)):
                    if prefix == sequence[i:i+len(prefix)]:
                        iplus1_count["{0}".format((", ").join(prefix))] += 1
                        break

            #drop
            for k in list(iplus1_count.keys()):
                if iplus1_count[k] < min_supp:
                    iplus1_count.pop(k)

        #output to proj_list
        proj_list.append(iplus1_count)

    return proj_list

#function to loop generate_proj with compute_and_drop
def mine_pattern(trace, min_supp):

    assert type(trace) == list
    assert type(min_supp) == int

    patts = {}
    docket = []
    doc_ind = 0

    #create initial bucket and projections
    init_buq, init_proj = create_init_prefix(trace, min_supp)

    docket.append((init_buq, init_proj))

    #loop until empty
    while(len(docket) > 0):

        curr_buq= docket[0][0] #current bucket
        curr_proj = docket[0][1] #current projections dict for entire bucket

        #compute and drop
        buckets = compute_and_drop(curr_proj, trace, min_supp)

        #for each bucket in buckets
        for bucket in buckets:
            #compute projection for bucket
            next_proj = generate_projections(bucket, trace)
            docket.append((bucket, next_proj))

            #add as patt
            for item in bucket:
                try:
                    patts[item] += bucket[item]
                except:
                    patts[item] = bucket[item]

        docket.pop(0) #finish item and move on

    return patts

def get_ranked_patterns (eventlog, top):

    assert type(eventlog) == list
    assert type(top) == int


    #get ITF counts for all traces
    itf = {}
    for trace in eventlog:
        pattern_count = mine_pattern(trace, 2)
        for k in list(pattern_count.keys()):
            try:
                itf[k] += pattern_count[k]
            except:
                itf[k] = pattern_count[k]

    #apply log(N/itf-count)
    for k in list(itf.keys()):
        itf[k] = np.log(len(itf)/itf[k])

    #calc total PF ITF for each trace patterns
    pf_itf = {}
    for trace in eventlog:
        tf = mine_pattern(trace, 2)
        for k in list(tf.keys()):
            pf_itf[k] = tf[k]*itf[k]

    pf_itf = {k: v for k, v in sorted(pf_itf.items(), key=lambda item: item[1], reverse=True)}

    ranked_patterns = []
    i = 0
    for k in list(pf_itf.keys()):
        ranked_patterns.append((k, pf_itf[k]))
        i += 1

        if i == top:
            break

    return ranked_patterns

#modify weights
def augment_weights(trace, trace_weights, rankings):

    assert type(trace_weights) == list #should be in the form [weight1, weight2, ..., weightN]
    assert type(rankings) == list #should be in the form [(pattern1, ranking_weight1), ... (patternN, ranking_weightN)]

    weights = trace_weights.copy()
    normalize_val = np.sum(weights)

    markers = {}
    for item in rankings:
        markers[item[0]] = []

    #search in trace for the top 10 rankings
    for item in rankings:
        pattern = item[0].split(", ")
        for i in range(0, len(trace), len(pattern)):
            if pattern == trace[i:i+len(pattern)]:
                markers[item[0]].append((item[1],i, i+len(pattern)-1)) #append (start, end)

    #drop empty
    for k in list(markers.keys()):
        if len(markers[k]) < 1:
            markers.pop(k)

    #use markers to augment weights
    for k in list(markers.keys()): #search all patterns
        for segment in markers[k]: #for each time pattern occurs
            rank_weight = segment[0]/normalize_val
            start = segment[1]
            end = segment[2]
            #augment
            for q in range(start, end+1):
                weights[q] += rank_weight

    return weights

#WEIGHTING FUNCTION
def calc_weights (trace, sub_convos):

    assert type(trace) == list #takes a list of traces, i.e: event log

    tw_id = {}
    center = len(trace)//2
    tf = {}
    idf = {}
    N = 0 #track number of sub_convos

    #instantiate tf and idf values
    for i in range(len(trace)):
        tf[i] = 0
        idf[i] = 0

    #find non-zero from sub_convos table to fill tf and idf table
    for key in list(sub_convos.keys()):
        if sub_convos[key][0] != 0:
            N += 1
            for entry in sub_convos[key][1]: #for each entry
                s = entry[0]
                e = entry[1]
                label1, label2 = key.split(" -> ")

                for i in range(s,e):
                    if trace[i] == label1 or trace[i] == label2:
                        tf[i] += 1
                    else:
                        idf[i] += 1

    #apply log onto idf vals
    for i in range(len(idf)):
        try:
            if idf[i] == 0:
                idf[i] = np.log2(N)
            else:
                idf[i] = np.log2(N/idf[i])
        except:
            idf[i] = np.log2(N)

    #calculate weights
    for i in range(len(trace)):
        if center == i:
            tw_id[i] = (tf[i]*idf[i]) #center
        else:
            tw_id[i] = (tf[i]*idf[i]) #w_i + tf-idf(a_i)

    #convert to list
    trace_weights = []
    for k in list(tw_id.keys()):
        trace_weights.append(tw_id[k])

    return trace_weights
