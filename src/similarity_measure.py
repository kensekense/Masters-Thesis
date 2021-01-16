import numpy as np
#STEPS for Substitution Score

#STEP 1: Define the symbols in the list of traces
def define_symbols (traces):
    assert type(traces) == list
    symbols = []
    for item in traces:
        symbols.append(set(item))
    x = symbols[0]
    for i in range(len(symbols)):
        x = x.union(symbols[i])

    return list(x)


#STEP 2: Define the set of all 3-grams in the logs and their frequencies
def three_grams (traces):
    assert type(traces) == list
    g3 = []
    g3_freq = {}
    for trace in traces:
        for i in range(len(trace)-2):
            g3.append(", ".join(list(trace[i:i+3])))
            try:
                g3_freq[", ".join(list(trace[i:i+3]))] += 1
            except:
                g3_freq[", ".join(list(trace[i:i+3]))] = 1
    return list(set(g3)), g3_freq


#STEP 3: Define the context for symbol a
def define_context(grams):

    assert type(grams) == list

    context = {}
    for gram in grams:
        x,a,y = gram.split(", ")
        try:
            context[a].append("{0}, {1}".format(x,y))
        except:
            context[a] = []
            context[a].append("{0}, {1}".format(x,y))

    #clear dups
    for k in list(context.keys()):
        context[k] = list(set(context[k]))

    return context


#STEP 4: define pairs of context
def context_pairs (context):

    assert type(context) == dict

    context_pairs = {}
    for a in list(context.keys()):
        for b in list(context.keys()):
            if a != b:
                context_pairs["{0}, {1}".format(a, b)] = list(set(context[a]).intersection(set(context[b])))

    return context_pairs


#STEP 5: define co-occurrence combinations
def define_cooccurrence(symbols, context_pairs, gram_freq):

    assert type(context_pairs) == dict
    assert type(gram_freq) == dict
    assert type(symbols) == list

    co_occur = {}
    for k in list(context_pairs.keys()):
        for item in context_pairs[k]:
            for a in symbols:
                for b in symbols:
                    x,y = item.split(", ")[0], item.split(", ")[1]
                    if a == b:
                        try:
                            n = gram_freq["{0}, {1}, {2}".format(x,a,y)]
                            co_occur["{0}, {1}({2}, {3})".format(x,y,a,b)] = (n*(n-1))/2
                        except:
                            co_occur["{0}, {1}({2}, {3})".format(x,y,a,b)] = 0.0

                    elif a != b:
                        try:
                            n_i = gram_freq["{0}, {1}, {2}".format(x,a,y)]
                            n_j = gram_freq["{0}, {1}, {2}".format(x,b,y)]
                            co_occur["{0}, {1}({2}, {3})".format(x,y,a,b)] = n_i*n_j
                        except:
                            co_occur["{0}, {1}({2}, {3})".format(x,y,a,b)] = 0.0

    return co_occur


#STEP 6: Define the count of co-occurrences for symbols a,b for all contexts
def co_occur_combos(symbols, con_pairs, co_occurs):
    assert type(symbols) == list
    assert type(con_pairs) == dict
    assert type(co_occurs) == dict

    co_occur_combos = {}
    for a in symbols:
        for b in symbols:
            total = 0.0
            for k in list(con_pairs.keys()):
                for item in con_pairs[k]:
                    total += co_occurs["{0}({1}, {2})".format(item,a,b)]
            co_occur_combos["{0}, {1}".format(a,b)] = total

    return co_occur_combos


#STEP 7: Define norm on the count of co-occur combos
def define_norm (co_combos):
    assert type(co_combos) == dict
    norm = 0.0
    for k in list(co_combos.keys()):
        norm += co_combos[k]

    return norm


#STEP 8: Define matrix M over A x A
def define_matrix (symbols, co_combos, norm):
    assert type(symbols) == list
    assert type(co_combos) == dict
    assert type(norm) == float

    mat_M = {}
    for a in symbols:
        for b in symbols:
            mat_M["{0}, {1}".format(a,b)] = co_combos["{0}, {1}".format(a,b)]/norm

    return mat_M


#STEP 9: Define the probability of occurrence
def prob_occur (symbols, mat_M):
    assert type(symbols) == list
    assert type(mat_M) == dict

    p = {}
    for a in symbols:
        total = 0
        for b in symbols:
            if a != b:
                total += mat_M["{0}, {1}".format(a,b)]
        total += mat_M["{0}, {1}".format(a,a)]
        p["{0}".format(a)] = total

    return p


#STEP 10: Define the expected values
def exp_val (symbols, prob):
    assert type(symbols) == list
    assert type(prob) == dict

    e_val = {}
    for a in symbols:
        for b in symbols:
            if a == b:
                e_val["{0}, {1}".format(a,b)] = prob["{0}".format(a)]**2
            else:
                e_val["{0}, {1}".format(a,b)] = 2*prob["{0}".format(a)]*prob["{0}".format(b)]

    return e_val


#STEP 11: Define the function for substitution scores
def sub_scores (traces):
    assert type(traces) == list

    symbols = define_symbols(traces)
    three_gs, three_gs_freq = three_grams(traces)
    cons = define_context(three_gs)
    con_pairs = context_pairs(cons)
    co_occurs = define_cooccurrence(symbols, con_pairs, three_gs_freq)
    co_combos = co_occur_combos(symbols, con_pairs, co_occurs)
    norm = define_norm(co_combos)
    matM = define_matrix(symbols, co_combos, norm)
    probs = prob_occur(symbols, matM)
    e_val = exp_val(symbols, probs)

    sub_costs = {}
    for a in symbols:
        for b in symbols:
            if a!=b:
                try:
                    sub_costs["{0}, {1}".format(a,b)] = np.log2(matM["{0}, {1}".format(a,b)]/e_val["{0}, {1}".format(a,b)])
                except:
                    sub_costs["{0}, {1}".format(a,b)] = -1000

    return sub_costs

#STEPS 1-3 are the same for Insertion Score

#STEP 4: Define occurence of 3-gram counts
def occ_count (symbols, cons, grams, gfreq):
    assert type(symbols) == list
    assert type(grams) == list
    assert type(cons) == dict

    o_counts = {}
    for a in list(cons.keys()):
        for pair in cons[a]:
            x = pair.split(", ")[0]
            y = pair.split(", ")[1]
            o_counts["{0}, {1}({2})".format(x,y,a)] = gfreq["{0}, {1}, {2}".format(x,a,y)]

    return o_counts


#STEP 5: define countRgivenL
def countRgL (symbols, ocounts):
    assert type(symbols) == list
    assert type(ocounts) == dict

    rgl_counts = {}

    for a in symbols:
        for x in symbols:
            #if a !=x:
            total = 0
            for k in list(ocounts.keys()):
                if k.split("(")[0].split(", ")[0] == x and k.split("(")[1] == "{0})".format(a):
                    total += ocounts[k]
            rgl_counts["{0}/{1}".format(a,x)] = total

    return rgl_counts


#STEP 6: define norm(a)
def rgl_norm (symbols, rgl_counts):
    assert type(symbols) == list
    assert type(rgl_counts) == dict

    rgl_norms = {}

    for a in symbols:
        total = 0
        for x in symbols:
            #if a !=x:
            total += rgl_counts["{0}/{1}".format(a,x)]
        rgl_norms["{0}".format(a)] = total

    return rgl_norms


#STEP 7: define the probability of all symbols
def rgl_prob (trace):
    assert type(trace) == list

    p = {}
    for item in trace:
        for a in item:
            try:
                p["{0}".format(a)] += 1
            except:
                p["{0}".format(a)] = 1

    tot_len = 0
    for item in trace:
        tot_len += len(item)

    for k in list(p.keys()):
        p[k] = p[k]/tot_len

    return p


#STEP 8: define rglNorm
def normed_counts (symbols, rgl, norms):
    assert type(symbols) == list
    assert type(rgl) == dict
    assert type(norms) == dict

    normed_rgls = {}

    for a in symbols:
        for b in symbols:
            normed_rgls["{0}/{1}".format(a,b)] = rgl["{0}/{1}".format(a,b)]/norms["{0}".format(a)]

    return normed_rgls


#STEP 9: Define the function for insertion score
def insert_scores (traces):
    assert type(traces) == list

    symbols = define_symbols(traces)
    grams, freq = three_grams(traces)
    cons = define_context(grams)
    oc = occ_count(symbols, cons, grams, freq)
    rgl = countRgL(symbols, oc)
    norms = rgl_norm(symbols, rgl)
    probs = rgl_prob(traces)
    norm_rgls = normed_counts(symbols ,rgl, norms)

    scores = {}
    for a in symbols:
        for b in symbols:
            scores["{0}/{1}".format(a,b)] = np.log2(norm_rgls["{0}/{1}".format(a,b)]/probs["{0}".format(a)]*probs["{0}".format(b)])

    #replace -inf
    for k in list(scores.keys()):
        if scores[k] == -np.inf:
            scores[k] = -1000

    return scores

#function definition for calculating similarity
def calc_similarity(trace1, trace2, sub_cost, ins_cost, probs):

    assert type(trace1) == type(trace2) == list

    #pad traces
    trace1 = ["_"] + trace1
    trace2 = ["_"] + trace2

    #set shorter one as tr1
    if len(trace1) > len(trace2):
        copy = trace1
        trace1 = trace2
        trace2 = copy

    M = len(trace1)
    N = len(trace2)
    sim_table = np.zeros((M,N)) #establish table
    s_score = sub_cost #get substitution score
    ins_score = ins_cost #get insertion score
    p = probs #get probabilities

    #fill table, horizontal -> vertical
    for i in range(M):
        for j in range(N):

            #original fill horizontal
            if i == 0:
                if j == 0: #first fill
                    sim_table[i][j] = 1000
                elif j == 1: #first insert
                    sim_table[i][j] = p["{0}".format(trace2[j])]
                else: #rest fill, base insert scores
                    sim_table[i][j] = ins_score["{0}/{1}".format(trace2[j], trace2[j-1])] + sim_table[i][j-1]

            #original fill vertical
            elif j == 0:
                if i == 0:#first fill
                    sim_table[i][j] = 1000
                elif i == 1:
                    sim_table[i][j] = p["{0}".format(trace1[i])]
                else: #rest fill, base is the opposite of insert scores
                    sim_table[i][j] = -1*ins_score["{0}/{1}".format(trace1[i], trace1[i-1])] + sim_table[i-1][j]

            elif trace1[i] == trace2[j]: #no changes
                sim_table[i][j] = sim_table[i-1][j-1]

            else: #substitution, insertion or deletion

                #determine the min
                op = np.argmax([sim_table[i-1][j], sim_table[i][j-1], sim_table[i-1][j-1]]) #in order, removal, insertion, substitution
                if op == 0:
                    sim_table[i][j] = -1 + sim_table[i-1][j]#-1*ins_score["{0}/{1}".format(trace2[j],trace1[i])] + sim_table[i-1][j] #removal
                elif op == 1:
                    sim_table[i][j] = ins_score["{0}/{1}".format(trace2[j],trace1[i])] + sim_table[i][j-1] #insertion
                elif op == 2:
                    sim_table[i][j] = s_score["{0}, {1}".format(trace1[i],trace2[j])] + sim_table[i-1][j-1] #substitution

    return sim_table[i][j] #final score
