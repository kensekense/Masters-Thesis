def after (trace, x):

    assert type(trace) == list
    assert type(x) == list

    #find first instance
    for pos in range(len(trace)):
        if trace[pos] in x:
            return [(pos,len(trace)-1)]

def before (trace, x):

    assert type(trace) == list
    assert type(x) == list

    #find first instance
    for pos in range(len(trace)):
        if trace[pos] in x:
            return [(0,pos)]

#return the sub-sequences for scope x UNTIL y
def until (trace, x, y):

    assert type(trace) == list
    assert type(x) == list
    assert type(y) == list

    #find all instances of until
    s = -1
    e = -1
    sol = []

    for pos in range(len(trace)):

        if trace[pos] in x and s == -1:
            s = pos

        if trace[pos] in y and s != -1:
            e = pos
            sol.append((s,e))
            s = -1
            e = -1

    #add if don't find anyways
    if s != -1 and pos == len(trace)-1:
        sol.append((s,len(trace)-1))

    return sol

#return the sub-sequences for scope x UNTIL y
def between (trace, x, y):

    assert type(trace) == list
    assert type(x) == list
    assert type(y) == list

    #find all instances of until
    s = -1
    e = -1
    sol = []

    for pos in range(len(trace)):

        if trace[pos] in x and s == -1:
            s = pos

        if trace[pos] in y and s != -1:
            e = pos
            sol.append((s,e))
            s = -1
            e = -1

    return sol


#PROPERTY FUNCTIONS

def absence (trace, scope, condition):

    assert type(trace) == list
    assert type(scope) == list
    assert type(condition) == list

    sol = []

    for subsequence in scope:
        for item in condition:
            if item in trace[subsequence[0]:subsequence[1]+1]:
                sol.append(False)
            else:
                sol.append(True)

    return sol

def existence (trace, scope, condition):

    assert type(trace) == list
    assert type(scope) == list
    assert type(condition) == list

    sol = []

    for subsequence in scope:
        for item in condition:
            if item in trace[subsequence[0]:subsequence[1]+1]:
                sol.append(True)
            else:
                sol.append(False)

    return sol

def universality (trace, scope, condition):

    assert type(trace) == list
    assert type(scope) == list
    assert type(condition) == list

    sol = []

    for subsequence in scope:
        flag = True
        for item in trace[subsequence[0]:subsequence[1]+1]:
            if item not in condition:
                sol.append(False)
                flag = False
                break
        if flag:
            sol.append(True)

    return sol

def response (trace, scope, condition1, condition2):

    assert type(trace) == list
    assert type(scope) == list
    assert type(condition1) == list
    assert type(condition2) == list

    sol = []

    for subsequence in scope:
        c1 = -1
        c2 = -1
        ind = 0
        for item in trace[subsequence[0]: subsequence[1]+1]:
            if item in condition1 and c1 == -1:
                c1 = ind
            if item in condition2 and c2 == -1:
                c2 = ind
            ind += 1
        if c1 > c2:
            sol.append(True)
        else:
            sol.append(False)

    return sol

def precedence (trace, scope, condition1, condition2):

    assert type(trace) == list
    assert type(scope) == list
    assert type(condition1) == list
    assert type(condition2) == list

    sol = []

    for subsequence in scope:
        c1 = -1
        c2 = -1
        ind = 0
        for item in trace[subsequence[0]: subsequence[1]+1]:
            if item in condition1 and c1 == -1:
                c1 = ind
            if item in condition2 and c2 == -1:
                c2 = ind
            ind += 1
        if (c1 < c2) or (c1 == -1 or c2 == -1):
            sol.append(True)
        else:
            sol.append(False)

    return sol


#UNTIL-N DEFINITION

#helper function to count the occurrences of violations to the label
def count_label (sub, labels):
    assert type(sub) == list
    assert type(labels) == list

    count = 0
    for item in sub:
        if item not in labels:
            count += 1
    return count

#until-N definition
def until_N (trace, x, y, N):

    assert type(trace) == list
    assert type(x) == list
    assert type(y) == list

    sol = []
    current = N
    s = -1
    e = -1
    for i in range(len(trace)):

        if (trace[i] in x) and e == -1 and s == -1: #finding first instance of x
            s = i

        if (s != -1) and (trace[i] not in x) and (trace[i] not in y): #started count and violates Until

            if current <= 0: #no more N to give
                s = -1
                e = -1
                current = N
                continue #search for next

            else: #more N to give, decrement
                current -= 1

        if s != -1 and (trace[i] in y): #found instance of y and x
            e = i
            sol.append((s,e, count_label(trace[s:e],x), e-s)) #append starting and ending index, with number of appearances of x
            s = -1
            e = -1

    return sol
