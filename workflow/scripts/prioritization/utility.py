def format_output(field):
    if field == None:
        return '.'
    else:
        return field

# deter
def find_all_lowercase(seq):
    matches = []
    for i in range(len(seq)):
        if seq[i].islower():
            matches.append(i)
    return matches
