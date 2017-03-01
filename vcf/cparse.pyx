from model import _Call

cdef _map(func, iterable, bad=['.', '']):
    '''``map``, but make bad values None.'''
    return [func(x) if x not in bad else None
            for x in iterable]

INTEGER = 'Integer'
FLOAT = 'Float'
NUMERIC = 'Numeric'

def _parse_filter(filt_str):
    '''Parse the FILTER field of a VCF entry into a Python list

    NOTE: this method has a python equivalent and care must be taken
    to keep the two methods equivalent
    '''
    if filt_str == '.':
        return None
    elif filt_str == 'PASS':
        return []
    else:
        return filt_str.split(';')

def parse_samples(
        list names, list samples, samp_fmt,
        list samp_fmt_types, list samp_fmt_nums, site):

    cdef char *name, *fmt, *entry_type, *sample
    cdef int i, j
    cdef list samp_data = []
    cdef dict sampdict
    cdef list sampvals
    n_samples = len(samples)
    n_formats = len(samp_fmt._fields)

    for i in range(n_samples):
        name = names[i]
        sample = samples[i]

        # parse the data for this sample
        sampdat = [None] * n_formats

        sampvals = sample.split(':')

        for j in range(n_formats):
            if j >= len(sampvals):
                break
            vals = sampvals[j]

            # short circuit the most common
            if samp_fmt._fields[j] == 'GT':
                sampdat[j] = vals
                continue
            # genotype filters are a special case
            elif samp_fmt._fields[j] == 'FT':
                sampdat[j] = _parse_filter(vals)
                continue
            elif not vals or vals == '.':
                sampdat[j] = None
                continue

            entry_type = samp_fmt_types[j]
            # TODO: entry_num is None for unbounded lists
            entry_num = samp_fmt_nums[j]

            # we don't need to split single entries
            if entry_num == 1:
                if entry_type == INTEGER:
                    try:
                        sampdat[j] = int(vals)
                    except ValueError:
                        sampdat[j] = float(vals)
                elif entry_type == FLOAT or entry_type == NUMERIC:
                    sampdat[j] = float(vals)
                else:
                    sampdat[j] = vals
                continue

            vals = vals.split(',')
            if entry_type == INTEGER:
                try:
                    sampdat[j] = _map(int, vals)
                except ValueError:
                    sampdat[j] = map(float, vals)
            elif entry_type == FLOAT or entry_type == NUMERIC:
                sampdat[j] = _map(float, vals)
            else:
                sampdat[j] = vals

        # create a call object
        call = _Call(site, name, samp_fmt(*sampdat))
        samp_data.append(call)

    return samp_data
