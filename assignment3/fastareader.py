def fastaread(fpath):
    '''
    :param fpath: File pach
    :return: list of headers / list of sequences of FASTA-file
    '''
    f = open(fpath, 'r')
    headers = []
    seqs = []
    temp = ''
    for line in f:
        if line:
            if line[0] == '>':
                headers.append(line)
                if temp:
                    seqs.append(temp)
                    temp = ''
            else:
                temp += line.strip()
    seqs.append(temp.upper())
    return headers, seqs
