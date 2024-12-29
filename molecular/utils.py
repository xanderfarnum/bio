def print_seq(seq):
    print(repr(str(seq)))


def get_subdata(obj):
    datafile = "/Users/alexanderfarnum/Documents/Code/Python/bio/molecular/data.json"
    from json import load
    data = load(open(datafile, 'r'))
    if obj.parameters.gene:
        data_dict = data[obj.parameters.genome][obj.parameters.gene]
    output = parse_data_dict(data_dict)


def parse_data_dict(data_dict,out=None):
    for k,v in data_dict.items():
        if isinstance(k,dict):
            _,out = parse_data_dict(v)
        else:
            data = v