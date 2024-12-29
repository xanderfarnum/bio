def print_seq(seq):
    print(repr(str(seq)))


def get_subdata(obj):
    datafile = "/Users/alexanderfarnum/Documents/Code/Python/bio/molecular/data.json"
    from json import load
    data = load(open(datafile, 'r'))
    if obj.parameters.gene:
        data_dict = data[obj.parameters.genome][obj.parameters.gene]
    return parse_dict(data_dict)


def parse_dict(_dict,out=[]):
    for k,v in _dict.items():       # only iterate through top level- change to recursive function to pull nested levels as necessary
        try:
            cls = get_class(mod_name='molecular.objects',
                            cls_name=k)
            out.append(cls(**v))
        except:
            raise(f'Error instantiating class: {k}')
    return out

def get_class(mod_name, cls_name):
    from importlib import import_module
    mod = import_module(name=mod_name)
    return getattr(mod, cls_name)

def list2attrs(obj,attr_list):
    for v in attr_list:
        setattr(obj,v.__class__.__name__.lower(),v)

def dict2attrs(obj,attr_dict): 
    for k,v in attr_dict.items():
        setattr(obj,k.lower(),v)   


class Query():
    def __init__(self,str,delim=None):
        self.str = [str]
        if delim is None:
            self.delim = [' ']
        else:
            self.delim = [' ' + delim + ' ']

    def append(self,str,delim=None):
        self.str.append(str)
        if delim is None:
            self.delim.append(self.delim[0])
        else:
            self.delim.append(' ' + delim + ' ')

    def construct(self):
        str = [f"{x}{y}" for x,y in zip(self.str[:-1],self.delim[:-1])]
        return "".join(str) + self.str[-1]