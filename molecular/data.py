from metatools import append_docstring
from tools import get_from_entrez


@append_docstring(func=get_from_entrez(db='pubmed'))
class Data():
    """
    Base class for creating objects containing data from local file or remote server.
    """
    path_config = '/Users/alexanderfarnum/Documents/Code/Python/bio/userconfig.json'
    def __init__(self,
                 source,
                 kwargs):

        from json import load
        self.config = load(open(self.path_config, 'r'))
        self.source = source
        self.parameters = Parameters(kwargs)

        match self.source:
            case 'remote':
                self._remote_source()
            case 'local':
                self._local_source()
    

    def get_data(self=None,
                 source='remote',
                 **kwargs):
        if self is None:
            self = Data(source,
                        kwargs)
        else:
            print(f'get_data method must be used prior to class object instantiation')
        return self

    def _remote_source(self):
        from Bio import Entrez
        Entrez.email = self.config['entrez']['email']

        search_term = []
        if self.species: 
            search_term.append(self.species+"[Orgn]")
        if self.gene: 
            search_term.apend(self.species+"[Gene]")
        search_term = ' AND '.join(search_term)

        stream = Entrez.esearch(db="nucleotide",
                                term=search_term,
                                idtype="acc")
        record = Entrez.read(stream)

    def _local_source(self):
        if self.parameters.gene:
            from molecular.gene import Gene
            data = Gene(obj=self)



        # print_seq(object.protein.seq)



class Parameters():
    def __init__(self,kwargs):
        self._set_params(kwargs)

    def _set_params(self,params):
        for k,v in params.items():
            setattr(self,k,v)