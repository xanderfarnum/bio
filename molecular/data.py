from metatools import append_docstring
from tools import get_from_entrez


@append_docstring(func=get_from_entrez(db='pubmed'))
class Data():
    """
    Base class for creating objects containing data from local file or remote server.
    """
    path_config = '/Users/alexanderfarnum/Documents/Code/Python/bio/userconfig.json'
    def __init__(self,config=path_config):
        from json import load
        self.config = load(open(config, 'r'))


    def get_data(self=None,source='remote',**kwargs):
        if self is None:
            construct = True
            if 'config' in kwargs:
                self = Data(config=kwargs.pop('config'))
            else:
                self = Data(config=Data.path_config)
        else:
            construct = False
        self.source = source
        self.parameters = Parameters(kwargs)
        self._get_from_source()
        if construct:
            return self

    def _get_from_source(self,**kwargs):
        """
        Pulls all relevant information based on specied user input (as arguments to self.get_data) and provided config
        #TODO: add kwarg specific functionalities
        """
        match self.source:
            case 'remote':
                from Bio import Entrez
                Entrez.email = self.config['entrez']['email']
            case 'local':
                pass
        if self.parameters.gene:
            from molecular.gene import Gene
            data = Gene(obj=self,source=self.source)


class Parameters():
    def __init__(self,kwargs):
        from molecular.utils import dict2attrs
        dict2attrs(obj=self,attr_dict=kwargs)