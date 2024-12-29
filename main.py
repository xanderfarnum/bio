
#TODO: create separate config files for loading data subsets based on intended functionality
config = '/Users/alexanderfarnum/Documents/Code/Python/bio/userconfig.json'

from molecular.data import Data

## Instantiate Data object via base constructor or get_data method

# data = Data(config)
data = Data.get_data(
    source='remote',
    species='Homo Sapien',
    gene='Pura',
    config=config
)

data.get_data(
    source='local',
    genome='GRCh37',
    species='Homo Sapien',
    gene='Pura'
)