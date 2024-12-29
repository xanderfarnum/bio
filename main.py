from molecular.data import Data

# rdata = Data.get_data(
#         source='remote',
#         species='Homo Sapien',
#         gene='Pura'
# )

lcl_data = Data.get_data(
        source='local',
        genome='GRCh37',
        species='Homo Sapien',
        gene='Pura'
)
