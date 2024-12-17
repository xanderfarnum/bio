import requests
import xml.etree.ElementTree as ET


source = requests.get('https://www.ncbi.nlm.nih.gov/gene/5813').json()

# payload = {'key1': 'value1', 'key2': ['value2', 'value3']}
resp = requests.get('https://www.ncbi.nlm.nih.gov/gene/5813')

root = ET.fromstring(resp.content)

temp = root.findall('td')

# Access elements and attributes
for element in root.findall('td'):  # Replace 'item' with the actual tag name
    print(element.tag, element.attrib)
    print(element.find('title').text)  # Access child elements

class Ncbi():
    def __init__(self):
        self