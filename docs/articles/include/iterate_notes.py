import xml.etree.cElementTree as ET
import os,re

from personal_package import readfiles,writefile

def folder_walk(folder):
    for root, dirs, files in os.walk(folder):
        return [[name,os.path.join(folder,name)] for name in files]

class n2c2xml():
    def __init__(self,path):
        self.path=path
        self.Tree=ET.parse(self.path)
        self.Root=self.Tree.getroot()
        self.TAGS = self.Root.find('TAGS')
    def GetText(self):
        return self.Root.find('TEXT').text
    def GetTAGS(self):
        return dict([[t.tag,t.attrib['met']] for t in self.TAGS])
def get_date(s):
    p = re.compile('\d+-\d+-\d+')
    res = p.search(s)
    if res:
        return res.group()
    else: return 'None'

if __name__ == '__main__':
    folder = PATH_TO_THE_DATA_FOLDER
    for line in folder_walk(folder):
        f,p = line
        data = n2c2xml(p)
        firstline = data.GetText().strip().split('\n')[0]
        date = get_date(firstline)