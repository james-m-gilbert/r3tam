# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 10:50:17 2025

@author: jgilbert
"""

grammar = r"""
    ?start: section*
    
    section: "[" SECTION_NAME "]" NEWLINE property*
    
    use: USE
    property: KEY "=" VALUE NEWLINE
    
    USE: /use/
    SECTION_NAME: /[a-zA-Z0-9_]+/
    KEY: /[a-zA-Z0-9_]+/
    VALUE: /.+/
    NEWLINE: /(\r?\n)+/
    %import common.WS
    %ignore WS
"""

from lark import Lark
from lark import Transformer



class ConfigTransformer(Transformer):
    def __init__(self):
        self.config = {}

    def section(self, items):
        section_name = str(items[0])
        self.config[section_name] = {}
        return section_name

    def property(self, items):
      key, value = map(str, items)
      self.config[self.current_section][key] = value

    def start(self, items):
      return self.config

    def SECTION_NAME(self, token):
        self.current_section = str(token)
        return token

    def KEY(self, token):
        return str(token)

    def VALUE(self, token):
        return str(token)

#%%    
    
# with open("config.lark", "r") as f:
#     grammar = f.read()


parser = Lark(grammar)


dts_fp = r'D:\02_Projects\CalSim\util\CalSim_Utilities\DSS_PostProc_DPPAT\DPPATWorkshop_20170721\DPPAT_20170721\ZacksInput\calsim2.dt'
#dts_fp = r''

with open(dts_fp, 'r') as f: #myconfig.conf", "r") as f:
    config_content = f.read()

tree = parser.parse(config_content)
#%%
config_transformer = ConfigTransformer()
config_data = config_transformer.transform(tree)