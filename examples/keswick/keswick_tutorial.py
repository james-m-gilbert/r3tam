# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 10:47:28 2024

@author: jgilbert
"""

from r3tam.longtemp import LongTemp as lt

input_fp = r'kwk_input.v20240729.yaml'

kwk = lt.initialize_longmod(input_fp, show_init=True)
