#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 14:34:59 2024

@author: yamkelamgwatyu
"""

#Calculating AT content

mydna = "ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"

a_count = mydna.count("A")
print (a_count)
t_count = mydna.count("T")
tot_len = len(mydna)

AT_content = (a_count + t_count)/tot_len *100

print(AT_content)

#reverse complement 

