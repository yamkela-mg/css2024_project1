#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 11:55:58 2024

@author: yamkelamgwatyu
"""

import numpy as np
from scipy.stats import chi2_contingency

ct = np.array([[28, 12], [16, 24]])

chi2_contingency(ct)


#titanic example

ct2 = np.array([[122, 203], [167, 118], [528, 178], [673, 212]])
chi2_contingency(ct2)

