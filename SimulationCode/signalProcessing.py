#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 14:35:58 2023

@author: bikashadhikari
"""
from scipy import signal
import statistics
class signalProcessing:
    
    def __init__(self, ref, b, a, TRANSOFF):
        self.ref = ref # reference signal
        self.b =b  # transfer function numberator
        self.a = a # transfer function denominator
        self.TRANSOFF = TRANSOFF
    
    @property
    def referenceFilter(self):
        filteredReference = signal.lfilter(self.b, self.a, self.ref) 
        return filteredReference[self.TRANSOFF:]
    
    def signalFilter(self, sig):
        sig = sig.squeeze()
        sig_len = sig.size 
        filteredSignal = signal.lfilter(self.b, self.a, sig)
        filteredSignal = filteredSignal[self.TRANSOFF:]
        errorWRTreference = self.referenceFilter[0:len(filteredSignal)]-filteredSignal
        varianceError = round(statistics.variance(errorWRTreference.squeeze()),8)
        return [filteredSignal, errorWRTreference, varianceError]
        # return filteredSignal
    
        
        
 