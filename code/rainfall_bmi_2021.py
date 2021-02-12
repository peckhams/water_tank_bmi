#      
# Copyright (c) 2009-2013, Scott D. Peckham
#
# Moved here from model_tests.py on 3/14/13.
# Modified for Python 3.0 on 2021-02-12
#
#----------------------------------------------------------------
#
# class rainfall()
#    __init__()
#    initialize()
#    update()
#    finalize()
#    -------------
#    read_from_file()
#    simulate()
#    print_data()
#
# make_rain_file()   # (function)
#
#----------------------------------------------------------------

import numpy

#----------------------------------------------------------------
#  Simulate rainfall or read data from a file.
#----------------------------------------------------------------
class rainfall:

    #------------------------------------------------------------
    def __init__(self, n_values=5):
        
        self.n_values   = n_values
        self.rates     = numpy.zeros([n_values]) + 50.0  # [mm/sec]
        self.durations = numpy.zeros([n_values]) + 1.0   # [hours]

    #   __init__()
    #------------------------------------------------------------
    def initialize(self):
        pass
    #   initialize()
    #------------------------------------------------------------
    def update(self):
        pass
    #   update()
    #------------------------------------------------------------    
    def finalize(self):
        pass
    #   finalize()
    #------------------------------------------------------------
    def read_from_file(self, file='rain_data.txt'):
        
        #------------------------------------------
        # Parts of this probably belong in methods
        # for a "text_file" class, e.g. n_lines
        #------------------------------------------
        try: f = open(file, 'r')
        except IOError as err:
            errno, strerror = err.args
            print('I/O error(%s): %s' % (errno, strerror))
            return
        
        n_lines = 0
        for line in f:
            #-------------------------------------
            # Note: len(line) == 1 for null lines
            #-------------------------------------
            if (len(line.strip()) > 0): n_lines = (n_lines + 1)
        self.n_values  = n_lines
        self.rates     = numpy.zeros([n_lines], dtype='d')
        self.durations = numpy.zeros([n_lines], dtype='d')
        
        f.seek(0)
        k = 0
        for line in f:
            words = line.split()
            if (len(words) > 1):
                self.rates[k]     = numpy.float64(words[0])
                self.durations[k] = numpy.float64(words[1])
                k = (k + 1)
        f.close()

    #   read_from_file()
    #------------------------------------------------------------    
    def simulate(self):
        
        nv = self.n_values
        self.rates = 10 * numpy.random.standard_exponential([nv])
        # self.rates = numpy.random.uniform([nv]) (only 1 value)
        self.durations = numpy.zeros([nv], dtype='d') + 1.0

    #   simulate()
    #------------------------------------------------------------
    def print_data(self):
        
        print('Here are the data members:')
        print('n_values  =', self.n_values)
        print('rates     =', self.rates)
        print('durations =', self.durations)

    #   print_data()
    #------------------------------------------------------------

#------------------------------------------
# Make a file of rainfall data for testing
#------------------------------------------
def make_rain_file(filename, n_vals):
    
    rate = 80.0
    file = open(filename, 'w')
    for k in range(0, n_vals):
        rate     = rate + 10.0
        duration = 1.0
        val_str = str(rate) + ' ' + str(duration) + '\n'
        file.write(val_str)
    file.close()

#   make_rain_file()
#-----------------------------------------------------------------


