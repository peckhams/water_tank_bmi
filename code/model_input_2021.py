#      
# Copyright (c) 2009-2013, Scott D. Peckham
#
# Moved here from model_tests.py on 3/14/13.
# Updated for Python 3.0 on 2021-02-12
#
#--------------------------------------------------------
# class input_file
#
#    __init__()
#    open()
#    check_format()
#    count_lines()
#    read_record()
#    read_value()
#    read_all()
#    close()
#
#--------------------------------------------------------
import numpy as np

#--------------------------------------------------------
#  Example utility class to read data from a text file.
#--------------------------------------------------------
class input_file:
    
    #------------------------------------------------------------
    def __init__(self, filename='data.txt'):
        
        self.filename = filename
        self.n_lines  = 0
        self.n_vals   = 0
        
        self.count_lines()   # (set self.n_lines)
        
        #------------------------------------------------
        # Check format of input data to set self.format
        # to 'numeric', 'key_value' or 'words'.
        #------------------------------------------------
        self.check_format()  # (set self.format string)

    #   __init__()    
    #------------------------------------------------------------
    def open(self):
        
        try:
            f = open(self.filename, 'r')
        except IOError as err:
            errno, strerror = err.args
            print('Could not find input file named:')
            print(self.filename)
            # print "I/O error(%s): %s" % (errno, strerror)
            print()
            return

        
        print('Reading file =', self.filename)   #######
        
        #------------------------------------------------------------
        #  The approach here is called "delegation" vs. inheritance.
        #  It would also be possible (in Python 2.2 and higher) for
        #  our input_file class to inherit from Python's built-in
        #  "file" type.  But that wouldn't necessarily be better.
        #  Type "Python, subclassing of built-in types" or
        #  "Python, how to subclass file" into a search engine for
        #  more information.

        #  If we subclass Python's "file" type, then how can we
        #  check for existence of the file, etc. ?
        #  Problem: file-object is a built-in type, not a class,
        #------------------------------------------------------------
        self.file = f   #  save the file object as an attribute

    #   open()
    #------------------------------------------------------------
    def check_format(self):
        
        self.open()
        record = self.read_record(format='words')
        n_words = len(record)
        
        #--------------------------------------
        # Is first value a number or keyword ?
        #--------------------------------------
        try:
            v = np.float64(record[0])
            format = 'numeric'
        except ValueError:
            format = 'unknown'
            if (n_words > 1):
                if (record[1] == "="):
                    format = 'key_value'

        self.format = format
        print('   data format =', format)
        print()
        self.close()

    #   check_format()
    #------------------------------------------------------------
    def count_lines(self):
        
        self.open()
        ## print 'Counting lines in file...'
        
        n_lines = 0
        n_vals  = 0
        for line in self.file:
            #-------------------------------------
            # Note: len(line) == 1 for null lines
            #-------------------------------------
            if (len(line.strip()) > 0):
                n_lines = (n_lines + 1)
                words   = line.split()
                n_words = len(words)
                n_vals  = max(n_vals, n_words)
        self.n_lines = n_lines
        self.n_vals  = n_vals
        #--------------------------------------
        # Initialize an array for reading data
        #--------------------------------------
        self.data = np.zeros([n_lines, n_vals], dtype='d')
        self.data = self.data - 9999.0
        print('   number of lines =', n_lines)
        print()
        self.close()

    #   count_lines()  
    #------------------------------------------------------------
    def read_record(self, format='not_set', dtype='double'):
        
        # Should this be named "next_record" ?
        if (format == 'not_set'): format=self.format
        
        line = self.file.readline()
        while (len(line.strip()) == 0):
            line = self.file.readline()
        words = line.split()
        n_words = len(words)
        
        if (format == 'numeric'):
            return list(map(np.float64, words))
        if (format == 'key_value'):
            key   = words[0]
            value = words[2]
            if (dtype == 'double'):  value = np.float64(value)
            if (dtype == 'integer'): value = np.int16(value)
            return [key, value]
        if (format == 'words'):
            return words

    #   read_record()
    #------------------------------------------------------------
    def read_value(self, dtype='double'):
        
        # Should this be named "next_value" ?
        record = self.read_record( dtype=dtype )
        if (self.format == 'numeric'):
            return record[0]
        if (self.format == 'key_value'):
            return record[1]

    #   read_value()
    #------------------------------------------------------------
    def read_all(self):
        
        #--------------------------------------------
        # Currently assumes all values are doubles,
        # but read_record() and read_value() do not.
        #--------------------------------------------
        if (self.n_lines == 0): self.count_lines()
        self.open()
        row = 0
        for line in self.file:
            record = self.read_record()
            if (self.format == 'numeric'):
                for col in range(0,len(record)):
                    self.data[row,col] = np.float64(record[col])
            if (self.format == 'key_value'):
                self.data[row,0] = np.float64(record[1])
            row = (row + 1)
        self.close()

    #   read_all()
    #------------------------------------------------------------
    def close(self):

        self.file.close()

    #   close()
    #------------------------------------------------------------        
        

