#      
# Copyright (c) 2009-2021, Scott D. Peckham
#
# Started from "model_tests.py" on 3/14/13 to create
# a BMI-compliant version for a hands-on example.
#
# Show how to inherit BMI functions from BMI_base.py.

# Previous version info.
# July 5-7, 2008
# Modified: July 23, 2008  (S. Peckham)
#    v = sqrt(2 * g * d))

# Modified: July 18, 2016  (S. Peckham)
#    v = sqrt(g * d / 2) * (A_top / sqrt(A_top^2 - A_out^2))
#    v ~ sqrt(g * d / 2), if A_top >> A_out
#    This equation is computed for draining case from dh/dt
#    (in the companion paper), after simplification.
#
# Modified: 2021-02-12 (S. Peckham)
#    Minor changes to run in Python 3.*
#
#-----------------------------------------------------------------
# 7/18/16.  Key equations
#
# Q_in  = R * A_top
# Q_out = v * A_out
# d(t)  = d0 * [1 - (t/td)]^2  # solution for draining; no rain)
#   Note that d(td) = 0.
#
# d'(t) = (-2 d0 / td) * [1 - (t/td)]
# V'(t) = A_top * d'(t) = Q_out

# t_drain = sqrt(2 * d0 / g) * sqrt( (A_top/A_out)^2 - 1 )
# t_fill  = (d_final / R)   # (time to fill from d=0)
#
# If t_fill = p * t_drain (and d_f = d0), then
#   (r_top/r_out)^2 = A_top/A_out = sqrt( df*g/(2 R^2 p^2) + 1 )
# If p = 1/2, then filling will overpower draining.
#
# Try these settings in CFG file:
#   dt         = 4000.0   [s]
#   n_steps    = 60
#   init_depth = 1.0     [m]
#   top_radius = 20.0    [m]
#   out_radius = 0.05    [m]
#
# Can also try these settings from the paper:
#   init_depth = 86.40 [cm] = 0.864 [m]
#   top_radius = 0.5 * (29.21)  [cm]
#   out_radius = 0.5 * (0.533, 0.668, 0.945, 1.087)  [cm]
#   td         = 1223, 767, 403, 288   [s]
#
#-----------------------------------------------------------------
# 9/26/13.
# The steady-state solution for this problem is:
#     Q_in = Q_out
#     R * A_top = v * A_out
#     d_f = (2/g) * ( R * A_top / A_out)^2
#
# For draining:  d'(t) = -(2 d0 / td) * (1 - t/td)
# For filling:   d'(t) = R
#
# This steady-state solution approached whether the initial
# water depth, d0, is greater than or less than d_f.
#
#-----------------------------------------------------------------
# The "water_tank" class defines a model of a cylindrical water
# tank that is open at the top (so rainfall can get in), with a
# small, circular outlet.
#
# The syntax to create an instance of the model is:
#     >>> tank_model = water_tank()
#  
# Once instantiated, we can call any of its methods or
# "member functions", most of which are BMI functions. e.g.
#     >>> tank_model.initialize(),
#     >>> tank_model.update()
#
#----------------------------------------------------------------
#
# test_tank_model()
#
# class water_tank()
#
#     get_attribute()
#     get_input_var_names()
#     get_output_var_names()
#     ------------------------
#     get_var_name()
#     get_var_units()
#     get_var_type()
#     get_var_rank()            # (CHECK)
#     -------------------
#     get_start_time()
#     get_end_time()
#     get_current_time()
#     get_time_step()
#     get_time_units()         # (can also use get_attribute() now)
#     -------------------
#     initialize()
#     update()
#     finalize()
#     run_model()    # (not BMI)
#     -------------
#     get_value()
#     set_value()
#
#     -------------------
#     Non-BMI functions
#     -------------------
#     read_cfg_file()
#     update_rain()
#     print_tank_data()
#     read_rain_data()
#
#----------------------------------------------------------------

# import BMI_base

import numpy as np  # (for: float64, int16, pi, sqrt, zeros)
import model_input_2021 as model_input
import time

#---------------------------------------------
# Create an instance of the water tank model
# and then call its "run_model()" method.
#---------------------------------------------
def test_tank_model():
    tank_model = water_tank()
    tank_model.run_model()
         
#-------------------------------
# Define the water_tank class.
#-------------------------------
# class water_tank( BMI_base.BMI_component ):   # (option to inherit)
class water_tank:

    #-------------------------------------------
    # Required, static attributes of the model
    #-------------------------------------------
    _att_map = {
        'model_name':         'Water_Tank_Model',
        'version':            '1.0',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'none',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        'time_units':         'seconds' }
    
    #----------------------------------------------
    # Input variable names (CSDMS Standard Names)
    #----------------------------------------------
    _input_var_names = [
        'atmosphere_water__liquid_equivalent_precipitation_rate',
        'atmosphere_water__precipitation_duration' ]

    #-----------------------------------------------
    # Output variable names (CSDMS Standard Names)
    #-----------------------------------------------
    _output_var_names = [
        'model__time_step',
        'tank_x-section__area',
        'tank_x-section__radius',
        'tank_outlet__area',
        'tank_outlet__radius',
        'tank_outlet_water__flow_speed',
        'tank_water__depth'
        'tank_water__initial_depth',
        'tank_water__volume' ]

    #------------------------------------------------------
    # Create a Python dictionary that maps CSDMS Standard
    # Names to the model's internal variable names.
    #------------------------------------------------------
    _var_name_map = {
        'atmosphere_water__liquid_equivalent_precipitation_rate': 'rain_rate',
        'atmosphere_water__precipitation_duration': 'rain_duration',
        #------------------------------------------------------------------------------
        'model__time_step':              'dt',
        'tank_x-section__area':          'top_area',
        'tank_x-section__radius':        'radius',
        'tank_outlet__area':             'out_area',
        'tank_outlet__radius':           'out_radius',
        'tank_outlet_water__flow_speed': 'out_speed',
        'tank_water__depth':             'depth',
        'tank_water__initial_depth':     'init_depth',
        'tank_water__volume':            'volume'}
        
    #------------------------------------------------------
    # Create a Python dictionary that maps CSDMS Standard
    # Names to the units of each model variable.
    #------------------------------------------------------
    _var_units_map = {
        'atmosphere_water__liquid_equivalent_precipitation_rate': 'm s-1',
        'atmosphere_water__precipitation_duration': 's',
        #--------------------------------------------------------------------
        'model__time_step':              's',
        'tank_x-section__area':          'm2',
        'tank_x-section__radius':        'm',
        'tank_outlet__area':             'm2',
        'tank_outlet__radius':           'm',
        'tank_outlet_water__flow_speed': 'm s-1',
        'tank_water__depth':             'm',
        'tank_water__initial_depth':     'm',
        'tank_water__volume':            'm3'}

    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    # BMI: Model Information Functions
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print('###################################################')
            print(' ERROR: Could not find attribute: ' + att_name)
            print('###################################################')
            print()

    #   get_attribute()    
    #-------------------------------------------------------------------
    def get_input_var_names(self):

        #--------------------------------------------------------
        # Note: These are currently variables needed from other
        #       components vs. those read from files or GUI.
        #--------------------------------------------------------   
        return self._input_var_names
    
    #   get_input_var_names()
    #-------------------------------------------------------------------
    def get_output_var_names(self):
 
        return self._output_var_names
    
    #   get_output_var_names()
    #-------------------------------------------------------------------
    # BMI: Variable Information Functions
    #-------------------------------------------------------------------
    def get_var_name(self, long_var_name):
            
        return self._var_name_map[ long_var_name ]

    #   get_var_name()
    #-------------------------------------------------------------------
    def get_var_units(self, long_var_name):

        return self._var_units_map[ long_var_name ]
   
    #   get_var_units()
    #-------------------------------------------------------------------
    def get_var_type(self, long_var_name):

        #---------------------------------------
        # So far, all vars have type "double".
        #---------------------------------------
        return 'float64'

        #--------------------------------------------------
        # A more general approach, with less maintenance.
        #--------------------------------------------------
        ## return str( self.get_value( long_var_name ).dtype )

        #-----------------------------------
        # Another approach (not ready yet)
        #-----------------------------------
##        dtype = getattr(self, var_name)
##        return str(dtype)
    
        #-------------------
        # Another approach
        #-------------------
##        var_name = self.get_var_name( long_var_name )
##
##        try:
##            exec( "dtype = self." + var_name + ".dtype" )
##        except:
##            dtype = 'unknown'
##        return str(dtype)       # (need str() here)
    
    #   get_var_type()
    #-------------------------------------------------------------------
    def get_var_rank(self, long_var_name):

        return np.int16(0)
    
    #   get_var_rank()
    #------------------------------------------------------------     
    def get_start_time( self ):
        
        return 0.0

    #   get_start_time()
    #------------------------------------------------------------ 
    def get_end_time( self ):

        return (self.n_steps * self.dt)
    
    #   get_end_time()
    #------------------------------------------------------------ 
    def get_current_time( self ):
        
        return self.time

    #   get_current_time()
    #------------------------------------------------------------ 
    def get_time_step( self ):
        
        return self.dt

    #   get_time_step()
    #------------------------------------------------------------ 
    def get_time_units( self ):
        
        return 'seconds'

        #--------------
        # Another way
        #--------------
        # units = self.get_attribute( 'time_units' )
        # return units
    
    #   get_time_units()
    #------------------------------------------------------------
    # BMI: Model Control Functions
    #------------------------------------------------------------ 
    def initialize( self, cfg_file=None ):

        self.SERIALIZABLE = True
        self.timer_start = time.time()

        #-------------------------------------------
        # Used in read_cfg_file(), so needed here.
        #-------------------------------------------
        self.g = 9.81    # [m/s^2]
        
        #--------------------------------------
        # Read tank settings from "tank_file"
        #--------------------------------------
        if (cfg_file == None):
            cfg_file = 'tank_info.cfg'
        self.cfg_file = cfg_file
        self.read_cfg_file()

        #-----------------------       
        # Initialize variables
        #-----------------------
        self.depth      = self.init_depth.copy()
        self.out_speed  = np.sqrt(self.g * self.depth / 2.0)
        self.volume     = self.depth * self.top_area  #[m^3]
        self.out_area   = np.pi * self.out_radius**2.0
        self.print_tank_data()

        #----------------------------
        # Initialize time variables
        #----------------------------
        self.time = np.float64(0)
        self.time_index = 0

        #-------------------------------------------------       
        # Use "input_file" class to create rain_file
        # object, then open "rain_file" to read data.
        # This will be used by the update_rain() method.
        #-------------------------------------------------
        if ('steady:' in self.rain_data_filename):
            #-------------------------------------------
            # New option to specify a steady rain rate
            # in mmph.
            #-------------------------------------------
            self.steady_rain = True
            words = self.rain_data_filename.split(':')
            R = np.float64( words[1] )
            self.rain_rate = R
            #-----------------------------------------
            # Print the predicted steady-state depth
            #-----------------------------------------
            R_mps   = R / (3600.0 * 1000.0)
            fac     = (1 / (2.0 * self.g))
            A_ratio = (self.top_area / self.out_area)
            d_steady = fac * (R_mps * A_ratio)**2
            # print 'Steady-state rainrate =', R, ' [mmph]'
            # print 'Predicted steady-state depth =', d_steady
        else:
            self.steady_rain = False
            if (self.SERIALIZABLE):
                # Can use dill to serialize.
                self.read_rain_data( self.rain_data_filename )
            else:
                self.rain_file = model_input.input_file(self.rain_data_filename)       
                self.rain_file.open()

    #   initialize()
    #------------------------------------------------------------       
    def update( self, dt=-1, REPORT=True ):

        #------------------------------------------------       
        # Read the next rainfall file data entry.
        # "rain_duration" is read, but ignored for now.
        #------------------------------------------------
        if not(self.steady_rain):
            self.update_rain()
        rain_rate_mps = self.rain_rate / (3600.0 * 1000.0)

        #--------------------------------------------------
        # Compute volume inflow rate from rain, Q_in
        # and volume outflow rate, Q_out, in [m^3 / sec].
        #--------------------------------------------------
        Q_in = rain_rate_mps * self.top_area
        if (self.depth > 0):
            Q_out = self.out_speed * self.out_area
        else:
            Q_out = 0.0
        dVol = (Q_in - Q_out) * self.dt

        #----------------------------
        # Store the state variables
        #----------------------------
        self.Q_out  = Q_out
        self.volume = max(self.volume + dVol, 0.0)
        self.depth  = (self.volume / self.top_area)
        self.out_speed = np.sqrt(2.0 * self.g * self.depth)
        ### self.out_speed = np.sqrt(self.g * self.depth / 2.0)
        
        #-------------------------
        # Optional status report
        #-------------------------
        if (REPORT):
            print('--------------------------------------')
            print('rain_rate =', self.rain_rate, ' [mmph]')
            print('depth     =', self.depth, '[meters]')
        
        #--------------------------------------
        # Write new depth to an output file ?
        #--------------------------------------

        #------------------------
        # Update the model time
        #------------------------
        if (dt == -1):
            dt = self.dt
        self.time += dt
        self.time_index += 1
        
    #   update()
    #------------------------------------------------------------       
##    def update_until( self, time ):
##        
##        #----------------------------------------------------
##        # Call update() method as many times as necessary
##        # in order to get to the requested time.  Note that
##        # we do not override the value of n_steps from
##        # the tank model's cfg_file.
##        #----------------------------------------------------
##        n_steps = np.int16(time / self.dt)
##        for k in xrange(1,n_steps+1):
##            self.update()
##
##    #   update_until()
    #------------------------------------------------------------    
    def finalize( self ):

        timer_stop = time.time()
        run_time = (timer_stop - self.timer_start)

        #-----------------------------------------        
        # Report simulation time with good units
        #-----------------------------------------
        sim_time  = self.time  # [secs]
        sim_units = ' [secs]'
        #----------------------------------
        if (sim_time > 3600 * 24):
            sim_time = sim_time / (3600.0 * 24)
            sim_units = ' [days]'
        elif (sim_time > 3600):
            sim_time = sim_time / 3600.0
            sim_units = ' [hrs]'       
        elif (sim_time > 60):
            sim_time = sim_time / 60.0
            sim_units = ' [mins]'

        print()
        print('Finished with water tank simulation.')
        print('Model run time =', run_time, ' [secs]')
        print('Simulated time =', sim_time, sim_units)
        print('Final depth    =', self.depth, ' [m]')
        print()

        #-------------------
        # Close input file
        #-------------------
        if not(self.SERIALIZABLE):
           if not(self.steady_rain):
              self.rain_file.close()

        #-----------------------        
        # Close output files ?
        #-----------------------
        
    #   finalize()
    #------------------------------------------------------------ 
    def run_model( self, cfg_file=None):

        #-------------------------------------------------------
        # Note: This is not a required BMI function, but gives
        #       an easy way to run the stand-alone model.
        #-------------------------------------------------------
        self.initialize( cfg_file=cfg_file )
        for k in range(1, self.n_steps+1):
            # print('k =', k)
            self.update()
        self.finalize()

    #   run_model()
    #------------------------------------------------------------
    # BMI: Variable Getters and Setters
    #------------------------------------------------------------       
    def get_value(self, long_var_name):

        var_name = self.get_var_name( long_var_name )

        try:
            return getattr(self, var_name)

            #----------------------------
            # This breaks the reference.
            #----------------------------
            ## return np.float64(result)
            
        except:
            print('ERROR in get_value() function')
            print('    for var_name =', var_name)
            print('    Returning 0.')
            return np.array(0, dtype='float64')

    #   get_value()
    #-------------------------------------------------------------------
    def set_value(self, long_var_name, value):

        #---------------------------------------------------------------
        # Notes: The "var_name" string cannot contain a ".". (5/17/12)
        #---------------------------------------------------------------
        # (2/7/13) We are now using 0D numpy arrays as a way to
        # produce "mutable scalars" that allow a component with a
        # reference to the scalar to see changes to its value.
        # But we can't apply np.float64() to the value as we did
        # before or it destroys the reference.
        # See BMI_base.initialize_scalar() for more information.
        #--------------------------------------------------------------- 
        var_name = self.get_var_name( long_var_name )
        setattr( self, var_name, value )
        
    #   set_value()
    #------------------------------------------------------------
    # Non-BMI functions that are only used internally.
    #------------------------------------------------------------
    def read_cfg_file( self, cfg_file=None ):
        
        #-----------------------------------
        # What if cfg_file doesn't exist ?
        #-----------------------------------
        if (cfg_file == None):
            cfg_file = self.cfg_file
      
        tank_file = model_input.input_file( cfg_file )
        tank_file.open()
        
        #------------------------------------------------
        # Read values from cfg_file and store in "self"
        #------------------------------------------------
        self.dt         = tank_file.read_value()
        self.n_steps    = tank_file.read_value( dtype='integer' )
        self.init_depth = tank_file.read_value()
        self.top_radius = tank_file.read_value()
        self.top_area   = np.pi * (self.top_radius)**2.0  # (7/18/16)
        # self.top_area   = tank_file.read_value()
        self.out_radius = tank_file.read_value()
        self.rain_data_filename = tank_file.read_value( dtype='string' )
        tank_file.close()
        
    #   read_cfg_file
    #------------------------------------------------------------
    def update_rain( self ):

        ## if (self.time_index < self.rain_file.n_lines):
        if (self.time_index < self.n_rain_lines):
            if (self.SERIALIZABLE):
                self.rain_rate     = self.rates[ self.time_index ]
                self.rain_duration = self.durations[ self.time_index ]
            else:
                record = self.rain_file.read_record()           
                self.rain_rate     = record[0]   ## (in mmph)
                self.rain_duration = record[1]   ## (in seconds)
            #---------------------------------------------------------
            ## print 'rain_rate =', self.rain_rate
            ## print 'duration  =', self.rain_duration
            ## print ' '
        else:
            self.rain_rate = 0.0
            self.rain_duration = self.dt

    #   update_rain()
    #------------------------------------------------------------
    def print_tank_data( self ):

        print('   dt         =', self.dt, '[sec]')
        print('   n_steps    =', self.n_steps)
        print('   init_depth =', self.init_depth, '[m]')
        print('   top_radius =', self.top_radius, '[m]')
        print('   top_area   =', self.top_area, '[m2]')
        print('   out_radius =', self.out_radius, '[m]')
        print('   out_speed  =', self.out_speed, '[m/s]')
        print('   depth      =', self.depth, '[m]')
        print('   volume     =', self.volume, '[m3]')
        print('   out_area   =', self.out_area, '[m2]')
        print('   rain_file  =', self.rain_data_filename)
        print()

    #   print_tank_data()
    #------------------------------------------------------------
    def read_rain_data(self, filename='rain_data.txt'):
        
        #------------------------------------------
        # Parts of this probably belong in methods
        # for a "text_file" class, e.g. n_lines
        #------------------------------------------
        try: f = open(filename, 'r')
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
        self.n_rain_lines  = n_lines
        self.rates     = np.zeros([n_lines], dtype='d')
        self.durations = np.zeros([n_lines], dtype='d')
        
        f.seek(0)
        k = 0
        for line in f:
            words = line.split()
            if (len(words) > 1):
                self.rates[k]     = np.float64(words[0])
                self.durations[k] = np.float64(words[1])
                k = (k + 1)
        f.close()

    #   read_rain_data()
    #------------------------------------------------------------
        
      
