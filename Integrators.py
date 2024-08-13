import numpy as np

G_const = 6.67e-11 # Grav constant [N*m^2/kg^2]
R_e = 6378.137 # Radius of the Earth in [m]
M_e = 5.972e24 # Mass of the Earth in [kg]

# RHS for the Keplerian problem, using dimensionless units
# Inputs:
# - pos[ 3 ] -> current moment position [R_e] dimensionless units
# - vel[ 3 ] -> current moment velocity [V_e] dimensionless units
# Ouputs:
# - rhs_pos[ 3 ] -> current moment RHS for position [-]
# - rhs_vel[ 3 ] -> current moment RHS for velocity [-]
# NOTE: ALL quantities plugged and returned to this function are already dimensionless and the equation of motion take that into account for R_e, V_e and T_e defined bellow
def Kepler_RHS( pos , vel ):

    #pos_norm = np.sqrt( pos[ 0 ]**2.0 + pos[ 1 ]**2.0 + pos[ 2 ]**2.0 )
    rhs_pos = vel
    rhs_vel = - 4.0*np.pi*np.pi*pos/( np.sqrt( np.dot( pos , pos ) )**3.0 )

    return rhs_pos, rhs_vel

# Keplerian orbital motion integrator using RK(4)
# Inputs:
# - pos_0[ 3 ] -> position in initial moment [ x , y , z ], [km]
# - vel_0[ 3 ] -> velocity in initial moment [ vx, vy, vz ], [km/sec]
# - t_params[ 3 ] -> [ t_i , t_f , Npoints ] 
# -- initial time (t_i) in [sec]
# -- final time (t_f) in [sec]
# -- Number of points Npoints [-]
# Outputs:
# - t_arr[ Npoints ] -> time array between t_i and t_f with Npoints in [sec]
# - pos[ Npoints ][ 3 ] -> position as function of time and [ X , Y , Z ] in [km]
# - vel[ Npoints ][ 3 ] -> velocity as function of time and [ VX, VY, VZ ] in [km/sec]
def Kepler_RK4( pos_0 , vel_0 , t_params ):

    # Compute dimensionfull quantities -> the ones we use for norming
    # We use R_e as spatial scale R_e is our dimensionful parameter in [km]
    T_e = 2.0*np.pi*np.sqrt( ( R_e*1000.0 )**3/( M_e*G_const ) ) # Characteristic period in [sec]
    V_e = R_e/T_e # Characteristic velocity in [km/sec]

    # Take out the initial and final time (1st and 2nd element of t_params)
    # NOTE: Renorm it to dimensionless
    t_i = t_params[ 0 ]/T_e # NOT [sec], but in units of T_e
    t_f = t_params[ 1 ]/T_e # NOT [sec], but in untis of T_e
    # Take out the number of points which is 3d element in the t_params
    Npoints = t_params[ 2 ]
    # Compute the time step dt 
    dt = ( t_f - t_i )/( Npoints - 1 ) # NOTE: Already dimensionless due to t_i and t_f
    # NOTE: Intervals are 1 less tha number of points -> | | | has 3 points, 2 interfvals

    # Initialize the vectors to hold the solution
    # First dimension gives you time cross-sections, 2nd dimension (index) gives you [ x , y , z ] 
    pos = np.zeros( ( Npoints , 3 ) )
    vel = np.zeros( ( Npoints , 3 ) )
    # For example, 10th point in time, Y position coordinate is pos[ 9 ][ 1 ], because you start from index 0!

    # Create the time array
    t_arr = np.linspace( t_i , t_f , Npoints )
    # Assign initial values for the vectors
    pos[ 0 ] = np.array( pos_0 )/R_e # Pass initial position value and renorm it in units of [R_e]
    vel[ 0 ] = np.array( vel_0 )/V_e # Pass initial velocity value and renorm it in units of [V_e]

    # Initializing the k vectors which will be used to hold the Runge-Kutta derivatives 
    k_pos = np.zeros( ( 4 , len( pos_0 ) ) )
    k_vel = np.zeros( ( 4 , len( vel_0 ) ) )
    # NOTE: Examples:
    # k_pos[ 0 ][ 0 ] -> is just k_1 for X coordinate
    # k_pos[ 1 ][ 0 ] -> is just k_2 for X coordinate
    # k_pos[ 2 ][ 0 ] -> is just k_3 for X coordinate
    # ...
    # k_pos[ 2 ][ 1 ] -> is just k_3 for Y coordinate
    # ...
    # k_pos[ 3 ][ 2 ] -> is just k_4 for Z coordinate

    # Doing the actual integration loop
    for i in range( 0 , Npoints - 1 ):
        # Compute RHS at current point f( t , x )
        rhs_pos, rhs_vel = Kepler_RHS( pos[ i ] , vel[ i ] ) # Pass current position and velocity in time step [ i ]
        # Assign to k1 for all quantities
        k_pos[ 0 ] = rhs_pos
        k_vel[ 0 ] = rhs_vel

        # Compute RHS at first half-step f( t + dt/2 , x + k1*dt/2 )
        rhs_pos, rhs_vel = Kepler_RHS( pos[ i ] + k_pos[ 0 ]*dt/2.0 , vel[ i ] + k_vel[ 0 ]*dt/2.0 ) # Pass current position and velocity in time step [ i ], incremented by k1/2 
        # Assign to k2 for all quantities
        k_pos[ 1 ] = rhs_pos
        k_vel[ 1 ] = rhs_vel

        # Compute RHS at second half-step f( t + dt/2 , x + k2*dt/2 )
        rhs_pos, rhs_vel = Kepler_RHS( pos[ i ] + k_pos[ 1 ]*dt/2.0 , vel[ i ] + k_vel[ 1 ]*dt/2.0 ) # Pass current position and velocity in time step [ i ], incremented by k2/2
        # Assign to k3 for all quantities
        k_pos[ 2 ] = rhs_pos
        k_vel[ 2 ] = rhs_vel

        # Compute RHS at next step f( t + dt , x + k3*dt )
        rhs_pos, rhs_vel = Kepler_RHS( pos[ i ] + k_pos[ 2 ]*dt , vel[ i ] + k_vel[ 2 ]*dt ) # Pass current position and velocity in time step [ i ], incremented by k3
        # Assign to k4 for all quantities
        k_pos[ 3 ] = rhs_pos
        k_vel[ 3 ] = rhs_vel

        # At this stage we have k1, k2, k3, k4 for each X, Y, Z, VX, VY, VZ computed (24 values), let's sum them up for each X, Y, Z, VX, VY, VZ to get the next step
        pos[ i + 1 ] = pos[ i ] + ( k_pos[ 0 ] + 2.0*k_pos[ 1 ] + 2.0*k_pos[ 2 ] + k_pos[ 3 ] )*dt/6.0
        vel[ i + 1 ] = vel[ i ] + ( k_vel[ 0 ] + 2.0*k_vel[ 1 ] + 2.0*k_vel[ 2 ] + k_vel[ 3 ] )*dt/6.0

    # Return after multiplying back to [km] and [km/sec] quantities with dimension
    return t_arr*T_e, pos*R_e , vel*V_e


# Right-Hand-Side function which allows us to abstract the forces
# for x' = rhs_x & v' = rhs_v -> you pass it x & v at some time and it gives you RHS
# Inputs: 
# - k spring const [N/m]
# - m spring mass [kg]
# - x position [m]
# - v velocity [m/sec]
# Outputs:
# - rhs_x [m/sec]
# - rhs_v [m/sec^2]
def rhs_osc( k , m , x , v ):

    rhs_x = v
    rhs_v = - k/m * x

    return rhs_x , rhs_v

# Integrate Harmonic oscillator through Euler's method
# Inputs:
# - k -> spring constant [N/m]
# - m -> mass [kg]
# - x_0 -> initial position [m]
# - v_0 -> initial velocity [m/sec]
# - t_params[ 3 ] -> [ t_i , t_f , Npoints ] 
# -- initial time (t_i) in [sec]
# -- final time (t_f) in [sec]
# -- Number of points Npoints [-]
# Outputs:
# - t_arr[ Npoints ] -> time array from t_i to t_f in [sec]
# - x_arr[ Npoints ] -> position solution in [m]
# - v_arr[ Npoints ] -> velocity solution in [m] 
def Euler_Oscillator( k , m , x_0 , v_0 , t_params ):

    # Take out the initial and final time (1st and 2nd element of t_params)
    t_i = t_params[ 0 ]
    t_f = t_params[ 1 ]
    # Take out the number of points which is 3d element in the t_params
    Npoints = t_params[ 2 ]
    # Compute the time step dt 
    dt = ( t_f - t_i )/( Npoints - 1 ) 
    # NOTE: Intervals are 1 less tha nnumber of points -> | | | has 3 points, 2 interfvals

    # Create empty arrays which will hold the solution
    x_arr = np.zeros( Npoints )
    v_arr = np.zeros( Npoints )
    # Create the time array
    t_arr = np.linspace( t_i , t_f , Npoints )

    # Assign the initial position and velocity
    x_arr[ 0 ] = x_0
    v_arr[ 0 ] = v_0

    # Integration loop -> go over each time step and extrapolate the next
    for i in range( 0 , Npoints - 1 ):
        # NOTE: THIS IS BAD, NEVER HARDCODE YOUR RHS
        # x_(i+1) = x_(i) + v_(i)*dt
        #x_arr[ i + 1 ] = x_arr[ i ] + v_arr[ i ]*dt
        # v_(i+1) = v_(i) + f_(i)*dt, where f_(i) = -k/m*x_(i) for oscillator
        #v_arr[ i + 1 ] = v_arr[ i ] + ( - k / m )*x_arr[ i ]*dt
        rhs_x, rhs_v = rhs_osc( m , k , x_arr[ i ] , v_arr[ i ] )
        x_arr[ i + 1 ] = x_arr[ i ] + rhs_x*dt
        v_arr[ i + 1 ] = v_arr[ i ] + rhs_v*dt

    # My x_arr and v_arr are populated with the solution - return them
    return t_arr , x_arr , v_arr


# Integrate Harmonic oscillator through Verlet's method
# Inputs:
# - k -> spring constant [N/m]
# - m -> mass [kg]
# - x_0 -> initial position [m]
# - v_0 -> initial velocity [m/sec]
# - t_params[ 3 ] -> [ t_i , t_f , Npoints ] 
# -- initial time (t_i) in [sec]
# -- final time (t_f) in [sec]
# -- Number of points Npoints [-]
# Outputs:
# - t_arr[ Npoints ] -> time array from t_i to t_f in [sec]
# - x_arr[ Npoints ] -> position solution in [m]
# - v_arr[ Npoints ] -> velocity solution in [m] 
def Verlet_Oscillator( k , m , x_0 , v_0 , t_params ):

    # Take out the initial and final time (1st and 2nd element of t_params)
    t_i = t_params[ 0 ]
    t_f = t_params[ 1 ]
    # Take out the number of points which is 3d element in the t_params
    Npoints = t_params[ 2 ]
    # Compute the time step dt 
    dt = ( t_f - t_i )/( Npoints - 1 ) 
    # NOTE: Intervals are 1 less tha nnumber of points -> | | | has 3 points, 2 interfvals

    # Create empty arrays which will hold the solution
    x_arr = np.zeros( Npoints )
    v_arr = np.zeros( Npoints )
    # Create the time array
    t_arr = np.linspace( t_i , t_f , Npoints )

    # Assign the initial position and velocity
    x_arr[ 0 ] = x_0
    v_arr[ 0 ] = v_0

    # Integration loop -> go over each time step and extrapolate the next
    for i in range( 0 , Npoints - 1 ):
        # Compute RHS in the first point (i)
        rhs_x, rhs_v = rhs_osc( m , k , x_arr[ i ] , v_arr[ i ] )
        # Compute half-velocity step -> intermediate step
        v_half = v_arr[ i ] + rhs_v*dt/2.0 # NOTE: dt/2 -> half-step
        # Compute RHS in the intermediate point (i+1/2)
        rhs_x, rhs_v = rhs_osc( m , k , x_arr[ i ] , v_half )
        # Compute full position step based on v_half
        x_arr[ i + 1 ] = x_arr[ i ] + rhs_x*dt
        # Compute RHS in the last point (i+1)
        rhs_x, rhs_v = rhs_osc( m , k , x_arr[ i + 1 ] , v_half )
        # Compute the full velocity step based on v_half and x[ i + 1 ]
        v_arr[ i + 1 ] = v_half + rhs_v*dt/2.0

    # My x_arr and v_arr are populated with the solution - return them
    return t_arr , x_arr , v_arr

# Integrate Harmonic oscillator through Runge-Kutta 4th order method
# Inputs:
# - k -> spring constant [N/m]
# - m -> mass [kg]
# - x_0 -> initial position [m]
# - v_0 -> initial velocity [m/sec]
# - t_params[ 3 ] -> [ t_i , t_f , Npoints ] 
# -- initial time (t_i) in [sec]
# -- final time (t_f) in [sec]
# -- Number of points Npoints [-]
# Outputs:
# - t_arr[ Npoints ] -> time array from t_i to t_f in [sec]
# - x_arr[ Npoints ] -> position solution in [m]
# - v_arr[ Npoints ] -> velocity solution in [m] 
def RK4_Oscillator( k , m , x_0 , v_0 , t_params ):

    # Take out the initial and final time (1st and 2nd element of t_params)
    t_i = t_params[ 0 ]
    t_f = t_params[ 1 ]
    # Take out the number of points which is 3d element in the t_params
    Npoints = t_params[ 2 ]
    # Compute the time step dt 
    dt = ( t_f - t_i )/( Npoints - 1 ) 
    # NOTE: Intervals are 1 less tha nnumber of points -> | | | has 3 points, 2 interfvals

    # Create empty arrays which will hold the solution
    x_arr = np.zeros( Npoints )
    v_arr = np.zeros( Npoints )
    # Create the time array
    t_arr = np.linspace( t_i , t_f , Npoints )

    # Assign the initial position and velocity
    x_arr[ 0 ] = x_0
    v_arr[ 0 ] = v_0

    # Integration loop -> go over each time step and extrapolate the next with RK(4)
    for i in range( 0 , Npoints - 1 ):
        # Compute RHS at first point (i) in preparation for k_1 computation
        rhs_x, rhs_v = rhs_osc( m , k , x_arr[ i ] , v_arr[ i ] )
        # Assign k_1 derivatives
        kx_1 = rhs_x
        kv_1 = rhs_v
        # Compute RHS at second point ( i+dt/2 , x+k_1/2 ) in preparation for k_2 computation
        rhs_x, rhs_v = rhs_osc( m , k , x_arr[ i ] + kx_1*dt/2.0 , v_arr[ i ] + kv_1*dt/2.0 )
        # Assign k_2 derivatives
        kx_2 = rhs_x
        kv_2 = rhs_v
        # Compute RHS at third point ( i+dt/2 , x+k_2/2 ) in preparation for k_3 computation
        rhs_x, rhs_v = rhs_osc( m , k , x_arr[ i ] + kx_2*dt/2.0 , v_arr[ i ] + kv_2*dt/2.0 )
        # Assign k_3 derivatives
        kx_3 = rhs_x
        kv_3 = rhs_v
        # Compute RHS at forth and last point ( i+dt , x+k_3 ) in preparation for final computation
        rhs_x, rhs_v = rhs_osc( m , k , x_arr[ i ] + kx_3*dt , v_arr[ i ] + kv_3*dt )
        # Assign k_4 derivatives
        kx_4 = rhs_x
        kv_4 = rhs_v

        # CONGRATS, WE HAVE OUR K VALUES, kx_1, kx_2, kx_3, kx_4 & kv_1, kv_2, kv_3, kv_4
        x_arr[ i + 1 ] = x_arr[ i ] + ( kx_1 + 2.0*kx_2 + 2.0*kx_3 + kx_4 )*dt/6.0
        v_arr[ i + 1 ] = v_arr[ i ] + ( kv_1 + 2.0*kv_2 + 2.0*kv_3 + kv_4 )*dt/6.0

    # My x_arr and v_arr are populated with the solution - return them
    return t_arr , x_arr , v_arr


# Compute the Right-Hand-Side (RHS) of ballistic trajectory in 2D
# Inputs:
# - g -> grav acc on surface: [m/sec^2]
# - beta -> drag coefficient [1/sec]
# - pos[ 2 ] -> [ x , y ] position in [m]
# - vel[ 2 ] -> [ vx , vy ] velocity in [m/sec]
# Output:
# - rhs_pos[ 2 ] -> the RHS of position -> [ f_x , f_y ]
# - rhs_vel[ 2 ] -> the RHS of velocity -> [ f_vx , f_vy ]
def rhs_ballistic( g , beta , pos , vel ):

    rhs_pos = [ vel[ 0 ] , # f_x = v_x
                vel[ 1 ] ] # f_y = v_y
    rhs_vel = [ - beta*vel[ 0 ] , # f_vx = - beta*v_x
                - g - beta*vel[ 1 ] ] # f_vy = - g - beta*v_y
    
    return rhs_pos , rhs_vel

# Integrate Ballistic Trajectory through Verlet's method
# Inputs:
# - g -> grav acceleration [m/sec^2]
# - beta -> drag coeff [1/sec]
# - pos_0[ 2 ] -> [ x_0 , y_0 ] initial position [m]
# - vel_0[ 2 ] -> [ vx_0 , vy_0 ] initial velocity [m/sec]
# - t_params[ 3 ] -> [ t_i , t_f , Npoints ] 
# -- initial time (t_i) in [sec]
# -- final time (t_f) in [sec]
# -- Number of points Npoints [-]
# Outputs:
# - t_arr[ Npoints ] -> time array from t_i to t_f in [sec]
# - pos_arr[ 2 ][ Npoints ] -> position solution in [m], first index 0 for x, 1 for y 
# - vel_arr[ 2 ][ Npoints ] -> velocity solution in [m], first index 0 for x, 1 for y
def Verlet_Ballistic( g , beta , pos_0 , vel_0 , t_params ):

    # Take out the initial and final time (1st and 2nd element of t_params)
    t_i = t_params[ 0 ]
    t_f = t_params[ 1 ]
    # Take out the number of points which is 3d element in the t_params
    Npoints = t_params[ 2 ]
    # Compute the time step dt 
    dt = ( t_f - t_i )/( Npoints - 1 ) 
    # NOTE: Intervals are 1 less tha nnumber of points -> | | | has 3 points, 2 interfvals

    # Create empty arrays which will hold the solution
    pos = np.zeros( ( 2 , Npoints ) )
    vel = np.zeros( ( 2 , Npoints ) )
    # Create the time array
    t_arr = np.linspace( t_i , t_f , Npoints )

    # Assign the initial position and velocity
    pos[ 0 ][ 0 ] = pos_0[ 0 ] # initial position in X
    pos[ 1 ][ 0 ] = pos_0[ 1 ] # initial position in Y
    vel[ 0 ][ 0 ] = vel_0[ 0 ] # initial velocity in X (Vx)
    vel[ 1 ][ 0 ] = vel_0[ 1 ] # initial velocity in Y (Vy)

    # Integration loop -> go over each time step and extrapolate the next
    for i in range( 0 , Npoints - 1 ):
        # Compute RHS in the first point (i)
        rhs_pos, rhs_vel = rhs_ballistic( g , beta , np.transpose( pos )[ i ] , np.transpose( vel )[ i ] )
        # Compute half-velocity step -> intermediate step # NOTE: dt/2 -> half-step
        vx_half = vel[ 0 ][ i ] + rhs_vel[ 0 ]*dt/2.0
        vy_half = vel[ 1 ][ i ] + rhs_vel[ 1 ]*dt/2.0
        # Compute RHS in the intermediate point (i+1/2)
        rhs_pos, rhs_vel = rhs_ballistic( g , beta , np.transpose( pos )[ i ] , [ vx_half , vy_half ] ) # NOTE: We are passing the half-point velocity
        # Compute full position step based on v_half
        pos[ 0 ][ i + 1 ] = pos[ 0 ][ i ] + rhs_pos[ 0 ]*dt
        pos[ 1 ][ i + 1 ] = pos[ 1 ][ i ] + rhs_pos[ 1 ]*dt
        # Compute RHS in the last point (i+1)
        rhs_pos, rhs_vel = rhs_ballistic( g , beta , np.transpose( pos )[ i + 1 ] , [ vx_half , vy_half ] ) # NOTE: We are passing the next moment in position
        # Compute the full velocity step based on v_half and x[ i + 1 ]
        vel[ 0 ][ i + 1 ] = vx_half + rhs_vel[ 0 ]*dt/2.0
        vel[ 1 ][ i + 1 ] = vy_half + rhs_vel[ 1 ]*dt/2.0
        # IF YOU PASS UNDERGROUND - STOP Y
        # Perfect plastic collision
        '''
        if pos[ 1 ][ i + 1 ] < 0.0:
            pos[ 1 ][ i + 1 ] = 0.0
            vel[ 1 ][ i + 1 ] = 0.0
        '''
        '''
        # IF YOU PASS UNDERGROUND - BOUNCE
        # Perfect ellastic collision
        if pos[ 1 ][ i + 1 ] < 0.0:
            pos[ 1 ][ i + 1 ] *= - 1.0
            vel[ 1 ][ i + 1 ] *= - 1.0
        '''
        # IF YOU PASS UNDERGROUND - BOUNCE BUT LOSE ENERGY
        # Realistic collision
        if pos[ 1 ][ i + 1 ] < 0.0:
            pos[ 1 ][ i + 1 ] *= - 1.0
            vel[ 1 ][ i + 1 ] *= - 0.5

    # My pos[ 2 ][ Npoints ] and vel[ 2 ][ Npoints ] are populated with the solution - return them
    return t_arr , pos , vel
