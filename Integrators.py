import numpy as np

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
