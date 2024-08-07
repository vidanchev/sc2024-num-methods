import numpy as np

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
        # x_(i+1) = x_(i) + v_(i)*dt
        x_arr[ i + 1 ] = x_arr[ i ] + v_arr[ i ]*dt
        # v_(i+1) = v_(i) + f_(i)*dt, where f_(i) = -k/m*x_(i) for oscillator
        v_arr[ i + 1 ] = v_arr[ i ] + ( - k / m )*x_arr[ i ]*dt

    # My x_arr and v_arr are populated with the solution - return them
    return t_arr , x_arr , v_arr