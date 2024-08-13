from Integrators import *
import matplotlib.pyplot as plt
from numpy import pi
import numpy as np

if __name__ == "__main__":
    
    pos_0 = [ 7000 , 0.000001 , -0.001608 ] # Initial position [km]
    vel_0 = [ 0.000002 , 1.310359 , 7.431412 ] # Initial velocity [km/sec]
    t_params = [ 0.0 , 100.0 , 10 ]
    Kepler_RK4( pos_0 , vel_0 , t_params )

    '''
    # Example parameters for ballistic
    g = 9.8 # [m/sec^2]
    beta = 0.01 # [1/sec] drag
    pos_0 = [ 0 , 0 ]
    alp = 30.0 # Alpha - angle to the horizon [deg]
    v_0 = 100.0 # Norm of the velocity [m/sec]
    # Compute velocity vector from trig
    vel_0 = [ v_0*np.cos( alp*np.pi/180.0 ) ,  # Converting Vx from V and Alpha
              v_0*np.sin( alp*np.pi/180.0 ) ]  # Converting Vy from V and Alpha
    Npoints = 1000 # Number of points
    t_params = [ 0.0 , # t_i 
                 100.0 , # t_f
                 Npoints ]

    t_arr, pos, vel = Verlet_Ballistic( g , beta , pos_0 , vel_0 , t_params )

    x_real = vel_0[ 0 ]*t_arr + pos_0[ 0 ]
    y_real = vel_0[ 1 ]*t_arr + pos_0[ 1 ] - 0.5*g*t_arr*t_arr

    fig, ax = plt.subplots()
    ax.scatter( pos[ 0 ] , pos[ 1 ] , color = "red" , label = r"Verlet sim" , linewidth = 4 )
    #ax.plot( x_real , y_real , color = "green" , linestyle = "solid" , label = r"real traj" , linewidth = 2 )

    ax.set_xlabel( r"X [m]" )
    ax.set_ylabel( r"Y [m]" )

    ax.legend( loc = "upper right" )
    plt.grid( True )
    plt.show()
    '''
    
    '''
    # Example parameters for which we know the solution
    k = 1.0 # Spring constant [N/kg]
    m = 1.0 # Mass [kg]
    x_0 = 1.0 # Initial position [m]
    v_0 = 0.0 # Initial speed [m/sec]
    Npoints = 1000 # Number of points
    t_params = [ 0.0 , # t_i 
                 4.0*np.pi , # t_f
                 Npoints ]
    
    # Obtain results with the Euler method 
    t_arr, x_eul, v_eul = Euler_Oscillator( k , m , x_0 , v_0 , t_params )
    # Obtain results with the Verlet method 
    t_arr, x_ver, v_ver = Verlet_Oscillator( k , m , x_0 , v_0 , t_params )
    # Obtain results with the RK(4) method 
    t_arr, x_rk, v_rk = RK4_Oscillator( k , m , x_0 , v_0 , t_params )
    # Real solution to compare with
    x_real = np.cos( np.sqrt( k/m )*t_arr )

    # 2D Plot of position vs time
    fig, ax = plt.subplots()
    ax.scatter( t_arr , x_eul , color = "red" , linestyle = "dashed" , label = r"Euler sim" , linewidth = 3 )
    ax.scatter( t_arr , x_ver , color = "green" , linestyle = "dashed" , label = r"Verlet sim" , linewidth = 3 )
    ax.scatter( t_arr , x_rk , color = "blue" , linestyle = "dashed" , label = r"RK(4) sim" , linewidth = 3 )
    ax.plot( t_arr , x_real , color = "black" , linestyle = "solid" , label = r"real X" , linewidth = 2 )

    ax.set_xlabel( r"Time [sec]" )
    ax.set_ylabel( r"X [m]" )

    ax.legend( loc = "upper right" )
    plt.grid( True )
    plt.show()

    # Create error arrays for Euler and Verlet to compare [%]
    err_eul = np.zeros( Npoints )
    err_ver = np.zeros( Npoints )
    err_rk = np.zeros( Npoints )

    for i in range( 0 , Npoints ):
        err_eul[ i ] = np.abs( ( x_real[ i ] - x_eul[ i ] )/x_real[ i ] )*100.0
        err_ver[ i ] = np.abs( ( x_real[ i ] - x_ver[ i ] )/x_real[ i ] )*100.0
        err_rk[ i ] = np.abs( ( x_real[ i ] - x_rk[ i ] )/x_real[ i ] )*100.0

    # 2D Plot example
    fig, ax = plt.subplots()
    ax.plot( t_arr , err_eul , color = "red" , linestyle = "dashed" , label = r"Euler err" , linewidth = 4 )
    ax.plot( t_arr , err_ver , color = "green" , linestyle = "dashed" , label = r"Verlet err" , linewidth = 4 )
    ax.plot( t_arr , err_rk , color = "blue" , linestyle = "dashed" , label = r"Verlet err" , linewidth = 4 )

    ax.set_xlabel( r"Time [sec]" )
    ax.set_ylabel( r"Error in [%]" )
    ax.set_yscale( "log" )

    ax.legend( loc = "upper right" )
    plt.grid( True )
    plt.show()
    print( "Max error for Euler = " + str( max( err_eul ) ) + " [%]" )
    print( "Max error for Verlet = " + str( max( err_ver ) ) + " [%]" )
    print( "Max error for RK(4) = " + str( max( err_rk ) ) + " [%]" )
    '''