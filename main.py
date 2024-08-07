from Integrators import *
import matplotlib.pyplot as plt
from numpy import pi
import numpy as np

if __name__ == "__main__":

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
    # Real solution to compare with
    x_real = np.cos( np.sqrt( k/m )*t_arr )

    # 2D Plot of position vs time
    fig, ax = plt.subplots()
    ax.scatter( t_arr , x_eul , color = "red" , linestyle = "dashed" , label = r"Euler sim" , linewidth = 4 )
    ax.scatter( t_arr , x_ver , color = "green" , linestyle = "dashed" , label = r"Verlet sim" , linewidth = 4 )
    ax.plot( t_arr , x_real , color = "blue" , linestyle = "solid" , label = r"real X" , linewidth = 2 )

    ax.set_xlabel( r"Time [sec]" )
    ax.set_ylabel( r"X [m]" )

    ax.legend( loc = "upper right" )
    plt.grid( True )
    plt.show()

    # Create error arrays for Euler and Verlet to compare [%]
    err_eul = np.zeros( Npoints )
    err_ver = np.zeros( Npoints )

    for i in range( 0 , Npoints ):
        err_eul[ i ] = np.abs( ( x_real[ i ] - x_eul[ i ] )/x_real[ i ] )*100.0
        err_ver[ i ] = np.abs( ( x_real[ i ] - x_ver[ i ] )/x_real[ i ] )*100.0

    # 2D Plot example
    fig, ax = plt.subplots()
    ax.plot( t_arr , err_eul , color = "red" , linestyle = "dashed" , label = r"Euler err" , linewidth = 4 )
    ax.plot( t_arr , err_ver , color = "green" , linestyle = "dashed" , label = r"Verlet err" , linewidth = 4 )


    ax.set_xlabel( r"Time [sec]" )
    ax.set_ylabel( r"Error in [%]" )

    ax.legend( loc = "upper right" )
    plt.grid( True )
    plt.show()
    print( "Max error for Euler = " + str( max( err_eul ) ) + " [%]" )
    print( "Max error for Verlet = " + str( max( err_ver ) ) + " [%]" )