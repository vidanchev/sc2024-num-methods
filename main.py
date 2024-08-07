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
    Npoints = 1000000 # Number of points
    t_params = [ 0.0 , # t_i 
                 4.0*np.pi , # t_f
                 Npoints ]

    t_arr, x_arr, v_arr = Euler_Oscillator( k , m , x_0 , v_0 , t_params )

    x_real = np.cos( np.sqrt( k/m )*t_arr )

    # 2D Plot of position vs time
    fig, ax = plt.subplots()
    ax.plot( t_arr , x_arr , color = "red" , linestyle = "dashed" , label = r"pos X" , linewidth = 4 )
    ax.plot( t_arr , x_real , color = "green" , linestyle = "solid" , label = r"pos X" , linewidth = 4 )

    ax.set_xlabel( r"Time [sec]" )
    ax.set_ylabel( r"X [m]" )

    ax.legend( loc = "upper right" )
    plt.grid( True )
    plt.show()

    err_arr = np.zeros( Npoints )
    for i in range( 0 , Npoints ):
        err_arr[ i ] = np.abs( ( x_real[ i ] - x_arr[ i ] )/x_real[ i ] )*100.0

    # 2D Plot example
    fig, ax = plt.subplots()
    ax.plot( t_arr , err_arr , color = "red" , linestyle = "dashed" , label = r"pos X" , linewidth = 4 )

    ax.set_xlabel( r"Time [sec]" )
    ax.set_ylabel( r"X [m]" )

    ax.legend( loc = "upper right" )
    plt.grid( True )
    plt.show()
    print( "Max error = " + str( max( err_arr ) ) + " [%]" )