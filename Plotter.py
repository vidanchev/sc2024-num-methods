import numpy as np 
import csv
import matplotlib.pyplot as plt
from Integrators import *

# Parse a .csv file with satellite time, position and velocity (PVT) for plotting
# Inputs:
# - filename of the .csv file where the results are saved
# NOTE: File structure is expected to be "Time (UTCG),x (km),y (km),z (km),vx (km/sec),vy (km/sec),vz (km/sec)"
# Output:
# It will return arrays of time, pos, vel as:
# - time[ Npoints ]: time in [sec] from simulation start
# - pos[ 3 ][ Npoints ]: 3D position in [km]
# - vel[ 3 ][ Npoints ]: 3D velocity in [km/s]
def parse_orbit_data( filename ):

    fp = open( filename , "r" )

    if fp.readable( ):
        data = csv.reader( fp )
        lst = [ ]
        for line in data:
            lst.append( line )
        ndata = len( lst ) - 1

        time = np.zeros( ndata )
        pos = np.zeros( ( 3 , ndata ) )
        vel = np.zeros( ( 3 , ndata ) )

        for i in range( 0 , ndata ):
            time[ i ] = float( lst[ i + 1 ][ 0 ] )
            for j in range( 0 , 3 ):
                pos[ j ][ i ] = float( lst[ i + 1 ][ j + 1 ] )
                vel[ j ][ i ] = float( lst[ i + 1 ][ j + 4 ] )
    else:
        print( "Unreadable data, something's wrong with the file " + filename )
    
    return time, pos, vel

if __name__ == "__main__":

    # File simulated from GMAT
    file_name = "Satellite_PVT_GMAT.csv"

    time , pos , vel = parse_orbit_data( file_name )

    # Run the simulation with the same initial data
    pos_0 = np.transpose( pos )[ 0 ] # Initial position [km]
    vel_0 = np.transpose( vel )[ 0 ] # Initial velocity [km/sec]
    t_params = [ time[ 0 ] , time[ -1 ] , 8641 ]
    t_RK, pos_RK, vel_RK = Kepler_RK4( pos_0 , vel_0 , t_params )

    # 3D Plot example
    # Data for a three-dimensional line
    ax = plt.axes( projection = "3d" )
    ax.plot3D( pos[ 0 ] , pos[ 1 ] , pos[ 2 ] , "green" )
    ax.plot3D( np.transpose( pos_RK )[ 0 ] , np.transpose( pos_RK )[ 1 ] , np.transpose( pos_RK )[ 2 ] , "red" )

    plt.show( )

    # 2D Plot example
    fig, ax = plt.subplots()
    ax.plot( time , pos[ 0 ] , color = "green" , linestyle = "dashed" , label = r"pos X STK" )
    ax.plot( time , pos[ 1 ] , color = "red" , linestyle = "dashed" , label = r"pos Y STK" )
    ax.plot( time , pos[ 2 ] , color = "blue" , linestyle = "dashed" , label = r"pos Z STK" )
    ax.scatter( t_RK , np.transpose( pos_RK )[ 0 ] , color = "green" , label = r"pos X RK" )
    ax.scatter( t_RK , np.transpose( pos_RK )[ 1 ] , color = "red" , label = r"pos Y RK" )
    ax.scatter( t_RK , np.transpose( pos_RK )[ 2 ] , color = "blue" , label = r"pos Z RK" )

    ax.set_xlabel( r"Time [sec]" )
    ax.set_ylabel( r"Position in [km]" )

    ax.legend( loc = "upper right" )

    #ax.set_ylim( 0.0 , pi/2.0 )
    plt.grid( True )

    #fig.savefig( "fig_name.pdf" , format = "pdf" )
    
    plt.show()

    # Error between us and STK
    pos_err = np.zeros( len( pos_RK ) )
    for i in range( 0 , len( pos_RK ) ):
        pos_err[ i ] = np.sqrt( ( pos_RK[ i ][ 0 ] - pos[ 0 ][ i ] )**2.0 + ( pos_RK[ i ][ 1 ] - pos[ 1 ][ i ] )**2.0 + ( pos_RK[ i ][ 2 ] - pos[ 2 ][ i ] )**2.0 )
        # We are computing sqrt( ( x_RK - x_R )^2 + ( y_RK - y_R )^2 + ( z_RK - z_R )^2 ) which is our error

    # 2D Plot of error
    fig, ax = plt.subplots()
    ax.plot( time , pos_err , color = "blue" , linestyle = "solid" , label = r"pos error" )

    ax.set_xlabel( r"Time [sec]" )
    ax.set_ylabel( r"Error in [km]" )

    ax.legend( loc = "upper right" )

    #ax.set_ylim( 0.0 , pi/2.0 )
    plt.grid( True )

    #fig.savefig( "fig_name.pdf" , format = "pdf" )
    
    plt.show()
