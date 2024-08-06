import numpy as np 
import csv
import matplotlib.pyplot as plt

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

    # 3D Plot example
    # Data for a three-dimensional line
    ax = plt.axes( projection = "3d" )
    ax.plot3D( pos[ 0 ] , pos[ 1 ] , pos[ 2 ] , "blue" )
    plt.show( )

    # 2D Plot example
    fig, ax = plt.subplots()
    ax.plot( time , pos[ 0 ] , color = "green" , linestyle = "solid" , label = r"pos X" )
    ax.plot( time , pos[ 1 ] , color = "red" , linestyle = "solid" , label = r"pos Y" )
    ax.plot( time , pos[ 2 ] , color = "blue" , linestyle = "solid" , label = r"pos Z" )

    ax.set_xlabel( r"Time [sec]" )
    ax.set_ylabel( r"Position in [km]" )

    ax.legend( loc = "upper right" )

    #ax.set_ylim( 0.0 , pi/2.0 )
    plt.grid( True )

    #fig.savefig( "fig_name.pdf" , format = "pdf" )
    
    plt.show()
