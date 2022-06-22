#include <stdio.h>
#include <math.h>
#include "mpi.h"

/* This example handles a m x n mesh, on p (always even) processors. */
#define m 18 // Assume m is row
#define n 24 // Assume n is column

int main(int argc, char **argv )
{
    int        rank, size, i, j, itcnt;
    int        j_first, j_last;
    MPI_Status status;
    double     diffnorm, gdiffnorm;
    double     start, finish;

    MPI_Init( &argc, &argv );
    start = MPI_Wtime();

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    
    if (size % 2 == 1) MPI_Abort( MPI_COMM_WORLD, 1 );
    if (n % (size / 2) != 0) MPI_Abort( MPI_COMM_WORLD, 1 );
    
    int n_row = m/2;
    int n_col = n/(size/2);
    
    double     xlocal[n_row+1][n_col+2];
    double     xnew[n_row+1][n_col+2];
    double     receive_col[n_row+1];
    double     receive_row[n_col+2];

    // Need to find appropriate first and last
	j_first = 1;
	j_last = n_col;
	if (rank == 0 || rank == 1) {
		j_first++;
	}
	if (rank == size - 1 || rank == size - 2) {
		j_last--;
	}

    /* Fill the data as specified */
    for (i=1; i<n_row; i++) {
		for (j=j_first - 1; j<=j_last + 1; j++) {
			xlocal[i][j] = rank;
		}
	}
    for (j=j_first - 1; j<=j_last + 1; j++) {
		xlocal[0][j] = -1;
		xlocal[n_row][j] = -1;
    }
    
    // Creating data type for column
    MPI_Datatype column_type;
	MPI_Type_vector(n_row+1, 1, n_col+2, MPI_DOUBLE, &column_type);
	MPI_Type_commit(&column_type);
	
	// Creating data type for row
	MPI_Datatype row_type;
	MPI_Type_vector(1, n_col+2, 0, MPI_DOUBLE, &row_type);
	MPI_Type_commit(&row_type);

    itcnt = 0;
    do {
		/* Send up unless I'm at the top, then receive from below */
		/* Note the use of xlocal[i] for &xlocal[i][0] */
		if (rank % 2 == 1) {
			MPI_Send(&xlocal[1][0], 1, row_type, rank - 1, 0, MPI_COMM_WORLD );
		}
		if (rank % 2 == 0) {
			MPI_Recv(&receive_row, n_col+2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status );
			
			for (int i = 0; i < n_col+2; i++) {
				xlocal[n_row][i] = receive_row[i];
			}
		}
		
		/* Send down unless I'm at the bottom */
		if (rank % 2 == 0) {
			MPI_Send(&xlocal[n_row-1][0], 1, row_type, rank + 1, 1, MPI_COMM_WORLD );
		}
		if (rank % 2 == 1) {
			MPI_Recv(&receive_row, n_col+2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status );
			
			for (int i = 0; i < n_col+2; i++) {
				xlocal[0][i] = receive_row[i];
			}
		}
		
		/* Send to right unless I'm Rank Size - 1 and Size - 2*/
		if (rank < size - 2) {
			MPI_Send(&xlocal[0][n_col], 1, column_type, rank + 2, 2, MPI_COMM_WORLD);
		}
		
		if (rank > 1) {
			MPI_Recv(&receive_col, n_row+1, MPI_DOUBLE, rank - 2, 2, MPI_COMM_WORLD, &status );
			
			for (int i = 0; i < n_row+1; i++) {
				xlocal[i][0] = receive_col[i];
			}
		}	
		
		/* Send to left unless I'm Rank 0 and 1*/
		if (rank > 1) {
			MPI_Send(&xlocal[0][1], 1, column_type, rank - 2, 3, MPI_COMM_WORLD);
		}
		
		if (rank < size - 2) {
			MPI_Recv(&receive_col, n_row+1, MPI_DOUBLE, rank + 2, 3, MPI_COMM_WORLD, &status );

			for (int i = 0; i < n_row+1; i++) {
				xlocal[i][n_col+1] = receive_col[i];
			}
		}
		
		/* Compute new values (but not on boundary) */
		itcnt ++;
		diffnorm = 0.0;
		for (i=1; i<n_row; i++) {
			for (j=j_first; j<=j_last; j++) {
				xnew[i][j] = (xlocal[i][j+1] + xlocal[i][j-1] + xlocal[i+1][j] + xlocal[i-1][j]) / 4.0;
				diffnorm += (xnew[i][j] - xlocal[i][j]) * (xnew[i][j] - xlocal[i][j]);
			}
		}
		/* Only transfer the interior points */
		for (i=1; i<n_row; i++) {
			for (j=j_first; j<=j_last; j++) {
				xlocal[i][j] = xnew[i][j];
			}
		}
			
		MPI_Allreduce( &diffnorm, &gdiffnorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
		
		gdiffnorm = sqrt( gdiffnorm );

		if (rank == 0) {
			printf( "At iteration %d, diff is %e\n", itcnt, gdiffnorm );
		}
    } while (gdiffnorm > 1.0e-2);
    
    finish = MPI_Wtime();
    
    printf("Rank: %d\nElapsed Time: %.4f seconds\n", rank, finish-start);
    
    MPI_Finalize( );
    return 0;
}

