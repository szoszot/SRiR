#include "mpi.h"
#include <iostream>

#include <math.h>
#include <time.h>
#include <cstdlib>


using namespace std;

#define BORN 1
#define DIES 0

struct Komorka
{
    int status, tstatus;
    int naGranicy;
    int wylosowana;
    int numer, tnumer;
}K;

MPI_Datatype MPI_KOM;

/* The Life function */
double life(int matrix_size, int ntimes, MPI_Comm comm)
{
    int      naGranicy = 0, grainCounter = 0;
    int      rank, size ;
    int      next, prev ;
    int      i, j, k;
    int      mysize, sum ;
    Komorka   **matrix, **temp, **addr ;
    double   slavetime, totaltime, starttime ;
    int      my_offset;

    /* Determine size and my rank in communicator */
    MPI_Comm_size(comm, &size) ;
    MPI_Comm_rank(comm, &rank) ;

    /* Set neighbors */
    if (rank == 0)
        prev = MPI_PROC_NULL;
    else
        prev = rank-1;
    if (rank == size - 1)
        next = MPI_PROC_NULL;
    else
        next = rank+1;

    /* Determine my part of the matrix */
    mysize = matrix_size/size + ((rank < (matrix_size % size)) ? 1 : 0 ) ;
    my_offset = rank * (matrix_size/size);
    if (rank > (matrix_size % size)) my_offset += (matrix_size % size);
    else                             my_offset += rank;

    matrix = new Komorka * [mysize + 2];
    temp= new Komorka * [mysize + 2];
    for (int i = 0; i<mysize + 2; i++) {
        matrix[i] = new Komorka [matrix_size+2];
        temp[i] = new Komorka [matrix_size+2];
    }

    /* Initialize the boundaries of the life matrix */
    for (j = 0; j < matrix_size+2; j++) {
        matrix[0][j].wylosowana = matrix[mysize+1][j].wylosowana  = temp[0][j].wylosowana  = temp[mysize+1][j].wylosowana  = 0 ;
        matrix[0][j].status     = matrix[mysize+1][j].status      = temp[0][j].status      = temp[mysize+1][j].status      = 0 ;
        matrix[0][j].naGranicy  = matrix[mysize+1][j].naGranicy   = temp[0][j].naGranicy   = temp[mysize+1][j].naGranicy   = 0 ;
        matrix[0][j].numer      = matrix[mysize+1][j].numer       = temp[0][j].numer       = temp[mysize+1][j].numer       = 0 ;
        matrix[0][j].tnumer      = matrix[mysize+1][j].tnumer       = temp[0][j].tnumer       = temp[mysize+1][j].tnumer       = 0 ;
        matrix[0][j].tstatus     = matrix[mysize+1][j].tstatus      = temp[0][j].tstatus      = temp[mysize+1][j].tstatus      = 0 ;
    }
    for (i = 0; i < mysize+2; i++) {
        matrix[i][0].wylosowana  = matrix[i][matrix_size+1].wylosowana  = temp[i][0].wylosowana  = temp[i][matrix_size+1].wylosowana  = 0 ;
        matrix[i][0].status      = matrix[i][matrix_size+1].status      = temp[i][0].status      = temp[i][matrix_size+1].status      = 0 ;
        matrix[i][0].naGranicy   = matrix[i][matrix_size+1].naGranicy   = temp[i][0].naGranicy   = temp[i][matrix_size+1].naGranicy   = 0 ;
        matrix[i][0].numer       = matrix[i][matrix_size+1].numer       = temp[i][0].numer       = temp[i][matrix_size+1].numer       = 0 ;
        matrix[i][0].tnumer       = matrix[i][matrix_size+1].tnumer       = temp[i][0].tnumer       = temp[i][matrix_size+1].tnumer       = 0 ;
        matrix[i][0].tstatus      = matrix[i][matrix_size+1].tstatus      = temp[i][0].tstatus      = temp[i][matrix_size+1].tstatus      = 0 ;
    }
    /* Initialize the life matrix */
#pragma omp parallel for private(i,j) shared (matrix)
    for (i = 1; i <= mysize; i++)  {
        srand((long)(1000^(i-1+mysize))) ;
        for (j = 1; j<= matrix_size; j++) {
            matrix[i][j].wylosowana = 0;
            matrix[i][j].naGranicy = 0;
            matrix[i][j].tnumer = 0;
            matrix[i][j].tstatus = 0;

            if (rand() < (RAND_MAX/20)) {
                matrix[i][j].numer = 1000 * rank + grainCounter++ ;
                matrix[i][j].status = 1;
            }
            else {
                matrix[i][j].numer = 0 ;
                matrix[i][j].status = 0;
            }
        }
    }

    //wypisz
    if(matrix_size < 11) {
        for (i = 1; i <= mysize; i++)  {
            for (j = 1; j<= matrix_size; j++) {
                cout<<matrix[i][j].numer<<"\t";
            }
            cout<<endl;
        }
    }
    /*ROZROST STRUKRURY*/
    int rs;
    int k_roz[6];
    int it = 0;
    int numer_koloru_tablica = 0;
    /*ROZROST - ZAPETLIC!*/
    while(it < matrix_size*matrix_size) {
        MPI_Request      req[4];
        MPI_Status       status[4];

        /* Send and receive boundary information */
        MPI_Isend(&matrix[1][0],matrix_size+2,MPI_KOM,prev,0,comm,req);
        MPI_Irecv(&matrix[0][0],matrix_size+2,MPI_KOM,prev,0,comm,req+1);
        MPI_Isend(&matrix[mysize][0],matrix_size+2,MPI_KOM,next,0,comm,req+2);
        MPI_Irecv(&matrix[mysize+1][0],matrix_size+2,MPI_KOM,next,0,comm,req+3);
        MPI_Waitall(4, req, status);
        for (i = 1; i <= mysize; i++)  {
            for (j = 1; j<= matrix_size; j++){
                if(matrix[i][j].status == 0 ) {
                    rs=rand() % 2;
                    if(rs){
                        //pierwsza sk_rozladowa numeru
                        k_roz[0] = matrix[i][j-1]		.numer;
                        k_roz[1] = matrix[i+1][j-1]		.numer;
                        k_roz[2] = matrix[i+1][j]		.numer;
                        k_roz[3] = matrix[i][j+1]		.numer;
                        k_roz[4] = matrix[i-1][j+1]		.numer;
                        k_roz[5] = matrix[i-1][j]		.numer;
                    }
                    //lewy gorny/prawy dolny pusty - j/w
                    else {
                        k_roz[0] = matrix[i-1][j-1]		.numer;
                        k_roz[1] = matrix[i][j-1]		.numer;
                        k_roz[2] = matrix[i+1][j]		.numer;
                        k_roz[3] = matrix[i+1][j+1]		.numer;
                        k_roz[4] = matrix[i][j+1]		.numer;
                        k_roz[5] = matrix[i-1][j]		.numer;
                    }

                    //wyszukanie koloru max

                    int temp_licznik1 = 0;
                    int temp_licznik2 = 0;

                    for ( int l = 0; l < 6; l++) {

                        for (int m = 0; m < 6; m++) {
                            if(k_roz[l] == k_roz[m]&& k_roz[l] != 0) temp_licznik1++;
                        }

                        if(temp_licznik1 > temp_licznik2) {
                            temp_licznik2 = temp_licznik1;
                            numer_koloru_tablica = l;
                        }
                        else if (temp_licznik2 == temp_licznik1) {
                            int r = rand() % 2;		// 0 lub 1
                            if(r) temp_licznik2 = temp_licznik1;
                            numer_koloru_tablica = l;
                        }
                        temp_licznik1 = 0;
                    }

                    //rozrost
                    if (matrix[i][j].status == 0 && temp_licznik2 > 0) {
                        matrix[i][j].tnumer = k_roz[numer_koloru_tablica];
                        matrix[i][j].tstatus = 1;
                    }
                    else if(temp_licznik2 ==0) {
                        matrix[i][j].tstatus = 0;
                    }

                }
            }
        }

        //update z rozrosrem

        for (i = 1; i <= mysize; i++)  {
            for (j = 1; j<= matrix_size; j++){
                if( matrix[i][j].status == 0) {
                    matrix[i][j].status   = matrix[i][j].tstatus;
                    matrix[i][j].numer = matrix[i][j].tnumer;
                }
                //zerowanie tempow
                matrix[i][j].tstatus=0;
                matrix[i][j].tnumer = 0;

            }
        }
        it++;
    }

    /*KONIEC ZAPETLANIA */
    /*USTAWIAMY GRANICE*/
#pragma omp parallel for private (i,j) shared (matrix)
    for( int i = 1; i <= mysize; i++ ){
        for( int j = 1; j<= matrix_size; j++ ){
            //ustawianie, ktora komorka jest na granicy - spr numery sasiadow
            if(		matrix[i][j].numer != matrix[i-1][j-1]	.numer
                ||	matrix[i][j].numer != matrix[i][j-1]	.numer
                ||	matrix[i][j].numer != matrix[i+1][j-1]	.numer
                ||	matrix[i][j].numer != matrix[i+1][j]	.numer
                ||	matrix[i][j].numer != matrix[i+1][j+1]	.numer
                ||	matrix[i][j].numer != matrix[i][j+1]	.numer
                ||	matrix[i][j].numer != matrix[i-1][j+1]	.numer
                ||	matrix[i][j].numer != matrix[i-1][j]	.numer
                )
                //+reszta numerow, ale wystarczy 1 skladowa
                {
                    matrix[i][j].naGranicy = 1;
                    naGranicy++;
                }
        }
    }
    cout<<"\nNa granicy: "<<naGranicy<<endl;

    cout<<"\nPo rozroscie: \n";
    if(matrix_size < 11) {
        for (i = 1; i <= mysize; i++)  {
            for (j = 1; j<= matrix_size; j++) {
                cout<<matrix[i][j].numer<<"\t";
            }
            cout<<endl;
        }
    }


    /* Play the game of life for given number of iterations */
    starttime = MPI_Wtime() ;
    int zliczacz = 0;
    for (k = 0; k < ntimes; k++) {

        MPI_Request      req[4];
        MPI_Status       status[4];

        /* Send and receive boundary information */
        MPI_Isend(&matrix[1][0],matrix_size+2,MPI_KOM,prev,0,comm,req);
        MPI_Irecv(&matrix[0][0],matrix_size+2,MPI_KOM,prev,0,comm,req+1);
        MPI_Isend(&matrix[mysize][0],matrix_size+2,MPI_KOM,next,0,comm,req+2);
        MPI_Irecv(&matrix[mysize+1][0],matrix_size+2,MPI_KOM,next,0,comm,req+3);
        MPI_Waitall(4, req, status);

        /*REKRYSTALIZACJA */
#pragma omp parallel shared (zliczacz,matrix)
        for(zliczacz = 0; zliczacz < naGranicy/2;){
            int numer_koloru_tablica = 0;
            int i = rand()%mysize+1;
            int j = rand()%matrix_size+1;
            int k0[8], k0T;
            if(matrix[i][j].wylosowana == 0 && matrix[i][j].naGranicy == 1 ) {
                k0[0] = matrix[i-1][j-1]	.numer;
                k0[1] = matrix[i][j-1]		.numer;
                k0[2] = matrix[i+1][j-1]	.numer;
                k0[3] = matrix[i+1][j]		.numer;
                k0[4] = matrix[i+1][j+1]	.numer;
                k0[5] = matrix[i][j+1]      .numer;
                k0[6] = matrix[i-1][j+1]	.numer;
                k0[7] = matrix[i-1][j]		.numer;

                int losowanie = rand()%7;
                k0T = k0[losowanie];

                int temp_licznik1 = 0;
                int temp_licznik2 = 0;
                int tyleSamo = 0;
                bool flaga = false;
                for ( int l = 0; l < 8; l++) {
                    for (int m = 0; m < 8; m++) {
                        if(k0[l] == k0[m] && k0[0] != 0) {
                            temp_licznik1++;
                        }
                        if(temp_licznik1 == 4 ) tyleSamo++;		//gdy 4 ma taki sam kolor
                    }
                    if(temp_licznik1 > temp_licznik2) {
                        temp_licznik2 = temp_licznik1;
                        numer_koloru_tablica = l;
                    }
                    temp_licznik1 = 0;
                    if (tyleSamo == 2 ) flaga = true;			//gdy mamy 2 kolory po 4 sasiadow
                    tyleSamo = 0;
                }

                if(numer_koloru_tablica == losowanie && !flaga) {
                    matrix[i][j].numer = k0T;

                    matrix[i][j].wylosowana = 1;
                    zliczacz++;

                }
                else if (flaga) {
                    int r = rand()%7;
                    matrix[i][j].numer = k0[r];

                    matrix[i][j].wylosowana = 1;
                    zliczacz++;
                }
            }
        }
        zliczacz = 0;

        //cout<<"\nhellol!\n";
        for (i = 0; i <= mysize; i++)  {
            for (j = 0; j<= matrix_size; j++) {
                matrix[i][j].wylosowana = 0;
            }
        }
    }

    /* Return the average time taken/processor */
    //wypisz
    cout<<"\nPo monte carlo: \n";
    if(matrix_size < 11) {
        for (i = 1; i <= mysize; i++)  {
            for (j = 1; j<= matrix_size; j++) {
                cout<<matrix[i][j].numer<<"\t";
            }
            cout<<endl;
        }
    }
    slavetime = MPI_Wtime() - starttime;
    MPI_Reduce (&slavetime, &totaltime, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    return (totaltime/(double)size);
}


int main(int argc, char *argv[])
{
    int rank, N = 8, iters = 100;
    double time ;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    /*DEFINICJA TYPU DANYCH*/

    int tab_dl_blokow[6]={1,1,1,1,1,1};
    MPI_Datatype tab_typow[6] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };
    MPI_Aint tab_odstepow[6], podstawa;
    MPI_Get_address(&K.status, &tab_odstepow[0]);
    MPI_Get_address(&K.naGranicy, &tab_odstepow[1]);
    MPI_Get_address(&K.wylosowana, &tab_odstepow[2]);
    MPI_Get_address(&K.numer, &tab_odstepow[3]);
    MPI_Get_address(&K.tnumer, &tab_odstepow[4]);
    MPI_Get_address(&K.tstatus, &tab_odstepow[5]);
    podstawa=tab_odstepow[0];
    for(int i=0; i<6; i++) tab_odstepow[i]-=podstawa;
    MPI_Type_struct(6, tab_dl_blokow, tab_odstepow, tab_typow, &MPI_KOM);
    MPI_Type_commit(&MPI_KOM);

    /* If I'm process 0, determine the matrix size and # of iterations */
    /* This relies on the MPI implementation properly flushing output
    that does not end in a newline.  MPI does not require this, though
    high-quality implementations will do this.
    */
/*#if !defined (SGI_MPI) && !defined (IBM_MPI)
    if ( rank == 0 ) {

        printf("Matrix Size : ") ;
        scanf("%d",&N) ;
        printf("Iterations : ") ;
        scanf("%d",&iters) ;
    }
#else
    N=20;
    iters=50;
#endif*/

    /* Broadcast the size and # of iterations to all processes */
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD) ;
    MPI_Bcast(&iters, 1, MPI_INT, 0, MPI_COMM_WORLD) ;

    /* Call the life routine */
    time = life ( N, iters, MPI_COMM_WORLD );

    /* Print the total time taken */
    if (rank == 0)
        printf("[%d] Czas wykonania: %f \n",rank,time);

    //MPE_Close_graphics(&graph);
    MPI_Finalize();
}
