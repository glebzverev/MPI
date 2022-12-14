#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <sys/time.h>

#define  Max(a,b) ((a)>(b)?(a):(b))
#define  N   (256)
float maxeps = 0.1e-7;
int itmax = 1000;
int i,j,k;
float w = 0.5;
float eps;
float A [N][N][N];
float B [N][N][N];
int I[N];

void printA(){
    printf("A---------------------A\n");
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            for(k=0; k<N; k++){
                printf("%f %d, %d, %d\n", A[i][j][k],i,j,k);
            }
        }
    }
    printf("A---------------------A\n");
}

MPI_Request request;

void init()
{ 
	for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            for(k=0; k<N; k++)
            {
                if(i==0 || i==N-1 || j==0 || j==N-1 || k==0 || k==N-1)     
                    A[i][j][k] = 0.;
                else 
                    A[i][j][k] = ( 4. + i + j + k) ;
            }
        }
    }
} 
void copy(int rank, int numprocess, int line);
void relax(int rank, int numprocs, int line);
void verify(int rank, int numprocs, int line); 

MPI_Request worker_send_requests[N-2];
MPI_Request root_rec_requests[N-2];
MPI_Status root_statuses[N-2];

int main(int an, char **as){

    int i, j, k, rank, numprocs;
    MPI_Init(NULL,NULL); //MPI Initialize
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); // Получить текущий номер процесса
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);// Получить количество процессов
    int line = N / (numprocs-1);
    int it;
    struct timeval start_time, stop_time, elapsed_time;

    // MPI_Barrier(MPI_COMM_WORLD);  
    
    MPI_Request root_send_requests[numprocs-1];
    MPI_Request worker_rec_requests[numprocs-1];
    MPI_Status worker_statuses[numprocs-1];

    if (rank == 0)
    {
        // Иницализируем матрицу в ROOT процессе
        init();
        // Устанавливаем время начала релаксации
        gettimeofday(&start_time, NULL);
        int worker_value = 0;
        MPI_Status worker_status;

        for(it=1; it < itmax; it++)
            {
            for(int INDEXER=1; INDEXER < numprocs; INDEXER++){
                // printf("MPI_Send %d %d\n", INDEXER, it);
                MPI_Isend( &A, N*N*N, MPI_INT, INDEXER, 0, MPI_COMM_WORLD, &root_send_requests[INDEXER-1]);
            }
                // printf("INDEXER %d %d\n", INDEXER, it);
            for(int INDEXER=1; INDEXER < numprocs; INDEXER++){

                int start =  (INDEXER == 1) ? 1:0;
                int finish = (INDEXER == numprocs-1)? 1:0;
                for( int i = start; i < line-finish; i++){
                    int base =(INDEXER-1)*line + i; 
                    float res[N][N];
                    // printf("base %d %d \n", base, INDEXER);
                    // MPI_Recv(&A[base], N*N, MPI_FLOAT, INDEXER, i, MPI_COMM_WORLD, &worker_status );
                    MPI_Irecv(&A[base], N*N, MPI_FLOAT, INDEXER, i, MPI_COMM_WORLD, &root_rec_requests[base-1] );
                    // printf("MPI_Recv %d INDEXER %d \n", base, INDEXER);
                } 
            }
            int flag[N-2];
            MPI_Testall(N - 2, root_rec_requests, flag, root_statuses);
            // printf("Tested\n");
        }

        printf("Relaxation SUCCESS! Amount of iterations: %d\n", it-1);
        // Отправлем процессам WORKER информацию, что матрица можно делать
        // итую итерацию релаксации
        worker_value = 1;
        
        // Релаксация закончена. Считаем затраченное время
        gettimeofday(&stop_time, NULL);
        timersub(&stop_time, &start_time, &elapsed_time);
        printf("ammount of iterations, %u\n", it-1);
        printf("eps: %f\n", eps);
        printf("time to work %f seconds\n", elapsed_time.tv_sec +
        elapsed_time.tv_usec/1000000.0);
        printf("Количество процессов %d \n",numprocs);
        int a = 10;
        MPI_Send(&a,1,MPI_INT,0,10,MPI_COMM_WORLD);
    }
    else{
        int root_value = 0;
        MPI_Status root_status;
        for(it=1; it < itmax; it++)
        {
            MPI_Irecv( &B, N*N*N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &worker_rec_requests[rank-1]);
            int flag;
            MPI_Wait( &worker_rec_requests[rank-1],  &worker_statuses[rank-1]);
            relax(rank, numprocs, line);
        }

    }
    MPI_Status status;
    int a;
    MPI_Recv(&a,1,MPI_INT,0,10,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    verify(rank, numprocs, line);
    // printf("Verifiyed");
    if (rank == 0){
    MPI_Abort(MPI_COMM_WORLD, 0);

    }
    return 0;
}

void relax(int rank, int numprocess, int line)
{
	int start, finish = (0,0);
    if (rank == 1)
        start = 1; 
    if(rank == numprocess-1)
        finish = 1;
    int base;
    for( int i = start; i < line-finish; i++){
        base =(rank-1)*line + i; 
        for(j = 1; j < N-1; j++){
            for(k=1+(i+j)%2; k < N-1; k+=2){ 
                float b;
                b = w*((B[base-1][j][k]+B[base+1][j][k]+B[base][j-1][k]
                +B[base][j+1][k] +B[base][j][k-1]+B[base][j][k+1] )/6. - B[(rank-1)*line+i][j][k]);
                eps =  Max(fabs(b),eps);
                B[base][j][k] = B[base][j][k] + b;
            }
        }
        MPI_Isend(&B[base], N*N, MPI_FLOAT,0, i,MPI_COMM_WORLD, &worker_send_requests[base-1]);
    }
}

void verify(int rank, int numprocs, int line){
    // printf("Verify\n");
    float a=0.;
    float s=0.;
    if (rank == 0)
    {
        for(i=0; i<=N-1; i++){
        for(j=0; j<=N-1; j++){
        for(k=0; k<=N-1; k++)
        {
        s=s+A[line*rank + i][j][k]*(line*rank+i+1)*(j+1)*(k+1)/(N*N*N);
        }
        }
        }
        printf("S = %f\n",s);
    }
    else{
        printf("RANK: %d", rank);
    }


    // {
    //     for (i=1;i<numprocs;i++)
    //     {
    //         MPI_Send(A,N*N*N,MPI_INT,i,0,MPI_COMM_WORLD);
    //     }
    //     for(i=0; i<=line-1; i++)
    //     for(j=0; j<=N-1; j++)
    //     for(k=0; k<=N-1; k++)
    //     {
    //     s=s+A[line*rank + i][j][k]*(line*rank+i+1)*(j+1)*(k+1)/(N*N*N);
    //     }

    //     for (i=1;i<numprocs;i++){
    //         MPI_Recv(&a,1,MPI_FLOAT,i,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //         s+=a;
    //     }
    //     printf("S = %f\n",s);
    // } else {
    //     MPI_Recv(A, N*N*N, MPI_FLOAT, 0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     for(i=0; i<=line-1; i++)
    //     for(j=0; j<=N-1; j++)
    //     for(k=0; k<=N-1; k++)
    //     {
    //        a=a+A[rank*line + i][j][k]*(rank*line + i+1)*(j+1)*(k+1)/(N*N*N);
    //     }
    //     MPI_Send(&a,1,MPI_FLOAT,0,2,MPI_COMM_WORLD);
    // }
    // printf("END Verify\n");

}


