#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv) {
    
    // Инициализируем MPI окружение
    MPI_Init(NULL, NULL);

    // Берём количество процессов (нитей)
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Берём rank выполняемого процесса
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Берём имя процесса
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // MPI окружение проинициализировано
    
    int value;
    int INDEXER;
    MPI_Status status;
    
    // Максимум между 25 и иаксимаьным количеством нитей 
    int MAX = (world_size <= 25) ? world_size : 25;
    if (world_rank == 0){
        // ROOT PROCESS
        printf("Hello world from ROOT processor %s, rank %d out of %d processors\n",
            processor_name, world_rank, world_size);
        value = 1;
        INDEXER = 1;

        while(INDEXER < MAX){
            // Отправляем сигнал процессу с индексом INDEXER на прохождение критической секции.
            MPI_Send( &value, 1, MPI_INT, INDEXER, 0, MPI_COMM_WORLD);
            // Ждём прохождения критической секции от процесса с индексом INDEXER.
            MPI_Recv( &value, 1, MPI_INT, INDEXER, 0, MPI_COMM_WORLD, &status );
            INDEXER++;
        }
    }
    else{
        printf("Hello world from WORKER processor %s, rank %d out of %d processors\n",
            processor_name, world_rank, world_size);
        MPI_Recv( &value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status );
        printf("Critical section for rank %d\n", world_rank);
        // Открываем файл на чтение
        FILE *fp = fopen("critical.txt", "rb+");
        if (fp) { 
            // Если существует, ошибка!
            printf("Exist\n");
            perror("Critical section racing!") ;
            return 1;
        }
        else{
            // Если не существует, создадим файл, sleep(0;10), удалим файл.
            srand(time(NULL));
            int r = rand() % 10;
            printf("Sleep: %d\n", r);
            FILE *fp = fopen("critical.txt", "wb+");
            sleep(r);
            remove("critical.txt");
        }
        // Конец критической секции

        // Шлём ROOT процессу сигнал об окончании работы с критической секцией
        MPI_Send( &value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);  
    MPI_Finalize();
}
