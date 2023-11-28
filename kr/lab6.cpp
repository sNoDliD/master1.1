#include <algorithm.h>
#include <conio.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static int ProcNum = 0;
static int ProcRank = -1;

void DataDistribution(double* pMatrix, double* pProcRows, int RowNum, int Size) {
    int* pSendNum;
    int* pSendInd;
    int RestRows = Size;

    pSendInd = new int[ProcNum];
    pSendNum = new int[ProcNum];

    RowNum = (Size - 2) / ProcNum + 2;
    pSendNum[0] = RowNum * Size;
    pSendInd[0] = 0;
    for (int i = 1; i < ProcNum; i++) {
        RestRows = RestRows - RowNum + 2;
        RowNum = (RestRows - 2) / (ProcNum - i) + 2;
        pSendNum[i] = (RowNum)*Size;
        pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1] - Size;
    }
    MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows, pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    delete[] pSendInd;
    delete[] pSendNum;
}

void ProcessTermination(double* pMatrix, double* pProcRows) {
    if (ProcRank == 0) delete[] pMatrix;
    delete[] pProcRows;
}

void PrintMatrix(double* pMatrix, int RowCount, int ColCount) {
    int i, j;
    for (int i = 0; i < RowCount; i++) {
        for (j = 0; j < ColCount; j++) printf("%7.4f ", pMatrix[i * ColCount + j]);
        printf("\n");
    }
}

double IterationCalculation(double* pProcRows, int Size, int RowNum) {
    int i, j;
    double dm, dmax, temp;
    dmax = 0;
    for (i = 1; i < RowNum - 1; i++)
        for (j = 1; j < Size - 1; j++) {
            temp = pProcRows[Size * i + j];
            pProcRows[Size * i + j] = 0.25 * (pProcRows[Size * i + j + 1] + pProcRows[Size * i + j - 1] +
                                              pProcRows[Size * (i + 1) + j] + pProcRows[Size * (i - 1) + j]);
            dm = fabs(pProcRows[Size * i + j] - temp);
            if (dmax < dm) dmax = dm;
        }
    return dmax;
}

void TestDistribution(double* pMatrix, double* pProcRows, int Size, int RowNum) {
    if (ProcRank == 0) {
        printf("Initial Matrix: \n");
        PrintMatrix(pMatrix, Size, Size);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < ProcNum; i++) {
        if (ProcRank == i) {
            printf("\nProcRank = %d \n", ProcRank);
            PrintMatrix(pProcRows, RowNum, Size);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
void DummyDataInitialization(double* pMatrix, int Size) {
    int i, j;
    double h = 1.0 / (Size - 1);
    for (i = 0; i < Size; i++) {
        for (j = 0; j < Size; j++)
            if ((i == 0) || (i == Size - 1) || (j == 0) || (j == Size - 1))
                pMatrix[i * Size + j] = 100;
            else
                pMatrix[i * Size + j] = 0;
    }
}

void ProcessInitialization(double*& pMatrix, double*& pProcRows, int& Size, int& RowNum, double& Eps) {
    int RestRows;
    if (ProcRank == 0) {
        do {
            printf("\nEnter the grid size: ");
            scanf("%d", &Size);
            if (Size <= 2) {
                printf("\n Size of grid must be greater than 2! \n");
            }
            if (Size < ProcNum) {
                printf(
                    "Size of grid must be greater than"
                    "the number of processes! \n ");
            }
        } while ((Size <= 2) || (Size < ProcNum));
        do {
            printf("\nEnter the required accuracy: ");
            scanf("%lf", &Eps);
            printf("\nChosen accuracy = %lf", Eps);
            if (Eps <= 0) printf("\nAccuracy must be greater than 0!\n");
        } while (Eps <= 0);
    }
    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    RestRows = Size;
    for (i = 0; i < ProcRank; i++) RestRows = RestRows - RestRows / (ProcNum - i);
    RowNum = (RestRows - 2) / (ProcNum - ProcRank) + 2
             pProcRows = new double[RowNum * Size];
    if (ProcRank == 0) {
        pMatrix = new double[Size * Size];
        DummyDataInitialization(pMatrix, Size);
    }
}

void ExchangeData(double* pProcRows, int Size, int RowNum) {
    MPI_Status status;
    int NextProcNum = (ProcRank == ProcNum - 1) ? MPI_PROC_NULL : ProcRank + 1;
    int PrevProcNum = (ProcRank == 0) ? MPI_PROC_NULL : ProcRank - 1;
    MPI_Sendrecv(pProcRows + Size * (RowNum - 2), Size, MPI_DOUBLE, NextProcNum, 4, pProcRows, Size, MPI_DOUBLE,
                 PrevProcNum, 4, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(pProcRows + Size, Size, MPI_DOUBLE, PrevProcNum, 5, pProcRows + (RowNum - 1) * Size, Size, MPI_DOUBLE,
                 NextProcNum, 5, MPI_COMM_WORLD, &status);
}

void ParallelResultCalculation(double* pProcRows, int Size, int RowNum, double Eps, int& Iterations) {
    double ProcDelta, Delta;
    Iterations = 0;
    do {
        Iterations++;
        ExchangeData(pProcRows, Size, RowNum);
        ProcDelta = IterationCalculation(pProcRows, Size, RowNum);
        MPI_Allreduce(&ProcDelta, &Delta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    } while (Delta > Eps);
}

void ResultCollection(double* pMatrix, double* pProcResult, int Size, int RowNum) {
    int* pReceiveNum;
    int* pReceiveInd;
    int RestRows = Size;
    int i;
    pReceiveNum = new int[ProcNum];
    pReceiveInd = new int[ProcNum];
    pReceiveInd[0] = 0;
    RowNum = (Size - 2) / ProcNum + 2;
    pReceiveNum[0] = RowNum * Size;
    for (i = 1; i < ProcNum; i++) {
        RestRows = RestRows - RowNum + 1;
        RowNum = (RestRows - 2) / (ProcNum - i) + 2;
        pReceiveNum[i] = RowNum * Size;
        pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1] - Size;
    }
    MPI_Allgatherv(pProcRows, pReceiveNum[ProcRank], MPI_DOUBLE, pMatrix, pReceiveNum, pReceiveInd, MPI_DOUBLE,
                   MPI_COMM_WORLD);
    delete[] pReceiveNum;
    delete[] pReceiveInd;
}

void SerialResultCalculation(double* pMatrixCopy, int Size, double Eps, int& Iter) {
    int i, j;
    double dm, dmax, temp;
    Iter = 0;
    do {
        dmax = 0;
        for (i = 1; i < Size - 1; i++)
            for (j = 1; j < Size - 1; j++) {
                temp = pMatrixCopy[Size * i + j];
                pMatrixCopy[Size * i + j] = 0.25 * (pMatrixCopy[Size * i + j + 1] + pMatrixCopy[Size * i + j - 1] +
                                                    pMatrixCopy[Size * (i + 1) + j] + pMatrixCopy[Size * (i - 1) + j]);
                dm = fabs(pMatrixCopy[Size * i + j] - temp);

                if (dmax < dm) dmax = dm;
            }
        Iter++;
    } while (dmax > Eps);
}

void CopyData(double* pMatrix, int Size, double* pSerialMatrix) {
    copy(pMatrix, pMatrix + Size, pSerialMatrix);
}

void TestResult(double* pMatrix, double* pSerialMatrix, int Size, double Eps) {
    int equal = 0;
    int Iter;
    if (ProcRank == 0) {
        SerialResultCalculation(pSerialMatrix, Size, Eps, Iter);
        for (int i = 0; i < Size * Size; i++) {
            if (fabs(pSerialMatrix[i] - pMatrix[i]) >= Eps) equal = 1;
            break;
        }
        if (equal == 1)
            printf(
                "The results of serial and parallel algorithms"
                "are NOT identical. Check your code.");
        else
            printf(
                "The results of serial and parallel algorithms"
                "are identical.");
    }
}

void RandowmDataInitialization(double* pMatrix, int Size) {
    int i, j;
    srand(unsigned(clock()));
    for (i = 0; i < Size; i++) {
        for (j = 0; j < Size; j++)
            if ((i == 0) || (i == Size - 1) || (j == 0) || (j == Size - 1))
                pMatrix[i * Size + j] = 100;
            else
                pMatrix[i * Size + j] = rand() / double(1000);
    }
}

void main(int argc, char* argv[]) {
    double* pMatrix;
    double* pProcRows;
    double* pSerialMatrix;
    int Size;
    int RowNum;
    double Eps;
    int Iterations;
    double currDelta, delta;
    setvbuf(stdout, 0, _IONBF, 0);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0) {
        printf("Parallel Gauss - Seidel algorithm \n");
        fflush(stdout);
    }
    ProcessInitialization(pMatrix, pProcRows, Size, RowNum, Eps);
    if (ProcRank == 0) {
        pSerialMatrix = new double[Size * Size];
        CopyData(pMatrix, Size, pSerialMatrix);
    }
    DataDistribution(pMatrix, pProcRows, Size, RowNum);
    ParallelResultCalculation(pProcRows, Size, RowNum, Eps, Iterations);
    ResultCollection(pProcRows, pMatrix, Size, RowNum);
    TestDistribution(pMatrix, pProcRows, Size, RowNum);

    printf("\n Iter %d \n", Iterations);
    printf("\nResult matrix: \n");
    if (ProcRank == 0) {
        PrintMatrix(pMatrix, Size, Size);
    }
    if (ProcRank == 0) delete[] pSerialMatrix;
    ProcessTermination(pMatrix, pProcRows);
    MPI_Finalize();
}
