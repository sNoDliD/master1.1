#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>

using namespace std;

int ProcRank;
int ProcNum;
const double InfinitiesPercent = 50.0;
const double RandomDataMultiplier = 10;

int Min(int A, int B) {
    int Result = (A < B) ? A : B;
    if ((A < 0) && (B >= 0)) Result = B;
    if ((B < 0) && (A >= 0)) Result = A;
    if ((A < 0) && (B < 0)) Result = -1;
    return Result;
}

void ProcessInitialization(int *&pMatrix, int *&pProcRows, int &Size, int &RowNum) {
    setvbuf(stdout, 0, _IONBF, 0);
    if (ProcRank == 0) {
        do {
            printf("Enter the number of vertices: ");
            scanf("%d", &Size);
            if (Size < ProcNum)
                printf(
                    "The number of vertices should be greater then"
                    "the number of processes\n");
        } while (Size < ProcNum);
        printf("Using the graph with %d vertices\n", Size);
    }
    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int RestRows = Size;
    for (int i = 0; i < ProcRank; i++) RestRows = RestRows - RestRows / (ProcNum - i);
    RowNum = RestRows / (ProcNum - ProcRank);
    pProcRows = new int[Size * RowNum];
    if (ProcRank == 0) {
        pMatrix = new int[Size * Size];
        DummyDataInitialization(pMatrix, Size);
    }
}

void ProcessTermination(int *pMatrix, int *pProcRows) {
    if (ProcRank == 0) delete[] pMatrix;
    delete[] pProcRows;
}

void DummyDataInitialization(int *pMatrix, int Size) {
    for (int i = 0; i < Size; i++)
        for (int j = i; j < Size; j++) {
            if (i == j)
                pMatrix[i * Size + j] = 0;
            else if (i == 0)
                pMatrix[i * Size + j] = j;
            else
                pMatrix[i * Size + j] = -1;
            pMatrix[j * Size + i] = pMatrix[i * Size + j];
        }
}

void RandomDataInitialization(int *pMatrix, int Size) {
    srand((unsigned)time(0));
    for (int i = 0; i < Size; i++)
        for (int j = 0; j < Size; j++)
            if (i != j) {
                if ((rand() % 100) < InfinitiesPercent)
                    pMatrix[i * Size + j] = -1;
                else
                    pMatrix[i * Size + j] = rand() + 1;
            } else
                pMatrix[i * Size + j] = 0;
}

void DataDistribution(int *pMatrix, int *pProcRows, int Size, int RowNum) {
    int *pSendNum;
    int *pSendInd;
    int RestRows = Size;
    pSendInd = new int[ProcNum];
    pSendNum = new int[ProcNum];
    RowNum = Size / ProcNum;
    pSendNum[0] = RowNum * Size;
    pSendInd[0] = 0;
    for (int i = 1; i < ProcNum; i++) {
        RestRows -= RowNum;
        RowNum = RestRows / (ProcNum - i);
        pSendNum[i] = RowNum * Size;
        pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
    }
    MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_INT, pProcRows, pSendNum[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);
    delete[] pSendNum;
    delete[] pSendInd;
}

void ResultCollection(int *pMatrix, int *pProcRows, int Size, int RowNum) {
    int *pReceiveNum;
    int *pReceiveInd;
    int RestRows = Size;

    pReceiveNum = new int[ProcNum];
    pReceiveInd = new int[ProcNum];
    RowNum = Size / ProcNum;
    pReceiveInd[0] = 0;
    pReceiveNum[0] = RowNum * Size;
    for (int i = 1; i < ProcNum; i++) {
        RestRows -= RowNum;
        RowNum = RestRows / (ProcNum - i);
        pReceiveNum[i] = RowNum * Size;
        pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
    }
    MPI_Gatherv(pProcRows, pReceiveNum[ProcRank], MPI_INT, pMatrix, pReceiveNum, pReceiveInd, MPI_INT, 0,
                MPI_COMM_WORLD);
    delete[] pReceiveNum;
    delete[] pReceiveInd;
}

void ParallelFloyd(int *pProcRows, int Size, int RowNum) {
    int *pRow = new int[Size];
    int t1, t2;
    for (int k = 0; k < Size; k++) {
        RowDistribution(pProcRows, Size, RowNum, k, pRow);
        for (int i = 0; i < RowNum; i++)
            for (int j = 0; j < Size; j++)
                if ((pProcRows[i * Size + k] != -1) && (pRow[j] != -1)) {
                    t1 = pProcRows[i * Size + j];
                    t2 = pProcRows[i * Size + k] + pRow[j];
                    pProcRows[i * Size + j] = Min(t1, t2);
                }
    }
    delete[] pRow;
}

void RowDistribution(int *pProcRows, int Size, int RowNum, int k, int *pRow) {
    int ProcRowRank;
    int ProcRowNum;
    int RestRows = Size;
    int Ind = 0;
    int Num = Size / ProcNum;
    for (ProcRowRank = 1; ProcRowRank < ProcNum + 1; ProcRowRank++) {
        if (k < Ind + Num) break;
        RestRows -= Num;
        Ind += Num;
        Num = RestRows / (ProcNum - ProcRowRank);
    }
    ProcRowRank = ProcRowRank - 1;
    ProcRowNum = k - Ind;
    if (ProcRowRank == ProcRank) copy(&pProcRows[ProcRowNum * Size], &pProcRows[(ProcRowNum + 1) * Size], pRow);
    MPI_Bcast(pRow, Size, MPI_INT, ProcRowRank, MPI_COMM_WORLD);
}

void ParallelPrintMatrix(int *pProcRows, int Size, int RowNum) {
    for (int i = 0; i < ProcNum; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (ProcRank == i) {
            printf("ProcRank = %d\n", ProcRank);
            fflush(stdout);
            printf("Proc rows:\n");
            fflush(stdout);
            PrintMatrix(pProcRows, RowNum, Size);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void TestDistribution(int *pMatrix, int *pProcRows, int Size, int RowNum) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (ProcRank == 0) {
        printf("Initial adjacency matrix:\n");
        PrintMatrix(pMatrix, Size, Size);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ParallelPrintMatrix(pProcRows, Size, RowNum);
}

void TestResult(int *pMatrix, int *pSerialMatrix, int Size) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (ProcRank == 0) {
        SerialFloyd(pSerialMatrix, Size);
        if (!CompareMatrices(pMatrix, pSerialMatrix, Size)) {
            printf(
                "Results of serial and parallel algorithms are "
                "NOT identical. Check your code\n");
        } else {
            printf(
                "Results of serial and parallel algorithms are "
                "identical\n");
        }
    }
}

int main(int argc, char *argv[]) {
    int *pMatrix;
    int Size;
    int *pProcRows;
    int RowNum;
    double start, finish;
    double duration = 0.0;
    int *pSerialMatrix = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    if (ProcRank == 0) printf("Parallel Floyd algorithm\n");

    ProcessInitialization(pMatrix, pProcRows, Size, RowNum);
    if (ProcRank == 0) {
        pSerialMatrix = new int[Size * Size];
        CopyMatrix(pMatrix, Size, pSerialMatrix);
    }
    start = MPI_Wtime();
    DataDistribution(pMatrix, pProcRows, Size, RowNum);
    ParallelFloyd(pProcRows, Size, RowNum);
    ResultCollection(pMatrix, pProcRows, Size, RowNum);
    finish = MPI_Wtime();

    duration = finish - start;
    if (ProcRank == 0) printf("Time of execution: %f\n", duration);
    if (ProcRank == 0) delete[] pSerialMatrix;
    ProcessTermination(pMatrix, pProcRows);
    MPI_Finalize();
    return 0;
}
