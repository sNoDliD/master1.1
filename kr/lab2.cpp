#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

int ProcNum = 0;
int ProcRank = 0;
int GridSize;


int GridCoords[2];
MPI_Comm GridComm;
MPI_Comm ColComm;
MPI_Comm RowComm;


void DummyDataInitialization (double* pAMatrix, double* pBMatrix,int Size) {
    int i, j; 
    for (i=0; i<Size; i++)
        for (j=0; j<Size; j++) {
            pAMatrix[i*Size+j] = 1;
            pBMatrix[i*Size+j] = 1;
        }
}
void RandomDataInitialization (double* pAMatrix, double* pBMatrix, int Size) {
    int i, j;
    srand(unsigned(clock()));
    for (i=0; i<Size; i++)
        for (j=0; j<Size; j++) {
            pAMatrix[i*Size+j] = rand()/double(1000);
            pBMatrix[i*Size+j] = rand()/double(1000);
        }
}

void PrintMatrix (double* pMatrix, int RowCount, int ColCount) {
    int i, j;
    for (i=0; i<RowCount; i++) {
        for (j=0; j<ColCount; j++)
            printf("%7.4f ", pMatrix[i*ColCount+j]);
            printf("\n");
        }
}

void SerialResultCalculation(double* pAMatrix, double* pBMatrix, double* pCMatrix, int Size) {
    int i, j, k;
    for (i=0; i<Size; i++) {
        for (j=0; j<Size; j++)
            for (k=0; k<Size; k++)
                pCMatrix[i*Size+j] += pAMatrix[i*Size+k]*pBMatrix[k*Size+j];
    }
}

void BlockMultiplication(double* pAblock, double* pBblock, double* pCblock, int Size) {
    SerialResultCalculation(pAblock, pBblock, pCblock, Size);
}

void CreateGridCommunicators() {
    int DimSize[2];
    int Periodic[2];
    int Subdims[2];
    DimSize[0] = GridSize;

    DimSize[1] = GridSize;
    Periodic[0] = 0;
    Periodic[1] = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 1, &GridComm);
    MPI_Cart_coords(GridComm, ProcRank, 2, GridCoords);
    Subdims[0] = 0;
    Subdims[1] = 1;
    MPI_Cart_sub(GridComm, Subdims, &RowComm);

    Subdims[0] = 1;
    Subdims[1] = 0;
    MPI_Cart_sub(GridComm, Subdims, &ColComm);
}

void ProcessInitialization (double* &pAMatrix, double* &pBMatrix, double* &pCMatrix,
    double* &pAblock, double* &pBblock, double* &pCblock,
    double* &pTemporaryAblock, int &Size, int &BlockSize ) {
    if (ProcRank == 0) {
        do {
            printf("\nEnter the size of matrices: ");
            scanf("%d", &Size);
            if (Size%GridSize != 0) {
                printf ("Size of matrices must be divisible by the grid size!\n");
            }
        }
        while (Size%GridSize != 0);
    }
    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    BlockSize = Size/GridSize;
    pAblock = new double [BlockSize*BlockSize];
    pBblock = new double [BlockSize*BlockSize];
    pCblock = new double [BlockSize*BlockSize];
    pTemporaryAblock = new double [BlockSize*BlockSize];
    for (int i=0; i<BlockSize*BlockSize; i++) {
        pCblock[i] = 0;
    }
    if (ProcRank == 0) {
        pAMatrix = new double [Size*Size];
        pBMatrix = new double [Size*Size];
        pCMatrix = new double [Size*Size];
        DummyDataInitialization(pAMatrix, pBMatrix, Size);
    }
}

void CheckerboardMatrixScatter(double* pMatrix, double* pMatrixBlock, int Size, int BlockSize) {
    double * MatrixRow = new double [BlockSize*Size];
    if (GridCoords[1] == 0) {
        MPI_Scatter(pMatrix, BlockSize*Size, MPI_DOUBLE, MatrixRow,
        BlockSize*Size, MPI_DOUBLE, 0, ColComm);
    }
    for (int i=0; i<BlockSize; i++) {
        MPI_Scatter(&MatrixRow[i*Size], BlockSize, MPI_DOUBLE,
            &(pMatrixBlock[i*BlockSize]), BlockSize, MPI_DOUBLE, 0, RowComm);
    }
    delete [] MatrixRow;
}

void DataDistribution(double* pAMatrix, double* pBMatrix, double* pMatrixAblock, double* pBblock, int Size, int BlockSize) {
    CheckerboardMatrixScatter(pAMatrix, pMatrixAblock, Size, BlockSize);
    CheckerboardMatrixScatter(pBMatrix, pBblock, Size, BlockSize);
}

void ResultCollection (double* pCMatrix, double* pCblock, int Size, int BlockSize) {
    double * pResultRow = new double [Size*BlockSize];
    for (int i=0; i<BlockSize; i++) {
        MPI_Gather( &pCblock[i*BlockSize], BlockSize, MPI_DOUBLE,
        &pResultRow[i*Size], BlockSize, MPI_DOUBLE, 0, RowComm);
    }
    if (GridCoords[1] == 0) {
        MPI_Gather(pResultRow, BlockSize*Size, MPI_DOUBLE, pCMatrix,
        BlockSize*Size, MPI_DOUBLE, 0, ColComm);
    }
    delete [] pResultRow;
}
void ABlockCommunication(int iter, double *pAblock, double* pMatrixAblock, int BlockSize) {
    int Pivot = (GridCoords[0] + iter) % GridSize;
    if (GridCoords[1] == Pivot) {
        for (int i=0; i<BlockSize*BlockSize; i++)
            pAblock[i] = pMatrixAblock[i];
    }

    MPI_Bcast(pAblock, BlockSize*BlockSize, MPI_DOUBLE, Pivot, RowComm);
}

void BblockCommunication (double *pBblock, int BlockSize) {
    MPI_Status Status;
    int NextProc = GridCoords[0] + 1;
    if ( GridCoords[0] == GridSize-1 ) NextProc = 0;
    int PrevProc = GridCoords[0] - 1;
    if ( GridCoords[0] == 0 ) PrevProc = GridSize-1;
    MPI_Sendrecv_replace( pBblock, BlockSize*BlockSize, MPI_DOUBLE,
    NextProc, 0, PrevProc, 0, ColComm, &Status);
}

void ParallelResultCalculation(double* pAblock, double* pMatrixAblock, double* pBblock, double* pCblock, int BlockSize) {
    for (int iter = 0; iter < GridSize; iter ++) {
        ABlockCommunication(iter, pAblock, pMatrixAblock, BlockSize);
        BlockMultiplication(pAblock, pBblock, pCblock, BlockSize);
        BblockCommunication(pBblock, BlockSize);
    }
}

void TestBlocks (double* pBlock, int BlockSize, char str[]) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (ProcRank == 0) {
        printf("%s \n", str);
    }
    for (int i=0; i<ProcNum; i++) {
        if (ProcRank == i) {
        printf ("ProcRank = %d \n", ProcRank);
            PrintMatrix(pBlock, BlockSize, BlockSize);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void TestResult(double* pAMatrix, double* pBMatrix, double* pCMatrix, int Size) {
    double* pSerialResult;
    double Accuracy = 1.e-6;
    int equal = 0;
    int i; // Loop variable
    if (ProcRank == 0) {
    pSerialResult = new double [Size*Size];
    for (i=0; i<Size*Size; i++) {
        pSerialResult[i] = 0;
    }
    BlockMultiplication(pAMatrix, pBMatrix, pSerialResult, Size);
    for (i=0; i<Size*Size; i++) {
        if (fabs(pSerialResult[i]-pCMatrix[i]) >= Accuracy)
            equal = 1;
    }
    if (equal == 1)
        printf("The results of serial and parallel algorithms are NOT"
            "identical. Check your code.");
    else
        printf("The results of serial and parallel algorithms are "
            "identical. ");
    }
}

void ProcessTermination (double* pAMatrix, double* pBMatrix, double* pCMatrix, double* pAblock,
    double* pBblock, double* pCblock, double* pMatrixAblock) {
    if (ProcRank == 0) {
        delete [] pAMatrix;
        delete [] pBMatrix;
        delete [] pCMatrix;
    }
    delete [] pAblock;
    delete [] pBblock;
    delete [] pCblock;
    delete [] pMatrixAblock;
}

void main(int argc, char* argv[]) {
    double* pAMatrix;
    double* pBMatrix;
    double* pCMatrix;
    int Size;
    int BlockSize;
    double *pAblock;
    double *pBblock;
    double *pCblock;
    double *pMatrixAblock;
    double Start, Finish, Duration;
    setvbuf(stdout, 0, _IONBF, 0);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    GridSize = sqrt((double)ProcNum);
    if (ProcNum != GridSize*GridSize) {
        if (ProcRank == 0) {
            printf ("Number of processes must be a perfect square \n");
        }
    }
    else {
        if (ProcRank == 0)
            printf("Parallel matrix multiplication program\n");
        CreateGridCommunicators();
        ProcessInitialization(pAMatrix, pBMatrix, pCMatrix, pAblock, pBblock,
        pCblock, pMatrixAblock, Size, BlockSize );
        DataDistribution(pAMatrix, pBMatrix, pMatrixAblock, pBblock, Size, BlockSize);
        ParallelResultCalculation(pAblock, pMatrixAblock, pBblock, pCblock, BlockSize);
        ResultCollection(pCMatrix, pCblock, Size, BlockSize);
        TestResult(pAMatrix, pBMatrix, pCMatrix, Size);
        ProcessTermination (pAMatrix, pBMatrix, pCMatrix, pAblock, pBblock, pCblock, pMatrixAblock);
    }
    MPI_Finalize();
}