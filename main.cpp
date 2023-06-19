/*****************************************
Program Title: Matrix Destruction
File Name: Software Math Assignment 2
Author: Matthew Molloy
Completion Date: 4/6/2022
Description: find the potential inverse of 1000 matrices using Gaussian elimination
*****************************************/
#include <iostream>
#include <vector> 
#include <iomanip>

using namespace std;

class MyMatrix
{
private:
    vector<double> elements;    // array containing data
    vector<double*> rows;       // array of pointers to start of each row
    unsigned ncols;             // number of columns

public:
    MyMatrix(unsigned nrows1, unsigned ncols1);
    bool gaussEliminate(MyMatrix og);
    double get(unsigned irow, unsigned jcol);
    void set(unsigned irow, unsigned jcol, double value);
    void interchange(unsigned irow, unsigned jrow);
    void multiply(unsigned irow, double x);
    void add(unsigned irow, unsigned jrow, double x);
    void fill();
    void print();
    void multiplyTest(MyMatrix matrix1, MyMatrix matrix2);
    // checks if diagonal values of reduced row echelon matrix all approx == 1
    bool invertibleTest(MyMatrix reducedMatrix);
};

int main()
{
    int howManyMatricesSir = 100; // will take approx 15seconds to output everything
    int invertible = 0;
    int notInvertible = 0;
    for (int i = 0; i < howManyMatricesSir; i++)
    {
        unsigned nrows = 13;
        unsigned ncols = 13;
        MyMatrix matrix(nrows, ncols);
        matrix.fill();
        cout << "the original matrix" << endl;
        matrix.print();
        if (matrix.gaussEliminate(matrix))
        {
            invertible++;
            continue;
        }
        notInvertible++;
    }
    // note: all tests resulted in 0 nonInvertibles, likely a result of a faulty fill()
    cout << "\n\n\n" << invertible << " Matrices were inverted\n" << notInvertible << " Have no inverse";
};


MyMatrix::MyMatrix(unsigned ncols, unsigned nrows)
{
    elements.resize(ncols * ncols);//resizing elements to nrows*ncols
    rows.resize(ncols);//resize rows to the number of rows
    for (unsigned i = 0; i < ncols; i++)
        rows[i] = &(elements[i * ncols]);
    
    this->ncols = ncols;
}

bool MyMatrix::gaussEliminate(MyMatrix A) {

    MyMatrix copy(ncols, ncols);    // copy of original matrix for multiplyTest()
    MyMatrix inverse(ncols, ncols); // another copy which will experience some serious elimination until reduced row echelon
    MyMatrix identity(ncols, ncols);// 13 by 13 identity matrix

    // filling in the copy, inverse, and identity matrices
    for (int row = 0; row < A.ncols; row++)
        for (int column = 0; column < A.ncols; column++)
        {
            copy.set(row, column, (A.get(row, column)));
            inverse.set(row, column, (A.get(row, column)));
            identity.set(row, column, 0.0);
            if (row == column)
                identity.set(row, column, 1);
        }

    // partial pivoting; ensure diagonal value is greater than all values below it in the column
    for (int row = 0; row < A.ncols; row++)
        for (int column = row + 1; column < A.ncols; column++)
            if (inverse.get(row, row) < inverse.get(column, row))
            {
                inverse.interchange(row, column);
                identity.interchange(row, column);
            }

    // forward elimination; multiplies pivot by its reciprocal, subtract(add) rows below it by multiples of of pivot
    for (int column = 0; column < A.ncols; column++)
    {
        double pivot = inverse.get(column, column);
        inverse.multiply(column, (1 / pivot));
        identity.multiply(column, (1 / pivot));

        for (int row = column + 1; row < A.ncols; row++)
        {
            double multiplyCoefficient = inverse.get(row, column);
            if (multiplyCoefficient != 0.0) // prevents the classic x / 0
            {
                inverse.add(row, column, -multiplyCoefficient);
                identity.add(row, column, -multiplyCoefficient);
            }
        }
    }

    // backward elimination, make all values above pivot == 0 using nice legal methods
    for (int column = A.ncols - 1; column > 0; column--)
    {
        for (int row = column - 1; row >= 0; row--)
        {
            double multiplyCoefficient = inverse.get(row, column);
            if (multiplyCoefficient != 0.0)
            {
                inverse.add(row, column, -multiplyCoefficient);
                identity.add(row, column, -multiplyCoefficient);
            }
        }
    }

    cout << "\nMatrix Reduced Form!\n";
    identity.print();

    // AA^-1
    multiplyTest(copy, identity);
    if (invertibleTest(inverse))
        return true;
    else
        return false;
}

// found a better way to access elements than lab9; doesn't complicate things with pointer math
double MyMatrix::get(unsigned irow, unsigned jcol)
{
    return this->rows[irow][jcol];
}

void MyMatrix::set(unsigned irow, unsigned jcol, double x)
{
    rows[irow][jcol] = x;
}

void MyMatrix::interchange(unsigned irow, unsigned jrow)
{
    double* ptrPlaceholder = rows[irow];
    rows[irow] = rows[jrow];
    rows[jrow] = ptrPlaceholder;
}

void MyMatrix::multiply(unsigned irow, double x)
{
    for (unsigned j = 0; j < ncols; j++)
    {
        rows[irow][j] = rows[irow][j] * x;
    }
}

void MyMatrix::add(unsigned irow, unsigned jrow, double x)
{
    for (unsigned i = 0; i < ncols; i++)
    {
        double value = get(irow, i) + (x * get(jrow, i));
        set(irow, i, value);
    }
}

void::MyMatrix::fill()
{
    double decimalCreator = 100;
    for (int i = 0; i < ncols; i++)
    {
        for (int j = 0; j < ncols; j++) {

            double random = (double)rand() / decimalCreator;
            if (i % 5 == 0 or j % 7 == 0 and j != 0) // creates some negative elements
                random *= -1;
            rows[i][j] = random;
            elements[0] = random;
        }
    }
}

void::MyMatrix::print()
{
    for (int i = 0; i < ncols; i++)
    {
        for (int j = 0; j < ncols; j++)
            cout << setw(10) << setprecision(3) << fixed << get(i, j);

        cout << endl;
    }
}

void::MyMatrix::multiplyTest(MyMatrix matrix1, MyMatrix matrix2) 
{
    MyMatrix id(ncols, ncols);

    for (int i = 0; i < ncols; i++)
        for (int j = 0; j < ncols; j++)
        {
            double sum = 0;
            for (int k = 0; k < ncols; k++)
            {
                sum += matrix1.get(i, k) * matrix2.get(k, j); 
                id.set(i, j, sum);                     
            }
        }

    cout << "\nResult of Inverse Multiplied by Identity\n";
    id.print();
}

bool::MyMatrix::invertibleTest(MyMatrix reducedMatrix)
{
    double determinator = 1;
    double acceptedError = 0.0001;
    double expectedResult = 1;
    for (int i = 0; i < ncols; i++)
        for (int j = 0; j < ncols; j++)
            if (i == j)
                determinator *= reducedMatrix.get(i, j);

    if (determinator >= determinator - acceptedError and determinator <= determinator + acceptedError)
        return true;
    return false;
}