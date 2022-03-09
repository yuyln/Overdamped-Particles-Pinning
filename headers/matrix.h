#ifndef __MATRIX
#define __MATRIX
#include <cstring>

template <typename T>
class Matrix
{
    T *m;
public:
    size_t nRows, nCols;
    Matrix(){};
    Matrix(size_t nRows_, size_t nCols_): nRows(nRows_), nCols(nCols_)
    {
        m = new T[nRows * nCols];
    }

    T & operator()(size_t i, size_t j) const noexcept { return m[i * nCols + j]; }
    T & operator()(size_t i) const noexcept { return m[i]; }
    void SetData(const T &data, size_t i, size_t j) noexcept { m[i * nCols + j] = data; }
    void SetData(const T** data) noexcept
    {
        for (size_t i = 0; i < nRows; ++i)
        {
            for (size_t j = 0; j < nCols; ++j)
            {
                m[i * nCols + j] = data[i][j];
            }
        }
    }
    void SetData(const T* data) noexcept
    {
        memcpy(m, data, sizeof(T) * nRows * nCols);
    }
};
#endif