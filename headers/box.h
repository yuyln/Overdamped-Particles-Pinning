#ifndef __BOX
#define __BOX
class Box
{
public:
    double lx, ly, x, y;
    size_t nIn, npar;
    size_t *indexes;
    Box(): indexes(nullptr) {}

    Box(const Box &o) = delete;
/*    {
        if (indexes)
            delete[] indexes;
        
        lx = o.lx; ly = o.ly;
        x = o.x; y = o.y;
        npar = o.npar;
        indexes = new size_t[o.npar];
        memcpy((void*)indexes, (void*)o.indexes, sizeof(size_t) * npar);
    }*/

    void operator=(const Box &o)
    {
        if (indexes)
            delete[] indexes;
        
        lx = o.lx; ly = o.ly;
        x = o.x; y = o.y;
        npar = o.npar;
        indexes = new size_t[o.npar];
        memcpy((void*)indexes, (void*)o.indexes, sizeof(size_t) * npar);
    }

    ~Box()
    {
        if (indexes)
            delete[] indexes;
    }

    Box(double Lx, double Ly, double X, double Y, size_t nPar): lx(Lx), ly(Ly), x(X), y(Y), npar(nPar)
    {
        indexes = new size_t[nPar];
    }

    template <typename T>
    void AttBox(T *pos, size_t n) noexcept
    {
        nIn = 0;
        for (size_t i = 0; i < n; ++i)
        {
            bool iX = (pos[i].x >= x) && (pos[i].x < (x + lx));
            bool iY = (pos[i].y >= y) && (pos[i].y < (y + ly));
            if (iX && iY)
            {
                indexes[nIn] = i;
                nIn++;
            }
        }
    }

    const size_t & GetIn() const noexcept { return nIn; }
    const size_t & GetIndex(const size_t &i) const noexcept { return indexes[i]; }
    const double & GetLx() const noexcept { return lx; }
    const double & GetLy() const noexcept { return ly; }
    const double & GetX() const noexcept { return x; }
    const double & GetY() const noexcept { return y; }
    void SetIn(size_t i) noexcept { nIn = i; }
    void SetIndex(size_t index, size_t value) { indexes[index] = value; }
};

#endif