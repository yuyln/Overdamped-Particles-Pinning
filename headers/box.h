#ifndef __BOX
#define __BOX
class Box
{
    double lx, ly, x, y;
    size_t nIn;
    size_t *indexes;
public:
    Box(){}

    Box(double Lx, double Ly, double X, double Y, size_t nPar): lx(Lx), ly(Ly), x(X), y(Y)
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
};

#endif