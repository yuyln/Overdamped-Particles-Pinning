#ifndef __LINE
#define __LINE
#include <point.h>
class LineSegment
{
public:
    double U0, R, R2, F0;
    double Vx, Vy;
    Point p0, p1;

    LineSegment() {}
    LineSegment(double x0, double y0, double x1, double y1, double R_, double U0_): 
    U0(U0_), R(R_), R2(R_ * R_), F0(2.0 * U0_ / (R_ * R_)), p0(x0, y0), p1(x1, y1)
    {
        Vx = p1.X() - p0.X();
        Vy = p1.Y() - p0.Y();
    }

    LineSegment(Point p0_, Point p1_, double R_, double U0_): 
    U0(U0_), R(R_), R2(R_ * R_), F0(2.0 * U0_ / (R_ * R_)), p0(p0_.X(), p0_.Y()), p1(p1_.X(), p1_.Y())
    {
        Vx = p1.X() - p0.X();
        Vy = p1.Y() - p0.Y();
    }

    const double &dirX() const noexcept { return Vx; }
    const double &dirY() const noexcept { return Vy; }
    void operator=(const LineSegment &o)
    {
        // Vx = o.Vx;
        // Vy = o.Vy;
        // p0.x = o.p0.X();
        // p0.y = o.p0.Y();
        // p1.x = o.p1.X();
        // p1.y = o.p1.Y();
        // U0 = o.U0;
        // R = o.R;
        // R2 = o.R2;
        // F0 = o.F0;
        memcpy((void*)this, (void*)&o, sizeof(LineSegment));
    }
    Point operator()(const double &t) noexcept
    {
        if (t >= 1.0)
        {
            return p1;
        }
        else if (t <= 0.0)
        {
            return p0;
        }
        else
        {
            return Point(p0.X() + Vx * t, p0.Y() + Vy * t);
        }
    }

    Vector DistanceVector(const double &x, const double &y)
    {
        double t = Vy * p0.Y() - Vy * y + Vx * p0.X() - Vx * x;
        t /= (Vy * Vy + Vx * Vx);
        Point ponL = (*this)(-t);
        return Vector(x - ponL.X(), y - ponL.Y());
    }

    double Distance(const double &x, const double &y)
    {
        Vector v = DistanceVector(x, y);
        return sqrt(v.X() * v.X() + v.Y() * v.Y());
    }

    static double Potential(const LineSegment *l, const double &x, const double &y, const Table &table) noexcept
    {
        Vector d = ((LineSegment*)l)->DistanceVector(x, y);
        double dx = d.X();
        double dy = d.Y();
        double d2 = dx * dx + dy * dy;
        return l->U0 * table(d2 / l->R2);
    }

    static void Force(const LineSegment *interact, const double &x, const double &y, const Table &table,
                       double *fx, double *fy) noexcept
    {
        double dx, dy;
        Vector dd = ((LineSegment*)interact)->DistanceVector(x, y);
        dx = dd.X();
        dy = dd.Y();
        if (dx > table.getMaxRange() || dy > table.getMaxRange())
        {
            *fx = 0.0;
            *fy = 0.0;
            return;
        }
        double d2 = dx * dx + dy * dy;
        double e = table(d2 / interact->R2);
        *fx = interact->F0 * e * dx;
        *fy = interact->F0 * e * dy;
    }
};

#endif