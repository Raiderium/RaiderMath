module raider.math.curve;

import raider.tools.array;
import raider.math;

/**
 * Animated float.
 *
 * A curve is a sequence of time-value pairs with slope handles.
 * 
 * Two successive points form a segment. The first point's out 
 * handle and the second point's in handle control the  
 * interpolation of the segment. They are constrained to point 
 * inwards, that is, right and left respectively, and also 
 * to not overlap on the time axis (the curve must not 'fold').
 * 
 * If a handle's length is zero, it changes the interpolation 
 * method. Two nonzero handles mean the curve is a cubic bezier.
 * Two zero handles mean a straight line. If one handle has 
 * nonzero length, the curve is linearly extrapolated along it
 * until meeting the time coordinate of the other point.
 * 
 * Extrapolation beyond the point sequence is linear, using the 
 * slope of the nearest handle or segment.
 * 
 * TODO Pre/post infinity looping.
 */
class Curve
{private:

    float t = 0.0;
    size_t i = 0; //closest point to the left of t
    bool loops;
    Array!Point points;

public:

    void add(Point point)
    {
        points.add!"a < b"(point);
    }

    /**
     * The value of the curve at the current time.
     */
    @property float value()
    {
        if(points.length == 0) return 0.0;
        if(points.length == 1)
        {
            Point p = points[0];
            if(t < p.t)
                return !p.hasIn ? p.v : p.v + (t - p.t)*p.slopeIn;
            if(t > p.t)
                return !p.hasOut ? p.v : p.v + (t - p.t)*p.slopeOut;
            return p.v;
        }

        //Get points defining segment
        Point a = points[i];
        Point b = points[i+1];

        float ab_time = b.t - a.t;
        float ab_value = b.v - a.v;
        float ab_slope = ab_time == 0.0 ? 0.0 : ab_value / ab_time;

        //Extrapolate left
        if(t < a.t)
        {
            //Linear extrapolation slope: First handle in, then 
            //handle out, then segment (if linear), then 0.
            float slope = 0.0;
            if(a.hasIn) slope = a.slopeIn;
            else if(a.hasOut) slope = a.slopeOut;
            else if(!b.hasIn) slope = ab_slope;

            return a.v + (t - a.t)*slope;
        }

        //Extrapolate right
        if(b.t < t)
        {
            float slope = 0.0;
            if(b.hasOut) slope = b.slopeOut;
            else if(b.hasIn) slope = b.slopeIn;
            else if(!a.hasOut) slope = ab_slope;

            return b.v + (t - b.t)*slope;
        }

        //Get normalized time
        float fac = (a.t == b.t) ? 0.0 : (t - a.t) / ab_time;

        //Linear interpolation
        if(!a.hasOut && !b.hasIn) return lerp(fac, a.v, b.v);

        //Linear from a
        if(!b.hasIn) return a.v + (t - a.t)*a.slopeOut;

        //Linear from b
        if(!a.hasOut) return b.v + (t - b.t)*b.slopeIn;

        //Bezier
        float x1 = a.handleOut[0] / ab_time;
        float x2 = (ab_time + b.handleIn[0]) / ab_time;

        float x_solve = solveMonotonicCubicBezier(fac, x1, x2);
        return bez(x_solve, a.v, a.v + a.handleOut[1], b.v + b.handleIn[1], b.v);
    }

    @property float time()
    {
        return t;
    } 

    @property void time(double value)
    {
        addTime(value - t);
    }

    void addTime(double dtime)
    {
        t += dtime;

        //Find the segment closest to the new time, and set i to 
        //the index of the leftmost point of the segment.
        if (dtime > 0.0) 
        {
            while(points[i+1].t < t)
            {
                if(i == points.length-1) break;
                i++;
            }
        }
        else
        {
            while(points[i].t > t) 
            {
                if(i == 0) break;
                i--;
            }
        }
    }

    /**
     * UUUUUUUUUUUUUHGHHGHHGHHHHGGGHGHGHGGGGGG.
     * 
     * This is a stupidly complex function that people forget to 
     * mention while writing tutorials for bezier animation.
     * 
     * An animation curve can't loop back on itself - that is,
     * at a given time, it can't have more than one solution.
     * That's the 'monotonic' part of this method.
     * 
     * For people who understand such things, this implementation 
     * uses a binary search via subdivision. An analytic solution 
     * provides a precise answer, but the math is facemelting.
     */
    static double solveMonotonicCubicBezier(double x, double p1, double p2) 
    { //p0 and p3 are implicitly 0.0 and 1.0
        assert(0.0 <= p1 && p1 <= p2 && p2 <= 1.0);
        assert(0.0 <= x && x <= 1.0);

        /* Instead of working in a 0 .. 1 floating point space, 
         * this algorithm converts everything to an integer space 
         * 0 .. int.max. This makes dividing by 2, which is 
         * important to the algorithm, 'faster'.
         * 
         * As with all poorly justified optimisations,
         * TODO profile to confirm this.
         */

        immutable uint tolerance = 100;
        uint loopMax = 32;

        //Note that a 64-bit double happily represents 32-bit 
        //int.max during the conversion.
        uint goal = cast(uint)(x*int.max);

        uint error;
        uint t = int.max/2; //Current binary split
        uint s = t/2;

        //Convert the bezier control points
        uint i0 = 0;
        uint i1 = cast(uint)(p1*int.max);
        uint i2 = cast(uint)(p2*int.max);
        uint i3 = int.max;

        //Binary search via de Casteljau subdivision
        do
        {//du hast
            //du hast mich

            //Control points are <= int.max so
            //halfway point computations fit in a uint.
            uint q0 = (i0 + i1) / 2; 
            uint q1 = (i1 + i2) / 2;
            uint q2 = (i2 + i3) / 2;

            uint r0 = (q0 + q1) / 2;
            uint r1 = (q1 + q2) / 2;

            uint result = (r0 + r1) / 2;

            if(goal == result) break;

            //Subdivide
            if(goal < result)
            {
                i1 = q0;
                i2 = r0;
                i3 = result;
                t -= s;
                error = result - goal;
            }
            else
            {
                i0 = result;
                i1 = r1;
                i2 = q2;
                t += s;
                error = goal - result;
            }

            s /= 2;
            --loopMax;
        }
        while(error > tolerance && loopMax);

        return cast(double)t / int.max;
    }
}

struct Point
{
    float t, v; //Time-value pair
    vec2f handleIn;
    vec2f handleOut;

    @property bool hasIn() { return handleIn[0] != 0.0 || handleIn[1] != 0.0; }
    @property bool hasOut() { return handleOut[0] != 0.0 || handleOut[1] != 0.0; }
    @property float slopeIn() { return handleIn[0] == 0.0 ? 0.0 : handleIn[1] / handleIn[0]; }
    @property float slopeOut() { return handleOut[0] == 0.0 ? 0.0 : handleOut[1] / handleOut[0]; }

    this(float time, float value, vec2f a, vec2f b)
    {
        t = time; v = value;
        this.handleIn = a; this.handleOut = b;
    }

    this(float time, float value)
    {
        t = time; v = value;
    }

    int opCmp(ref const Point p) const
    {
        if(t < p.t) return -1;
        if(t > p.t) return 1;
        return 0;
    }
}

unittest
{
    auto c = new Curve;
    c.add(Point(0.0, 1.0, vec2f(-1.0, 0.0), vec2f(2, 4)));
    c.add(Point(8.0, -1.0, vec2f(-4, -2), vec2f(1.0, 0.0)));
    c.time = 4.0;
    import std.stdio : write;
    foreach(x; 0..41)
    {
        c.time = (x / 40.0)*8.0;
        write("(", c.time, ", ", c.value, "), ");
    }
}

/* TODO: Profile this alternative method.

public float newtonRaphson(float x, CubicBezier b)
{
    // Finds the y value of a cubic Bezier curve at a value of x,
    // assuming y is a function of x. Pretty errors ensue if it isn't.

    // The x component of a cubic Bezier curve as a function of t is
    // x = a0*t^3 + a1*t^2 + a2*t + a3
    float a0 = -b.points[0].x + 3*b.points[1].x - 3*b.points[2].x + b.points[3].x;
    float a1 = 3*b.points[0].x - 6*b.points[1].x + 3*b.points[2].x;
    float a2 = -3*b.points[0].x + 3*b.points[1].x;
    float a3 = b.points[0].x - x;

    float t = 0.5f;
    float y;

    // Approximate t using the Newton-Raphson method
    for (int i=0; i<NRiterations; i++) {
        t = t - (a0*pow(t, 3)+a1*pow(t, 2)+a2*t+a3)/(3*a0*pow(t, 2)+2*a1*t+a2);
    }

    // Find the y component, using t
    y = b.points[0].y*pow((1-t), 3) + 3*b.points[1].y*pow((1-t), 2)*t
    + 3*b.points[2].y*(1-t)*t*t + b.points[3].y*pow(t, 3);
    return y;
} */
