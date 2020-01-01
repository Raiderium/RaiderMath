module raider.math.misc;

import std.traits : isFloatingPoint, isIntegral;
import raider.math.vec;

public import std.math : PI, fmin, fmax, fabs;

T fmin3(T)(T a, T b, T c)
{
    return a < b ? (a < c ? a : c) : (b < c ? b : c);
}

T clamp(T)(T f, T min, T max) 
{
    return f < min ? min : (f > max ? max : f);
}

bool eqZero(T)(T f) if(isFloatingPoint!T)
{
    return fabs(f) < T.epsilon;
}

//Inexact equality check
bool eq(T, T margin = T.epsilon)(T a, T b) if(isFloatingPoint!T)
{
    T ab = fabs(a - b);
    if(ab < margin) return true;
    a = fabs(a); b = fabs(b);
    return ab < margin * (b > a ? b : a);
}

/**
 * Float to normalised integer.
 * Converts any float type in the range [0, 1] to a normalised 
 * unsigned integer in range [0, T.max]. Input is clamped.
 */
T ftni(T, F)(F f) if(isIntegral!T && isUnsigned!(T) && isFloatingPoint!(F))
{
    return cast(T)(clamp(f, 0.0, 1.0) * cast(F)T.max);
}

/**
 * Normalised integer to float.
 * Converts a normalised unsigned integer in the 
 * range [0, T.max] to any float type in the range [0, 1].
 */
F nitf(F, T)(T t) if(isIntegral!T && isUnsigned!(T) && isFloatingPoint!(F))
{
    return cast(F)(t) / cast(F)(T.max);
}

///Upper multiple of
private auto umo(T)(T v, T b) if(isIntegral!T)
{
    assert(b);
    auto r = v % b;
    return r ? v + b - r : v;
    
    /* Explicitly branchless variant for purely academic purposes:
     * auto r = v % b;
     * auto s = -T(!r);
     * return ((v + b - r) & ~s) | (v & s);
     */
}

unittest
{
    assert(umo(0, 10) == 0);
    assert(5.umo(4) == 8);
    assert(umo(9, 10) == 10);
    assert(umo(10, 10) == 10);
    assert(umo(11, 10) == 20);
    assert(umo(8, 1) == 8);
}

///Upper power of 2
auto upo2(T)(T v) if(isIntegral!T)
{
    size_t p = 1;
    while(p <= x) p <<= 1;
    return p;
}

///Is power of 2
bool ispo2(T)(T v)
{
    return v && !(v & (v - 1));
}

///Linear interpolation
T lerp(T)(double t, T a, T b)
{
    return a + (b - a)*t;
}

///Spherical linear interpolation
vec3 slerp(F)(F t, Vec!(3, F) a, Vec!(3, F) b)
{
    F dot = clamp(a.dot(b), -1.0, 1.0);
    F theta = acos(dot)*t;
    vec3 r = b - a*dot;
    return a*cast(F)cos(theta) + r*cast(F)sin(theta);
}

///Normalised linear interpolation
vec3 nlerp(F)(F t, Vec!(3, F) a, Vec!(3, F) b)
{
    auto result = lerp(t, a, b);
    result.normalize();
    return result;
}

///Cubic bezier interpolation
T bez(T)(double t, T p0, T p1, T p2, T p3)
{
    //De Casteljau's algorithm
    T l0, l1, l2;
    
    l0 = lerp(t, p0, p1);
    l1 = lerp(t, p1, p2);
    l2 = lerp(t, p2, p3);
    
    l0 = lerp(t, l0, l1);
    l1 = lerp(t, l1, l2);
    
    l0 = lerp(t, l0, l1);
    
    return l0;

    /* For the sake of curiosity, this is the polynomial form
     * A really curious person might profile this against de casteljau

    T c0, c1, c2, c3;
    
    c0 = p0;
    c1 = 3*(p1 - p0);
    c2 = 3*(p0 - 2*p1 + p2);
    c3 = p3 - p0 + 3*(p1 - p2);

    T r = c0;
    r+=t*c1; t*=t;
    r+=t*c2; t*=t;
    r+=t*c3;
    return r;
    */
}

///Use forward differencing to evaluate multiple points on a cubic bezier.
void bez(T)(ref T[] result, T p0, T p1, T p2, T p3)
{
    double t = 1.0 / result.length;
    double t2 = t*t;
    double t3 = t2*t;

    //Find polynomial coefficients
    T c1 = 3.0*(p1 - p0);
    T c2 = 3.0*(p0 - 2.0*p1 + p2);
    T c3 = p3 - p0 + 3.0*(p1 - p2);

    //Find initial forward differences
    T c3t3 = c3*t3;
    T c2t2 = c2*t2;
    T c3t36 = c3t3*6.0;

    T f0 = p0;
    T f1 = c3t3 + c2t2 + c1*t;
    T f2 = c3t36 + 2.0*c2t2;
    T f3 = c3t36;

    //Create curvy goodness
    foreach(ref T r; result)
    {
        r = f0;
        f0 += f1;
        f1 += f2;
        f2 += f3;
    }
}

///Cubic hermite interpolation
double her4(double t, double p0, double m0, double p1, double m1)
{
    double t12, t22, t32;
    double t13, t23, t33;

    t12 = t*t;              //1*t^2
    t22 = t12 + t12;        //2*t^2
    t32 = t22 + t12;        //3*t^2;

    t13 = t12 * t;          //1*t^3
    t23 = t13 + t13;        //2*t^3
    t33 = t23 + t13;        //3*t^3;
    
    //Basis functions
    double b0, b1, b2, b3;
    
    b0 = t23 - t32 + 1;     //For p0
    b1 = t13 - t22 + t;     //For m0
    b2 = t32 - t23;         //For p1
    b3 = t13 - t12;         //For m1

    return b0*p0 + b1*m0 + b2*p1 + b3*m1;
}
