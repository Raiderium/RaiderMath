module raider.math.vec;

import std.conv;
import std.math;

import raider.math.mat;
static import raider.math.misc;

alias vec2 = Vec!(2, double);
alias vec3 = Vec!(3, double);
alias vec4 = Vec!(4, double);
alias vec2f = Vec!(2, float);
alias vec3f = Vec!(3, float);
alias vec4f = Vec!(4, float);
alias vec2i = Vec!(2, int);
alias vec3i = Vec!(3, int);
alias vec4i = Vec!(4, int);
alias vec2u = Vec!(2, uint);
alias vec3u = Vec!(3, uint);
alias vec4u = Vec!(4, uint);

//TODO Fixed-points.
//This would be nice:
// vec3q!(1,15)
//Common Q-types
// vec3q1_15
// vec4q16

package template isVecType(T)
{
    const isVecType = is(T == double) || is(T == float) || is(T == int) || is(T == uint);
}

struct Vec(int _D, _F) if( 2 <= _D && _D <= 4 && isVecType!_F)
{
    static const int D = _D;
    alias _F F;

    alias Vec!(D, F) Vec_;

    F[D] f = 0;

    this(T)(in T s)    if(isVecType!T) { f[] = cast(F)s; }
    this(T)(in T* o)   if(isVecType!T) { static foreach(x; 0..D) f[x] = cast(F) o[x]; }
    this(T)(in T[] o)  if(isVecType!T) { static foreach(x; 0..D) if(x < o.length) f[x] = cast(F)o[x]; }
    //this(T)(in T[D] o) if(isVecType!T) { mixin(fer!("f[x] = cast(F) o[x];")); }
    this(T)(in Vec!(D, T) o)           { static foreach(x; 0..D) f[x] = cast(F) o[x]; }
    void opAssign(T)(auto ref in Vec!(D, T) o) if(!is(T == F)) { static foreach(x; 0..D) f[x] = cast(F) o[x]; }
    void opAssign(T)(T o) if(isVecType!T) { static foreach(x; 0..D) f[x] = cast(F)(o); }

    //TODO Why do method overloads need the same template parameters?
    static if(D == 2) this(T=int)(in F x, in F y)                 { f = [x, y]; }
    static if(D == 3) this(T=int)(in F x, in F y, in F z)         { f = [x, y, z]; }
    static if(D == 4) this(T=int)(in F x, in F y, in F z, in F w) { f = [x, y, z, w]; }

    const F opIndex(const size_t i) { return f[i]; }
    ref   F opIndex(const size_t i) { return f[i]; }

    F[] opSlice() { return f[]; }
    F[] opSlice(size_t x, size_t y) { return f[x..y]; }

    void opSliceAssign(in F[] t) { f[] = t[]; }
    void opSliceAssign(in F[] t, size_t x, size_t y) { f[x..y] = t[]; }

    @property F* ptr()
    {
        return f.ptr;
    }

    //this op= vec
    void opOpAssign(string op, T)(auto ref in Vec!(D, T) o)
    {
        static foreach(x; 0..D) mixin("f[x] = cast(F)(f[x]" ~ op ~ "o.f[x]);");
    }

    //this op= scalar
    void opOpAssign(string op, T)(T o) if(isVecType!T)
    {
        static foreach(x; 0..D) mixin("f[x] = cast(F)(f[x]" ~ op ~ "o);");
    }

    static if(is(F == float) || is(F == double))
    {
        bool isZero()
        {
            static foreach(x; 0..D) if(fabs(f[x]) > F.epsilon) return false;
            return true;
        }

        bool eq()(auto ref in Vec_ o)
        {
            static foreach(x; 0..D) if(!raider.math.misc.eq(f[x], o[x])) return false;
            return true;
        }
    }

    static immutable Vec_ zero = Vec_();

    //op this
    Vec_ opUnary(string s)() if(s == "-")
    {
        Vec_ result = void;
        static foreach(x; 0..D) result[x] = -f[x];
        return result;
    }

    static if(isMatType!F)
    {
        //vec = this * matrix
        Vec_ opBinary(string op, T)(auto ref in Mat!(D, T) o) if(op == "*")
        {
            Vec_ result;
            static foreach(r; 0..D)
                static foreach(c; 0..D)
                    result[c] += f[r]*o[r][c];
            return result;
        }
    }

    //vec = this op vec
    Vec!(D, typeof(T + F)) opBinary(string op, T)(auto ref in Vec!(D, T)o) const
    {
        typeof(return) r = void;
        static foreach(x; 0..D) mixin("r[x] = f[x]" ~ op ~ " o[x];");
        return r;
    }

    //vec = this op scalar
    Vec!(D, typeof(T + F)) opBinary(string op, T)(T o) const if(isVecType!T)
    {
        typeof(return) r = void;
        static foreach(x; 0..D) mixin("r[x] = f[x]" ~ op ~ " o;");
        return r;
    }

    ///Returns a string containing the alias of the vec instantiation followed by its elements inside parentheses.
    //eg vec3f(0.0, 1.2, 9.8)
    string toString() const
    {
        string result = "vec" ~ to!string(D);
        static if(is(F == double)) result ~= "(";
        static if(is(F == float)) result ~= "f(";
        static if(is(F == uint)) result ~= "u(";
        static if(is(F == int)) result ~= "i(";

        foreach(uint x, F ff; f)
        {
            result ~= to!string(ff);
            if(x != D-1) result ~= ", ";
        }
        result ~= ")";
        return result;
    }

    ///Returns a copy with the absolute value of each element.
    Vec_ abs() const
    {
        Vec_ r = void;
        static foreach(x; 0..D) r[x] = std.math.abs(f[x]);
        return r;
    }

    ///Angle 'twixt this and another vec with same dimensions.
    double angle(T)(auto ref in Vec!(D, T) o) const
    {
        Vec!(D, double) v = o;
        Vec!(D, double) t = this;
        v.normalize();
        t.normalize();
        return acos(t.dot(v));
    }

    ///Set length to 1.
    void normalize()
    {
        length = 1.0;
    }

    ///Clamp length between a minimum and maximum value.
    void clampLength(double min, double max)
    {
        double length2 = dot(this, this);
        if(length2 != 0.0)
        {
            if(length2 > max*2) this *= max / sqrt(length2);
            if(length2 < min*2) this *= min / sqrt(length2);
        }
    }

    alias clampLength clampMagnitude;

    ///Length of vector.
    @property double length()
    {
        return sqrt(dot(this, this));
    }

    ///Squared length of vector.
    @property double lengthSquared() const
    {
        return dot(this, this);
    }

    ///ditto
    @property void length(double value)
    {
        double length2 = dot(this, this);
        if(length2 != 0.0) this *= (value / sqrt(length2));
    }

    alias length magnitude;
    
    @property double manhattan() const
    {
        double mh;
        static foreach(x; 0..D) mh += std.math.abs(f[x]);
        return mh;
    }
    
    static if(D == 2 && !is(T == uint))
    {
        ///Rotate a signed vec2 clockwise.
        void roll()
        {
            F temp = f[0];
            f[0] = f[1];
            f[1] = -temp;
        }
    }
}

/**
 * Dot product
 */
double dot(int D, T)(auto ref in Vec!(D, T) a, auto ref in Vec!(D, T) b)
{
    double r = 0.0;
    static foreach(x; 0..D) r += a[x] * b[x];
    return r;
}

/**
 * Cross product.
 *
 * With your right hand, point index finger forwards,
 * middle finger to the left, and thumb upwards.
 * thumb = cross(index, middle)
 *
 * In other words, this is a right cross. Powie!
 * Only available on vec3 and vec3f.
 */
Vec!(3, F) cross(F)(auto ref in Vec!(3, F) a, auto ref in Vec!(3, F) b)
    if(is(F == float) || is(F == double))
{
    Vec!(3, F) r = void;
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
    return r;
}

unittest
{
    //UUUUAAAAAGH I HAVE NO IDEA WHAT I'M DOING
    //TODO Proper unittests.

    bool eq(T)(T x, T y) { return feqrel(x, y) > 20; }

    //Constructors
    vec3 a3d = vec3();                   assert(eq(a3d[0], 0.0) && eq(a3d[1], 0.0));
    vec2  a2d = vec2(-2.4, 0);           assert(eq(a2d[0], -2.4) && eq(a2d[1], 0.0));
    vec3i a3i = vec3i(-2, 0, 3);         assert(a3i[0] == -2 && a3i[1] == 0 && a3i[2] == 3);
    
    vec4f a4f = vec4f(-2.6, 0, 4, 2.3);  assert(eq(a4f[0], -2.6f) &&
                                                eq(a4f[1], 0.0f) &&
                                                eq(a4f[2], 4.0f) &&
                                                eq(a4f[3], 2.3f));
    //Casting constructors
    vec3 a = vec3(1.0, -2.3, 0.0);
    vec3i b = a;                         assert(b[0] == 1 && b[1] == -2 && b[2] == 0);

    //Casts
    b = a;                               assert(b[0] == 1 && b[1] == -2 && b[2] == 0);

    //[] reference assignment
    b[2] = 2;                            assert(b[2] == 2);

    //Binary assignment operators
    a *= 2.0;                            assert(eq(a[0], 2.0) && eq(a[1], -4.6) && eq(a[2], 0.0));
    a /= b;                              assert(eq(a[0], 2.0) && eq(a[1], 2.3) && eq(a[2], 0.0));

    //Binary operators
    b = a + a - a*4;                     assert(b[0] == -4 && b[1] == -4 && b[2] == 0);

    //abs
    a = vec3(-2.0, -0.0, 1.0);
    a = a.abs();                         assert(eq(a[0], 2.0) && eq(a[1], 0.0) && eq(a[2], 1.0));

    //Scalar assignment
    a = 4.0;                             assert(eq(a[0], 4.0) && eq(a[1], 4.0) && eq(a[2], 4.0));
}

//In quat:
/*
Vec_ opBinary(string op, F)(auto ref in Vec v) if(op == "*")
{
    F x = f[1] * v[2] - f[2] * v[1];
    F y = f[2] * v[0] - f[0] * v[2];
    F z = f[0] * v[1] - f[1] * v[0];
    x *= 2; y *= 2; z *= 2; //shifts with fixed point
    Vec_ r = void;
    r[0] = v[0] + f[3] * x + f[1] * z - f[2] * y;
    r[1] = v[1] + f[3] * y + f[2] * x - f[0] * z;
    r[2] = v[2] + f[3] * z + f[0] * y - f[1] * x;
    return r;
}*/
