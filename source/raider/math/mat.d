module raider.math.mat;

import std.math;
import std.conv;
import raider.math.vec;

alias mat2 = Mat!(2, double);
alias mat3 = Mat!(3, double);
alias mat4 = Mat!(4, double);
alias mat2f = Mat!(2, float);
alias mat3f = Mat!(3, float);
alias mat4f = Mat!(4, float);

package template isMatType(T)
{
    const isMatType = is(T == double) || is(T == float);
}

struct Mat(int _D, _F) if(2 <= _D && _D <= 4 && isMatType!_F)
{
    static const int D = _D;
    alias F = _F;

    alias Vec_ = Vec!(D, F);
    alias Mat_ = Mat!(D, F);

    union
    {
        F[D*D] f = 0;
        Vec_[D] v;
    }

    static if(D == 2)
    {
        this(T=int)(Vec_ v0, Vec_ v1)
        { v[] = [v0, v1]; }

        this(T=int)(F f0, F f1, 
                    F f2, F f3)
        { f[] = [f0, f1, f2, f3]; }

        static enum Mat_ identity = Mat_(1,0,
                                         0,1);
    }

    static if(D == 3)
    {
        this(T=int)(Vec_ v0, Vec_ v1, Vec_ v2)
        { v[] = [v0, v1, v2]; }

        this(T=int)(F f0, F f1, F f2, 
                    F f3, F f4, F f5, 
                    F f6, F f7, F f8)
        { f[] = [f0, f1, f2, f3, f4, f5, f6, f7, f8]; }

        static enum Mat_ identity = Mat_(1,0,0,
                                         0,1,0,
                                         0,0,1);
    }

    static if(D == 4)
    {
        this(T=int)(Vec_ v0, Vec_ v1, Vec_ v2, Vec_ v3, Vec_ v4)
        { v[] = [v0, v1, v2, v3, v4]; }

        this(T=int)(F f0,  F f1,  F f2,  F f3,  
                    F f4,  F f5,  F f6,  F f7,  
                    F f8,  F f9,  F f10, F f11,
                    F f12, F f13, F f14, F f15)
        { f[] = [f0,  f1,  f2,  f3, 
                 f4,  f5,  f6,  f7,
                 f8,  f9,  f10, f11,
                 f12, f13, f14, f15]; }

        static enum Mat_ identity = Mat_(1,0,0,0,
                                         0,1,0,0,
                                         0,0,1,0,
                                         0,0,0,1);

        //Construct affine orientation and translation matrix
        this(T=int)(auto ref in Mat!(3, T) o, Vec!(3, T) p)
        {
            f[0..3]   = o.f[0..3]; f[3]  = 0.0;
            f[4..7]   = o.f[3..6]; f[7]  = 0.0;
            f[8..11]  = o.f[6..9]; f[11] = 0.0;
            f[12..15] =   p[0..3]; f[15] = 1.0;
        }
    }

    //Construct submatrix, cutting out row r and column c.
    static if(D < 4)
    {
        this(T=int)(auto ref in Mat!(D+1, F) o, uint r, uint c)
        {
            assert(r < D+1);
            assert(c < D+1);

            static foreach(ri; 0..D)
                static foreach(ci; 0..D)
                    v[ri][ci] = o[ri+(ri >= r)][ci+(ci >= c)];
        }
    }

    this(T)(const T* o) if(isMatType!T)
    {
        static foreach(x; 0..D*D) f[x] = cast(F) o[x];
    }

    this(T)(auto ref in Mat!(D, T) o)
    {
        static foreach(x; 0..D*D) f[x] = cast(F) o.f[x];
    }

    void opAssign(T)(auto ref in Mat!(D, T) o)
    {
        static foreach(x; 0..D*D) f[x] = cast(F) o.f[x];
    }

    string toString() const
    {
        string result = "mat" ~ to!string(D);
        static if(is(F == double)) result ~= "(";
        static if(is(F == float)) result ~= "f(";

        foreach(uint x, F ff; f)
        {
            result ~= to!string(ff);
            if(x != (D*D)-1) result ~= ", ";
        }
        result ~= ")";
        return result;
    }

    /* Checks for diagonal values of exactly 1.0.
     * Only use on matrices known to be orthogonal.
     * Does NOT check off-diagonal values are 0.0.
     */
    bool isIdentityFast()
    {
        static foreach(x; 0..D) if(v[x][x] != 1.0) return false;
        return true;
    }

    const Vec_ opIndex(const size_t i) { return v[i]; }
    ref   Vec_ opIndex(const size_t i) { return v[i]; }

    @property F* ptr() { return f.ptr; }

    //matrix = this * matrix
    Mat_ opBinary(string op)(auto ref in Mat_ o) if(op == "*")
    {
        Mat_ result;
        static foreach(i; 0..D)
            static foreach(j; 0..D)
                static foreach(k; 0..D)
                    result[i][j] += v[i][k] * o[k][j];

        /*
        static if(D == 4) //TODO Profile static foreach vs explicit unroll for 4x4. 
        {
            result[0][0] = v[0][0] * o[0][0] + v[0][1] * o[1][0] + v[0][2] * o[2][0] + v[0][3] * o[3][0];
            result[1][0] = v[1][0] * o[0][0] + v[1][1] * o[1][0] + v[1][2] * o[2][0] + v[1][3] * o[3][0];
            result[2][0] = v[2][0] * o[0][0] + v[2][1] * o[1][0] + v[2][2] * o[2][0] + v[2][3] * o[3][0];
            result[3][0] = v[3][0] * o[0][0] + v[3][1] * o[1][0] + v[3][2] * o[2][0] + v[3][3] * o[3][0];

            result[0][1] = v[0][0] * o[0][1] + v[0][1] * o[1][1] + v[0][2] * o[2][1] + v[0][3] * o[3][1];
            result[1][1] = v[1][0] * o[0][1] + v[1][1] * o[1][1] + v[1][2] * o[2][1] + v[1][3] * o[3][1];
            result[2][1] = v[2][0] * o[0][1] + v[2][1] * o[1][1] + v[2][2] * o[2][1] + v[2][3] * o[3][1];
            result[3][1] = v[3][0] * o[0][1] + v[3][1] * o[1][1] + v[3][2] * o[2][1] + v[3][3] * o[3][1];

            result[0][2] = v[0][0] * o[0][2] + v[0][1] * o[1][2] + v[0][2] * o[2][2] + v[0][3] * o[3][2];
            result[1][2] = v[1][0] * o[0][2] + v[1][1] * o[1][2] + v[1][2] * o[2][2] + v[1][3] * o[3][2];
            result[2][2] = v[2][0] * o[0][2] + v[2][1] * o[1][2] + v[2][2] * o[2][2] + v[2][3] * o[3][2];
            result[3][2] = v[3][0] * o[0][2] + v[3][1] * o[1][2] + v[3][2] * o[2][2] + v[3][3] * o[3][2];

            result[0][3] = v[0][0] * o[0][3] + v[0][1] * o[1][3] + v[0][2] * o[2][3] + v[0][3] * o[3][3];
            result[1][3] = v[1][0] * o[0][3] + v[1][1] * o[1][3] + v[1][2] * o[2][3] + v[1][3] * o[3][3];
            result[2][3] = v[2][0] * o[0][3] + v[2][1] * o[1][3] + v[2][2] * o[2][3] + v[2][3] * o[3][3];
            result[3][3] = v[3][0] * o[0][3] + v[3][1] * o[1][3] + v[3][2] * o[2][3] + v[3][3] * o[3][3];
        }
        */

        return result;
    }

    //vec = this * vec
    Vec_ opBinary(string op)(auto ref in Vec_ o) if(op == "*")
    {
        //TODO Profile these different levels of unrolling.
        //I would love to trust the compiler, but trust
        //must be earned in this specific situation.

        Vec_ result;
        static foreach(r; 0..D)
            static foreach(c; 0..D)
                result[r] += o[c]*v[r][c];
        return result;

        /*
        Vec_ result;
        foreach(r; 0..D)
            foreach(c; 0..D)
                result[r] += o[c]*v[r][c];
        return result;
        */

        /*
        static if(D == 2)
            return Vec_(
                v[0][0] * o[0] + v[0][1] * o[1], 
                v[1][0] * o[0] + v[1][1] * o[1]);

        static if(D == 3)
            return Vec_(
                v[0][0] * o[0] + v[0][1] * o[1] + v[0][2] * o[2],
                v[1][0] * o[0] + v[1][1] * o[1] + v[1][2] * o[2],
                v[2][0] * o[0] + v[2][1] * o[1] + v[2][2] * o[2]);

        static if(D == 4)
            return Vec_(
                v[0][0] * o[0] + v[0][1] * o[1] + v[0][2] * o[2] + v[0][3] * o[3],
                v[1][0] * o[0] + v[1][1] * o[1] + v[1][2] * o[2] + v[1][3] * o[3],
                v[2][0] * o[0] + v[2][1] * o[1] + v[2][2] * o[2] + v[2][3] * o[3],
                v[3][0] * o[0] + v[3][1] * o[1] + v[3][2] * o[2] + v[3][3] * o[3]);
        */
        //TODO Test ODE's mysterious aliasing avoidance magic

        /*
        static if(D == 2)
            return Vec_(dot(v[0], o), dot(v[1], o));
        static if(D == 3)
            return Vec_(dot(v[0], o), dot(v[1], o), dot(v[2], o));
        static if(D == 4)
            return Vec_(dot(v[0], o), dot(v[1], o), dot(v[2], o), dot(v[3], o));
        */

        /*
        Vec_ result = void;
        foreach(r; 0..D) result[r] = dot(v[r], o);
        return result;
        */


    }

    Mat_ scaled()(auto ref in Vec_ o)
    {
        Mat_ result = void;
        static foreach(c; 0..D) result[c] = o * v[c];
        return r;
    }

    //matrix = this * scalar
    Mat_ opBinary(string op)(F o) if(op == "*")
    {
        Mat_ r = void;
        r.f[0..D*D] = f[0..D*D]*o;
        return r;
    }

    Mat_ transposed()
    {
        Mat_ t = void;
        static foreach(r; 0..D)
            static foreach(c; 0..D)
                t[c][r] = v[r][c];

        return t;
    }

    /**
     * Compute the adjoint matrix (that is, the transposed cofactor matrix).
     * 
     * Inefficient for 4x4 inversion. Use invert().
     */
    Mat_ adjoint()
    {
        //There is no 1x1 submatrix constructor. Implement 2x2 adjoint directly.
        static if(D == 2) return Mat_(f[3], -f[2], f[1], -f[0]);
        else
        {
            Mat_ adj = void;

            static foreach(r; 0..D)
                static foreach(c; 0..D)
                    adj[c][r] = Mat!(D-1, F)(this, r, c).determinant * ((r+c)%2 ? -1 : 1);

            return adj;
        }
    }

    F determinant()
    {
        static if(D == 2) return f[0]*f[3] - f[1]*f[2];
        else
        {
            F r = 0.0;
            static foreach(c; 0..D)
                r += v[0][c] * Mat!(D-1, F)(this, 0, c).determinant * (c%2 ? -1 : 1);
            return r;
        }
    }

    static if(D == 2)
    {
        void invert()
        {
            F det = determinant;
            if(det == 0) return;

            f[1] /= -det;
            f[2] /= -det;
            F temp = f[3]/det;
            f[3] = f[0]/det;
            f[0] = temp;
        }

        Mat_ rotation(F angle)
        {
            angle *= PI/180;
            F c = cos(angle);
            F s = sin(angle);

            return Mat_(c, s, -s, c);
        }
    }

    static if(D == 3)
    {
        void invert()
        {
            //m stands for minor
            F m0 = f[4]*f[8] - f[5]*f[7];
            F m1 = f[8]*f[3] - f[5]*f[6];
            F m2 = f[3]*f[7] - f[4]*f[6];
            F m3 = f[1]*f[8] - f[2]*f[7];
            F m4 = f[0]*f[8] - f[2]*f[6];
            F m5 = f[0]*f[7] - f[1]*f[6];
            F m6 = f[1]*f[5] - f[2]*f[4];
            F m7 = f[0]*f[5] - f[2]*f[3];
            F m8 = f[0]*f[4] - f[1]*f[3];

            F det = 
                f[0] * m0 -
                f[1] * m1 +
                f[2] * m2;

            if(det == 0.0) return;

            double invdet = 1.0/det;

            //Notice this adds the sign to produce the cofactor matrix, 
            //transposes to produce the adjoint matrix,
            //and divides by the determinant to produce the inverse matrix.
            f[0] =  m0 * invdet;
            f[1] = -m3 * invdet;
            f[2] =  m6 * invdet;
            f[3] = -m1 * invdet;
            f[4] =  m4 * invdet;
            f[5] = -m7 * invdet;
            f[6] =  m2 * invdet;
            f[7] = -m5 * invdet;
            f[8] =  m8 * invdet;
        }

        /**
         * Build a rotation matrix.
         * 
         * This accepts either a scalar angle and integer axis (012 = XYZ)
         * or a rotation 'vector', where the direction is the axis, and the 
         * magnitude is the angle. All angles are in radians.
         */
        static Mat_ rotation(double angle, int axis)
        {
            assert(0 <= axis && axis < D);
            //angle *= PI/180;
            //why oh why did we ever choose pi
            //pls embrace tau; the pi is a lie

            F c = cos(angle);
            F s = sin(angle);

            switch(axis)
            {
                case 0: return Mat_( 1, 0, 0, 
                                     0, c, s, 
                                     0,-s, c);
                case 1: return Mat_( c, 0,-s, 
                                     0, 1, 0, 
                                     s, 0, c);
                case 2: return Mat_( c, s, 0, 
                                    -s, c, 0, 
                                     0, 0, 1);
                default: assert(0);
            }
        }

        ///ditto
        static Mat_ rotation(Vec_ r)
        {
            //Adapted from 
            //http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToMatrix/

            F a = r.length;
            if(a != 0.0) r /= a;
            F c = cos(a);
            F s = sin(a);
            F t = 1.0 - c;

            Mat_ m = void;
            m[0][0] = c + r[0] * r[0] * t;
            m[1][1] = c + r[1] * r[1] * t;
            m[2][2] = c + r[2] * r[2] * t;

            F u = r[0] * r[1]*t;
            F v = r[2] * s;
            m[1][0] = u + v;
            m[0][1] = u - v;
            u = r[0] * r[2] * t;
            v = r[1] * t;
            m[2][0] = u - v;
            m[0][2] = u + v;
            u = r[1] * r[2] * t;
            v = r[0] * s;
            m[2][1] = u + v;
            m[1][2] = u - v;

            return m;
        }
    }

    static if(D == 4)
    {
        void invert()
        {
            /*
             * This is an implementation of the Laplace Expansion Theorem which is described here:
             * http://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf
             * 
             * Instead of finding the 4x4 determinant through the determinants of four 3x3 matrices
             * and four 1x1 matrices (expansion by columns or rows), we find the determinants of 
             * twelve 2x2 matrices (expansion by 2 columns or 2 rows). This requires fewer operations.
             * 
             * And, as it happens, all those 2x2 determinants may be reused while computing the 
             * cofactor matrix, by using particular submatrices.
             * 
             * Finding determinant by...
             * row expansion = 40 multiplications, 23 additions
             * 2 row expansion = 30 mult, 17 add
             * 
             * Finding cofactors by...
             * naive expansion = 160 mult, 80 add
             * reusing determinants = 48 mult, 32 add
             * 
             * Finding inverse by...
             * row expansion & naive expansion = 216 mult, 103 add
             * 2 row expansion & reusing determinants = 94 mult, 49 add
             */

            //s stands for submatrix
            F s0 = f[0] * f[5] - f[4] * f[1];
            F s1 = f[0] * f[6] - f[4] * f[2];
            F s2 = f[0] * f[7] - f[4] * f[3];
            F s3 = f[1] * f[6] - f[5] * f[2];
            F s4 = f[1] * f[7] - f[5] * f[3];
            F s5 = f[2] * f[7] - f[6] * f[3];

            //c stands for complementary submatrix
            F c5 = f[10] * f[15] - f[14] * f[11];
            F c4 = f[9]  * f[15] - f[13] * f[11];
            F c3 = f[9]  * f[14] - f[13] * f[10];
            F c2 = f[8]  * f[15] - f[12] * f[11];
            F c1 = f[8]  * f[14] - f[12] * f[10];
            F c0 = f[8]  * f[13] - f[12] * f[9];

            F det = s0 * c5 - s1 * c4 + s2 * c3 +
                    s3 * c2 - s4 * c1 + s5 * c0;

            if(det == 0.0) return;

            F invdet = 1.0/det;

            Mat_ i = void;

            i.f[0]  = ( f[5]  * c5 - f[6]  * c4 + f[7]  * c3) * invdet;
            i.f[1]  = (-f[1]  * c5 + f[2]  * c4 - f[3]  * c3) * invdet;
            i.f[2]  = ( f[13] * s5 - f[14] * s4 + f[15] * s3) * invdet;
            i.f[3]  = (-f[9]  * s5 + f[10] * s4 - f[11] * s3) * invdet;

            i.f[4]  = (-f[4]  * c5 + f[6]  * c2 - f[7]  * c1) * invdet;
            i.f[5]  = ( f[0]  * c5 - f[2]  * c2 + f[3]  * c1) * invdet;
            i.f[6]  = (-f[12] * s5 + f[14] * s2 - f[15] * s1) * invdet;
            i.f[7]  = ( f[8]  * s5 - f[10] * s2 + f[11] * s1) * invdet;

            i.f[8]  = ( f[4]  * c4 - f[5]  * c2 + f[7]  * c0) * invdet;
            i.f[9]  = (-f[0]  * c4 + f[1]  * c2 - f[3]  * c0) * invdet;
            i.f[10] = ( f[12] * s4 - f[13] * s2 + f[15] * s0) * invdet;
            i.f[11] = (-f[8]  * s4 + f[9]  * s2 - f[11] * s0) * invdet;

            i.f[12] = (-f[4]  * c3 + f[5]  * c1 - f[6]  * c0) * invdet;
            i.f[13] = ( f[0]  * c3 - f[1]  * c1 + f[2]  * c0) * invdet;
            i.f[14] = (-f[12] * s3 + f[13] * s1 - f[14] * s0) * invdet;
            i.f[15] = ( f[8]  * s3 - f[9]  * s1 + f[10] * s0) * invdet;

            this = i;
        }

        void invertAffine()
        {
            //Invert rotation and scale
            Mat!(3, F) r = Mat!(3, F)(this, 3, 3);
            r.invert;

            //Invert translation
            Vec!(3, F) t = Vec!(3, F)(-f[12], -f[13], -f[14]);
            t = t * r;

            this = Mat_(r, t);
        }
    }
}

unittest
{
    mat3f m3fa;
    mat3f m3fb;
    auto r = m3fa * m3fb;
}
