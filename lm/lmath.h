// lm.h : Header file for your target.

#pragma once

#include <cmath>

namespace lm
{
	struct Vec2
	{
		inline Vec2() : Vec2(0.f, 0.f) {}
		inline Vec2(float _x, float _y) : x(_x), y(_y) {}

		inline float& operator[](int i) { return a[i]; }
		inline const float& operator[](int i) const { return a[i]; }

		inline Vec2 operator- () const { return Vec2(-x, -y); }

		inline Vec2 operator+ (const Vec2& v) const { return Vec2(x + v.x, y + v.y); }
		inline Vec2 operator- (const Vec2& v) const { return Vec2(x - v.x, y - v.y); }
		inline Vec2 operator* (const Vec2& v) const { return Vec2(x * v.x, y * v.y); }
		inline Vec2 operator/ (const Vec2& v) const { return Vec2(x / v.x, y / v.y); }

		inline Vec2 operator* (float v) const { return Vec2(x * v, y * v); }
		inline Vec2 operator/ (float v) const { return Vec2(x / v, y / v); }

		inline Vec2 operator*= (float v) { x *= v; y *= v; return *this; }
		inline Vec2 operator+= (const Vec2& v) { return Vec2(x += v.x, y += v.y); }
		inline Vec2 operator-= (const Vec2& v) { return Vec2(x -= v.x, y -= v.y); }
		inline Vec2 operator/= (float v) { x /= v; y /= v; return *this; }

		union
		{
			float a[2];
			struct { float x, y; };
		};
	};

	struct Vec3
	{
		inline Vec3() : Vec3(0.f, 0.f, 0.f) {}
		inline Vec3(float _x, float _y, float _z) { x = _x; y = _y; z = _z; }

		inline float& operator[](int i) { return a[i]; }
		inline const float& operator[](int i) const { return a[i]; }

		inline Vec3 operator- () const { return Vec3(-x, -y, -z); }

		inline Vec3 operator+ (const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
		inline Vec3 operator- (const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
		inline Vec3 operator* (const Vec3& v) const { return Vec3(x * v.x, y * v.y, z * v.z); }
		inline Vec3 operator/ (const Vec3& v) const { return Vec3(x / v.x, y / v.y, z / v.z); }
		inline Vec3 operator* (float v) const { return Vec3(x * v, y * v, z * v); }
		inline Vec3 operator/ (float v) const { return Vec3(x / v, y / v, z / v); }

		inline Vec3 operator+= (const Vec3& v) { return Vec3(x += v.x, y += v.y, z += v.z); }
		inline Vec3 operator-= (const Vec3& v) { return Vec3(x -= v.x, y -= v.y, z -= v.z); }
		inline Vec3 operator*= (float v) { x *= v; y *= v; z *= v; return *this; }
		inline Vec3 operator/= (float v) { x /= v; y /= v; z /= v; return *this; }
        
        inline Vec3 Left() const { return Vec3(y, -x, z); }
        inline Vec3 Right() const { return -Left(); }

		union
		{
			struct { float a[3]; };
			struct { float x, y, z; };
			struct { Vec2 xy; } v2;
		};
	};

	struct Vec4
	{
		inline Vec4() : Vec4(0.f, 0.f, 0.f, 0.f) {}
		inline Vec4(const Vec3& v, float _w) : v3{v} { w = _w; }
		inline Vec4(float _x, float _y, float _z, float _w) { x = _x; y = _y; z = _z; w = _w; }

		inline float& operator[](int i) { return a[i]; }
		inline const float& operator[](int i) const { return a[i]; }

		inline Vec4 operator+= (Vec4 v) { return Vec4(x += v.x, y += v.y, z += v.z, w += v.w); }
		inline Vec4 operator-= (Vec4 v) { return Vec4(x -= v.x, y -= v.y, z -= v.z, w -= v.w); }

		inline Vec4 operator+ (const Vec4& v) const { return Vec4(x + v.x, y + v.y, z + v.z, w + v.w); }
		inline Vec4 operator- (const Vec4& v) const { return Vec4(x - v.x, y - v.y, z - v.z, w - v.w); }
		inline Vec4 operator* (const Vec4& v) const { return Vec4(x * v.x, y * v.y, z * v.z, w * v.w); }
		inline Vec4 operator/ (const Vec4& v) const { return Vec4(x / v.x, y / v.y, z / v.z, w / v.w); }
		inline Vec4 operator* (float v) const { return Vec4(x * v, y * v, z * v, w * v); }
		inline Vec4 operator/ (float v) const { return Vec4(x / v, y / v, z / v, w / v); }

		union
		{
			struct { float a[4]; };
			struct { Vec2 xy; } v2;
			struct { Vec3 xyz; } v3;
			struct { float x, y, z, w; };
		};
	};

	inline float Dot(const Vec2& v1, const Vec2& v2) { return v1.x * v2.x + v1.y * v2.y; }
	inline float Dot(const Vec3& v1, const Vec3& v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }
	inline float Dot(const Vec4& v1, const Vec4& v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w; }

	inline float Length(const Vec2& v) { return sqrtf(Dot(v, v)); }
	inline float Length(const Vec3& v) { return sqrtf(Dot(v, v)); }
	inline float Length(const Vec4& v) { return sqrtf(Dot(v, v)); }

	inline Vec3 Cross(const Vec3& v1, const Vec3& v2) {
		return Vec3(
			v1.y * v2.z - v1.z * v2.y,
			v1.z * v2.x - v1.x * v2.z,
			v1.x * v2.y - v1.y * v2.x
		);
	}

	inline Vec2 Normalize(const Vec2& v)
	{
		float freq = 1.f / Length(v);
		return Vec2(v.x * freq, v.y * freq);
	}

	inline Vec3 Normalize(const Vec3& v)
	{
		float freq = 1.f / Length(v);
		return Vec3(v.x * freq, v.y * freq, v.z * freq);
	}

	inline Vec4 Normalize(const Vec4& v)
	{
		float freq = 1.f / Length(v);
		return Vec4(v.x * freq, v.y * freq, v.z * freq, v.w * freq);
	}

	inline float Distance(const Vec2& v1, const Vec2& v2)
	{
		float dx = v2.x - v1.x, dy = v2.y - v1.y;
		return sqrtf((dx * dx) + (dy * dy));
	}

	inline float Distance(const Vec3& v1, const Vec3& v2)
	{
		float dx = v2.x - v1.x, dy = v2.y - v1.y, dz = v2.z - v1.z;
		return sqrtf((dx * dx) + (dy * dy) + (dz * dz));
	}

	inline float Distance(const Vec4& v1, const Vec4& v2)
	{
		float dx = v2.x - v1.x, dy = v2.y - v1.y, dz = v2.z - v1.z, dw = v2.w - v1.w;
		return sqrtf((dx * dx) + (dy * dy) + (dz * dz) + (dw * dw));
	}

	struct Matrix4
	{
		inline Matrix4() {}
		inline Matrix4(const Vec4& v1, const Vec4& v2, const Vec4& v3, const Vec4& v4) :
			v{v1, v2, v3, v4} {}

		inline Matrix4(float _m11, float  _m12, float  _m13, float _m14,
			float  _m21, float  _m22, float  _m23, float  _m24,
			float  _m31, float  _m32, float  _m33, float  _m34,
			float  _m41, float  _m42, float  _m43, float  _m44)
		{
			m11 = _m11; m12 = _m12; m13 = _m13; m14 = _m14;
			m21 = _m21; m22 = _m22; m23 = _m23; m24 = _m24;
			m31 = _m31; m32 = _m32; m33 = _m33; m34 = _m34;
			m41 = _m41; m42 = _m42; m43 = _m43; m44 = _m44;
		}

		inline Matrix4 operator*(const Matrix4& m) const
		{
			return Matrix4(
				Dot(v.v1, Vec4(m.m11, m.m21, m.m31, m.m41)),
				Dot(v.v1, Vec4(m.m12, m.m22, m.m32, m.m42)),
				Dot(v.v1, Vec4(m.m13, m.m23, m.m33, m.m43)),
				Dot(v.v1, Vec4(m.m14, m.m24, m.m34, m.m44)),

				Dot(v.v2, Vec4(m.m11, m.m21, m.m31, m.m41)),
				Dot(v.v2, Vec4(m.m12, m.m22, m.m32, m.m42)),
				Dot(v.v2, Vec4(m.m13, m.m23, m.m33, m.m43)),
				Dot(v.v2, Vec4(m.m14, m.m24, m.m34, m.m44)),

				Dot(v.v3, Vec4(m.m11, m.m21, m.m31, m.m41)),
				Dot(v.v3, Vec4(m.m12, m.m22, m.m32, m.m42)),
				Dot(v.v3, Vec4(m.m13, m.m23, m.m33, m.m43)),
				Dot(v.v3, Vec4(m.m14, m.m24, m.m34, m.m44)),

				Dot(v.v4, Vec4(m.m11, m.m21, m.m31, m.m41)),
				Dot(v.v4, Vec4(m.m12, m.m22, m.m32, m.m42)),
				Dot(v.v4, Vec4(m.m13, m.m23, m.m33, m.m43)),
				Dot(v.v4, Vec4(m.m14, m.m24, m.m34, m.m44))
			);
		}
        
        inline Matrix4 Transpose() const
        {
            Matrix4 m;
            for(int i = 0; i < 4; i++)
                for(int j = 0; j < 4; j++)
                    m.m[i][j] = this->m[j][i];
            
            return m;
        }

		inline static Matrix4 Identity()
		{
			return Matrix4(
				1.f, 0.f, 0.f, 0.f,
				0.f, 1.f, 0.f, 0.f,
				0.f, 0.f, 1.f, 0.f,
				0.f, 0.f, 0.f, 1.f
			);
		}

		inline static Matrix4 Translation(const Vec3& v)
		{
			return Matrix4(
				1.f, 0.f, 0.f, 0.f,
				0.f, 1.f, 0.f, 0.f,
				0.f, 0.f, 1.f, 0.f,
				v.x, v.y, v.z, 1.f
			);
		}

		inline static Matrix4 Ortho(float w, float h, float zn, float zf)
		{
			return Matrix4(
				2.f / w, 0.f, 0.f,
				0.f, 0.f, 2.f / h, 0.f, 0.f,
				0.f, 0.f, 1.f / (zf - zn), 0.f,
				0.f, 0.f, zn / (zn - zf), 1.f
			);
		}

		inline static Matrix4 Ortho2(float w, float h, float zn, float zf)
		{
			return Matrix4(
				1.f / w, 0, 0,
				0, 0, 1.f / h, 0, 0,
				0, 0, 1.f / (zf - zn), 0.f,
				0.f, 0.f, zn / (zn - zf), 1.f
			);
		}

		inline static Matrix4 LookAtLH(const Vec3& eye, const Vec3& at, const Vec3& up)
		{
			Vec3 zaxis = Normalize(at - eye);
			Vec3 xaxis = Normalize(Cross(up, zaxis));
			Vec3 yaxis = Cross(zaxis, xaxis);

			return Matrix4(
				Vec4(xaxis.x, yaxis.x, zaxis.x, 0.f),
				Vec4(xaxis.y, yaxis.y, zaxis.y, 0.f),
				Vec4(xaxis.z, yaxis.z, zaxis.z, 0.f),
				Vec4(-Dot(xaxis, eye), -Dot(yaxis, eye), -Dot(zaxis, eye), 1.f)
			);
		}

		inline static Matrix4 LookToLH(const Vec3& eye, const Vec3& dir, const Vec3& up)
		{
			Vec3 zaxis = Normalize(dir);
			Vec3 xaxis = Normalize(Cross(up, zaxis));
			Vec3 yaxis = Cross(zaxis, xaxis);

			return Matrix4(
				Vec4(xaxis.x, yaxis.x, zaxis.x, 0.f),
				Vec4(xaxis.y, yaxis.y, zaxis.y, 0.f),
				Vec4(xaxis.z, yaxis.z, zaxis.z, 0.f),
				Vec4(-Dot(xaxis, eye), -Dot(yaxis, eye), -Dot(zaxis, eye), 1.f)
			);
		}

		inline static Matrix4 PerspectiveFovLH(float fovy, float aspect, float zn, float zf)
		{
			const float fovd2 = fovy / 2.f;
			const float yScale = cosf(fovd2) / sinf(fovd2); // cot
			const float xScale = yScale / aspect;

			return Matrix4(
				xScale, 0.f, 0.f, 0.f,
				0.f, yScale, 0.f, 0.f,
				0.f, 0.f, zf / (zf - zn), 1.f,
				0.f, 0.f, -zn * zf / (zf - zn), 1.f
			);
		}

		inline static Matrix4 RotationAxis(const Vec3& v, float angle)
		{
			const float sangle = sinf(angle);
			const float cangle = cosf(angle);
			const float cdiff = 1.f - cangle;

			Matrix4 r;
			r.m[0][0] = cdiff * v.x * v.x + cangle;
			r.m[1][0] = cdiff * v.x * v.y - sangle * v.z;
			r.m[2][0] = cdiff * v.x * v.z + sangle * v.y;
			r.m[3][0] = 0.0f;
			r.m[0][1] = cdiff * v.y * v.x + sangle * v.z;
			r.m[1][1] = cdiff * v.y * v.y + cangle;
			r.m[2][1] = cdiff * v.y * v.z - sangle * v.x;
			r.m[3][1] = 0.0f;
			r.m[0][2] = cdiff * v.z * v.x - sangle * v.y;
			r.m[1][2] = cdiff * v.z * v.y + sangle * v.x;
			r.m[2][2] = cdiff * v.z * v.z + cangle;
			r.m[3][2] = 0.0f;
			r.m[0][3] = 0.0f;
			r.m[1][3] = 0.0f;
			r.m[2][3] = 0.0f;
			r.m[3][3] = 1.0f;

			return r;
		}

		union
		{
			struct { float m[4][4]; };
			struct { Vec4 v1, v2, v3, v4; } v;
			struct
			{
				float m11, m12, m13, m14,
					m21, m22, m23, m24,
					m31, m32, m33, m34,
					m41, m42, m43, m44;
			};
		};
	};

	inline Matrix4 Mul(const Matrix4& m1, const Matrix4& m2) { return m1 * m2; }
	inline Vec3 Mul(const Vec3& v, const Matrix4& m)
	{
		return Vec3(
			m.m[0][0] * v.x + m.m[1][0] * v.y + m.m[2][0] * v.z + m.m[3][0],
			m.m[0][1] * v.x + m.m[1][1] * v.y + m.m[2][1] * v.z + m.m[3][1],
			m.m[0][2] * v.x + m.m[1][2] * v.y + m.m[2][2] * v.z + m.m[3][2]
		);
	}

	inline Vec3 Abs(const Vec3& v) { return Vec3(fabsf(v.x), fabsf(v.y), fabsf(v.z)); }
}

