#ifndef __OUR_GL_H__
#define __OUR_GL_H__

#include "tgaimage.h"
#include "geometry.h"

extern Matrix ModelView;
extern Matrix Viewport;
extern Matrix Projection;
const float depth = 2000.f;

void viewport(int x, int y, int w, int h);
void projection(float coeff=0.f); // coeff = -1/c
void lookat(Vec3f eye, Vec3f center, Vec3f up);

struct IShader {
    virtual ~IShader();
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    virtual bool fragment(Vec3f bar, TGAColor &color) = 0;
};

struct ZShader {
    virtual ~ZShader();
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor& color) = 0;
};

//void triangle(Vec4f *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer);
void triangle(Vec4f* pts, IShader& shader, TGAImage& image, float* zbuffer);
void triangle(mat<4, 3, float>& clipc, ZShader& shader, TGAImage& image, float* zbuffer);

#endif //__OUR_GL_H__

