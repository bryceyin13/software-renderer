#define _USE_MATH_DEFINES
#include <vector>
#include <iostream>
#include <limits>
#include <math.h>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "our_gl.h"

Model* model = NULL;
float* shadowbuffer = NULL;
const int width = 800;
const int height = 800;

Vec3f light_dir(1,1,1);
Vec3f       eye(1,1,3);
Vec3f    center(0,0,0);
Vec3f        up(0,1,0);

struct GouraudShader : public IShader {
    Vec3f varying_intensity; // written by vertex shader, read by fragment shader
    mat<2, 3, float> varying_uv;

    virtual Vec4f vertex(int iface, int nthvert) {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));// get uv of vertex
        Vec4f gl_Vertex = embed<4, 3>(model->vert(iface, nthvert)); // read the vertex from .obj file 
        // embd--> transform points represented in 3D coordinates to homogeneous coordinates
        varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity
        return Viewport * Projection * ModelView * gl_Vertex; // return screen coordinates
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        float intensity = varying_intensity*bar;   // interpolate intensity for the current pixel
        Vec2f uv = varying_uv * bar;
        color = model->diffuse(uv)*intensity; // well duh
        return false;                              // no, we do not discard this pixel
    }
};

struct MyShader :public IShader {
    Matrix Trans_M;
    Matrix Trans_MIT;
    mat<2, 3, float> varying_uv;

    virtual Vec4f vertex(int iface, int nthvert)
    {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        Vec4f gl_Vertex = embed<4, 3>(model->vert(iface, nthvert));
        Trans_M = Projection * ModelView;
        Trans_MIT = Trans_M.invert_transpose();
        return Viewport * Projection * ModelView * gl_Vertex;
    }

    virtual bool fragment(Vec3f barycentric, TGAColor& color)
    {
        Vec2f uv = varying_uv * barycentric;
        Vec3f normal = proj<3, 4>(Trans_MIT * embed<4, 3>(model->normal(uv))).normalize();
        Vec3f light = proj<3, 4>(Trans_M * embed<4, 3>(light_dir)).normalize();
        float intensity = std::max(0.f, normal * light);
        color = model->diffuse(uv);
        return false;
    }
};

struct PhongShader :public IShader {
    Matrix Trans_M;
    Matrix Trans_MIT;
    mat<2, 3, float> varying_uv;

    virtual Vec4f vertex(int iface, int nthvert)
    {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        Vec4f gl_Vertex = embed<4, 3>(model->vert(iface, nthvert));
        Trans_M = Projection * ModelView;
        Trans_MIT = Trans_M.invert_transpose();
        return Viewport * Projection * ModelView * gl_Vertex;
    }

    virtual bool fragment(Vec3f barycentric, TGAColor& color)
    {
        Vec2f uv = varying_uv * barycentric;
        Vec3f normal = proj<3, 4>(Trans_MIT * embed<4, 3>(model->normal(uv))).normalize();
        Vec3f light = proj<3, 4>(Trans_M * embed<4, 3>(light_dir)).normalize();
        Vec3f reflect = normal * (normal * light * 2.f) - light;
        float p = model->specular(uv);
        TGAColor c = model->diffuse(uv);
        float ka = 0.05;
        float kd = 0.9;
        float ks = 0.05;
        float Ia = 10;
        for (int i = 0; i < 3; i++)
        {
            color[i] = std::min(c[i] * (ka * Ia + kd * std::max(0.f, normal * light) + ks * std::pow(std::max(0.f, normal * reflect), p)), 255.f);
        }
        return false;
    }
};

struct DepthShader :public IShader{
    mat<3, 3, float>varying_coord;
    //DepthShader() : varying_coord() {}
    virtual Vec4f vertex(int iface, int nthvert)
    {
        Vec4f gl_Vertex = embed<4, 3>(model->vert(iface, nthvert));
        //compute 3D triangle coordinates
        varying_coord.set_col(nthvert, proj<3, 4>(gl_Vertex / gl_Vertex[3]));
        return Viewport * Projection * ModelView * gl_Vertex;
    }
    virtual bool fragment(Vec3f barycentric, TGAColor& color)
    {
        Vec3f p = varying_coord * barycentric;
        color = TGAColor(255, 255, 255) * (p.z / depth);
        return false;
    }
};

struct ShadowPhongShader :public IShader {
    Matrix Trans_M;
    Matrix Trans_MIT;
    Matrix Trans_Shadow;
    mat<3, 3, float> varying_coord;
    mat<2, 3, float> varying_uv;

    ShadowPhongShader(Matrix M, Matrix MIT, Matrix MS) : Trans_M(M), Trans_MIT(MIT), Trans_Shadow(MS), varying_uv(), varying_coord() {}

    virtual Vec4f vertex(int iface, int nthvert)
    {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        Vec4f gl_Vertex = Viewport * Projection * ModelView * embed<4, 3>(model->vert(iface, nthvert));
        // gl_vertex in screen space, (x, y, depth, homogenerous coord dim 4)
        varying_coord.set_col(nthvert, proj<3, 4>(gl_Vertex / gl_Vertex[3]));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f barycentric, TGAColor& color)
    {
        // compute point coord in shadow buffer and its shadow value
        Vec4f sb_p = Trans_Shadow * embed<4>(varying_coord * barycentric); // corresponding point in the shadow buffer
        sb_p = sb_p / sb_p[3];
        int idx = int(sb_p[0]) + int(sb_p[1]) * width; // index in the shadowbuffer array
        float shadow = .3 + .7 * (shadowbuffer[idx] < sb_p[2] + 43.34); // magic coeff to avoid z-fighting
        // phong shading
        Vec2f uv = varying_uv * barycentric;
        Vec3f normal = proj<3, 4>(Trans_MIT * embed<4, 3>(model->normal(uv))).normalize();
        Vec3f light = proj<3, 4>(Trans_M * embed<4, 3>(light_dir)).normalize();
        Vec3f reflect = (normal * (normal * light * 2.f) - light).normalize();
        Vec3f half = (light + reflect).normalize();
        //Vec3f reflect = (normal *  2.f - light).normalize();
        float p = model->specular(uv);
        TGAColor c = model->diffuse(uv);
        float ka = 0.05;
        float kd = 0.9;
        float ks = 0.05;
        float Ia = 10;
        float ambient = ka * Ia;
        float diffuse = kd * std::max(0.f, normal * light);
        float specular = ks * std::pow(std::max(0.f, normal * half), p);
        for (int i = 0; i < 3; i++)
        {
            color[i] = std::min(c[i] * (ambient + diffuse + specular) * shadow, 255.f);
        }
        return false;
    }
};

struct ZTShader : public ZShader {
    mat<4, 3, float> varying_tri;

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = Projection * ModelView * embed<4>(model->vert(iface, nthvert));
        varying_tri.set_col(nthvert, gl_Vertex);
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor& color) {
        color = TGAColor(0, 0, 0);
        return false;
    }
};

float max_elevation_angle(float* zbuffer, Vec2f p, Vec2f dir) {
    float maxangle = 0;
    for (float t = 0.; t < 30.; t += 1.) {
        Vec2f cur = p + dir * t;
        if (cur.x >= width || cur.y >= height || cur.x < 0 || cur.y < 0) return maxangle;

        float distance = (p - cur).norm();
        if (distance < 1.f) continue;
        float elevation = zbuffer[int(cur.x) + int(cur.y) * width] - zbuffer[int(p.x) + int(p.y) * width];
        maxangle = std::max(maxangle, atanf(elevation / distance));
    }
    return maxangle;
}

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    float* zbuffer = new float[width * height];
    shadowbuffer = new float[width * height];
    for (int i = width * height; --i; ) {
        zbuffer[i] = shadowbuffer[i] = -std::numeric_limits<float>::max();
    }
    {
        // first time rendering: generate shadow map
        TGAImage depth(width, height, TGAImage::RGB);
        lookat(light_dir, center, up);
        viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
        projection(0);
        DepthShader depthshader;
        for (int i = 0; i < model->nfaces(); i++)
        {
            Vec4f screen_coords[3];
            for (int j = 0; j < 3; j++)
            {
                screen_coords[j] = depthshader.vertex(i, j);
            }
            triangle(screen_coords, depthshader, depth, shadowbuffer);
        }
    }
    {
        // second time rendering: generate image
        Matrix M = Viewport * Projection * ModelView; //transform from World Coord to Shadow Map
        lookat(eye, center, up);
        viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
        projection(-1.f / (eye - center).norm());
        light_dir.normalize();

        TGAImage image(width, height, TGAImage::RGB);
        //TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

        GouraudShader shader;
        MyShader myshader;
        PhongShader phongshader;
        ShadowPhongShader shadowphongshader(Projection * ModelView, (Projection * ModelView).invert_transpose(), M * (Viewport * Projection * ModelView).invert());

        for (int i = 0; i < model->nfaces(); i++) {
            Vec4f screen_coords[3];
            for (int j = 0; j < 3; j++) {
                screen_coords[j] = shadowphongshader.vertex(i, j);
            }
            triangle(screen_coords, shadowphongshader, image, zbuffer);
        }

        //image.flip_vertically(); // to place the origin in the bottom left corner of the image
        //image.write_tga_file("output.tga");
        //zbuffer.write_tga_file("zbuffer.tga");
        
        // generate SSAO and add it to image
        TGAImage frame(width, height, TGAImage::RGB);
        lookat(eye, center, up);
        viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
        projection(-1.f / (eye - center).norm());

        ZTShader zshader;
        for (int i = 0; i < model->nfaces(); i++) {
            for (int j = 0; j < 3; j++) {
                zshader.vertex(i, j);
            }
            triangle(zshader.varying_tri, zshader, frame, zbuffer);
        }

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (zbuffer[x + y * width] < -1e5) continue;
                float total = 0;
                for (float a = 0; a < M_PI * 2 - 1e-4; a += M_PI / 4) {
                    total += M_PI / 2 - max_elevation_angle(zbuffer, Vec2f(x, y), Vec2f(cos(a), sin(a)));
                }
                total /= (M_PI / 2) * 8;
                //std::cout << total << " ";
                //total = pow(total, 100);
                TGAColor color = image.get(x, y);
                TGAColor c;
                /*for (int i = 0; i < 3; i++)
                {
                    c[i] = std::min(color[i] + total * 255, 255.f);
                    c[i] = color[i] * total;
                }*/
                if (total > 0.6)c = color;
                else c = color * 0.3;
                frame.set(x, y, c);
            }
        }
        frame.flip_vertically();
        frame.write_tga_file("SSAO2.tga");
    }
    delete model;
    delete[] zbuffer;
    delete[] shadowbuffer;
    return 0;
}
