#pragma once

#include <ai.h>
#include "alUtil.h"
#include <cassert>
#include <map>
#include "stats.h"

#define SSS_MAX_SAMPLES 16
#define SSS_MAX_RADIUS 25.0f
#define SSS_ALBEDO_LUT_SZ 256
#define SSS_MAX_PROFILES 9

struct ScatteringProfileDirectional
{
    ScatteringProfileDirectional(){}
    ScatteringProfileDirectional(float Rd, float scale);
    ScatteringProfileDirectional(float sigma_s, float sigma_a, float g);

    const static float eta;
    const static float C_phi;
    const static float C_phi_inv;
    const static float C_E;
    const static float _3C2;
    const static float A;
    const static float _albedo_lut[SSS_ALBEDO_LUT_SZ];

    // numerically calculate Rd from alpha_prime according to the better dipole model.
    // see http://graphics.pixar.com/library/TexturingBetterDipole/paper.pdf
    // {
    float computeRd(float alpha_prime_c)
    {
       float _4A = (1.0f + _3C2) / C_phi;
       float sigma_tr_D = sqrtf((1.0f-alpha_prime_c) * (2.0f-alpha_prime_c) / 3.0f);
       float ex = expf(-_4A * sigma_tr_D);
       return 0.5f * SQR(alpha_prime_c)
                   * expf(-sqrtf(3.0f*(1.0f - alpha_prime_c)/(2.0f-alpha_prime_c))) 
                   * (C_E * (1.0f+ex) + C_phi/sigma_tr_D * (1.0f-ex));
    }

    float computeAlphaPrime(float rd)
    {
       int i, niter = 50;
       float x0 = 0, x1 = 1;
       float xmid, fmid;

       for (i=0; i < niter; ++i)
       {
          xmid = 0.5f * (x0+x1);
          fmid = computeRd(xmid);
          fmid < rd ? x0 = xmid : x1 = xmid;
       }

       return xmid;
    }
    // }

    float de;
    float safe_radius;
    float sigma_t;
    float sigma_t_prime;
    float sigma_tr;
    float D;
    float zr;
    float alpha_prime;
    float albedo;
    float pdf;
};

struct DiffusionSample
{
    AtVector P;     //< sampled point
    AtVector N;     //< normal at sample
    AtVector Ng;    //< geometric normal at sample
    AtRGB R;        //< diffusion result at sample
    AtRGB Rd;       //< convolved diffusion result
    AtVector S;     //< vector from shading point to sample
    float r;        //< distance from shading point to sample
    float b;        //< bounce attenuation factor
    AtRGB lightGroupsDirect[8];
    AtRGB lightGroupsIndirect[8];
};

struct DirectionalMessageData
{
    int sss_depth;
    float maxdist;
    AtVector Po;
    AtVector No;
    // AtVector U;
    // AtVector V;
    AtVector wo;
    ScatteringProfileDirectional sp[9];
    AtRGB weights[9];
    int numComponents;
    bool directional;
    AtShaderGlobals* sg;
    DiffusionSample samples[SSS_MAX_SAMPLES];
    AtRGB* deepGroupPtr;
    AtNode* shader_orig;
};

inline float diffusionSampleDistance(float u1, float sigma)
{
    return -logf(1.0f - std::min(0.999999f,u1)) / sigma;
}

inline float diffusionSampleDisk(float u1, float u2, float sigma, float& dx, float& dy, float& r)
{
    r = -logf(1.0f - std::min(0.999999f,u1)) / sigma;
    float phi = u2 * AI_PITIMES2;
    dx = r * cosf(phi);
    dy = r * sinf(phi);

    // return pdf
    return sigma * AiFastExp(-sigma * r);
}

inline float diffusionPdf(float r, float sigma)
{
    return sigma * AiFastExp(-sigma * r);
}

inline float dipoleProfileRd(float r, float sigma_tr, float zr, float zv)
{
    float dr = sqrtf(r*r + zr*zr);
    float dv = sqrtf(r*r + zv*zv);
    dr = std::max(dr, zr);
    dv = std::max(dv, zv);

    float inv_dr3 = 1.0f / (dr*dr*dr);
    float inv_dv3 = 1.0f / (dv*dv*dv);

    float sigma_tr_dr = sigma_tr * dr;
    float sigma_tr_dv = sigma_tr * dv;

    return zr*(sigma_tr_dr+1.0f) * AiFastExp(-sigma_tr_dr) * inv_dr3 + zv*(sigma_tr_dv+1.0f) * AiFastExp(-sigma_tr_dv) * inv_dv3;
}

inline AtRGB dipoleProfileRd(float r, const AtRGB& sigma_tr, const AtRGB& zr, const AtRGB& zv)
{
    return rgb(
        dipoleProfileRd(r, sigma_tr.r, zr.r, zv.r),
        dipoleProfileRd(r, sigma_tr.g, zr.g, zv.g),
        dipoleProfileRd(r, sigma_tr.b, zr.b, zv.b)
    );
}

#define DD_SINGLE_PRECISION 
#ifdef DD_SINGLE_PRECISION
// Directional dipole profile evaluation
// see: http://www.ci.i.u-tokyo.ac.jp/~hachisuka/dirpole.pdf
// and: http://www.ci.i.u-tokyo.ac.jp/~hachisuka/dirpole.cpp
inline float Sp_d(const AtVector x, const AtVector w, const float r, const AtVector n, const float sigma_tr, const float D, const float Cp_norm, 
                    const float Cp, const float Ce) 
{
    // evaluate the profile
    const float s_tr_r = sigma_tr * r;
    const float s_tr_r_one = 1.0f + s_tr_r;
    const float x_dot_w = AiV3Dot(x, w);
    const float r_sqr = r * r;

    const float t0 = Cp_norm * (1.0f / (4.0f * AI_PI * AI_PI)) * expf(-s_tr_r) / (r * r_sqr);
    const float t1 = r_sqr / D + 3.0f * s_tr_r_one * x_dot_w;
    const float t2 = 3.0f * D * s_tr_r_one * AiV3Dot(w, n);
    const float t3 = (s_tr_r_one + 3.0f * D * (3.0f * s_tr_r_one + s_tr_r * s_tr_r) / r_sqr * x_dot_w) * AiV3Dot(x, n);

    const float Sp = t0 * (Cp * t1 - Ce * (t2 - t3));
    assert(AiIsFinite(Sp));
    return std::max(Sp, 0.0f);
    return Sp;
}
#else
inline double Sp_d(const AtVector x, const AtVector w, const double r, const AtVector n, const double sigma_tr, const double D, const double Cp_norm, 
                    const double Cp, const double Ce) 
{
    // evaluate the profile
    const double s_tr_r = sigma_tr * r;
    const double s_tr_r_one = 1.0 + s_tr_r;
    const double x_dot_w = AiV3Dot(x, w);
    const double r_sqr = r * r;
    // if (r_sqr < 1.e-5) return 1.0;

    const double t0 = Cp_norm * (1.0 / (4.0 * AI_PI * AI_PI)) * exp(-s_tr_r) / (r * r_sqr);
    const double t1 = r_sqr / D + 3.0 * s_tr_r_one * x_dot_w;
    const double t2 = 3.0 * D * s_tr_r_one * AiV3Dot(w, n);
    const double t3 = (s_tr_r_one + 3.0 * D * (3.0 * s_tr_r_one + s_tr_r * s_tr_r) / r_sqr * x_dot_w) * AiV3Dot(x, n);

    const double Sp = t0 * (Cp * t1 - Ce * (t2 - t3));
    assert(AiIsFinite(Sp));
    return std::max(Sp, 0.0);
    return Sp;
}
#endif


inline float directionalDipole(AtPoint xi, AtVector ni, AtPoint xo, AtVector no, AtVector wi, AtVector wo, ScatteringProfileDirectional& sp)
{
    assert(sp.D > 0.0f);
    // distance
    AtVector xoxi = xo - xi;
    float r = AiV3Length(xoxi);
    if (r == 0.0f) return 0.0f;

    // modified normal
    AtVector ni_s = AiV3Cross(AiV3Normalize(xoxi), AiV3Normalize(AiV3Cross(ni, xoxi)));

    // directions of ray sources
    float nnt = 1.0f / sp.eta;
    float ddn = -AiV3Dot(wi, ni);
    AtVector wr = AiV3Normalize(wi * -nnt - ni * (ddn * nnt + sqrtf(1.0f - nnt*nnt * (1.0f - ddn*ddn))));
    AtVector wv = wr - ni_s * (2.0f * AiV3Dot(wr, ni_s));

    // distance to real sources
    const float cos_beta = -sqrtf((r * r - AiV3Dot(xoxi, wr) * AiV3Dot(xoxi, wr)) / (r * r + sp.de * sp.de));
    float dr;
    const float mu0 = -AiV3Dot(ni, wr);
    if (mu0 > 0.0) 
    {
        dr = sqrtf((sp.D * mu0) * ((sp.D * mu0) - sp.de * cos_beta * 2.0) + r * r);
    } 
    else 
    {
        dr = sqrtf(1.0f / (3.0f * sp.sigma_t * 3.0f * sp.sigma_t) + r * r);
    }

    AtVector xoxv = xo - (xi + ni_s * (2.0f * sp.A * sp.de));
    const float dv = AiV3Length(xoxv);

    const float real = Sp_d(xoxi, wr, dr, no, sp.sigma_tr, sp.D, sp.C_phi_inv, sp.C_phi, sp.C_E);
    const float virt = Sp_d(xoxv, wv, dv, no, sp.sigma_tr, sp.D, sp.C_phi_inv, sp.C_phi, sp.C_E);
    const float result = real - virt;
    assert(AiIsFinite(result));
    return std::max(0.0f, result); 
}

#if 0
inline AtRGB integrateDirectional(const ScatteringParamsDirectional& sp, float rmax, int steps)
{
    
	float rstep = 1.0f / float(steps);

    float ns = 0.0f;
    AtPoint Po(0.0f, 0.0f, 0.0f);
    AtVector up(0.0f, 1.0f, 0.0f);
    
	AtRGB result = AI_RGB_BLACK;
	/*
    float sigma = minh(sp.sigma_tr);
    for (float r = rstep/2; r < 1.0f; r += rstep)
    {
        int component;
        float u = drand48();
        if (u < 1.0f/3.0f) component = 0;
        else if (u < 2.03 / 3.0f) component = 1;
        else component = 2;

        float rs = diffusionSampleDistance(r, sp.sigma_tr[component]);
        if (rs > rmax) continue;

        float pdf = (diffusionPdf(rs, sp.sigma_tr[0]) + diffusionPdf(rs, sp.sigma_tr[1]) + diffusionPdf(rs, sp.sigma_tr[2])) / 3.0f;

        AtPoint Pi(rs, 0.0f, 0.0f);
        result += directionalDipole(Pi, up, Po, up, up, up, sp) * rs / pdf;

        ns++;
    }

    result /= ns;
	*/
    return result;
}

#define SSS_INT_HEMI_SAMPLES 30
inline AtRGB integrateDirectionalHemi(const ScatteringParamsDirectional& sp, float rmax, int steps)
{
    float rstep = 1.0f / float(steps);

    float ns = 0.0f;
    AtPoint Po(0.0f, 0.0f, 0.0f);
    AtVector up(0.0f, 1.0f, 0.0f);
    AtRGB result = AI_RGB_BLACK;
	/*
    float sigma = minh(sp.sigma_tr);
    float hemi_step = 1.0f / SSS_INT_HEMI_SAMPLES;
    for (float r = rstep/2; r < 1.0f; r += rstep)
    {
        int component;
        float u = drand48();
        if (u < 1.0f/3.0f) component = 0;
        else if (u < 2.03 / 3.0f) component = 1;
        else component = 2;

        float rs = diffusionSampleDistance(r, sp.sigma_tr[component]);
        if (rs > rmax) continue;

        float pdf = (diffusionPdf(rs, sp.sigma_tr[0]) + diffusionPdf(rs, sp.sigma_tr[1]) + diffusionPdf(rs, sp.sigma_tr[2])) / 3.0f;

        AtPoint Pi(rs, 0.0f, 0.0f);

        for (float u1 = hemi_step/2; u1 < 1.0f; u1 += hemi_step)
        {
            for (float u2 = hemi_step/2; u2 < 1.0f; u2 += hemi_step)
            {
                AtVector wi = cosineSampleHemisphere(u1, u2);

                result += directionalDipole(Pi, up, Po, up, wi, up, sp) * rs / pdf;

                ns++;
            }
        }
    }

    result /= ns;
	*/
    return result;
}
#endif

void alsIrradiateSample(AtShaderGlobals* sg, DirectionalMessageData* dmd, AtSampler* diffuse_sampler, 
                        AtVector U, AtVector V, std::map<const AtNode*, int>& lightGroupMap, AtRGB path_throughput,
                        const AtString& trace_set, bool trace_set_enabled, bool trace_set_inclusive, float sssMix);
AtRGB alsDiffusion(AtShaderGlobals* sg, DirectionalMessageData* dmd, AtSampler* sss_sampler, 
                   bool directional, int numComponents, 
                   AtRGB& result_direct, AtRGB& result_indirect, AtRGB* lightGroupsDirect, AtRGB* deepGroupsSss,
                   AtRGB* deepGroupPtr,
                   const AtString& trace_set, bool trace_set_enabled, bool trace_set_inclusive);
