/// @file alGaborNoise.cpp
/// Improved Gabor Noise based on:
/// [1] Improving Gabor Noise - Lagae et al. - 2011
/// [2] OSL https://github.com/imageworks/OpenShadingLanguage/blob/master/src/liboslexec/gabornoise.cpp

#define _USE_MATH_DEFINES

#include "Remap.h"
#include <ai.h>
#include <cmath>
#include <cstring>

AI_SHADER_NODE_EXPORT_METHODS(alGaborNoise);

static const float GABOR_TRUNCATE = 0.02f;
static const float GABOR_FREQUENCY = 2.0f;

enum GaborAnisotropy
{
    GP_ISOTROPIC,
    GP_ANISOTROPIC,
    GP_HYBRID
};

static const char* gaborAnisotropyNames[] = 
{
    "isotropic",
    "anisotropic",
    "hybrid",
    NULL
};

enum NoiseSpaceEnum
{
    NS_WORLD = 0,
    NS_OBJECT,
    NS_PREF,
    NS_UV
};

static const char* noiseSpaceNames[] =
{
        "world",
        "object",
        "Pref",
        "UV",
        NULL
};

struct GaborParams
{
    AtVector omega;
    int anisotropic;
    bool do_filter;
    float a;
    float weight;
    AtVector N;
    float filter[2][2];
    float local[3][3];
    float det_filter;
    float bandwidth;
    bool periodic;
    AtVector period;
    float lambda;
    float sqrt_lambda_inv;
    float radius, radius2, radius3, radius_inv;

    GaborParams(AtVector aniso_direction, int aniso, bool df, float bw, float imp) :
        omega(aniso_direction),
        anisotropic(aniso),
        do_filter(df),
        weight(1.0f),
        bandwidth(clamp(bw, 0.01f, 100.0f)),
        periodic(false)
    {
        float two_to_bw = powf(2.0f, bandwidth);
        static const float SQRT_PI_OVER_LN2 = sqrtf(AI_PI / 0.693147180559945309417);
        a = GABOR_FREQUENCY * ((two_to_bw - 1.0f) / (two_to_bw + 1.0f)) * SQRT_PI_OVER_LN2;

        // calculate maximum extent
        radius = sqrtf(-logf(GABOR_TRUNCATE) / AI_PI) / a;
        radius2 = SQR(radius);
        radius3 = radius * radius2;
        radius_inv = 1.0f / radius;
        // lambda is the impulse density
        float impulses = clamp(imp, 1.0f, 32.0f);
        lambda = impulses / (1.33333333f * AI_PI * radius3);
        sqrt_lambda_inv = 1.0f / sqrtf(lambda);
    }
};

/// gabor kernel. cosine weighted by a gaussian, augmented with a phase as per [1]
/// @param weight magnitude of the pulse
/// @param omega orientation of the harmonic
/// @param phi phase of the harmonic
/// @param bandwidth width of the gaussian envelope
/// @param x the position being sampled
inline float gaborKernel(float weight, const AtVector& omega, float phi, float bandwidth, const AtVector& x)
{
    float g = fast_exp(-AI_PI * SQR(bandwidth) * AiV3Dot(x, x));
    float h = cosf(AI_PITIMES2 * AiV3Dot(omega, x) + phi);
    return weight * g * h;
}

inline void sliceGaborKernel3d(float d, float w, float a, const AtVector& omega, float phi, float& w_s, AtVector& omega_s, float& phi_s)
{
    // [1] eq. (6)
    w_s = w * fast_exp(-AI_PI * SQR(a) * SQR(d));
    omega_s.x = omega.x;
    omega_s.y = omega.y;
    phi_s = phi - AI_PITIMES2 * d * omega.z;
}

/// choose an omega and phi value for a particular gabor impulse based on the user parameters
void gaborSample(GaborParams& gp, const AtVector& x_c, LCG& rng, AtVector& omega, float& phi)
{
    if (gp.anisotropic == GP_ANISOTROPIC)
    {
        omega = gp.omega;
    }
    else if (gp.anisotropic == GP_ISOTROPIC)
    {
        float omega_t = AI_PITIMES2 * rng();
        float cos_omega_p = lerp(-1.0f, 1.0f, rng());
        float sin_omega_p = sqrtf(std::max(0.0f, 1.0f - SQR(cos_omega_p)));
        float sin_omega_t, cos_omega_t;
        sincosf_(omega_t, &sin_omega_t, &cos_omega_t);
        AiV3Create(omega, cos_omega_t*sin_omega_p, sin_omega_t*sin_omega_p, cos_omega_p);
        AiV3Normalize(omega);
    }
    else // hybrid
    {
        float omega_r = AiV3Length(gp.omega);
        float omega_t = AI_PITIMES2 * rng();
        float sin_omega_t, cos_omega_t;
        sincosf_(omega_t, &sin_omega_t, &cos_omega_t);
        AtVector omega_tt;
        AiV3Create(omega_tt, cos_omega_t, sin_omega_t, 0.0f);
        omega = omega_r * omega_tt;
    }
    phi = AI_PITIMES2 * rng();
}

/// Evaluate the summed contribution of all gabor impulses in the cell whose corner is c_i. x_c_i is the vector from the point we are
/// evaluating to the cell corner
float gaborCell(GaborParams& gp, const AtVector& c_i, const AtVector& x_c_i, int seed=0)
{
    //LCG rng(gp.periodic ? wrap(c_i, gp.period) : c_i, seed);
    LCG rng(c_i, seed);
    int n_impules = rng.poisson(gp.lambda * gp.radius3);
    float sum = 0;
    for (int i=0; i < n_impules; ++i)
    {
        float z_rng = rng(); float y_rng = rng(); float x_rng = rng();
        AtVector x_i_c; AiV3Create(x_i_c, x_rng, y_rng, z_rng);
        AtVector x_k_i = gp.radius * (x_c_i - x_i_c);
        float phi_i;
        AtVector omega_i;
        gaborSample(gp, c_i, rng, omega_i, phi_i);
        if (AiV3Dot(x_k_i, x_k_i) < gp.radius2)
        {
            // just do unfiltered for now...
            sum += gaborKernel(gp.weight, omega_i, phi_i, gp.a, x_k_i);
        }
    }

    return sum;
}

float gaborGrid(GaborParams& gp, const AtVector& x_g, int seed=0)
{
    AtVector floor_x_g = floor(x_g);
    AtVector x_c = x_g - floor_x_g;
    float sum = 0;

    for (int k = -1; k <= 1; ++k)
    {
        for (int j = -1; j <= 1; ++j)
        {
            for (int i = -1; i <= 1; ++i)
            {
                AtVector c;
                AiV3Create(c, i, j, k);
                AtVector c_i = floor_x_g + c;
                AtVector x_c_i = x_c - c;

                sum += gaborCell(gp, c_i, x_c_i, seed);
            }
        }
    }

    return sum * gp.sqrt_lambda_inv;
}

inline float gaborEvaluate(GaborParams& gp, const AtVector& x, int seed=0)
{
    AtVector x_g = x * gp.radius_inv;
    return gaborGrid(gp, x_g, seed);
}

inline float gabor(const AtVector& P, GaborParams& gp)
{
    if (gp.do_filter)
    {
        // gaborSetupFilter(...);
    }

    /// TODO: if we're happy with the parameters being constant, the scale calculation can be moved into the param constructor
    float result = gaborEvaluate(gp, P);
    float gaborVariance = 1.0f / (4.0f * AI_SQRT2 * gp.a * gp.a * gp.a);
    float scale = 1.0f / (3.0f * sqrtf(gaborVariance));
    scale *= 0.5f;

    return result * scale;
}

enum alGaborNoiseParams
{
    p_space,
    p_frequency,
    p_anisotropy,
    p_anisotropyDirection,
    p_filter,
    p_bandwidth,
    p_impulses,
    p_turbulent,
    REMAP_FLOAT_PARAM_ENUM,

    p_color1,
    p_color2,

    p_P,

    p_end_of_params    
};

node_parameters
{
    AiParameterEnum("space", NS_WORLD, noiseSpaceNames);
    AiParameterFlt("frequency", 1.0f);
    AiParameterEnum("anisotropy", GP_ISOTROPIC, gaborAnisotropyNames);
    AiParameterVec("anisotropyDirection", 0.0f, 1.0f, 0.0f);
    AiParameterBool("filter", false);
    AiParameterFlt("bandwidth", 1.0f);
    AiParameterFlt("impulses", 8.0f);
    AiParameterBool("turbulent", false);
    REMAP_FLOAT_PARAM_DECLARE;
    AiParameterRGB("color1", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("color2", 1.0f, 1.0f, 1.0f);

    AiParameterPnt("P", 0.0f, 0.0f, 0.0f);
}

node_loader
{
    if (i>0) return 0;
    node->methods     = alGaborNoise;
    node->output_type = AI_TYPE_RGB;
    node->name        = "alGaborNoise";
    node->node_type   = AI_NODE_SHADER;
    ::strcpy(node->version, AI_VERSION);
    return true;
}

struct ShaderData
{
    int space;
    int anisotropy;
    bool filter;
    bool turbulent;
    int impulses;
};

node_initialize
{
    ShaderData* data = new ShaderData;
    AiNodeSetLocalData(node, data);
}

node_finish
{
    ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
    delete data;
}

node_update
{
    ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
    data->anisotropy = params[p_anisotropy].INT;
    data->filter = params[p_filter].BOOL;
    data->impulses = params[p_impulses].INT;
    data->space = params[p_space].INT;
    data->turbulent = params[p_turbulent].BOOL;

}

shader_evaluate
{
    ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
    AtRGB color1 = AiShaderEvalParamRGB(p_color1);
    AtRGB color2 = AiShaderEvalParamRGB(p_color2);
    float frequency = AiShaderEvalParamFlt(p_frequency);
    GaborParams gp(AiShaderEvalParamVec(p_anisotropyDirection), data->anisotropy, 
                                            data->filter, AiShaderEvalParamFlt(p_bandwidth), data->impulses);

    // choose what space we want to calculate in
    AtPoint P;
    if (AiNodeIsLinked(node, "P"))
    {
        AtPoint Pin = AiShaderEvalParamPnt(p_P);
        P = Pin;
    }
    else
    {
        switch (data->space)
        {
        case NS_OBJECT:
            P = sg->Po;
            break;
        case NS_UV:
            P.x = sg->u;
            P.y = sg->v;
            P.z = 0.0f;
            break;
        case NS_PREF:
            if (!AiUDataGetPnt("Pref", &P))
                P = sg->Po;
            break;
        default:
            P = sg->P;
            break;
        }
    }

    P *= frequency;
    
    float result = gabor(P, gp);
    if (data->turbulent) result = fabsf(result);

    RemapFloat r = REMAP_FLOAT_CREATE;
    result = r.remap(result);

    sg->out.RGB() = lerp(color1, color2, result);
}
