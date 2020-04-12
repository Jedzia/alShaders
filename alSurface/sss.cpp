#include "sss.h"
#include "aovs.h"

const float ScatteringProfileDirectional::eta(1.3f);
const float ScatteringProfileDirectional::C_phi(0.175626f);
const float ScatteringProfileDirectional::C_phi_inv(1.03668f);
const float ScatteringProfileDirectional::C_E(0.27735f);
const float ScatteringProfileDirectional::_3C2(0.0611156f);
const float ScatteringProfileDirectional::A(2.05736f);
const float ScatteringProfileDirectional::_albedo_lut[SSS_ALBEDO_LUT_SZ] = {0.004660, 0.005024, 0.005285, 0.005503, 0.005697, 0.005875, 0.006041, 0.006197, 0.006345, 0.006488, 0.006625, 0.006758, 0.006887, 0.007013, 0.007136, 0.007256, 0.007374, 0.007490, 0.007604, 0.007716, 0.007826, 0.007935, 0.008043, 0.008149, 0.008255, 0.008358, 0.008461, 0.008564, 0.008665, 0.008765, 0.008865, 0.008964, 0.009063, 0.009160, 0.009257, 0.009354, 0.009450, 0.009546, 0.009641, 0.009735, 0.009829, 0.009923, 0.010016, 0.010110, 0.010203, 0.010295, 0.010387, 0.010479, 0.010571, 0.010662, 0.010754, 0.010845, 0.010935, 0.011026, 0.011116, 0.011207, 0.011297, 0.011387, 0.011477, 0.011566, 0.011656, 0.011746, 0.011835, 0.011925, 0.012014, 0.012103, 0.012192, 0.012281, 0.012370, 0.012459, 0.012548, 0.012637, 0.012726, 0.012815, 0.012904, 0.012993, 0.013082, 0.013171, 0.013260, 0.013349, 0.013438, 0.013527, 0.013616, 0.013705, 0.013794, 0.013883, 0.013972, 0.014061, 0.014150, 0.014240, 0.014329, 0.014418, 0.014508, 0.014597, 0.014687, 0.014777, 0.014866, 0.014956, 0.015046, 0.015136, 0.015226, 0.015317, 0.015407, 0.015497, 0.015588, 0.015678, 0.015769, 0.015860, 0.015951, 0.016042, 0.016133, 0.016224, 0.016315, 0.016407, 0.016498, 0.016590, 0.016682, 0.016774, 0.016866, 0.016958, 0.017051, 0.017143, 0.017236, 0.017329, 0.017421, 0.017514, 0.017608, 0.017701, 0.017794, 0.017888, 0.017982, 0.018076, 0.018170, 0.018264, 0.018359, 0.018453, 0.018548, 0.018643, 0.018738, 0.018833, 0.018928, 0.019024, 0.019120, 0.019215, 0.019311, 0.019408, 0.019504, 0.019601, 0.019698, 0.019795, 0.019892, 0.019990, 0.020087, 0.020185, 0.020283, 0.020381, 0.020479, 0.020577, 0.020676, 0.020774, 0.020873, 0.020973, 0.021072, 0.021172, 0.021272, 0.021372, 0.021472, 0.021573, 0.021673, 0.021774, 0.021875, 0.021976, 0.022078, 0.022180, 0.022282, 0.022383, 0.022486, 0.022588, 0.022690, 0.022794, 0.022897, 0.023001, 0.023104, 0.023208, 0.023311, 0.023417, 0.023521, 0.023626, 0.023730, 0.023837, 0.023942, 0.024047, 0.024153, 0.024260, 0.024366, 0.024472, 0.024579, 0.024687, 0.024794, 0.024901, 0.025010, 0.025118, 0.025226, 0.025335, 0.025444, 0.025554, 0.025663, 0.025773, 0.025883, 0.025992, 0.026104, 0.026214, 0.026326, 0.026436, 0.026548, 0.026661, 0.026772, 0.026886, 0.026997, 0.027111, 0.027225, 0.027338, 0.027452, 0.027567, 0.027682, 0.027795, 0.027911, 0.028027, 0.028143, 0.028259, 0.028376, 0.028491, 0.028609, 0.028726, 0.028844, 0.028962, 0.029082, 0.029201, 0.029320, 0.029439, 0.029559, 0.029678, 0.029800, 0.029921, 0.030041, 0.030164, 0.030286, 0.030409, 0.030531, 0.030655, 0.030778, 0.030903, 0.031028, 0.031151, 0.031277, 0.031403};

ScatteringProfileDirectional::ScatteringProfileDirectional(float Rd, float scale)
{
    // get temp reduced scattering coefficients from Rd
    const float ap = computeAlphaPrime(Rd * 0.439f);
    const float str = 1.0f / Rd;
    const float stp = str / sqrtf(3.0f * (1.0f - ap));

    // calculate actual scattering coefficients from the temps
    float sigma_s_prime = stp * ap;
    float sigma_a = stp - sigma_s_prime;

    // factor of 2.5 is eyeballed to roughly match the look of the cubic
    sigma_s_prime *= scale * AI_PI * 2;
    sigma_a *= scale * AI_PI * 2;

    const float sigma_s = sigma_s_prime;
    // safe_radius = 0.5f / sigma_s;
    safe_radius = 0.0f;

    sigma_t_prime = sigma_s_prime + sigma_a;
    sigma_t = sigma_s + sigma_a;

    alpha_prime = sigma_s_prime / sigma_t_prime;

    D = (2*sigma_t_prime) / (3*SQR(sigma_t_prime));
    sigma_tr = sqrt(sigma_a / D);
    de = 2.131 * D / sqrt(alpha_prime);
    zr = 1.0f / sigma_t_prime;

    assert(AiIsFinite(D));
    assert(D > 0.0f);
    assert(AiIsFinite(de));
    assert(AiIsFinite(sigma_tr));
    assert(AiIsFinite(sigma_t_prime));
    assert(AiIsFinite(alpha_prime));
    assert(AiIsFinite(zr));

    const float maxdist = zr * SSS_MAX_RADIUS;

    // grab the precalculated albedo for this Rd
    const int idx = int(Rd * (SSS_ALBEDO_LUT_SZ-1));
    albedo = _albedo_lut[idx];
    assert(albedo > AI_EPSILON);
    assert(albedo < 1.0f - AI_EPSILON);
    // printf("Rd: %f - albedo: %f\n", Rd, albedo);
}

ScatteringProfileDirectional::ScatteringProfileDirectional(float sigma_s, float sigma_a, float g)
{
    float sigma_s_prime = sigma_s * (1.0f - g);
    sigma_t_prime = sigma_s_prime + sigma_a;
    sigma_t = sigma_s + sigma_a;

    safe_radius = 0.5f / sigma_s;

    alpha_prime = sigma_s_prime / sigma_t_prime;

    D = (2*sigma_t_prime) / (3*SQR(sigma_t_prime));
    sigma_tr = sqrt(sigma_a / D);
    de = 2.131 * D / sqrt(alpha_prime);
    zr = 1.0f / sigma_t_prime;

    assert(AiIsFinite(D));
    assert(AiIsFinite(de));
    assert(AiIsFinite(sigma_tr));
    assert(AiIsFinite(sigma_t_prime));
    assert(AiIsFinite(alpha_prime));
    assert(AiIsFinite(zr));

    const float maxdist = zr * SSS_MAX_RADIUS;
    albedo = 1.0f;
    // printf("FUCK\n");
}

void matchNormals(const AtVector Nref, AtVector& Nmatch)
{
    if (AiV3Dot(Nref, Nmatch) < 0.0f)
    {
        Nmatch = -Nmatch;
    }
}

void alsIrradiateSample(AtShaderGlobals* sg, DirectionalMessageData* dmd, AtSampler* diffuse_sampler, 
                        AtVector U, AtVector V, std::map<AtNode*, int>& lightGroupMap, AtRGB path_throughput,
                        const char* trace_set, bool trace_set_enabled, bool trace_set_inclusive, float sssMix)
{
    if (dmd->sss_depth >= SSS_MAX_SAMPLES) return;
    
    void* orig_op;
    AiStateGetMsgPtr("als_sss_op", &orig_op);
    
    int als_context = ALS_CONTEXT_NONE;
    AiStateGetMsgInt("als_context", &als_context);

    AtRay ray;
    AtScrSample scrs;


    // There are a few contexts in which we can be called here:
    // 1) Intersecting same object with same shader
    // 2) Intersecting same object with different shader
    // 3) Intersecting different object with same shader
    // 4) Intersecting different object with different shader
    // 5) Called from an indirect ray from another shader while tracing for SSS
    // 2 and 4 could possibly be because we're running from a layer shader. if we try to evaluate in that context we will simply crash
    // 5 will cause fireflies so we need to detect this and return

    // just return if the shader that called us was not alSurface
    if (sg->psg && AiNodeGetNodeEntry(sg->shader) != AiNodeGetNodeEntry(sg->psg->shader))
    {
        return;
    }

    // if we're in a layered context and running a different shader from the one that's tracing for sss, just return to let the other shader handle it
    if (als_context == ALS_CONTEXT_LAYER && (dmd->shader_orig != sg->shader))
    {
        return;
    }

    DiffusionSample& samp = dmd->samples[dmd->sss_depth];
    samp.S = sg->P - dmd->Po;
    samp.r = AiV3Length(samp.S);
    dmd->maxdist -= sg->Rl;
    samp.Rd = AI_RGB_BLACK;

    // if we're not using trace sets to explicitly define what objects we want to trace against, then assume we only want to trace against ourselves
    // or if we're not doing sss in this shader, just continue the ray
    if ((!trace_set_enabled && orig_op != sg->Op) || sssMix == 0.0f)
    {
        if (dmd->sss_depth < SSS_MAX_SAMPLES && dmd->maxdist > 0.0f)
        {
            AiMakeRay(&ray, AI_RAY_SUBSURFACE, &sg->P, &sg->Rd, dmd->maxdist, sg);
            AiTrace(&ray, &scrs);
        }
        return;
    }

    AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);
    
    // assert(sg->P != dmd->Po);
    // put normals the right way round
    AtVector Nref = sg->N;
    matchNormals(Nref, sg->Nf);
    matchNormals(Nref, sg->Ng);
    matchNormals(Nref, sg->Ngf);
    matchNormals(Nref, sg->Ns);

    //void* brdf_data = AiOrenNayarMISCreateData(sg, 0.0f);
    AiLightsPrepare(sg);
    AtRGB result_direct = AI_RGB_BLACK;
    uint32_t old_fi = sg->fi;

    AtRGB Rnond = AI_RGB_BLACK;
    bool directional = dmd->directional;
    
    if (!directional)
    {
        for (int c = 0; c < dmd->numComponents; ++c)
        {
            Rnond += directionalDipole(sg->P, sg->N, dmd->Po, dmd->No, sg->N, dmd->No, dmd->sp[c]) * dmd->weights[c];
        }
    }

    assert(samp.lightGroupsDirect[0] == AI_RGB_BLACK);
    while (AiLightsGetSample(sg))
    {
        // get the group assigned to this light from the hash table using the light's pointer
        int lightGroup = lightGroupMap[sg->Lp];
        float diffuse_strength = AiLightGetDiffuse(sg->Lp);

        // Does not using MIS here get us anything?
        // AtRGB L = AiEvaluateLightSample(sg, brdf_data, AiOrenNayarMISSample, AiOrenNayarMISBRDF, AiOrenNayarMISPDF);
        AtRGB L = sg->Li * MAX(AiV3Dot(sg->Ld, sg->N), 0.0f) * AI_ONEOVERPI * sg->we * diffuse_strength;
        if (AiColorIsZero(L)) continue;
        if (directional)
        {
            AtRGB R = AI_RGB_BLACK;
            for (int c=0; c < dmd->numComponents; ++c)
            {
                if (samp.r < dmd->sp[c].safe_radius)
                    R += directionalDipole(sg->P, sg->N, dmd->Po, dmd->No, sg->N, dmd->No, dmd->sp[c]) * dmd->weights[c];
                else
                    R += directionalDipole(sg->P, sg->N, dmd->Po, dmd->No, sg->Ld, dmd->wo, dmd->sp[c]) * dmd->weights[c]; 
            }
            L *= R;
            result_direct += L;
            if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
            {
                samp.lightGroupsDirect[lightGroup] += L;
            }
            assert(AiIsFinite(result_direct));
        }
        else
        {
            L *= Rnond;
            result_direct += L;
            if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
            {
                samp.lightGroupsDirect[lightGroup] += L;
            }
            assert(AiIsFinite(result_direct));
        }
    }

    AtRGB result_indirect = AI_RGB_BLACK;
    AtSamplerIterator* sampit = AiSamplerIterator(diffuse_sampler, sg);
    float samples[2];
    
    AiMakeRay(&ray, AI_RAY_DIFFUSE, &sg->P, &sg->N, AI_BIG, sg);
    while (AiSamplerGetSample(sampit, samples))
    {
        float stheta = sqrtf(float(samples[0]));
        float phi = float(AI_PITIMES2 * samples[1]);
        ray.dir.x = stheta * cosf(phi);
        ray.dir.y = stheta * sinf(phi);
        ray.dir.z = sqrtf(1.0f - float(samples[0]));
        AiV3RotateToFrame(ray.dir, U, V, sg->N);

        bool hit = AiTrace(&ray, &scrs);

        AtRGB f;
        if (directional)
        {
            AtRGB R = AI_RGB_BLACK;
            for (int c = 0; c < dmd->numComponents; ++c)
            {
                R += directionalDipole(sg->P, sg->N, dmd->Po, dmd->No, ray.dir, dmd->wo, dmd->sp[c]) * dmd->weights[c];
            }
            f = R;
            result_indirect += scrs.color * R;;
            assert(AiIsFinite(result_indirect));
        }
        else
        {
            f = Rnond;
            result_indirect += scrs.color * Rnond;
            assert(AiIsFinite(result_indirect));
        }

        if (hit && !AiColorIsZero(f) && dmd->deepGroupPtr)
        {
            for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
            {
                samp.lightGroupsIndirect[i] += dmd->deepGroupPtr[i] * f;
            }
        }
    }
    

    // TODO: this is guaranteed to be 1 in every case, right?
    // result_indirect *= AiSamplerGetSampleInvCount(sampit);
    samp.Rd = result_direct + result_indirect;

    assert(AiIsFinite(samp.Rd));

   
    samp.N = sg->N;
    samp.Ng = sg->Ng;
    samp.P = sg->P;

    dmd->sss_depth++;

    AiStateSetMsgInt("als_raytype", ALS_RAY_SSS);
    
    if (dmd->sss_depth < SSS_MAX_SAMPLES && dmd->maxdist > 0.0f)
    {
        
        if (trace_set_enabled)
        {
            AiShaderGlobalsSetTraceSet(sg, trace_set, trace_set_inclusive);
        }
        AiMakeRay(&ray, AI_RAY_SUBSURFACE, &sg->P, &sg->Rd, dmd->maxdist, sg);
        AiTrace(&ray, &scrs);
    }
}

AtRGB alsDiffusion(AtShaderGlobals* sg, DirectionalMessageData* dmd, AtSampler* sss_sampler, 
                   bool directional, int numComponents,
                   AtRGB& result_direct, AtRGB& result_indirect, AtRGB* lightGroupsDirect, AtRGB* lightGroupsIndirect,
                   AtRGB* deepGroupPtr,
                   const char* trace_set, bool trace_set_enabled, bool trace_set_inclusive)

{
    AtVector U, V;
    AiBuildLocalFrameShirley(&U, &V, &sg->Ng);

    numComponents = std::min(numComponents, SSS_MAX_PROFILES);
    float l = 0.0f;
    float inv_pdf_sum = 0.0f;
    float comp_pdf[9];
    float comp_cdf[9+1];
    comp_cdf[0] = 0.0f;
    int last_nonzero = numComponents;
    for (int i=0; i < numComponents; ++i)
    {
        // dmd->sp[i] = ScatteringProfileDirectional(Rd[i], sssDensityScale/radii[i]);
        float w = maxh(dmd->weights[i]);
        float pdf = dmd->sp[i].alpha_prime;
        comp_pdf[i] = pdf * w;
        comp_cdf[i+1] = comp_cdf[i] + comp_pdf[i];
        inv_pdf_sum += comp_pdf[i];

        if (w > 0.0f)
        {
            // track the last non-zero weight we encounter so that we can ignore completely missing lobes
            last_nonzero = i+1;  
            // track the largest mean free path so we can set that to be our maximum raytracing distance
            l = std::max(l, dmd->sp[i].zr);
        } 
    }

    // set the number of components to be the number of non-zero-weight components
    numComponents = std::min(numComponents, last_nonzero);
    // set the maximum raytracing distance to be some multiple of the largest mean-free path
    // the choice of SSS_MAX_RADIUS is a quality/speed tradeoff. The default of 25 seems to work well for most cases
    const float R_max = l * SSS_MAX_RADIUS;

    // normalize the PDF and CDF
    inv_pdf_sum = 1.0f / inv_pdf_sum;
    for (int i=0; i < numComponents; ++i)
    {
        comp_pdf[i] *= inv_pdf_sum;
        comp_cdf[i+1] *= inv_pdf_sum;
    }

    // trick Arnold into thinking we're shooting from a different face than we actually are so he doesn't ignore intersections
    uint32_t old_fi = sg->fi;
    sg->fi = UINT_MAX;

    AtRGB result_sss = AI_RGB_BLACK;
    
    AtRGB Rd_sum = AI_RGB_BLACK;
    int samplesTaken = 0;
    float samples[2];
    AtRay wi_ray;
    AtScrSample scrs;
    AtSamplerIterator* sampit = AiSamplerIterator(sss_sampler, sg);
    dmd->wo = -sg->Rd;
    dmd->numComponents = numComponents;
    dmd->directional = directional;
    dmd->deepGroupPtr = deepGroupPtr;
    dmd->shader_orig = sg->shader;
    AtVector axes[3] = {U, V, sg->Ng};
    double sss_samples_c  = 0;
    if (trace_set_enabled)
    {
        AiShaderGlobalsSetTraceSet(sg, trace_set, trace_set_inclusive);
    }
    while (AiSamplerGetSample(sampit, samples))
    {
        float dx = 0.0f, dy = 0.0f;

        AtVector Wsss, Usss, Vsss, Usss_1, Vsss_1, Usss_2, Vsss_2;
        float c_axis = 1.0f, c_axis_1, c_axis_2;
        float axes_c[3] = {0.25f, 0.25f, 0.5f};
        int chosen_axis;
        if (samples[0] < 0.5f)
        {
            chosen_axis = 2;
            samples[0] *= 2.0f;
            c_axis = 0.5f;
            c_axis_1 = 0.25f;
            c_axis_2 = 0.25f;
            Wsss = sg->Ng;
            
            Usss = U;
            Vsss = V;

            Usss_1 = sg->Ng;
            Vsss_1 = V;

            Usss_2 = U;
            Vsss_2 = sg->Ng;
        }
        else if (samples[0] < 0.75f)
        {
            chosen_axis = 0;
            samples[0] = (samples[0] - 0.5f) * 4.0f;
            c_axis = 0.25f;
            c_axis_1 = 0.5f;
            c_axis_2 = 0.25f;
            Wsss = U;
            
            Usss = sg->Ng;
            Vsss = V;

            Usss_1 = U;
            Vsss_1 = V;

            Usss_2 = U;
            Vsss_2 = sg->Ng;
        }
        else
        {
            chosen_axis = 1;
            samples[0] = (1.0f-samples[0])* 4.0f;
            c_axis = 0.25f;
            c_axis_1 = 0.25f;
            c_axis_2 = 0.5f;
            Wsss = V;
            
            Usss = U;
            Vsss = sg->Ng;

            Usss_1 = sg->Ng;
            Vsss_1 = V;

            Usss_2 = U;
            Vsss_2 = V;
        }

        float r_disk[3] = {0.0f, 0.0f, 0.0f};
        for (int i=0; i < numComponents; ++i)
        {
            if (samples[1] < comp_cdf[i+1])
            {
                samples[1] -= comp_cdf[i];
                samples[1] /= comp_pdf[i];
                diffusionSampleDisk(samples[0], samples[1], dmd->sp[i].sigma_tr, dx, dy, r_disk[0]);
                break;
            }
        }

        // if the sampled distance is greater than we'll consider anyway, don't bother tracing at all.
        if (r_disk[0] > R_max) continue;

        AtVector dir = -Wsss;
        float dz = R_max;
        AtPoint origin = sg->P + Wsss*(dz*.25) + Usss * dx + Vsss * dy;
        float maxdist = R_max ;//* 2.0f;

        AiMakeRay(&wi_ray, AI_RAY_SUBSURFACE, &origin, &dir, maxdist, sg);

        AiStateSetMsgInt("als_raytype", ALS_RAY_SSS);
        AiStateSetMsgPtr("als_sss_op", sg->Op);
        dmd->sss_depth = 0;
        dmd->maxdist = R_max;
        dmd->Po = sg->P;
        dmd->No = sg->N;
        float geom[3];
        
        memset(dmd->samples, 0, sizeof(DiffusionSample)*SSS_MAX_SAMPLES);
        dmd->sss_depth = 0;
        if (AiTrace(&wi_ray, &scrs))
        {            
            for (int i=0; i < dmd->sss_depth; ++i)
            {                
                if (AiColorIsZero(dmd->samples[i].Rd)) continue;

                geom[0] = fabsf(AiV3Dot(dmd->samples[i].Ng, Wsss));
                geom[1] = fabsf(AiV3Dot(dmd->samples[i].Ng, Usss));
                geom[2] = fabsf(AiV3Dot(dmd->samples[i].Ng, Vsss));

                float r_u_1 = AiV3Dot(dmd->samples[i].S, Usss_1);
                float r_v_1 = AiV3Dot(dmd->samples[i].S, Vsss_1);
                r_disk[1] = sqrtf(SQR(r_u_1)+SQR(r_v_1));

                float r_u_2 = AiV3Dot(dmd->samples[i].S, Usss_2);
                float r_v_2 = AiV3Dot(dmd->samples[i].S, Vsss_2);
                r_disk[2] = sqrtf(SQR(r_u_2)+SQR(r_v_2));

                float pdf_sum = 0.0f;

                for (int c=0; c < numComponents; ++c)
                {
                    pdf_sum += diffusionPdf(r_disk[0], dmd->sp[c].sigma_tr) * comp_pdf[c] * geom[0] * c_axis;
                    pdf_sum += diffusionPdf(r_disk[1], dmd->sp[c].sigma_tr) * comp_pdf[c] * geom[1] * c_axis_1; 
                    pdf_sum += diffusionPdf(r_disk[2], dmd->sp[c].sigma_tr) * comp_pdf[c] * geom[2] * c_axis_2;   
                }

                if (pdf_sum == 0.0f) continue;

                float f = r_disk[0] / pdf_sum;
                result_sss += dmd->samples[i].Rd * f;

                for (int g=0; g < NUM_LIGHT_GROUPS; ++g)
                {
                    lightGroupsDirect[g] += dmd->samples[i].lightGroupsDirect[g] * f;
                    assert(AiIsFinite(lightGroupsDirect[g]));
                    lightGroupsIndirect[g] += dmd->samples[i].lightGroupsIndirect[g] * f;
                    assert(AiIsFinite(lightGroupsIndirect[g]));
                }
                assert(AiIsFinite(result_sss));
            }
            
        }
    }
    assert(AiIsFinite(result_sss));
    float w = AiSamplerGetSampleInvCount(sampit);
    result_sss *= w;
    assert(AiIsFinite(result_sss));
    AtRGB norm_factor = AI_RGB_BLACK;
    for (int c=0; c < numComponents; ++c)
    {
        norm_factor += dmd->sp[c].albedo * dmd->weights[c];
    }

    result_sss /= norm_factor;
    for (int g=0; g < NUM_LIGHT_GROUPS; ++g)
    {
        lightGroupsDirect[g] *= rgb(w) / norm_factor;
        lightGroupsIndirect[g] *= rgb(w) / norm_factor;
    }

    assert(AiIsFinite(result_sss));
    sg->fi = old_fi;

    if (trace_set_enabled)
    {
        AiShaderGlobalsUnsetTraceSet(sg);
    }
    // Optimization hack: do a regular indirect diffuse and colour it like the subsurface instead of allowing sss rays to
    // continue for another diffuse bounce.
    /*
    if (!sssDoIndirect)
    {
        AtRGB result_sss_indirect = AiIndirectDiffuse(&sg->N, sg); 
        if (!data->sss_normalize) result_sss_indirect *= sp.albedo;
        result_sss += result_sss_indirect;
    }
    */
    AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);

    return result_sss;
}
