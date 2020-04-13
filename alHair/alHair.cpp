
// Hair shader based on 
// [1] ISHair: Importance Sampling for Hair Scattering by Ou et al. 2012
// http://www.cs.dartmouth.edu/~ouj/site/Research/Entries/2012/6/21_ISHair__Importance_Sampling_for_Hair_Scattering.html
// [2] Dual Scattering Approximation For Fast Multiple Scattering in Hair by Zinke et al. 2008


#include <ai.h>
#include "alUtil.h"
#include "stats.h"
#include "scattering.h"
#include <vector>
#include <algorithm>
#include <map>
#include <cassert>
#include "cryptomatte/cryptomatte.h"

//#define DEBUG_LUTS
#ifdef DEBUG_LUTS
#include "exr.h"
#endif

#define NUM_LIGHT_GROUPS 16

AI_SHADER_NODE_EXPORT_METHODS(alHair);

enum alHairParams
{
    p_dyeColor,
    p_melanin,
    p_specularShift,
    p_specularWidth,
    p_extraSamplesDiffuse,
    p_extraSamplesGlossy,
    p_diffuseStrength,
    p_diffuseColor,
    p_diffuseScatteringMode,
    p_specular1Strength,
    p_specular1Color,
    p_specular2Strength,
    p_specular2Color,
    p_glintStrength,
    p_glintRolloff,
    p_transmissionStrength,
    p_transmissionColor,
    p_transmissionRolloff,
    p_opacity,
    p_dualDepth,
    p_densityFront,
    p_densityBack,
    p_specular1WidthScale,
    p_specular2WidthScale,
    p_transmissionWidthScale,
    p_specular1Shift,
    p_specular2Shift,
    p_transmissionShift,
    p_diffuseIndirectStrength,
    p_glossyIndirectStrength,

    p_randomTangent,
    p_randomHue,
    p_randomSaturation,
    p_randomMelanin,

    p_uparam,
    p_vparam,

    p_id1,
    p_id2,
    p_id3,
    p_id4,
    p_id5,
    p_id6,
    p_id7,
    p_id8,

    p_aov_diffuse_color,
    p_aov_direct_diffuse,
    p_aov_indirect_diffuse,
    p_aov_direct_local,
    p_aov_indirect_local,
    p_aov_direct_global,
    p_aov_indirect_global,
    p_aov_direct_specular,
    p_aov_indirect_specular,
    p_aov_direct_specular_2,
    p_aov_indirect_specular_2,
    p_aov_direct_glint,
    p_aov_indirect_glint,
    p_aov_direct_transmission,
    p_aov_indirect_transmission,
    p_aov_depth,
    p_aov_light_group_1,
    p_aov_light_group_2,
    p_aov_light_group_3,
    p_aov_light_group_4,
    p_aov_light_group_5,
    p_aov_light_group_6,
    p_aov_light_group_7,
    p_aov_light_group_8,
    p_aov_light_group_9,
    p_aov_light_group_10,
    p_aov_light_group_11,
    p_aov_light_group_12,
    p_aov_light_group_13,
    p_aov_light_group_14,
    p_aov_light_group_15,
    p_aov_light_group_16,
    p_aov_id_1,
    p_aov_id_2,
    p_aov_id_3,
    p_aov_id_4,
    p_aov_id_5,
    p_aov_id_6,
    p_aov_id_7,
    p_aov_id_8,
    p_crypto_asset_override,
    p_crypto_object_override,
    p_crypto_material_override,
    p_aov_crypto_asset,
    p_aov_crypto_object,
    p_aov_crypto_material
};

enum AovIndices
{
    k_diffuse_color=0,
    k_direct_diffuse,
    k_indirect_diffuse,
    k_direct_local,
    k_indirect_local,
    k_direct_global,
    k_indirect_global,
    k_direct_specular,
    k_indirect_specular,
    k_direct_specular_2,
    k_indirect_specular_2,
    k_direct_glint,
    k_indirect_glint,
    k_direct_transmission,
    k_indirect_transmission,
    k_depth,
    k_light_group_1,
    k_light_group_2,
    k_light_group_3,
    k_light_group_4,
    k_light_group_5,
    k_light_group_6,
    k_light_group_7,
    k_light_group_8,
    k_light_group_9,
    k_light_group_10,
    k_light_group_11,
    k_light_group_12,
    k_light_group_13,
    k_light_group_14,
    k_light_group_15,
    k_light_group_16,
    k_id_1,
    k_id_2,
    k_id_3,
    k_id_4,
    k_id_5,
    k_id_6,
    k_id_7,
    k_id_8
};

#define NUM_ID_AOVS 8
static const char* id_names[NUM_ID_AOVS] =
{
    "id_1",
    "id_2",
    "id_3",
    "id_4",
    "id_5",
    "id_6",
    "id_7",
    "id_8",
};

enum ScatteringMode
{
    KK = 0,
    DUAL
};

const char* ScatteringModeNames[] =
{
    "kajiya-kay",
    "dual-scattering",
    NULL
};

node_parameters
{
    AiParameterRGB("dyeColor", 1.0f, 1.0f, 1.0f);
    AiParameterFlt("melanin", 0.15f);
    AiParameterFlt("specularShift", 4.0f);
    AiParameterFlt("specularWidth", 3.0f);
    AiParameterInt("extraSamplesDiffuse", 0);
    AiParameterInt("extraSamplesGlossy", 0);
    AiParameterFlt("diffuseStrength", 1.0f);
    AiParameterRGB("diffuseColor", 1.0f, 1.0f, 1.0f);
    AiParameterEnum("diffuseScatteringMode", DUAL, ScatteringModeNames);
    AiParameterFlt("specular1Strength", 1.0f);
    AiParameterRGB("specular1Color", 1.0f, 1.0f, 1.0f);
    AiParameterFlt("specular2Strength", 1.0f);
    AiParameterRGB("specular2Color",1.0f, 1.0f, 1.0f);
    AiParameterFlt("glintStrength", 2.5f);
    AiParameterFlt("glintRolloff", 5.0f);
    AiParameterFlt("transmissionStrength", 1.0f);
    AiParameterRGB("transmissionColor", 1.0f, 1.0f, 1.0f);
    AiParameterFlt("transmissionRolloff", 10.0f);
    AiParameterRGB("opacity", 1.0f, 1.0f, 1.0f);
    AiParameterInt("dualDepth", 0);
    AiParameterFlt("diffuseForward", 0.7f);
    AiParameterFlt("diffuseBack", 0.7f);
    AiParameterFlt("specular1WidthScale", 1.0f);
    AiParameterFlt("specular2WidthScale", 1.0f);
    AiParameterFlt("transmissionWidthScale", 1.0f);
    AiParameterFlt("specular1Shift", 0.0f);
    AiParameterFlt("specular2Shift", 0.0f);
    AiParameterFlt("transmissionShift", 0.0f);
    AiParameterFlt("diffuseIndirectStrength", 0.0f);
    AiParameterFlt("glossyIndirectStrength", 0.0f);

    AiParameterFlt("randomTangent", 0.0f);
    AiParameterFlt("randomHue", 0.0f);
    AiParameterFlt("randomSaturation", 0.0f);
    AiParameterFlt("randomMelanin", 0.08f);

    AiParameterStr("uparam", "uparamcoord");
    AiParameterStr("vparam", "vparamcoord");

    AiParameterRGB("id1", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id2", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id3", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id4", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id5", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id6", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id7", 0.0f, 0.0f, 0.0f);
    AiParameterRGB("id8", 0.0f, 0.0f, 0.0f);

    AiParameterStr("aov_diffuse_color", "diffuse_color");
    AiParameterStr("aov_direct_diffuse", "direct_diffuse");
    AiParameterStr("aov_indirect_diffuse", "indirect_diffuse");
    AiParameterStr("aov_direct_local", "direct_local");
    AiParameterStr("aov_indirect_local", "indirect_local");
    AiParameterStr("aov_direct_global", "direct_global");
    AiParameterStr("aov_indirect_global", "indirect_global");
    AiParameterStr("aov_direct_specular", "direct_specular");
    AiParameterStr("aov_indirect_specular", "indirect_specular");
    AiParameterStr("aov_direct_specular_2", "direct_specular_2");
    AiParameterStr("aov_indirect_specular_2", "indirect_specular_2");
    AiParameterStr("aov_direct_glint", "direct_glint");
    AiParameterStr("aov_indirect_glint", "indirect_glint");
    AiParameterStr("aov_direct_transmission", "direct_transmission");
    AiParameterStr("aov_indirect_transmission", "indirect_transmission");
    AiParameterStr("aov_depth", "depth");
    AiParameterStr("aov_light_group_1", "light_group_1");
    AiParameterStr("aov_light_group_2", "light_group_2");
    AiParameterStr("aov_light_group_3", "light_group_3");
    AiParameterStr("aov_light_group_4", "light_group_4");
    AiParameterStr("aov_light_group_5", "light_group_5");
    AiParameterStr("aov_light_group_6", "light_group_6");
    AiParameterStr("aov_light_group_7", "light_group_7");
    AiParameterStr("aov_light_group_8", "light_group_8");
    AiParameterStr("aov_light_group_9", "light_group_9");
    AiParameterStr("aov_light_group_10", "light_group_10");
    AiParameterStr("aov_light_group_11", "light_group_11");
    AiParameterStr("aov_light_group_12", "light_group_12");
    AiParameterStr("aov_light_group_13", "light_group_13");
    AiParameterStr("aov_light_group_14", "light_group_14");
    AiParameterStr("aov_light_group_15", "light_group_15");
    AiParameterStr("aov_light_group_16", "light_group_16");
    AiParameterStr("aov_id_1", "id_1");
    AiParameterStr("aov_id_2", "id_2");
    AiParameterStr("aov_id_3", "id_3");
    AiParameterStr("aov_id_4", "id_4");
    AiParameterStr("aov_id_5", "id_5");
    AiParameterStr("aov_id_6", "id_6");
    AiParameterStr("aov_id_7", "id_7");
    AiParameterStr("aov_id_8", "id_8");
    AiParameterStr("crypto_asset_override", "");
    AiParameterStr("crypto_object_override", "");
    AiParameterStr("crypto_material_override", "");
    AiParameterStr("aov_crypto_asset", "crypto_asset");
    AiParameterStr("aov_crypto_object", "crypto_object");
    AiParameterStr("aov_crypto_material", "crypto_material");
}

#define B_WIDTH_SCALE 2.0f
#define F_WIDTH_SCALE 4.0f

struct LutgenThreadData
{
   DualScattering* ds;
   ScatteringParams* sp;
   int start;
   int end;
};

unsigned int LutgenThreadFunc(void* data)
{
   LutgenThreadData* td = (LutgenThreadData*)data;
   ScatteringParams& sp = *td->sp;

   for (int idx=td->start; idx < td->end; ++idx)
   {
      float a = float(idx) / float(DS_MASTER_LUT_SZ);
      float sigma = -logf(std::max(a, 0.0001f));
      td->ds->_luts[idx] = new ScatteringLut(sp.ior, sp.alpha_R, sp.alpha_TT, sp.alpha_TRT, sp.beta_R2, sp.beta_TT2,
                                               sp.beta_TRT2, sp.gamma_TT, sp.gamma_g, sp.phi_g, sigma, sp.shape);
   }

   return 0;
}

struct HairBsdf
{
    struct ShaderData
    {
        ShaderData()
        : sampler_glossy(NULL), sampler_diffuse(NULL), ds(NULL)
        {
            _last_sp.shape = -1.0f;
            _last_sp.ior = 0.0f;
        }

        ~ShaderData()
        {
            AiSamplerDestroy(sampler_glossy);
            AiSamplerDestroy(sampler_diffuse);
            delete ds;
        }

        void update(/*AtParamValue* params*/)
        {
            AtNode *options   = AiUniverseGetOptions();
            int glossy_samples = std::max(0, AiNodeGetInt(options, "GI_glossy_samples") +AiNodeGetInt(options, "extraSamplesGlossy"));
            int diffuse_samples = std::max(0, AiNodeGetInt(options, "GI_diffuse_samples") +AiNodeGetInt(options, "extraSamplesDiffuse"));

            static const uint32_t seed =
                    static_cast<uint32_t>(AiNodeEntryGetNameAtString(AiNodeGetNodeEntry(options)).hash());
            AiSamplerDestroy(sampler_glossy);
            sampler_glossy = AiSampler(seed, glossy_samples, 2);
            AiSamplerDestroy(sampler_diffuse);
            sampler_diffuse = AiSampler(seed, diffuse_samples, 2);

            alpha_R = -AiNodeGetFlt(options, "specularShift") * AI_DTOR;
            beta_R =AiNodeGetFlt(options, "specularWidth") * AI_DTOR;

            beta_TT = beta_R * 0.5f;

            alpha_TT = -alpha_R * 0.5f;

            beta_TRT = beta_R * 2.0f;
            alpha_TRT = -alpha_R * 1.5f;

            beta_R2 = beta_R*beta_R;
            beta_TT2 = beta_TT*beta_TT;
            beta_TRT2 = beta_TRT*beta_TRT;

            gamma_TT =AiNodeGetFlt(options, "transmissionRolloff") * AI_DTOR;
            gamma_g =AiNodeGetFlt(options, "glintRolloff") * AI_DTOR;

            glintStrength =AiNodeGetFlt(options, "glintStrength");

            dual_depth =AiNodeGetInt(options, "dualDepth");

            sampleLobesIndividually = false;

            uparam =AiNodeGetStr(options, "uparam");
            vparam =AiNodeGetStr(options, "vparam");

            aovs.clear();
            aovs.push_back(AiNodeGetStr(options, "aov_diffuse_color"));
            aovs.push_back(AiNodeGetStr(options, "aov_direct_diffuse"));
            aovs.push_back(AiNodeGetStr(options, "aov_indirect_diffuse"));
            aovs.push_back(AiNodeGetStr(options, "aov_direct_local"));
            aovs.push_back(AiNodeGetStr(options, "aov_indirect_local"));
            aovs.push_back(AiNodeGetStr(options, "aov_direct_global"));
            aovs.push_back(AiNodeGetStr(options, "aov_indirect_global"));
            aovs.push_back(AiNodeGetStr(options, "aov_direct_specular"));
            aovs.push_back(AiNodeGetStr(options, "aov_indirect_specular"));
            aovs.push_back(AiNodeGetStr(options, "aov_direct_specular_2"));
            aovs.push_back(AiNodeGetStr(options, "aov_indirect_specular_2"));
            aovs.push_back(AiNodeGetStr(options, "aov_direct_glint"));
            aovs.push_back(AiNodeGetStr(options, "aov_indirect_glint"));
            aovs.push_back(AiNodeGetStr(options, "aov_direct_transmission"));
            aovs.push_back(AiNodeGetStr(options, "aov_indirect_transmission"));
            aovs.push_back(AiNodeGetStr(options, "aov_depth"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_1"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_2"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_3"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_4"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_5"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_6"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_7"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_8"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_9"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_10"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_11"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_12"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_13"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_14"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_15"));
            aovs.push_back(AiNodeGetStr(options, "aov_light_group_16"));
            aovs.push_back(AiNodeGetStr(options, "aov_id_1"));
            aovs.push_back(AiNodeGetStr(options, "aov_id_2"));
            aovs.push_back(AiNodeGetStr(options, "aov_id_3"));
            aovs.push_back(AiNodeGetStr(options, "aov_id_4"));
            aovs.push_back(AiNodeGetStr(options, "aov_id_5"));
            aovs.push_back(AiNodeGetStr(options, "aov_id_6"));
            aovs.push_back(AiNodeGetStr(options, "aov_id_7"));
            aovs.push_back(AiNodeGetStr(options, "aov_id_8"));
            for (size_t i=0; i < aovs.size(); ++i)
                AiAOVRegister(aovs[i].c_str(), AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);

            // Only update dual scattering LUTs if we're actually using dual scattering
            if (AiNodeGetInt(options, "diffuseScatteringMode") == DUAL)
            {
                delete ds;
                ds = new DualScattering;

                ScatteringParams sp;
                sp.alpha_R = alpha_R;
                sp.beta_R = beta_R;
                sp.beta_R2 = beta_R2;
                sp.alpha_TT = alpha_TT;
                sp.beta_TT = beta_TT;
                sp.beta_TT2 = beta_TT2;
                sp.alpha_TRT = alpha_TRT;
                sp.beta_TRT = beta_TRT;
                sp.beta_TRT2 = beta_TRT2;
                sp.gamma_TT = gamma_TT;
                sp.gamma_g = gamma_g;
                sp.ior = 1.55f;
                sp.phi_g = 35.0f;
                sp.shape = 0.0f;
#if 0
                if (sp != _last_sp)
                {
                    int num_threads = AiNodeGetInt(options, "threads");
                    LutgenThreadData* thread_data = new LutgenThreadData[num_threads];
                    void** threads = new void*[num_threads];

                    int chunk_size = (DS_MASTER_LUT_SZ / num_threads)+1;
                    int start = 0;
                    float time = AiMsgUtilGetElapsedTime();
                    for (int i=0; i < num_threads; ++i)
                    {
                    thread_data[i].ds = ds;
                    thread_data[i].sp = &sp;
                    thread_data[i].start = start;
                    thread_data[i].end = std::min(start + chunk_size, DS_MASTER_LUT_SZ);
                    start += chunk_size;
                    threads[i] = AiThreadCreate(LutgenThreadFunc, &thread_data[i], AI_PRIORITY_HIGH);
                    }

                    for (int i=0; i < num_threads; ++i)
                    {
                    AiThreadWait(threads[i]);
                    }

                    delete[] thread_data;

                    time = AiMsgUtilGetElapsedTime() - time;
                    AiMsgInfo("[alHair] lutgen time on %d threads: %.2fs", num_threads, time/1000.0f);
                    _last_sp = sp;
                }
#endif
            }
        }

        AtSampler* sampler_glossy;
        AtSampler* sampler_diffuse;

        int dual_depth;

        float beta_R;       //< R width
        float beta_R2;
        float alpha_R;      //< R shift
        float beta_TT;      //< TT width
        float beta_TT2;
        float alpha_TT;     //< TT shift
        float beta_TRT;     //< TRT width
        float beta_TRT2;
        float alpha_TRT;    //< TRT shift
        float gamma_TT;     //< TT rolloff
        float gamma_g;      //< g rolloff
        float phi_g;        //< g separation

        float glintStrength;

        bool sampleLobesIndividually;

        DualScattering* ds;

        std::vector<AtString> aovs;
        std::map<const AtNode*, int> lightGroups;

        AtString uparam, vparam;

        CryptomatteData* cryptomatte;
        ScatteringParams _last_sp;
    };

    HairBsdf(AtNode* n, AtShaderGlobals* sg, ShaderData* d) :
    numBlendHairs(2), density_front(0.7f), density_back(0.7f), node(n), data(d), _sg(sg)
    {
        depth = sg->bounces;

        result_R_direct = AI_RGB_BLACK;
        result_R_indirect = AI_RGB_BLACK;
        result_TT_direct = AI_RGB_BLACK;
        result_TT_indirect = AI_RGB_BLACK;
        result_TRT_direct = AI_RGB_BLACK;
        result_TRT_indirect = AI_RGB_BLACK;
        result_TRTg_direct = AI_RGB_BLACK;
        result_TRTg_indirect = AI_RGB_BLACK;
        result_Pg_direct = AI_RGB_BLACK;
        result_Pl_direct = AI_RGB_BLACK;
        result_Pg_indirect = AI_RGB_BLACK;
        result_Pl_indirect = AI_RGB_BLACK;
        result_id1 = AI_RGB_BLACK;
        result_id2 = AI_RGB_BLACK;
        result_id3 = AI_RGB_BLACK;
        result_id4 = AI_RGB_BLACK;
        result_id5 = AI_RGB_BLACK;
        result_id6 = AI_RGB_BLACK;
        result_id7 = AI_RGB_BLACK;
        result_id8 = AI_RGB_BLACK;


        memset(lightGroupsDirect, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
    }

    /// Parameter evaluation. This should be called after opacity() and before anything else.
    inline void evaluateParameters(AtShaderGlobals* sg, ShaderData* data)
    {
        // save uv coord 
        float u = sg->u;
        float v = sg->v;

        // replace uvs
        float s, t;
        static const AtString maya_ramp_uv_override("maya_ramp_uv_override");
        if (AiUDataGetFlt(data->uparam, s) && AiUDataGetFlt(data->vparam, t))
        {
            // set the hair u, v coords so other shaders (e.g. ramp) can use them
        AtPoint2 hair_uv(sg->u, sg->v);
#if AI_VERSION_MINOR_NUM >= 6
        AiStateSetMsgVec2(maya_ramp_uv_override, hair_uv);
#else
        AiStateSetMsgVec2(AtString("maya_ramp_uv_override"), hair_uv);
#endif

        sg->u = s;
        sg->v = t;

        }


        // Get a random value per curve
        uint32_t curve_idi = 0;
        int curve_idii = 0;
        float curve_id = 0.0f;
        cn = 1.0f;
        AtVector cv = aivec(0.0f);
        float randomTangent = AiShaderEvalParamFlt(p_randomTangent) * 0.3f;
        if (AiUDataGetUInt(AtString("curve_id"), curve_idi))
        {
            AtPoint2 p; p.x = float(curve_idi); p.y = 0.0f;
            cn = AiCellNoise2(p);

            AtPoint p2 = aivec(float(curve_idi)+17.0f, 0.0f, 0.0f);
            cv = (AiVCellNoise3(p2)*2.0f - aivec(1.0f));
        }
        else if (AiUDataGetFlt(AtString("curve_id"), curve_id))
        {
            AtPoint2 p; p.x = curve_id; p.y = 0.0f;
            cn = AiCellNoise2(p);

            AtPoint p2 = aivec(curve_id+17.0f, 0.0f, 0.0f);
            cv = (AiVCellNoise3(p2)*2.0f - aivec(1.0f));

        }
        else if (AiUDataGetInt(AtString("curve_id"), curve_idii))
        {
            AtPoint2 p; p.x = float(curve_idii); p.y = 0.0f;
            cn = AiCellNoise2(p);

            AtPoint p2 = aivec(float(curve_idii)+17.0f, 0.0f, 0.0f);
            cv = (AiVCellNoise3(p2)*2.0f - aivec(1.0f));
        }
        else
        {
            cn = 0.5f;
        }

        // Get a local coordinate frame based on the hair fibre direction
        U = AiV3Normalize(sg->dPdv + cv*randomTangent);
        V = AiV3Cross(U, sg->N);
        W = AiV3Cross(V, U);

        // Get the spherical angles of the exitant direction relative to the hair fibre
        wo = -sg->Rd;
        theta_r = AI_PIOVER2 - sphericalTheta(wo, U);
        phi_r = sphericalPhi(wo, V, W);

        sp.beta_R = data->beta_R;
        sp.alpha_R = data->alpha_R;
        sp.alpha_R_offset = -AiShaderEvalParamFlt(p_specular1Shift) * AI_DTOR;

        sp.beta_TT = data->beta_TT;
        sp.alpha_TT = data->alpha_TT;
        sp.alpha_TT_offset = -AiShaderEvalParamFlt(p_transmissionShift) * AI_DTOR;

        sp.beta_TRT = data->beta_TRT;
        sp.alpha_TRT = data->alpha_TRT;
        sp.alpha_TRT_offset = -AiShaderEvalParamFlt(p_specular2Shift) * AI_DTOR;

        sp.beta_R2 = data->beta_R2;
        sp.beta_TT2 = data->beta_TT2;
        sp.beta_TRT2 = data->beta_TRT2;

        sp.gamma_TT = data->gamma_TT;
        sp.gamma_g = data->gamma_g;

        float glintAngleMin = 0.0f * AI_DTOR;
        float glintAngleMax = 360.0f * AI_DTOR;

        sp.phi_g = lerp(glintAngleMin, glintAngleMax, cn);

        sp.ior = 1.55;

        AB(theta_r, sp.alpha_R, sp.beta_R, A_R, B_R);
        AB(theta_r, sp.alpha_TT, sp.beta_TT, A_TT, B_TT);
        AB(theta_r, sp.alpha_TRT, sp.beta_TRT, A_TRT, B_TRT);
        AB(theta_r, sp.alpha_R, sp.beta_R*B_WIDTH_SCALE, A_b, B_b);
        AB(theta_r, sp.alpha_TT, sp.beta_TT*F_WIDTH_SCALE, A_f, B_f);

        AtRGB dyeColor = clamp(AiShaderEvalParamRGB(p_dyeColor), AI_RGB_BLACK, AI_RGB_WHITE);
        float melanin = AiShaderEvalParamFlt(p_melanin);
        opacity = AiShaderEvalParamFlt(p_opacity);

        float randomHue = AiShaderEvalParamFlt(p_randomHue) * 0.1f;
        float randomSaturation = AiShaderEvalParamFlt(p_randomSaturation);
        float randomMelanin = AiShaderEvalParamFlt(p_randomMelanin);
        melanin = AiClamp(melanin+randomMelanin*cv.z, 0.0f, 1.0f);

        if (randomHue != 0.0f || randomSaturation != 0.0f)
        {
            dyeColor = rgb2hsv(dyeColor);
            dyeColor.r += cv.x * randomHue;
            // dyeColor.r = clamp(-1.0f, 1.0f, dyeColor.r);
            dyeColor.g += cv.y * randomSaturation;
            dyeColor.g = clamp(0.0f, 1.0f, dyeColor.g);
            dyeColor = hsv2rgb(dyeColor);
        }

        float m = AiMax(powf(melanin, 2.0f)*33.0f, 1.0e-2f);
        hairColor = exp(m * -rgb(0.187f, 0.4f, 1.05f));
        hairColor *= dyeColor;

        sp.hairColor = hairColor;

        sp.absorption = -log(hairColor);
        sp.dabsorption = -log(hairColor);

        diffuseColor = AiShaderEvalParamRGB(p_diffuseColor) * AiShaderEvalParamFlt(p_diffuseStrength);
        diffuseScatteringMode = AiShaderEvalParamInt(p_diffuseScatteringMode);
        specular1Color = AiShaderEvalParamRGB(p_specular1Color) * AiShaderEvalParamFlt(p_specular1Strength);
        specular2Color = AiShaderEvalParamRGB(p_specular2Color) * AiShaderEvalParamFlt(p_specular2Strength);
        transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) * AiShaderEvalParamFlt(p_transmissionStrength);
        glintStrength = AiShaderEvalParamFlt(p_glintStrength);

        specular1WidthScale = AiShaderEvalParamFlt(p_specular1WidthScale);
        specular2WidthScale = AiShaderEvalParamFlt(p_specular2WidthScale);
        transmissionWidthScale = AiShaderEvalParamFlt(p_transmissionWidthScale);

        diffuseIndirectStrength = AiShaderEvalParamFlt(p_diffuseIndirectStrength);
        glossyIndirectStrength = AiShaderEvalParamFlt(p_glossyIndirectStrength);

        density_front = AiShaderEvalParamFlt(p_densityFront);
        density_back = AiShaderEvalParamFlt(p_densityBack);

        if (density_back*maxh(diffuseColor) > IMPORTANCE_EPS && density_front*maxh(diffuseColor) > IMPORTANCE_EPS)
            dualImportant = true;
        else
            dualImportant = false;

        do_diffuse = do_glossy = true;
        if (sg->bounces_diffuse > 0)
        {
            do_glossy = false;
        }

        hairColor *= diffuseColor;

        doLightGroups = true;
        if (sg->Rt & AI_RAY_CAMERA)
        {
            // if this is a camera ray allocate the group storage
            deepGroupPtr = (AtRGB*)AiShaderGlobalsQuickAlloc(sg, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
            memset(deepGroupPtr, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
            AiStateSetMsgPtr(AtString("als_deepGroupPtr"), deepGroupPtr);
        }
        else
        {
            // secondary ray hit - get the pointer from the state
            // if the pointer hasn't been set we're being called from a BSDF that doesn't have deep group support
            // so don't try and do it or we'll be in (crashy) trouble!
            if (!AiStateGetMsgPtr(AtString("als_deepGroupPtr"), (void**)&deepGroupPtr)) deepGroupPtr = NULL;
        }
        if (deepGroupPtr)
        {
            memset(deepGroupsSpecular, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
            memset(deepGroupsSpecular2, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
            memset(deepGroupsDiffuse, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
            memset(deepGroupsTransmission, 0, sizeof(AtRGB)*NUM_LIGHT_GROUPS);
        }

#if AI_VERSION_MINOR_NUM >= 8
        AiStateUnsetMsgPnt2(maya_ramp_uv_override);
#endif
        // restore saved uv coords for hygiene
        sg->u = u;
        sg->v = v;
    }

    inline void AB(float theta_r, float alpha, float beta, float& A, float& B)
    {
        A = atanf((PIOVER4 + theta_r*0.5f - alpha) / beta);
        B = atanf((-PIOVER4 + theta_r*0.5f - alpha) / beta);
    }

    inline float sampleLong(double u, float theta_r, float alpha, float beta, float A, float B)
    {
        float t = 2.0f*beta * tanf(u*(A-B) + B) + 2.0f*alpha - theta_r;
        return t;//clamp( -0.4999f * AI_PI, 0.4999f * AI_PI, t); //< TODO: not handling the t == +/- AI_PIOVER2 case doesn't seem to hurt?
    }

    inline AtVector sample_R(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_R - sp.alpha_R_offset, sp.beta_R*specular1WidthScale, A_R, B_R);
        float phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f));
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_R(const SctGeo& geo)
    {
        AtRGB result = rgb(bsdfR(sp.beta_R2*SQR(specular1WidthScale), sp.alpha_R - sp.alpha_R_offset, geo.theta_h, geo.cosphi2)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVERPI;
        assert(isValidColor(result));
        return result;
    }

    inline float pdf_R(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_R - sp.alpha_R_offset;
        float pdf_theta = (1.0f / (2.0f*std::max(geo.cos_theta_i, 1e-7f)*(A_R-B_R))) * (sp.beta_R*specular1WidthScale / (t*t + sp.beta_R2*SQR(specular1WidthScale)));
        float pdf_phi = geo.cosphi2*0.25f;
        float pdf = pdf_theta * pdf_phi;
        assert(isPositiveReal(pdf));
        return pdf;
    }

    inline AtVector sample_TRT(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_TRT - sp.alpha_TRT_offset, sp.beta_TRT*specular2WidthScale, A_R, B_R);
        float phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f));
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_TRT(const SctGeo& geo)
    {
        AtRGB result = rgb(bsdfR(sp.beta_TRT2*SQR(specular2WidthScale), sp.alpha_TRT - sp.alpha_TRT_offset, geo.theta_h, geo.cosphi2)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVERPI;
        assert(isValidColor(result));
        return result;
    }

    inline float pdf_TRT(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_TRT - sp.alpha_TRT_offset;
        float pdf_theta = (1.0f / (2.0f*std::max(geo.cos_theta_i, 1e-7f)*(A_R-B_R))) * (sp.beta_TRT*specular2WidthScale / (t*t + sp.beta_TRT2*SQR(specular2WidthScale)));
        float pdf_phi = geo.cosphi2*0.25f;
        float pdf = pdf_theta * pdf_phi;
        assert(isPositiveReal(pdf));
        return pdf;
    }

    inline AtVector sample_TRTg(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_TRT - sp.alpha_TRT_offset, sp.beta_TRT*specular2WidthScale, A_TRT, B_TRT);
        float sign;
        if (u2 < 0.5f)
        {
            sign = 1.0f;
            u2 = 2.0f * u2;
        }
        else
        {
            sign = -1.0f;
            u2 = 2.0f * (1.0f-u2);
        }
        Cg = atanf((AI_PIOVER2 - sp.phi_g)/sp.gamma_g);
        Dg = atanf(-sp.phi_g/sp.gamma_g);
        float phi = sp.gamma_g * tanf(u2*(Cg-Dg)+Dg) + sp.phi_g;
        phi *= sign;
        float phi_i = phi_r + phi;

        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_TRTg(const SctGeo& geo)
    {
        AtRGB result = rgb(bsdfg(sp.beta_TRT2*SQR(specular2WidthScale), sp.alpha_TRT - sp.alpha_TRT_offset, geo.theta_h, sp.gamma_g, geo.phi, sp.phi_g)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVERPI;
        assert(isValidColor(result));
        return result;
    }

    inline float pdf_TRTg(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_TRT - sp.alpha_TRT_offset;
        float pdf_theta = (1.0f / (2.0f*std::max(geo.cos_theta_i, 1e-7f)*(A_TRT-B_TRT))) * (sp.beta_TRT*specular2WidthScale / (t*t + sp.beta_TRT2*SQR(specular2WidthScale)));
        float p = fabsf(geo.phi) - sp.phi_g;
        float Cg = atanf((AI_PIOVER2 - sp.phi_g)/sp.gamma_g);
        float Dg = atanf(-sp.phi_g/sp.gamma_g);
        float pdf_phi = (1.0f / (2.0f * (Cg-Dg))) * (sp.gamma_g / (p*p + sp.gamma_g*sp.gamma_g));
        float pdf = pdf_theta * pdf_phi;
        assert(isPositiveReal(pdf));
        return pdf;
    }

    inline AtVector sample_TT(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_TT - sp.alpha_TT_offset, sp.beta_TT*transmissionWidthScale, A_TT, B_TT);
        C_TT = 2.0f * atanf(AI_PI/sp.gamma_TT);
        float phi = sp.gamma_TT * tanf(C_TT * (u2-0.5f)) + AI_PI;
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline AtRGB bsdf_TT(const SctGeo& geo)
    {
        AtRGB result = rgb(bsdfTT(sp.beta_TT2*SQR(transmissionWidthScale), sp.alpha_TT - sp.alpha_TT_offset, geo.theta_h, sp.gamma_TT, geo.phi_d)) * geo.cos_theta_i * geo.inv_cos_theta_d2 * AI_ONEOVERPI;
        assert(isValidColor(result));
        return result;
    }

    inline float pdf_TT(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_TT - sp.alpha_TT_offset;
        float pdf_theta = (1.0f / (2.0f*std::max(geo.cos_theta_i, 1e-7f)*(A_TT-B_TT))) * (sp.beta_TT*transmissionWidthScale / (t*t + sp.beta_TT2*SQR(transmissionWidthScale)));
        float p = geo.phi-AI_PI;
        float C_TT = 2.0f * atanf(AI_PI/sp.gamma_TT);
        float pdf_phi = (1.0f / C_TT) * (sp.gamma_TT / (p*p + sp.gamma_TT*sp.gamma_TT));

        float pdf = pdf_theta * pdf_phi;
        assert(isPositiveReal(pdf));
        return pdf;
    }

    inline AtVector sample_Sb(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_R, sp.beta_R*B_WIDTH_SCALE, A_b, B_b);
        float phi = 2.0f * asinf(clamp(2.0f*u2 - 1.0f, -1.0f, 1.0f));
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline float pdf_Sb(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_R*B_WIDTH_SCALE;
        float pdf_theta = (1.0f / (2.0f*std::max(geo.cos_theta_i, 1e-7f)*(A_b-B_b))) * (sp.beta_R*B_WIDTH_SCALE / (t*t + sp.beta_R2*B_WIDTH_SCALE*B_WIDTH_SCALE));
        float pdf_phi = geo.cosphi2*0.25f;
        float pdf = pdf_theta * pdf_phi;
        assert(isPositiveReal(pdf));
        return pdf;
    }

    inline AtVector sample_Sf(float u1, float u2)
    {
        float theta_i = sampleLong(u1, theta_r, sp.alpha_TT, sp.beta_TT*F_WIDTH_SCALE, A_f, B_f);
        C_TT = 2.0f * atanf(AI_PI/(sp.gamma_TT*F_WIDTH_SCALE));
        float phi = sp.gamma_TT*F_WIDTH_SCALE * tanf(C_TT * (u2-0.5f)) + AI_PI;
        float phi_i = phi_r - phi;
        AtVector wi;
        sphericalDirection(theta_i, phi_i, V, W, U, wi);
        return wi;
    }

    inline float pdf_Sf(const SctGeo& geo)
    {
        float t = geo.theta_h-sp.alpha_TT;
        float pdf_theta = (1.0f / (2.0f*std::max(geo.cos_theta_i, 1e-7f)*(A_f-B_f))) * (sp.beta_TT*F_WIDTH_SCALE / (t*t + sp.beta_TT2*F_WIDTH_SCALE*F_WIDTH_SCALE));
        float p = geo.phi-AI_PI;
        float C_TT = 2.0f * atanf(AI_PI/(sp.gamma_TT*F_WIDTH_SCALE));
        float pdf_phi = (1.0f / C_TT) * (sp.gamma_TT*F_WIDTH_SCALE / (p*p + sp.gamma_TT*sp.gamma_TT*F_WIDTH_SCALE*F_WIDTH_SCALE));

        float pdf = pdf_theta * pdf_phi;
        assert(isPositiveReal(pdf));
        return pdf;
    }

    static AtVector Hair_Sample_R(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        return hb->sample_R(u1, u2);
    }

    static AtVector Hair_Sample_TT(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        return hb->sample_TT(u1, u2);
    }

    static AtVector Hair_Sample_TRT(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        return hb->sample_TRT(u1, u2);
    }

    static AtVector Hair_Sample_TRTg(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        return hb->sample_TRTg(u1, u2);
    }

    static AtRGB Hair_Bsdf_R(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        return hb->bsdf_R(geo) * kfr[0] * hb->specular1Color;
    }
    static AtRGB Hair_Bsdf_TT(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        return hb->bsdf_TT(geo) * kfr[1] * hb->transmissionColor;
    }
    static AtRGB Hair_Bsdf_TRT(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        return hb->bsdf_TRT(geo) * kfr[2] * hb->specular2Color;
    }
    static AtRGB Hair_Bsdf_TRTg(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);
        return hb->bsdf_TRTg(geo) * kfr[2] * hb->specular2Color * hb->glintStrength;
    }

    static AtRGB HairGlossyBsdf(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        AtRGB result = AI_RGB_BLACK;

        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);

        AtRGB kfr[3];
        hairAttenuation(hb->sp.ior, geo.theta_d, geo.phi_d, hb->sp.absorption, kfr);

        AtRGB f_R = hb->do_glossy ? hb->bsdf_R(geo) * kfr[0] * hb->specular1Color : AI_RGB_BLACK;
        AtRGB f_TT = hb->do_glossy ? hb->bsdf_TT(geo) * kfr[1] * hb->transmissionColor : AI_RGB_BLACK;
        AtRGB f_TRT = hb->do_glossy ? hb->bsdf_TRT(geo) * kfr[2] * hb->specular2Color : AI_RGB_BLACK;
        AtRGB f_TRTg = hb->do_glossy ? hb->bsdf_TRTg(geo) * kfr[2] * hb->specular2Color * hb->glintStrength : AI_RGB_BLACK;

        assert(isValidColor(f_R));
        assert(isValidColor(f_TT));
        assert(isValidColor(f_TRT));
        assert(isValidColor(f_TRTg));

        return f_R + f_TT + f_TRT + f_TRTg;
    }

    static AtVector HairGlossySample(const void* brdf_data, float u1, float u2)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        hb->is_bsdf_sample = true;
        if (u1 < 0.5f && u2 < 0.5f)
        {
            u1 *= 2.0f;
            u2 *= 2.0f;
            return hb->sample_R(u1, u2);
        }
        else if (u1 < 0.5f && u2 >= 0.5f)
        {
            u1 *= 2.0f;
            u2 = 2.0f * (1.0f-u2);
            return hb->sample_TT(u1, u2);
        }
        else if (u1 >= 0.5f && u2 < 0.5f)
        {
            u1 = 2.0f * (1.0f-u1);
            u2 *= 2.0f;
            return hb->sample_TRT(u1, u2);
        }
        else
        {
            u2 = 2.0f * (1.0f-u1);
            u2 = 2.0f * (1.0f-u2);
            return hb->sample_TRTg(u1, u2);
        }


    }

    static float HairGlossyPdf(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;

        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);

        float p = (hb->pdf_R(geo) + hb->pdf_TT(geo) + hb->pdf_TRT(geo) + hb->pdf_TRTg(geo))*0.25f;
        assert(isPositiveReal(p));

        return p;
    }

    static float Hair_Pdf_R(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return hb->pdf_R(geo);
    }

    static float Hair_Pdf_TT(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return hb->pdf_TT(geo);
    }

    static float Hair_Pdf_TRT(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return hb->pdf_TRT(geo);
    }

    static float Hair_Pdf_TRTg(const void* brdf_data, const AtVector* wi)
    {
        HairBsdf* hb = (HairBsdf*)brdf_data;
        SctGeo geo(*wi, hb->theta_r, hb->phi_r, hb->U, hb->V, hb->W);
        return hb->pdf_TRTg(geo);
    }

    /// Integrate the direct illumination for all diffuse and glossy lobes
    inline void integrateDirectMis(AtShaderGlobals* sg)
    {
        // Tell Arnold we want the full sphere for lighting.
        sg->fhemi = false;

        if (do_glossy)
        {

            //sg->skip_shadow = true;
            AiLightsPrepare(sg);

            AtLightSample ls;
            while (AiLightsGetSample(sg, ls))
            {
                // per-light specular and diffuse strength multipliers
                float specular_strength = AiLightGetSpecular(ls.Lp);
                if (specular_strength < IMPORTANCE_EPS) continue;

                // get the group assigned to this light from the hash table using the light's pointer
                // ToDoJed: Fix for Porting->v6, is the "const AtNode*" of the map appropriate?
                int lightGroup = data->lightGroups[static_cast<const AtNode*>(ls.Lp)];

                is_bsdf_sample = false;
                pdf_bb = pdf_lb = 0.0f;
                AtRGB f = AiEvaluateLightSample(sg, this, HairGlossySample, HairGlossyBsdf, HairGlossyPdf);

                // mesh lights return infinite PDFs?
                if (!AiIsFinite(pdf_ll) || !AiIsFinite(pdf_lb)) continue;

                result_R_direct += f * specular_strength;
                assert(isValidColor(result_R_direct));
                assert(isValidColor(result_TT_direct));
                assert(isValidColor(result_TRT_direct));
                assert(isValidColor(result_TRTg_direct));

                //// if the light is assigned a valid group number, add this sample's contribution to that light group
                if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                {
                    lightGroupsDirect[lightGroup] += f * specular_strength;
                }
            }

        }

        if (do_diffuse)
        {
            AiLightsPrepare(sg);
            AtLightSample ls;
            while (AiLightsGetSample(sg, ls))
            {
                // per-light specular and diffuse strength multipliers
                float diffuse_strength = AiLightGetDiffuse(ls.Lp);
                if (diffuse_strength < IMPORTANCE_EPS) continue;

                // get the group assigned to this light from the hash table using the light's pointer
                int lightGroup = data->lightGroups[ls.Lp];

                float cos_theta_i = AiV3Dot(U, ls.Ld);
                float sin_theta_i = sqrtf(AiMax(1.0f - SQR(cos_theta_i), 0.0f));

                AtRGB L = ls.Li * JedPortGetShaderWeight(sg) * diffuse_strength * sin_theta_i * hairColor * (1.0f / (4*AI_PI));
                assert(isValidColor(L));
                result_Pl_direct += L;

                // if the light is assigned a valid group number, add this sample's contribution to that light group
                if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                {
                    lightGroupsDirect[lightGroup] += L;
                }
            }
        }

        sg->fhemi = true;
    }

    inline void integrateDirectDualMis(AtShaderGlobals* sg)
    {
        float als_hairNumIntersections = 0;
        AtRGB T_f = AI_RGB_BLACK;
        AtRGB als_hairOpacity = AI_RGB_WHITE;
        AtRGB sigma_f = AI_RGB_BLACK;
        bool old_hemi = sg->fhemi;
        sg->fhemi = false;
        AiLightsPrepare(sg);
        AtRay ray;
        AtScrSample scrs;
        float directFraction;
        AtRGB occlusion;
        AtRGB F_direct;
        AtRGB F_scatter;

        AtRGB kfr[3];

        if (do_glossy)
        {

            AiLightsPrepare(sg);

            AtLightSample ls;
            while (AiLightsGetSample(sg, ls))
            {
                // per-light specular and diffuse strength multipliers
                float specular_strength = AiLightGetSpecular(ls.Lp);
                if (specular_strength < IMPORTANCE_EPS) continue;

                // get the group assigned to this light from the hash table using the light's pointer
                int lightGroup = data->lightGroups[ls.Lp];

                is_bsdf_sample = false;
                pdf_bb = pdf_lb = 0.0f;
                AtRGB f = AiEvaluateLightSample(sg, this, HairGlossySample, HairGlossyBsdf, HairGlossyPdf);

                // mesh lights return infinite PDFs?
                if (!AiIsFinite(pdf_ll) || !AiIsFinite(pdf_lb)) continue;

                result_R_direct += f * specular_strength;
                assert(isValidColor(result_R_direct));
                assert(isValidColor(result_TT_direct));
                assert(isValidColor(result_TRT_direct));
                assert(isValidColor(result_TRTg_direct));

                //// if the light is assigned a valid group number, add this sample's contribution to that light group
                if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                {
                    lightGroupsDirect[lightGroup] += f * specular_strength;
                }
            }

        }

        if (do_diffuse)
        {
            bool old_skipshadow = sg->skip_shadow;
            sg->skip_shadow = true;
            AiLightsPrepare(sg);
            AtLightSample ls;
            while (AiLightsGetSample(sg, ls))
            {
                // per-light specular and diffuse strength multipliers
                float diffuse_strength = AiLightGetDiffuse(ls.Lp);
                if (diffuse_strength < IMPORTANCE_EPS) continue;

                // get the group assigned to this light from the hash table using the light's pointer
                int lightGroup = data->lightGroups[ls.Lp];
                SctGeo geo_l(ls.Ld, theta_r, phi_r, U, V, W);

                AiStateSetMsgInt(AtString("als_raytype"), ALS_RAY_DUAL);
                AiStateSetMsgFlt(AtString("als_hairNumIntersections"), 0);
                AiStateSetMsgRGB(AtString("als_T_f"), AI_RGB_WHITE);
                AiStateSetMsgRGB(AtString("als_hairOpacity"), AI_RGB_WHITE);
                AiStateSetMsgRGB(AtString("als_sigma_bar_f"), AI_RGB_BLACK);
                ray = AiMakeRay(AI_RAY_SHADOW, sg->P, &ls.Ld, ls.Ldist, sg);
                AiTrace(ray, JedPortGetAiTraceWeight(sg), scrs);
                AiStateGetMsgFlt(AtString("als_hairNumIntersections"), &als_hairNumIntersections);
                AiStateGetMsgRGB(AtString("als_T_f"), &T_f);
                AiStateGetMsgRGB(AtString("als_hairOpacity"), &als_hairOpacity);
                AiStateGetMsgRGB(AtString("als_sigma_bar_f"), &sigma_f);

                directFraction = 1.0f - std::min(als_hairNumIntersections, float(numBlendHairs))/float(numBlendHairs);
                occlusion = AI_RGB_WHITE - scrs.opacity;

                T_f = pow(data->ds->forward_attenuation(sp, geo_l.theta_i), als_hairNumIntersections);
                assert(AiIsFinite(T_f));
                F_direct = directFraction * density_back*data->ds->f_back_direct(sp, geo_l); //< do f_s_direct separately for AOVs
                F_scatter = (1.0f-directFraction) * T_f * density_front * (data->ds->f_s_scatter(sp, geo_l, sigma_f) + AI_PI * density_back * data->ds->f_back_scatter(sp, geo_l, sigma_f));

                assert(isValidColor(F_direct));
                assert(isValidColor(F_scatter));

                AtRGB L = ls.Li * JedPortGetShaderWeight(sg) * occlusion * diffuse_strength;
                AtRGB f_Pl = L * F_direct * geo_l.cos_theta_i * AI_ONEOVERPI * diffuseColor;
                AtRGB f_Pg = L * F_scatter * geo_l.cos_theta_i * AI_ONEOVERPI * diffuseColor;
                result_Pl_direct += f_Pl;
                result_Pg_direct += f_Pg;

                // if the light is assigned a valid group number, add this sample's contribution to that light group
                if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
                {
                    lightGroupsDirect[lightGroup] += f_Pl + f_Pg;
                }
            }
            sg->skip_shadow = old_skipshadow;
        }

        AiStateSetMsgInt(AtString("als_raytype"), ALS_RAY_UNDEFINED);


        sg->fhemi = old_hemi;

    }

    /// Integrate the indirect illumination for all diffuse and glossy lobes
    inline void integrateIndirect(AtShaderGlobals* sg)
    {
        if (do_glossy && glossyIndirectStrength > 0.0f)
        {
            wi_ray = AiMakeRay(AI_RAY_SPECULAR_REFLECT, sg->P, NULL, AI_BIG, sg);

            sampit = AiSamplerIterator(data->sampler_glossy, sg);
            AiStateSetMsgInt(AtString("als_raytype"), ALS_RAY_HAIR);
            while(AiSamplerGetSample(sampit, samples))
            {
                wi_ray.dir = HairGlossySample(this, samples[0], samples[1]);

                AtScrSample scrs;

                // trace our ray
                bool hit = AiTrace(wi_ray, JedPortGetAiTraceWeight(sg), scrs);
                // calculate result
                const float p = HairGlossyPdf(this, &wi_ray.dir);
                const AtRGB f_R = Hair_Bsdf_R(this, &wi_ray.dir) / p;
                const AtRGB f_TT = Hair_Bsdf_TT(this, &wi_ray.dir) / p;
                const AtRGB f_TRT = Hair_Bsdf_TRT(this, &wi_ray.dir) / p;
                const AtRGB f_TRTg = Hair_Bsdf_TRTg(this, &wi_ray.dir) / p;
                result_R_indirect += scrs.color * (f_R + f_TT + f_TRT + f_TRTg) * AI_ONEOVERPI;

                if (deepGroupPtr && hit)
                {
                    for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
                    {
                        deepGroupsSpecular[i] += deepGroupPtr[i] * (f_R + f_TT + f_TRT + f_TRTg) * AI_ONEOVERPI;
                    }
                }
            }
            AiStateSetMsgInt(AtString("als_raytype"), ALS_RAY_UNDEFINED);
            float weight = AiSamplerGetSampleInvCount(sampit) * glossyIndirectStrength;
            result_R_indirect *= weight;

            if (deepGroupPtr)
            {
                for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
                {
                    deepGroupsSpecular[i] *= weight;
                    //deepGroupsSpecular2[i] *= weight;
                    //deepGroupsTransmission[i] *= weight;
                }
            }
        }

        if (diffuseIndirectStrength > 0.0f)
        {
            // result_Pl_indirect = AiIndirectDiffuse(&wo, sg) * hairColor * diffuseIndirectStrength;
            AtRay wi_ray;
            AtVector wi;
            AtScrSample scrs;
            wi_ray = AiMakeRay(AI_RAY_ALL_DIFFUSE | AI_RAY_VOLUME, sg->P, NULL, AI_BIG, sg);
            sampit = AiSamplerIterator(data->sampler_diffuse, sg);
            result_Pl_indirect = AI_RGB_BLACK;
            while (AiSamplerGetSample(sampit, samples))
            {
                // cosine hemisphere sampling as O-N sampling does not work outside of a light loop
                float stheta = sqrtf(float(samples[0]));
                float phi = float(AI_PITIMES2 * samples[1]);
                wi.x = stheta * cosf(phi);
                wi.y = stheta * sinf(phi);
                wi.z = sqrtf(1.0f - float(samples[0]));
                AiV3RotateToFrame(wi, U, V, sg->Nf);

                float cos_theta = AiV3Dot(wi, sg->Nf);
                if (cos_theta <= 0.0f) continue;

                // trace the ray
                wi_ray.dir = wi;

                bool hit = AiTrace(wi_ray, JedPortGetAiTraceWeight(sg), scrs);

                result_Pl_indirect += scrs.color * hairColor * diffuseIndirectStrength;
                assert(isValidColor(result_Pl_indirect));
                if (deepGroupPtr && hit)
                {
                    for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
                    {
                        deepGroupsDiffuse[i] += deepGroupPtr[i] * hairColor * diffuseIndirectStrength;
                    }
                }

            }
            result_Pl_indirect *= AiSamplerGetSampleInvCount(sampit);
            if (deepGroupPtr)
            {
                for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
                {
                    deepGroupsDiffuse[i] *= AiSamplerGetSampleInvCount(sampit);
                }
            }
        }
    }

    inline void integrateIndirectDual(AtShaderGlobals* sg)
    {
        wi_ray = AiMakeRay(AI_RAY_SPECULAR_REFLECT, sg->P, NULL, AI_BIG, sg);
        float weight;

        bool do_glossy = true;
        if (sg->Rt & AI_RAY_ALL_DIFFUSE | AI_RAY_VOLUME || sg->bounces > 0) do_glossy = false;

        if (do_glossy && glossyIndirectStrength > 0.0f)
        {
            AiStateSetMsgInt(AtString("als_raytype"), ALS_RAY_HAIR);
            sampit = AiSamplerIterator(data->sampler_glossy, sg);
            while(AiSamplerGetSample(sampit, samples))
            {
                wi_ray.dir = HairGlossySample(this, samples[0], samples[1]);

                AtScrSample scrs;

                // trace our ray
                AiStateSetMsgFlt(AtString("alsPreviousRoughness"), 1.0f);
                // trace our ray
                bool hit = AiTrace(wi_ray, JedPortGetAiTraceWeight(sg), scrs);
                // calculate result
                const float p = HairGlossyPdf(this, &wi_ray.dir);
                const AtRGB f_R = Hair_Bsdf_R(this, &wi_ray.dir) / p;
                const AtRGB f_TT = Hair_Bsdf_TT(this, &wi_ray.dir) / p;
                const AtRGB f_TRT = Hair_Bsdf_TRT(this, &wi_ray.dir) / p;
                const AtRGB f_TRTg = Hair_Bsdf_TRTg(this, &wi_ray.dir) / p;
                result_R_indirect += scrs.color * (f_R + f_TT + f_TRT + f_TRTg) * AI_ONEOVERPI;
                //result_TT_indirect += scrs.color * f_TT;
                //result_TRT_indirect += scrs.color * f_TRT;
                //result_TRTg_indirect += scrs.color * f_TRTg;

                if (deepGroupPtr && hit)
                {
                    for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
                    {
                        deepGroupsSpecular[i] += deepGroupPtr[i] * (f_R + f_TT + f_TRT + f_TRTg) * AI_ONEOVERPI;
                        //deepGroupsSpecular[i] += deepGroupPtr[i] * f_R;
                        //deepGroupsSpecular2[i] += deepGroupPtr[i] * (f_TRT + f_TRTg);
                        //deepGroupsTransmission[i] += deepGroupPtr[i] * f_TT;
                    }
                }
            }
            AiStateSetMsgInt(AtString("als_raytype"), ALS_RAY_UNDEFINED);
            weight = AiSamplerGetSampleInvCount(sampit) * glossyIndirectStrength;
            result_R_indirect *= weight;

            if (deepGroupPtr)
            {
                for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
                {
                    deepGroupsSpecular[i] *= weight;
                    //deepGroupsSpecular2[i] *= weight;
                    //deepGroupsTransmission[i] *= weight;
                }
            }
        }

        if (diffuseIndirectStrength > 0.0f)
        {
            float als_hairNumIntersections = 0;
            AtRGB T_f = AI_RGB_BLACK;
            AtRGB als_hairOpacity = AI_RGB_WHITE;
            AtRGB sigma_f = AI_RGB_BLACK;
            sampit = AiSamplerIterator(data->sampler_diffuse, sg);
            AiStateSetMsgInt(AtString("als_raytype"), ALS_RAY_DUAL);
            while (AiSamplerGetSample(sampit, samples))
            {
                if (samples[0] < 0.5f)
                {
                    samples[0] *= 2.0f;
                    wi_ray.dir = sample_Sb(samples[0], samples[1]);
                }
                else
                {
                    samples[0] = 1.0f - (2.0f * samples[0]);
                    wi_ray.dir = sample_Sf(samples[0], samples[1]);
                }

                SctGeo geo(wi_ray.dir, theta_r, phi_r, U, V, W);

                // float p = 1.0f / ((pdf_Sf(geo) + pdf_Sb(geo)) * 0.5f);

                AiStateSetMsgFlt(AtString("als_hairNumIntersections"), 0);
                AiStateSetMsgRGB(AtString("als_T_f"), AI_RGB_WHITE);
                AiStateSetMsgRGB(AtString("als_hairOpacity"), AI_RGB_WHITE);
                AiStateSetMsgRGB(AtString("als_sigma_bar_f"), AI_RGB_BLACK);
                AiStateSetMsgFlt(AtString("alsPreviousRoughness"), 1.0f);
                bool als_hitHair = false;
                AiStateSetMsgBool(AtString("als_hitHair"), als_hitHair);
                bool hit = AiTrace(wi_ray, JedPortGetAiTraceWeight(sg), scrs);
                AiStateGetMsgFlt(AtString("als_hairNumIntersections"), &als_hairNumIntersections);
                AiStateGetMsgRGB(AtString("als_T_f"), &T_f);
                AiStateGetMsgRGB(AtString("als_hairOpacity"), &als_hairOpacity);
                AiStateGetMsgRGB(AtString("als_sigma_bar_f"), &sigma_f);

                float directFraction = 1.0f - std::min(als_hairNumIntersections, float(numBlendHairs))/float(numBlendHairs);
                // AtRGB occlusion = AI_RGB_WHITE - scrs.opacity;

                T_f = pow(data->ds->forward_attenuation(sp, geo.theta_i), als_hairNumIntersections);
                assert(AiIsFinite(T_f));
                AtRGB F_direct = directFraction * density_back*data->ds->f_back_direct(sp, geo); //< do f_s_direct separately for AOVs
                AtRGB F_scatter = (1.0f-directFraction) * T_f * density_front *
                                    (data->ds->f_s_scatter(sp, geo, sigma_f) + AI_PI * density_back * data->ds->f_back_scatter(sp, geo, sigma_f));

                const AtRGB f_Pl =  F_direct * geo.cos_theta_i * AI_ONEOVERPI;
                const AtRGB f_Pg =  F_scatter * geo.cos_theta_i * AI_ONEOVERPI;
                result_Pl_indirect += scrs.color * f_Pl;
                result_Pg_indirect += scrs.color * f_Pg;

                // only accumulate deep groups if we hit something that wasn't a hair
                if (deepGroupPtr && (hit && !als_hitHair))
                {
                    for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
                    {
                        deepGroupsDiffuse[i] += deepGroupPtr[i] * (f_Pl + f_Pg);
                    }
                }

            }
            AiStateSetMsgInt(AtString("als_raytype"), ALS_RAY_UNDEFINED);
            weight = AiSamplerGetSampleInvCount(sampit) * diffuseIndirectStrength;
            result_Pl_indirect *= weight * diffuseColor;
            result_Pg_indirect *= weight * diffuseColor;

            if (deepGroupPtr)
            {
                for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
                {
                    deepGroupsDiffuse[i] *= weight * diffuseColor;
                }
            }
        }
    }

    inline void writeResult(AtShaderGlobals* sg)
    {
        if (deepGroupPtr)
        {
            for (int i=0; i < NUM_LIGHT_GROUPS; ++i)
            {
                deepGroupPtr[i] = lightGroupsDirect[i] + deepGroupsDiffuse[i] + deepGroupsSpecular[i] ;/*
                                + deepGroupsSpecular2[i] + deepGroupsTransmission[i];*/
            }
        }

        if (sg->Rt & AI_RAY_CAMERA)
        {
            AiAOVSetRGB(sg, data->aovs[k_diffuse_color], hairColor);
            AiAOVSetRGB(sg, data->aovs[k_direct_diffuse], result_Pg_direct + result_Pl_direct);
            AiAOVSetRGB(sg, data->aovs[k_indirect_diffuse], result_Pg_indirect + result_Pl_indirect);
            AiAOVSetRGB(sg, data->aovs[k_direct_specular], result_R_direct);
            AiAOVSetRGB(sg, data->aovs[k_indirect_specular], result_R_indirect);
            AiAOVSetRGB(sg, data->aovs[k_direct_specular_2], result_TRT_direct + result_TRTg_direct);
            AiAOVSetRGB(sg, data->aovs[k_indirect_specular_2], result_TRT_indirect + result_TRTg_indirect);
            AiAOVSetRGB(sg, data->aovs[k_direct_glint], result_TRTg_direct);
            AiAOVSetRGB(sg, data->aovs[k_indirect_glint], result_TRTg_indirect);
            AiAOVSetRGB(sg, data->aovs[k_direct_transmission], result_TT_direct);
            AiAOVSetRGB(sg, data->aovs[k_indirect_transmission], result_TT_indirect);
            AiAOVSetRGB(sg, data->aovs[k_direct_global], result_Pg_direct);
            AiAOVSetRGB(sg, data->aovs[k_direct_local], result_Pl_direct);
            AiAOVSetRGB(sg, data->aovs[k_indirect_global], result_Pg_indirect);
            AiAOVSetRGB(sg, data->aovs[k_indirect_local], result_Pl_indirect);

            // write depth AOV
            AtRGB depth = AtRGB(float(sg->Rl), AiV3Dot(sg->Nf, wo), sg->P.y);
            AiAOVSetRGB(sg, data->aovs[k_depth], depth);

            assert(isValidColor(result_Pg_direct));
            assert(isValidColor(result_Pg_indirect));
            assert(isValidColor(result_Pl_direct));
            assert(isValidColor(result_Pl_indirect));
            assert(isValidColor(result_R_direct));
            assert(isValidColor(result_R_indirect));
            assert(isValidColor(result_TT_direct));
            assert(isValidColor(result_TT_indirect));
            assert(isValidColor(result_TRT_direct));
            assert(isValidColor(result_TRT_indirect));
            assert(isValidColor(result_TRTg_direct));
            assert(isValidColor(result_TRTg_indirect));

            for (int i=0; i < NUM_ID_AOVS; ++i)
            {
                AtRGB tmp;
                // check if output is enabled first in case we have an expensive network upstream
                if (AiAOVEnabled(data->aovs[k_id_1+i], AI_TYPE_RGB))
                {
                    tmp = AiShaderEvalParamRGB(p_id1 + i);

                    // check if we're overriding it with a per-object id
                    if (AiNodeLookUpUserParameter(sg->Op, id_names[i]))
                    {
                        tmp = AiNodeGetRGB(sg->Op, id_names[i]);
                    }


                    AiAOVSetRGB(sg, data->aovs[k_id_1+i], tmp);
                }
            }

            for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
            {
                AiAOVSetRGB(sg, data->aovs[k_light_group_1+i], deepGroupPtr[i]);
            }
        }

        sg->out.RGB() =   result_R_direct +
                        result_R_indirect +
                        result_TT_direct +
                        result_TT_indirect +
                        result_TRT_direct +
                        result_TRT_indirect +
                        result_TRTg_direct +
                        result_TRTg_indirect +
                        result_Pg_direct +
                        result_Pl_direct +
                        result_Pg_indirect +
                        result_Pl_indirect;
        // ToDoJed: Fix for Porting->v6, opacity needs fix, out.CLOSURE()
        // sg->out_opacity = AI_RGB_WHITE;
    }

    AtVector U, V, W;   //< local coordinate frame
    float theta_r;      //< exitant spherical theta
    float phi_r;        //< existant spherical phi
    float cn;           //< random value per curve in [0,1)

    ScatteringParams sp; //< varying scattering parameters

    float diffuseStrength;
    AtRGB diffuseColor;
    int diffuseScatteringMode;
    AtRGB hairColor;
    AtRGB specular1Color;
    AtRGB specular2Color;
    AtRGB transmissionColor;
    float glintStrength;
    float diffuseIndirectStrength;
    float glossyIndirectStrength;

    float A_R;
    float B_R;
    float A_TT;
    float B_TT;
    float C_TT;
    float A_TRT;
    float B_TRT;
    float Cg;
    float Dg;

    float A_f;
    float A_b;
    float B_f;
    float B_b;

    // dual-scattering parameters
    int numBlendHairs;
    float density_front;
    float density_back;
    float colorNorm;

    AtRGB result_R_direct;
    AtRGB result_R_indirect;
    AtRGB result_TT_direct;
    AtRGB result_TT_indirect;
    AtRGB result_TRT_direct;
    AtRGB result_TRT_indirect;
    AtRGB result_TRTg_direct;
    AtRGB result_TRTg_indirect;

    // dual-scattering aovs
    AtRGB result_Pg_direct;
    AtRGB result_Pl_direct;
    AtRGB result_Pg_indirect;
    AtRGB result_Pl_indirect;

    // IDs
    AtRGB result_id1;
    AtRGB result_id2;
    AtRGB result_id3;
    AtRGB result_id4;
    AtRGB result_id5;
    AtRGB result_id6;
    AtRGB result_id7;
    AtRGB result_id8;

    bool do_diffuse;
    bool do_glossy;

    int depth;

    AtNode* node;
    ShaderData* data;
    AtRay wi_ray;
    AtScrSample scrs;
#if AI_VERSION_MAJOR_NUM > 0 || AI_VERSION_ARCH_NUM > 4
    float samples[2];
#else
    double samples[2];
#endif
    AtSamplerIterator* sampit;

    AtVector wo;

    AtShaderGlobals* _sg;

    float specular1WidthScale;
    float specular2WidthScale;
    float transmissionWidthScale;

    bool dualImportant;

    // MIS temporaries
    // {
    bool is_bsdf_sample;
    AtRGB L_l;
    AtRGB f_R_l;
    AtRGB f_TT_l;
    AtRGB f_TRT_l;
    AtRGB f_TRTg_l;
    float pdf_ll;
    float pdf_lb;
    float d_l;
    AtVector w_l;
    AtRGB L_b;
    AtRGB f_R_b;
    AtRGB f_TT_b;
    AtRGB f_TRT_b;
    AtRGB f_TRTg_b;
    float pdf_bb;
    float pdf_bl;
    float d_b;
    AtVector w_b;
    // }

    bool doLightGroups;
    AtRGB lightGroupsDirect[NUM_LIGHT_GROUPS];
    AtRGB* deepGroupPtr;
    AtRGB deepGroupsSpecular[NUM_LIGHT_GROUPS];
    AtRGB deepGroupsSpecular2[NUM_LIGHT_GROUPS];
    AtRGB deepGroupsDiffuse[NUM_LIGHT_GROUPS];
    AtRGB deepGroupsTransmission[NUM_LIGHT_GROUPS];

    AtRGB opacity;
};

node_loader
{
   if (i>0) return 0;
   node->methods     = alHair;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alHair";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return true;
}

node_initialize
{
    HairBsdf::ShaderData* data = new HairBsdf::ShaderData;
    data->cryptomatte = new CryptomatteData();
    AiNodeSetLocalData(node, data);
}

node_finish
{
    if (AiNodeGetLocalData(node))
    {
        HairBsdf::ShaderData* data = (HairBsdf::ShaderData*)AiNodeGetLocalData(node);
        if (data->cryptomatte)
            delete data->cryptomatte;
        delete data;
    }
}

node_update
{
    HairBsdf::ShaderData* data = (HairBsdf::ShaderData*)AiNodeGetLocalData(node);
    data->update();

    AiAOVRegister("diffuse_color", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_specular", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_specular", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_specular_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_specular_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_transmission", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_transmission", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_glint", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_global", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_global", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("direct_local", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("indirect_local", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_1", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_3", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_4", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_5", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_6", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_7", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("id_8", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_1", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_2", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_3", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_4", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_5", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_6", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_7", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_8", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_9", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_10", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_11", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_12", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_13", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_14", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_15", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);
    AiAOVRegister("light_group_16", AI_TYPE_RGB, AI_AOV_BLEND_OPACITY);

    data->cryptomatte->setup_all(AiNodeGetStr(node, "aov_crypto_asset"),
        AiNodeGetStr(node, "aov_crypto_object"), AiNodeGetStr(node, "aov_crypto_material"));
    // Get all the light nodes in the scene and try and find their light group parameter
    // we'll store this based on the light pointer for fast access during rendering
    AtNodeIterator* it = AiUniverseGetNodeIterator(AI_NODE_LIGHT);
    while (!AiNodeIteratorFinished(it))
    {
        AtNode* light = AiNodeIteratorGetNext(it);
        if (AiNodeLookUpUserParameter(light, "lightGroup"))
            data->lightGroups[light] = AiNodeGetInt(light, "lightGroup") - 1;
        else
            data->lightGroups[light] = -1;
    }
}

shader_evaluate
{

    // Get shader data
    HairBsdf::ShaderData* data = (HairBsdf::ShaderData*)AiNodeGetLocalData(node);

    // Create HairBsdf object 
    HairBsdf hb(node, sg, data);
    // Get parameters
    hb.evaluateParameters(sg, data);

    float geo_opacity = 1.0f;
    if (AiUDataGetFlt(AtString("geo_opacity"), geo_opacity))
    {
        hb.opacity *= geo_opacity;
    }

    const bool transp_early_out = AiShaderGlobalsApplyOpacity(sg, hb.opacity);
    data->cryptomatte->do_cryptomattes(sg, node, p_crypto_asset_override, p_crypto_object_override, p_crypto_material_override);
    if (transp_early_out)
      return;

    // ToDoJed: Fix for Porting->v6, opacity needs fix, out.CLOSURE()
    // hb.opacity = sg->out_opacity;

    AiStateSetMsgBool(AtString("als_hitHair"), true);

    float als_hairNumIntersections = 0;
    AtRGB als_sigma_bar_f = AI_RGB_BLACK;
    bool do_dual = false;
    if (sg->bounces >= data->dual_depth) do_dual = true;

    int als_raytype = ALS_RAY_UNDEFINED;
    AiStateGetMsgInt(AtString("als_raytype"), &als_raytype);

    if (do_dual && als_raytype == ALS_RAY_DUAL)
    {
        if (AiStateGetMsgFlt(AtString("als_hairNumIntersections"), &als_hairNumIntersections)
            && AiStateGetMsgRGB(AtString("als_sigma_bar_f"), &als_sigma_bar_f))
        {
            float theta_i = AI_PIOVER2 - sphericalTheta(sg->Rd, hb.U);

            als_sigma_bar_f += hb.sp.beta_R2 + hb.sp.beta_TRT2 + hb.sp.beta_TT2;
            AiStateSetMsgRGB(AtString("als_sigma_bar_f"), als_sigma_bar_f);

            als_hairNumIntersections+=minh(hb.opacity);
            AiStateSetMsgFlt(AtString("als_hairNumIntersections"), als_hairNumIntersections);

            // ToDoJed: Fix for Porting->v6, opacity needs fix, out.CLOSURE()
            // sg->out_opacity = AI_RGB_BLACK;
        }
        else
        {
            // ToDoJed: Fix for Porting->v6, opacity needs fix, out.CLOSURE()
            // sg->out_opacity = AI_RGB_WHITE;
        }

        return; // early out
    }

    // early-out regardless if we're in a shadow ray, or if opacity is zero
#if AI_VERSION_MINOR_NUM >= 2
    if (sg->Rt & AI_RAY_SHADOW || AiColorIsSmall(hb.opacity) || AiShaderGlobalsIsObjectMatte(sg)) return;
#else
    if (sg->Rt & AI_RAY_SHADOW || AiColorIsSmall(hb.opacity)) return;
#endif
    // early out if we're a hair-hair glossy ray and the ray depth says we should be calculating dual scattering only
    if (sg->Rr_gloss > data->dual_depth && als_raytype == ALS_RAY_HAIR)
    {
        sg->out.RGB() = AI_RGB_BLACK;
        return;
    }

    // calculate scattering explicitly up to the dual depth cutoff.
    // in other words, do a brute force path trace for x=dual-depth bounces, then fall back to dual scattering for the rest.
    if (do_dual && hb.dualImportant && hb.diffuseScatteringMode == DUAL)
    {
        hb.integrateDirectDualMis(sg);
        hb.integrateIndirectDual(sg);
    }
    else
    {
        hb.integrateDirectMis(sg);
        hb.integrateIndirect(sg);
    }

    // Write shader result
    hb.writeResult(sg);
    // ToDoJed: Fix for Porting->v6, opacity needs fix, out.CLOSURE()
    // sg->out_opacity = hb.opacity;
}


