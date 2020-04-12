#include <iostream>
#include <ai_sampler.h>
#include <map>
#include <cassert>
#include <cstdlib>
#include <vector>

#include "alUtil.h"
#include "MIS.h"
#include "alSurface.h"
#include "aovs.h"
#include "tea.h"
#include "fresnel.h"
#include "microfacet.h"
#include "sss.h"

#define RR_BOUNCES

AI_SHADER_NODE_EXPORT_METHODS(alSurfaceMtd)

/*
#define GlossyMISBRDF AiCookTorranceMISBRDF
#define GlossyMISPDF AiCookTorranceMISPDF
#define GlossyMISSample AiCookTorranceMISSample

#define GlossyMISBRDF_wrap AiCookTorranceMISBRDF_wrap
#define GlossyMISPDF_wrap AiCookTorranceMISPDF_wrap
#define GlossyMISSample_wrap AiCookTorranceMISSample_wrap
*/

#define GlossyMISBRDF AiMicrofacetMISBRDF
#define GlossyMISPDF AiMicrofacetMISPDF
#define GlossyMISSample AiMicrofacetMISSample

#define GlossyMISBRDF_wrap AiMicrofacetMISBRDF_wrap
#define GlossyMISPDF_wrap AiMicrofacetMISPDF_wrap
#define GlossyMISSample_wrap AiMicrofacetMISSample_wrap
#define GlossyMISCreateData AiMicrofacetMISCreateData

#define NUM_ID_AOVS 8
static const char* id_names[NUM_ID_AOVS] = {
    "id_1", "id_2", "id_3", "id_4", "id_5", "id_6", "id_7", "id_8",
};

inline void flipNormals(AtShaderGlobals* sg)
{
   sg->Nf = -sg->Nf;
   sg->Ngf = -sg->Ngf;
}

static const char* fresnel_mode_names[] = {"dielectric", "metallic", NULL};

enum FresnelModes
{
   FM_DIELECTRIC = 0,
   FM_METALLIC
};

static const char* distribution_names[] = {"beckmann", "ggx", NULL};

enum alSurfaceParams
{
   // diffuse
   p_diffuseStrength = 0,
   p_diffuseColor,
   p_diffuseRoughness,

   // backlight
   p_backlightStrength,
   p_backlightColor,
   p_backlightIndirectStrength,

   p_emissionStrength,
   p_emissionColor,

   // sss
   p_sssMix,
   p_sssMode,
   p_sssRadius,
   p_sssWeight1,
   p_sssRadiusColor,
   p_sssRadius2,
   p_sssWeight2,
   p_sssRadiusColor2,
   p_sssRadius3,
   p_sssWeight3,
   p_sssRadiusColor3,
   p_sssDensityScale,
   p_sssTraceSet,
   p_sssExtraSamples,

   p_ssInScatteringStrength,
   p_ssAttenuationColor,
   p_ssSpecifyCoefficients,
   p_ssScattering,
   p_ssAbsorption,
   p_ssDensityScale,
   p_ssDirection,

   p_diffuseExtraSamples,
   p_diffuseIndirectStrength,
   p_diffuseIndirectClamp,
   p_diffuseNormal,

   // specular
   p_specular1Strength,
   p_specular1Color,
   p_specular1Roughness,
   p_specular1Anisotropy,
   p_specular1Rotation,
   p_specular1FresnelMode,
   p_specular1Ior,
   p_specular1Reflectivity,
   p_specular1EdgeTint,
   p_specular1RoughnessDepthScale,
   p_specular1ExtraSamples,
   p_specular1Normal,
   p_specular1IndirectStrength,
   p_specular1IndirectClamp,
   p_specular1CausticPaths,
   p_specular1InternalDirect,
   p_specular1Distribution,

   p_specular2Strength,
   p_specular2Color,
   p_specular2Roughness,
   p_specular2Anisotropy,
   p_specular2Rotation,
   p_specular2FresnelMode,
   p_specular2Ior,
   p_specular2Reflectivity,
   p_specular2EdgeTint,
   p_specular2RoughnessDepthScale,
   p_specular2ExtraSamples,
   p_specular2Normal,
   p_specular2IndirectStrength,
   p_specular2IndirectClamp,
   p_specular2CausticPaths,
   p_specular2InternalDirect,
   p_specular2Distribution,

   // transmission
   p_transmissionStrength,
   p_transmissionColor,
   p_transmissionLinkToSpecular1,
   p_transmissionRoughness,
   p_transmissionIor,
   p_transmissionRoughnessDepthScale,
   p_transmissionEnableCaustics,
   p_transmissionExtraSamples,
   p_transmissionClamp,
   p_transmissionDoDirect,
   p_transmissionNormal,
   p_transmissionCausticPaths,

   p_id1,
   p_id2,
   p_id3,
   p_id4,
   p_id5,
   p_id6,
   p_id7,
   p_id8,

   p_aiEnableMatte,
   p_aiMatteColor,
   p_aiMatteColorA,

   p_aov_diffuse_color,
   p_aov_direct_diffuse,
   p_aov_direct_diffuse_raw,
   p_aov_indirect_diffuse,
   p_aov_indirect_diffuse_raw,
   p_aov_direct_backlight,
   p_aov_indirect_backlight,
   p_aov_direct_specular,
   p_aov_indirect_specular,
   p_aov_direct_specular_2,
   p_aov_indirect_specular_2,
   p_aov_single_scatter,
   p_aov_sss,
   p_aov_refraction,
   p_aov_emission,
   p_aov_uv,
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
   p_aov_crypto_asset,
   p_aov_crypto_object,
   p_aov_crypto_material,

   p_aov_shadow_group_1,
   p_aov_shadow_group_2,
   p_aov_shadow_group_3,
   p_aov_shadow_group_4,
   p_aov_shadow_group_5,
   p_aov_shadow_group_6,
   p_aov_shadow_group_7,
   p_aov_shadow_group_8,
   p_aov_shadow_group_9,
   p_aov_shadow_group_10,
   p_aov_shadow_group_11,
   p_aov_shadow_group_12,
   p_aov_shadow_group_13,
   p_aov_shadow_group_14,
   p_aov_shadow_group_15,
   p_aov_shadow_group_16,

   p_standardAovs,
   p_transmitAovs,
   p_rrTransmissionDepth,

   p_opacity,

   p_rr,

   p_trace_set_all,
   p_trace_set_shadows,
   p_trace_set_diffuse,
   p_trace_set_backlight,
   p_trace_set_specular1,
   p_trace_set_specular2,
   p_trace_set_transmission,

   p_debug,

   p_aov_diffuse_color_clamp,
   p_aov_direct_diffuse_clamp,
   p_aov_direct_diffuse_raw_clamp,
   p_aov_indirect_diffuse_clamp,
   p_aov_indirect_diffuse_raw_clamp,
   p_aov_direct_backlight_clamp,
   p_aov_indirect_backlight_clamp,
   p_aov_direct_specular_clamp,
   p_aov_indirect_specular_clamp,
   p_aov_direct_specular_2_clamp,
   p_aov_indirect_specular_2_clamp,
   p_aov_single_scatter_clamp,
   p_aov_sss_clamp,
   p_aov_refraction_clamp,
   p_aov_emission_clamp,
   p_aov_light_group_1_clamp,
   p_aov_light_group_2_clamp,
   p_aov_light_group_3_clamp,
   p_aov_light_group_4_clamp,
   p_aov_light_group_5_clamp,
   p_aov_light_group_6_clamp,
   p_aov_light_group_7_clamp,
   p_aov_light_group_8_clamp,

   p_aov_light_group_9_clamp,
   p_aov_light_group_10_clamp,
   p_aov_light_group_11_clamp,
   p_aov_light_group_12_clamp,
   p_aov_light_group_13_clamp,
   p_aov_light_group_14_clamp,
   p_aov_light_group_15_clamp,
   p_aov_light_group_16_clamp,

   p_crypto_asset_override,
   p_crypto_object_override,
   p_crypto_material_override,

   p_bump
};

enum SssMode
{
   SSSMODE_CUBIC = 0,
   SSSMODE_DIFFUSION = 1,
   SSSMODE_DIRECTIONAL = 2,
   SSSMODE_EMPIRICAL = 3
};

const char* SssModeNames[] = {"cubic", "diffusion", "directional", "empirical",
                              NULL};

enum Debug
{
   DB_OFF = 0,
   DB_diffuseStrength,
   DB_diffuseColor,
   DB_diffuseRoughness,
   DB_backlightStrength,
   DB_backlightColor,
   DB_sssMix,
   DB_sssWeight1,
   DB_sssColor1,
   DB_sssWeight2,
   DB_sssColor2,
   DB_sssWeight3,
   DB_sssColor3,
   DB_specular1Strength,
   DB_specular1Color,
   DB_specular1Roughness,
   DB_specular2Strength,
   DB_specular2Color,
   DB_specular2Roughness,
   DB_transmissionStrength,
   DB_transmissionColor,
   DB_transmissionRoughness,
   DB_END
};

const char* debugNames[] = {
    "off", "diffuse-strength", "diffuse-color", "diffuse-roughness",
    "backlight-strength", "backlight-color", "sss-mix", "sss-weight-1",
    "sss-color-1", "sss-weight-2", "sss-color-2", "sss-weight-3", "sss-color-3",
    "specular-1-strength", "specular-1-color", "specular-1-roughness",
    "specular-2-strength", "specular-2-color", "specular-2-roughness",
    "transmission-strength", "transmission-color", "transmission-roughness",
    NULL};

node_parameters
{
   AiParameterFLT("diffuseStrength", 1.0f);
   AiParameterRGB("diffuseColor", 0.5f, 0.5f, 0.5f);
   AiParameterFLT("diffuseRoughness", 0.0f);

   AiParameterFLT("backlightStrength", 0.0f);
   AiParameterRGB("backlightColor", 0.18f, 0.18f, 0.18f);
   AiParameterFLT("backlightIndirectStrength", 1.0f);

   AiParameterFLT("emissionStrength", 0.0f);
   AiParameterRGB("emissionColor", 1.0f, 1.0f, 1.0f);

   AiParameterFLT("sssMix", 0.0f);
   AiParameterEnum("sssMode", SSSMODE_CUBIC, SssModeNames);
   AiParameterFLT("sssRadius", 1.5f);
   AiParameterFLT("sssWeight1", 1.0f);
   AiParameterRGB("sssRadiusColor", .439, .156, .078);
   AiParameterFLT("sssRadius2", 4.0f);
   AiParameterFLT("sssWeight2", 0.0f);
   AiParameterRGB("sssRadiusColor2", .439, .08, .018);
   AiParameterFLT("sssRadius3", .75f);
   AiParameterFLT("sssWeight3", .0f);
   AiParameterRGB("sssRadiusColor3", .523, .637, .667);
   AiParameterFLT("sssDensityScale", 1.0f);
   AiParameterSTR("sssTraceSet", "");
   AiParameterINT("sssExtraSamples", 0);

   AiParameterFLT("ssInScatteringStrength", 0.0f);
   AiParameterRGB("ssAttenuationColor", 1.0f, 1.0f, 1.0f);
   AiParameterBOOL("ssSpecifyCoefficients", false);
   AiParameterRGB("ssScattering", 1.0f, 1.0f, 1.0f);
   AiParameterRGB("ssAbsorption", 1.0f, 1.0f, 1.0f);
   AiParameterFLT("ssDensityScale", 1.0f);
   AiParameterFLT("ssDirection", 0.0f);

   AiParameterINT("diffuseExtraSamples", 0);
   AiParameterFLT("diffuseIndirectStrength", 1.0f);
   AiParameterFLT("diffuseIndirectClamp", 0.0f);
   AiParameterVec("diffuseNormal", 0, 0, 0);

   AiParameterFLT("specular1Strength", 1.0f);
   AiParameterRGB("specular1Color", 1.0f, 1.0f, 1.0f);
   AiParameterFLT("specular1Roughness", 0.3f);
   AiParameterFLT("specular1Anisotropy", 0.5f);
   AiParameterFLT("specular1Rotation", 0.0f);
   AiParameterENUM("specular1FresnelMode", FM_DIELECTRIC, fresnel_mode_names);
   AiParameterFLT("specular1Ior", 1.4f);
   AiParameterRGB("specular1Reflectivity", 0.548, .549, .570);
   AiParameterRGB("specular1EdgeTint", 0.579, .598, .620);
   AiParameterFLT("specular1RoughnessDepthScale", 1.0f);
   AiParameterINT("specular1ExtraSamples", 0);
   AiParameterVec("specular1Normal", 0, 0, 0);
   AiParameterFLT("specular1IndirectStrength", 1.0f);
   AiParameterFLT("specular1IndirectClamp", 0.0f);
   AiParameterBOOL("specular1CausticPaths", false);
   AiParameterBOOL("specular1InternalDirect", true);
   AiParameterENUM("specular1Distribution", AI_MICROFACET_BECKMANN,
                   distribution_names);

   AiParameterFLT("specular2Strength", 0.0f);
   AiParameterRGB("specular2Color", 1.0f, 1.0f, 1.0f);
   AiParameterFLT("specular2Roughness", 0.5f);
   AiParameterFLT("specular2Anisotropy", 0.5f);
   AiParameterFLT("specular2Rotation", 0.0f);
   AiParameterENUM("specular2FresnelMode", FM_DIELECTRIC, fresnel_mode_names);
   AiParameterFLT("specular2Ior", 1.4f);
   AiParameterRGB("specular2Reflectivity", 0.548, .549, .570);
   AiParameterRGB("specular2EdgeTint", 0.579, .598, .620);
   AiParameterFLT("specular2RoughnessDepthScale", 1.0f);
   AiParameterINT("specular2ExtraSamples", 0);
   AiParameterVec("specular2Normal", 0, 0, 0);
   AiParameterFLT("specular2IndirectStrength", 1.0f);
   AiParameterFLT("specular2IndirectClamp", 0.0f);
   AiParameterBOOL("specular2CausticPaths", false);
   AiParameterBOOL("specular2InternalDirect", true);
   AiParameterENUM("specular2Distribution", AI_MICROFACET_BECKMANN,
                   distribution_names);

   AiParameterFLT("transmissionStrength", 0.0f);
   AiParameterRGB("transmissionColor", 1.0f, 1.0f, 1.0f);
   AiParameterBOOL("transmissionLinkToSpecular1", true);
   AiParameterFLT("transmissionRoughness", 0.f);
   AiParameterFLT("transmissionIor", 1.4f);
   AiParameterFLT("transmissionRoughnessDepthScale", 1.0f);
   AiParameterBOOL("transmissionEnableCaustics", true);
   AiParameterINT("transmissionExtraSamples", 0);
   AiParameterFLT("transmissionClamp", 0.0f);
   AiParameterBOOL("transmissionDoDirect", false);
   AiParameterVec("transmissionNormal", 0, 0, 0);
   AiParameterBOOL("transmissionCausticPaths", false);

   AiParameterRGB("id1", 0.0f, 0.0f, 0.0f);
   AiParameterRGB("id2", 0.0f, 0.0f, 0.0f);
   AiParameterRGB("id3", 0.0f, 0.0f, 0.0f);
   AiParameterRGB("id4", 0.0f, 0.0f, 0.0f);
   AiParameterRGB("id5", 0.0f, 0.0f, 0.0f);
   AiParameterRGB("id6", 0.0f, 0.0f, 0.0f);
   AiParameterRGB("id7", 0.0f, 0.0f, 0.0f);
   AiParameterRGB("id8", 0.0f, 0.0f, 0.0f);

   AiParameterBOOL("aiEnableMatte", false);
   AiParameterRGB("aiMatteColor", 0.0f, 0.0f, 0.0f);
   AiParameterFlt("aiMatteColorA", 0.0f);

   AiParameterStr("aov_diffuse_color", "diffuse_color");
   AiParameterStr("aov_direct_diffuse", "direct_diffuse");
   AiParameterStr("aov_direct_diffuse_raw", "direct_diffuse_raw");
   AiParameterStr("aov_indirect_diffuse", "indirect_diffuse");
   AiParameterStr("aov_indirect_diffuse_raw", "indirect_diffuse_raw");
   AiParameterStr("aov_direct_backlight", "direct_backlight");
   AiParameterStr("aov_indirect_backlight", "indirect_backlight");
   AiParameterStr("aov_direct_specular", "direct_specular");
   AiParameterStr("aov_indirect_specular", "indirect_specular");
   AiParameterStr("aov_direct_specular_2", "direct_specular_2");
   AiParameterStr("aov_indirect_specular_2", "indirect_specular_2");
   AiParameterStr("aov_single_scatter", "single_scatter");
   AiParameterStr("aov_sss", "sss");
   AiParameterStr("aov_refraction", "refraction");
   AiParameterStr("aov_emission", "emission");
   AiParameterStr("aov_uv", "uv");
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
   AiParameterStr("aov_crypto_asset", "crypto_asset");
   AiParameterStr("aov_crypto_object", "crypto_object");
   AiParameterStr("aov_crypto_material", "crypto_material");
   AiParameterStr("aov_shadow_group_1", "shadow_group_1");
   AiParameterStr("aov_shadow_group_2", "shadow_group_2");
   AiParameterStr("aov_shadow_group_3", "shadow_group_3");
   AiParameterStr("aov_shadow_group_4", "shadow_group_4");
   AiParameterStr("aov_shadow_group_5", "shadow_group_5");
   AiParameterStr("aov_shadow_group_6", "shadow_group_6");
   AiParameterStr("aov_shadow_group_7", "shadow_group_7");
   AiParameterStr("aov_shadow_group_8", "shadow_group_8");

   AiParameterStr("aov_shadow_group_9", "shadow_group_9");
   AiParameterStr("aov_shadow_group_10", "shadow_group_10");
   AiParameterStr("aov_shadow_group_11", "shadow_group_11");
   AiParameterStr("aov_shadow_group_12", "shadow_group_12");
   AiParameterStr("aov_shadow_group_13", "shadow_group_13");
   AiParameterStr("aov_shadow_group_14", "shadow_group_14");
   AiParameterStr("aov_shadow_group_15", "shadow_group_15");
   AiParameterStr("aov_shadow_group_16", "shadow_group_16");

   AiParameterBool("standardCompatibleAOVs", false);
   AiParameterBool("transmitAovs", false);

   AiParameterInt("rrTransmissionDepth", 1);

   AiParameterRGB("opacity", 1.0f, 1.0f, 1.0f);

   AiParameterBool("rr", true);

   AiParameterSTR("traceSetAll", "");
   AiParameterSTR("traceSetShadows", "");
   AiParameterSTR("traceSetDiffuse", "");
   AiParameterSTR("traceSetBacklight", "");
   AiParameterSTR("traceSetSpecular1", "");
   AiParameterSTR("traceSetSpecular2", "");
   AiParameterSTR("traceSetTransmission", "");

   AiParameterENUM("debug", DB_OFF, debugNames);

   AiParameterFlt("aov_diffuse_color_clamp", 0.0f);
   AiParameterFlt("aov_direct_diffuse_clamp", 0.0f);
   AiParameterFlt("aov_direct_diffuse_raw_clamp", 0.0f);
   AiParameterFlt("aov_indirect_diffuse_clamp", 0.0f);
   AiParameterFlt("aov_indirect_diffuse_raw_clamp", 0.0f);
   AiParameterFlt("aov_direct_backlight_clamp", 0.0f);
   AiParameterFlt("aov_indirect_backlight_clamp", 0.0f);
   AiParameterFlt("aov_direct_specular_clamp", 0.0f);
   AiParameterFlt("aov_indirect_specular_clamp", 0.0f);
   AiParameterFlt("aov_direct_specular_2_clamp", 0.0f);
   AiParameterFlt("aov_indirect_specular_2_clamp", 0.0f);
   AiParameterFlt("aov_single_scatter_clamp", 0.0f);
   AiParameterFlt("aov_sss_clamp", 0.0f);
   AiParameterFlt("aov_refraction_clamp", 0.0f);
   AiParameterFlt("aov_emission_clamp", 0.0f);
   AiParameterFlt("aov_light_group_1_clamp", 0.0f);
   AiParameterFlt("aov_light_group_2_clamp", 0.0f);
   AiParameterFlt("aov_light_group_3_clamp", 0.0f);
   AiParameterFlt("aov_light_group_4_clamp", 0.0f);
   AiParameterFlt("aov_light_group_5_clamp", 0.0f);
   AiParameterFlt("aov_light_group_6_clamp", 0.0f);
   AiParameterFlt("aov_light_group_7_clamp", 0.0f);
   AiParameterFlt("aov_light_group_8_clamp", 0.0f);
   AiParameterFlt("aov_light_group_9_clamp", 0.0f);
   AiParameterFlt("aov_light_group_10_clamp", 0.0f);
   AiParameterFlt("aov_light_group_11_clamp", 0.0f);
   AiParameterFlt("aov_light_group_12_clamp", 0.0f);
   AiParameterFlt("aov_light_group_13_clamp", 0.0f);
   AiParameterFlt("aov_light_group_14_clamp", 0.0f);
   AiParameterFlt("aov_light_group_15_clamp", 0.0f);
   AiParameterFlt("aov_light_group_16_clamp", 0.0f);

   AiParameterStr("crypto_asset_override", "");
   AiParameterStr("crypto_object_override", "");
   AiParameterStr("crypto_material_override", "");
}

#ifdef MSVC
#define _CRT_SECURE_NO_WARNINGS 1
#endif
node_loader
{
   if (i > 0) return 0;
   node->methods = alSurfaceMtd;
   node->output_type = AI_TYPE_RGB;
   node->name = "alSurface";
   node->node_type = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return true;
}

node_initialize
{
   ShaderData* data = new ShaderData();
   AiNodeSetLocalData(node, data);
   data->diffuse_sampler = NULL;
   data->sss_sampler = NULL;
   data->glossy_sampler = NULL;
   data->glossy2_sampler = NULL;
   data->refraction_sampler = NULL;
   data->backlight_sampler = NULL;
   data->perm_table = NULL;
   data->perm_table_diffuse = NULL;
   data->perm_table_spec1 = NULL;
   data->perm_table_spec2 = NULL;
   data->perm_table_backlight = NULL;
   data->perm_table_sss = NULL;

   data->cryptomatte = new CryptomatteData();
};

node_finish
{
   if (AiNodeGetLocalData(node))
   {
      ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
      AiSamplerDestroy(data->diffuse_sampler);
      AiSamplerDestroy(data->sss_sampler);
      AiSamplerDestroy(data->glossy_sampler);
      AiSamplerDestroy(data->glossy2_sampler);
      AiSamplerDestroy(data->refraction_sampler);
      AiSamplerDestroy(data->backlight_sampler);
      delete[] data -> perm_table;

      delete[] data -> perm_table_diffuse;
      delete[] data -> perm_table_spec1;
      delete[] data -> perm_table_spec2;
      delete[] data -> perm_table_backlight;
      delete[] data -> perm_table_sss;

      if (data->cryptomatte)
         delete data->cryptomatte;

      AiNodeSetLocalData(node, NULL);
      delete data;
   }
}

/// Generate a randomly permted [0,n) sequence
inline void permute(int* perm, int n)
{
   int i;
   for (i = 0; i < n; ++i)
   {
      perm[i] = i;
   }

   while (--i)
   {
      int tmp = perm[i];
      int rindex = rand0n(i);
      perm[i] = perm[rindex];
      perm[rindex] = tmp;
   }
}

node_update
{
   ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);

   // set up AOVs
   REGISTER_AOVS

   data->standardAovs = params[p_standardAovs].BOOL;
   data->transmitAovs = params[p_transmitAovs].BOOL;

   // store some options we'll reuse later
   AtNode* options = AiUniverseGetOptions();
   data->GI_diffuse_depth = AiNodeGetInt(options, "GI_diffuse_depth");
   data->GI_reflection_depth = AiNodeGetInt(options, "GI_reflection_depth");
   data->GI_refraction_depth = AiNodeGetInt(options, "GI_refraction_depth");
   data->GI_glossy_depth = AiNodeGetInt(options, "GI_glossy_depth");
   data->GI_glossy_samples = AiNodeGetInt(options, "GI_glossy_samples") +
                             params[p_specular1ExtraSamples].INT;
   data->glossy_samples2 = SQR(data->GI_glossy_samples);
   data->GI_glossy2_samples = AiNodeGetInt(options, "GI_glossy_samples") +
                              params[p_specular2ExtraSamples].INT;
   data->glossy2_samples2 = SQR(data->GI_glossy2_samples);
   data->GI_diffuse_samples = AiNodeGetInt(options, "GI_diffuse_samples") +
                              params[p_diffuseExtraSamples].INT;
   data->diffuse_samples2 = SQR(data->GI_diffuse_samples);
   data->GI_refraction_samples =
       AiNodeGetInt(options, "GI_refraction_samples") +
       params[p_transmissionExtraSamples].INT;
   data->refraction_samples2 = SQR(data->GI_refraction_samples);
   data->sss_bssrdf_samples = AiNodeGetInt(options, "sss_bssrdf_samples") +
                              params[p_sssExtraSamples].INT;
   data->sss_bssrdf_samples2 = SQR(data->sss_bssrdf_samples);

   // setup samples
   AiSamplerDestroy(data->diffuse_sampler);
   AiSamplerDestroy(data->sss_sampler);
   AiSamplerDestroy(data->glossy_sampler);
   AiSamplerDestroy(data->glossy2_sampler);
   AiSamplerDestroy(data->refraction_sampler);
   AiSamplerDestroy(data->backlight_sampler);
   data->diffuse_sampler = AiSampler(data->GI_diffuse_samples, 2);
   data->sss_sampler = AiSampler(data->sss_bssrdf_samples, 2);
   data->glossy_sampler = AiSampler(data->GI_glossy_samples, 2);
   data->glossy2_sampler = AiSampler(data->GI_glossy2_samples, 2);
   data->refraction_sampler = AiSampler(data->GI_refraction_samples, 2);
   data->backlight_sampler = AiSampler(data->GI_diffuse_samples, 2);

   // Get all the light nodes in the scene and try and find their light group
   // parameter
   // we'll store this based on the light pointer for fast access during
   // rendering
   AtNodeIterator* it = AiUniverseGetNodeIterator(AI_NODE_LIGHT);

   while (!AiNodeIteratorFinished(it))
   {
      AtNode* light = AiNodeIteratorGetNext(it);
      if (AiNodeLookUpUserParameter(light, "lightGroup"))
         data->lightGroups[light] = AiNodeGetInt(light, "lightGroup") - 1;
      else
         data->lightGroups[light] = -1;

      data->shadowDensities[light] = AiNodeGetFlt(light, "shadow_density");
   }
   AiNodeIteratorDestroy(it);

   // check whether the normal parameters are connected or not
   data->specular1NormalConnected = AiNodeIsLinked(node, "specular1Normal");
   data->specular2NormalConnected = AiNodeIsLinked(node, "specular2Normal");
   data->diffuseNormalConnected = AiNodeIsLinked(node, "diffuseNormal");
   data->transmissionNormalConnected =
       AiNodeIsLinked(node, "transmissionNormal");

   data->rrTransmissionDepth = params[p_rrTransmissionDepth].INT;

   data->specular1IndirectClamp = params[p_specular1IndirectClamp].FLT;
   if (data->specular1IndirectClamp == 0.0f)
      data->specular1IndirectClamp = AI_INFINITE;
   data->specular2IndirectClamp = params[p_specular2IndirectClamp].FLT;
   if (data->specular2IndirectClamp == 0.0f)
      data->specular2IndirectClamp = AI_INFINITE;
   data->transmissionClamp = params[p_transmissionClamp].FLT;
   if (data->transmissionClamp == 0.0f) data->transmissionClamp = AI_INFINITE;
   data->diffuseIndirectClamp = params[p_diffuseIndirectClamp].FLT;
   if (data->diffuseIndirectClamp == 0.0f)
      data->diffuseIndirectClamp = AI_INFINITE;

   data->transmissionDoDirect = params[p_transmissionDoDirect].BOOL;

   // Set up info for RR
   data->do_rr = params[p_rr].BOOL;
   data->AA_samples = SQR(std::max(1, AiNodeGetInt(options, "AA_samples")));
   data->AA_samples_inv = 1.0f / float(data->AA_samples);

   data->total_depth = AiNodeGetInt(options, "GI_total_depth") + 1;
   delete[] data -> perm_table;
   data->perm_table = new int[data->AA_samples * data->total_depth];
   delete[] data -> perm_table_diffuse;
   data->perm_table_diffuse =
       new int[data->AA_samples * data->diffuse_samples2 * data->total_depth];
   delete[] data -> perm_table_spec1;
   data->perm_table_spec1 =
       new int[data->AA_samples * data->glossy_samples2 * data->total_depth];
   delete[] data -> perm_table_spec2;
   data->perm_table_spec2 =
       new int[data->AA_samples * data->glossy2_samples2 * data->total_depth];
   delete[] data -> perm_table_backlight;
   data->perm_table_backlight =
       new int[data->AA_samples * data->diffuse_samples2 * data->total_depth];
   delete[] data -> perm_table_sss;
   data->perm_table_sss = new int[data->AA_samples * data->sss_bssrdf_samples2 *
                                  data->total_depth];
   // permute uses rand() to generate the random number stream so seed it first
   // so we get a determistic sequence between renders
   srand(RAND_STREAM_ALSURFACE_RR_PERMUTE);
   // generate the permutation table for RR;
   for (int d = 0; d < data->total_depth; ++d)
   {
      permute(&(data->perm_table[d * data->AA_samples]), data->AA_samples);
   }

   srand(RAND_STREAM_ALSURFACE_RR_DIFF_PERMUTE);
   // generate the permutation table for rr_diffuse
   if (data->diffuse_samples2 > 0)
   {
      for (int d = 0; d < data->total_depth; ++d)
      {
         permute(&(data->perm_table_diffuse[d * data->AA_samples *
                                            data->diffuse_samples2]),
                 data->AA_samples * data->diffuse_samples2);
      }
   }

   srand(RAND_STREAM_ALSURFACE_RR_SPEC1_PERMUTE);
   // generate the permutation table for rr_spec1
   if (data->glossy_samples2 > 0)
   {
      for (int d = 0; d < data->total_depth; ++d)
      {
         permute(&(data->perm_table_spec1[d * data->AA_samples *
                                          data->glossy_samples2]),
                 data->AA_samples * data->glossy_samples2);
      }
   }

   srand(RAND_STREAM_ALSURFACE_RR_SPEC2_PERMUTE);
   // generate the permutation table for rr_spec2
   if (data->glossy2_samples2 > 0)
   {
      for (int d = 0; d < data->total_depth; ++d)
      {
         permute(&(data->perm_table_spec2[d * data->AA_samples *
                                          data->glossy2_samples2]),
                 data->AA_samples * data->glossy2_samples2);
      }
   }

   srand(RAND_STREAM_ALSURFACE_RR_BACKLIGHT_PERMUTE);
   // generate the permutation table for rr_backlight
   if (data->diffuse_samples2 > 0)
   {
      for (int d = 0; d < data->total_depth; ++d)
      {
         permute(&(data->perm_table_backlight[d * data->AA_samples *
                                              data->diffuse_samples2]),
                 data->AA_samples * data->diffuse_samples2);
      }
   }

   srand(RAND_STREAM_ALSURFACE_RR_SSS_PERMUTE);
   // generate the permutation table for rr_backlight
   if (data->sss_bssrdf_samples2 > 0)
   {
      for (int d = 0; d < data->total_depth; ++d)
      {
         permute(&(data->perm_table_sss[d * data->AA_samples *
                                        data->sss_bssrdf_samples2]),
                 data->AA_samples * data->sss_bssrdf_samples2);
      }
   }

   data->xres = AiNodeGetInt(options, "xres");

   // fresnel
   data->specular1FresnelMode = params[p_specular1FresnelMode].INT;
   data->specular2FresnelMode = params[p_specular2FresnelMode].INT;

   // Trace sets setup
   data->trace_set_all_enabled = false;
   data->trace_set_shadows_enabled = false;
   data->trace_set_diffuse_enabled = false;
   data->trace_set_specular1_enabled = false;
   data->trace_set_specular2_enabled = false;
   data->trace_set_transmission_enabled = false;

   if (strlen(params[p_trace_set_all].STR))
   {
      std::string tmp(params[p_trace_set_all].STR);
      data->trace_set_all_enabled = true;
      if (tmp[0] == '-')
      {
         data->trace_set_all_inclusive = false;
         data->trace_set_all = tmp.substr(1);
      }
      else
      {
         data->trace_set_all_inclusive = true;
         data->trace_set_all = tmp;
      }
   }

   if (strlen(params[p_trace_set_shadows].STR))
   {
      std::string tmp(params[p_trace_set_shadows].STR);
      data->trace_set_shadows_enabled = true;
      if (tmp[0] == '-')
      {
         data->trace_set_shadows_inclusive = false;
         data->trace_set_shadows = tmp.substr(1);
      }
      else
      {
         data->trace_set_shadows_inclusive = true;
         data->trace_set_shadows = tmp;
      }
   }

   if (strlen(params[p_trace_set_diffuse].STR))
   {
      std::string tmp(params[p_trace_set_diffuse].STR);
      data->trace_set_diffuse_enabled = true;
      if (tmp[0] == '-')
      {
         data->trace_set_diffuse_inclusive = false;
         data->trace_set_diffuse = tmp.substr(1);
      }
      else
      {
         data->trace_set_diffuse_inclusive = true;
         data->trace_set_diffuse = tmp;
      }
   }

   if (strlen(params[p_trace_set_backlight].STR))
   {
      std::string tmp(params[p_trace_set_backlight].STR);
      data->trace_set_backlight_enabled = true;
      if (tmp[0] == '-')
      {
         data->trace_set_backlight_inclusive = false;
         data->trace_set_backlight = tmp.substr(1);
      }
      else
      {
         data->trace_set_backlight_inclusive = true;
         data->trace_set_backlight = tmp;
      }
   }

   if (strlen(params[p_trace_set_specular1].STR))
   {
      std::string tmp(params[p_trace_set_specular1].STR);
      data->trace_set_specular1_enabled = true;
      if (tmp[0] == '-')
      {
         data->trace_set_specular1_inclusive = false;
         data->trace_set_specular1 = tmp.substr(1);
      }
      else
      {
         data->trace_set_specular1_inclusive = true;
         data->trace_set_specular1 = tmp;
      }
   }

   if (strlen(params[p_trace_set_specular2].STR))
   {
      std::string tmp(params[p_trace_set_specular2].STR);
      data->trace_set_specular2_enabled = true;
      if (tmp[0] == '-')
      {
         data->trace_set_specular2_inclusive = false;
         data->trace_set_specular2 = tmp.substr(1);
      }
      else
      {
         data->trace_set_specular2_inclusive = true;
         data->trace_set_specular2 = tmp;
      }
   }

   if (strlen(params[p_trace_set_transmission].STR))
   {
      std::string tmp(params[p_trace_set_transmission].STR);
      data->trace_set_transmission_enabled = true;
      if (tmp[0] == '-')
      {
         data->trace_set_transmission_inclusive = false;
         data->trace_set_transmission = tmp.substr(1);
      }
      else
      {
         data->trace_set_transmission_inclusive = true;
         data->trace_set_transmission = tmp;
      }
   }

   if (strlen(params[p_sssTraceSet].STR))
   {
      std::string tmp(params[p_sssTraceSet].STR);
      data->trace_set_sss_enabled = true;
      if (tmp[0] == '-')
      {
         data->trace_set_sss_inclusive = false;
         data->trace_set_sss = tmp.substr(1);
      }
      else
      {
         data->trace_set_sss_inclusive = true;
         data->trace_set_sss = tmp;
      }
   }

   // check if we're connected to an alCel shader
   // TODO: no easy way to do this... we'll want to traverse all shading nodes
   // and
   // see if they're alCel and their input is connected to us
   // for now just assume we're connected
   data->cel_connected = true;

   data->sssMode = params[p_sssMode].INT;

   data->debug = params[p_debug].INT;

   // sample clamping values
   data->aov_diffuse_color_clamp = params[p_aov_diffuse_color_clamp].FLT;
   if (data->aov_diffuse_color_clamp == 0.0f)
      data->aov_diffuse_color_clamp = AI_INFINITE;
   data->aov_direct_diffuse_clamp = params[p_aov_direct_diffuse_clamp].FLT;
   if (data->aov_direct_diffuse_clamp == 0.0f)
      data->aov_direct_diffuse_clamp = AI_INFINITE;
   data->aov_direct_diffuse_raw_clamp =
       params[p_aov_direct_diffuse_raw_clamp].FLT;
   if (data->aov_direct_diffuse_raw_clamp == 0.0f)
      data->aov_direct_diffuse_raw_clamp = AI_INFINITE;
   data->aov_indirect_diffuse_clamp = params[p_aov_indirect_diffuse_clamp].FLT;
   if (data->aov_indirect_diffuse_clamp == 0.0f)
      data->aov_indirect_diffuse_clamp = AI_INFINITE;
   data->aov_indirect_diffuse_raw_clamp =
       params[p_aov_indirect_diffuse_raw_clamp].FLT;
   if (data->aov_indirect_diffuse_raw_clamp == 0.0f)
      data->aov_indirect_diffuse_raw_clamp = AI_INFINITE;
   data->aov_direct_backlight_clamp = params[p_aov_direct_backlight_clamp].FLT;
   if (data->aov_direct_backlight_clamp == 0.0f)
      data->aov_direct_backlight_clamp = AI_INFINITE;
   data->aov_indirect_backlight_clamp =
       params[p_aov_indirect_backlight_clamp].FLT;
   if (data->aov_indirect_backlight_clamp == 0.0f)
      data->aov_indirect_backlight_clamp = AI_INFINITE;
   data->aov_direct_specular_clamp = params[p_aov_direct_specular_clamp].FLT;
   if (data->aov_direct_specular_clamp == 0.0f)
      data->aov_direct_specular_clamp = AI_INFINITE;
   data->aov_indirect_specular_clamp =
       params[p_aov_indirect_specular_clamp].FLT;
   if (data->aov_indirect_specular_clamp == 0.0f)
      data->aov_indirect_specular_clamp = AI_INFINITE;
   data->aov_direct_specular_2_clamp =
       params[p_aov_direct_specular_2_clamp].FLT;
   if (data->aov_direct_specular_2_clamp == 0.0f)
      data->aov_direct_specular_2_clamp = AI_INFINITE;
   data->aov_indirect_specular_2_clamp =
       params[p_aov_indirect_specular_2_clamp].FLT;
   if (data->aov_indirect_specular_2_clamp == 0.0f)
      data->aov_indirect_specular_2_clamp = AI_INFINITE;
   data->aov_single_scatter_clamp = params[p_aov_single_scatter_clamp].FLT;
   if (data->aov_single_scatter_clamp == 0.0f)
      data->aov_single_scatter_clamp = AI_INFINITE;
   data->aov_sss_clamp = params[p_aov_sss_clamp].FLT;
   if (data->aov_sss_clamp == 0.0f) data->aov_sss_clamp = AI_INFINITE;
   data->aov_refraction_clamp = params[p_aov_refraction_clamp].FLT;
   if (data->aov_refraction_clamp == 0.0f)
      data->aov_refraction_clamp = AI_INFINITE;
   data->aov_emission_clamp = params[p_aov_emission_clamp].FLT;
   if (data->aov_emission_clamp == 0.0f) data->aov_emission_clamp = AI_INFINITE;
   data->aov_light_group_clamp[0] = params[p_aov_light_group_1_clamp].FLT;
   if (data->aov_light_group_clamp[0] == 0.0f)
      data->aov_light_group_clamp[0] = AI_INFINITE;
   data->aov_light_group_clamp[1] = params[p_aov_light_group_2_clamp].FLT;
   if (data->aov_light_group_clamp[1] == 0.0f)
      data->aov_light_group_clamp[1] = AI_INFINITE;
   data->aov_light_group_clamp[2] = params[p_aov_light_group_3_clamp].FLT;
   if (data->aov_light_group_clamp[2] == 0.0f)
      data->aov_light_group_clamp[2] = AI_INFINITE;
   data->aov_light_group_clamp[3] = params[p_aov_light_group_4_clamp].FLT;
   if (data->aov_light_group_clamp[3] == 0.0f)
      data->aov_light_group_clamp[3] = AI_INFINITE;
   data->aov_light_group_clamp[4] = params[p_aov_light_group_5_clamp].FLT;
   if (data->aov_light_group_clamp[4] == 0.0f)
      data->aov_light_group_clamp[4] = AI_INFINITE;
   data->aov_light_group_clamp[5] = params[p_aov_light_group_6_clamp].FLT;
   if (data->aov_light_group_clamp[5] == 0.0f)
      data->aov_light_group_clamp[5] = AI_INFINITE;
   data->aov_light_group_clamp[6] = params[p_aov_light_group_7_clamp].FLT;
   if (data->aov_light_group_clamp[6] == 0.0f)
      data->aov_light_group_clamp[6] = AI_INFINITE;
   data->aov_light_group_clamp[7] = params[p_aov_light_group_8_clamp].FLT;
   if (data->aov_light_group_clamp[7] == 0.0f)
      data->aov_light_group_clamp[7] = AI_INFINITE;

   data->aov_light_group_clamp[8] = params[p_aov_light_group_9_clamp].FLT;
   if (data->aov_light_group_clamp[8] == 0.0f)
      data->aov_light_group_clamp[8] = AI_INFINITE;
   data->aov_light_group_clamp[9] = params[p_aov_light_group_10_clamp].FLT;
   if (data->aov_light_group_clamp[9] == 0.0f)
      data->aov_light_group_clamp[9] = AI_INFINITE;
   data->aov_light_group_clamp[10] = params[p_aov_light_group_11_clamp].FLT;
   if (data->aov_light_group_clamp[10] == 0.0f)
      data->aov_light_group_clamp[10] = AI_INFINITE;
   data->aov_light_group_clamp[11] = params[p_aov_light_group_12_clamp].FLT;
   if (data->aov_light_group_clamp[11] == 0.0f)
      data->aov_light_group_clamp[11] = AI_INFINITE;
   data->aov_light_group_clamp[12] = params[p_aov_light_group_13_clamp].FLT;
   if (data->aov_light_group_clamp[12] == 0.0f)
      data->aov_light_group_clamp[12] = AI_INFINITE;
   data->aov_light_group_clamp[13] = params[p_aov_light_group_14_clamp].FLT;
   if (data->aov_light_group_clamp[13] == 0.0f)
      data->aov_light_group_clamp[13] = AI_INFINITE;
   data->aov_light_group_clamp[14] = params[p_aov_light_group_15_clamp].FLT;
   if (data->aov_light_group_clamp[14] == 0.0f)
      data->aov_light_group_clamp[14] = AI_INFINITE;
   data->aov_light_group_clamp[15] = params[p_aov_light_group_16_clamp].FLT;
   if (data->aov_light_group_clamp[15] == 0.0f)
      data->aov_light_group_clamp[15] = AI_INFINITE;

   data->cryptomatte->setup_all(AiNodeGetStr(node, "aov_crypto_asset"), 
      AiNodeGetStr(node, "aov_crypto_object"), AiNodeGetStr(node, "aov_crypto_material"));

   // caustic flags
   data->specular1CausticPaths = params[p_specular1CausticPaths].BOOL;
   data->specular2CausticPaths = params[p_specular2CausticPaths].BOOL;
   data->transmissionCausticPaths = params[p_transmissionCausticPaths].BOOL;
   if (data->transmissionCausticPaths) data->transmissionDoDirect = true;
};

shader_evaluate
{
   ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);

   // do debug mode if requested
   if (data->debug > DB_OFF && data->debug < DB_END)
   {
      AtRGB result = AI_RGB_BLACK;

      switch (data->debug)
      {
      case DB_diffuseStrength:
         result = AiColor(AiShaderEvalParamFlt(p_diffuseStrength));
         break;
      case DB_diffuseColor:
         result = AiShaderEvalParamRGB(p_diffuseColor);
         break;
      case DB_diffuseRoughness:
         result = AiColor(AiShaderEvalParamFlt(p_diffuseRoughness));
         break;
      case DB_backlightStrength:
         result = AiColor(AiShaderEvalParamFlt(p_backlightStrength));
         break;
      case DB_backlightColor:
         result = AiShaderEvalParamRGB(p_backlightColor);
         break;
      case DB_sssMix:
         result = AiColor(AiShaderEvalParamFlt(p_sssMix));
         break;
      case DB_sssWeight1:
         result = AiColor(AiShaderEvalParamFlt(p_sssWeight1));
         break;
      case DB_sssColor1:
         result = AiShaderEvalParamRGB(p_sssRadiusColor);
         break;
      case DB_sssWeight2:
         result = AiColor(AiShaderEvalParamFlt(p_sssWeight2));
         break;
      case DB_sssColor2:
         result = AiShaderEvalParamRGB(p_sssRadiusColor2);
         break;
      case DB_sssWeight3:
         result = AiColor(AiShaderEvalParamFlt(p_sssWeight3));
         break;
      case DB_sssColor3:
         result = AiShaderEvalParamRGB(p_sssRadiusColor3);
         break;
      case DB_specular1Strength:
         result = AiColor(AiShaderEvalParamFlt(p_specular1Strength));
         break;
      case DB_specular1Color:
         result = AiShaderEvalParamRGB(p_specular1Color);
         break;
      case DB_specular1Roughness:
         result = AiColor(AiShaderEvalParamFlt(p_specular1Roughness));
         break;
      case DB_specular2Strength:
         result = AiColor(AiShaderEvalParamFlt(p_specular2Strength));
         break;
      case DB_specular2Color:
         result = AiShaderEvalParamRGB(p_specular2Color);
         break;
      case DB_specular2Roughness:
         result = AiColor(AiShaderEvalParamFlt(p_specular2Roughness));
         break;
      case DB_transmissionStrength:
         result = AiColor(AiShaderEvalParamFlt(p_transmissionStrength));
         break;
      case DB_transmissionColor:
         result = AiShaderEvalParamRGB(p_transmissionColor);
         break;
      case DB_transmissionRoughness:
         result = AiColor(AiShaderEvalParamFlt(p_transmissionRoughness));
         break;
      }

      sg->out.RGB = result;
      sg->out_opacity = AI_RGB_WHITE;
      return;
   }

   float ior, ior2;
   float eta, eta2;
   ior = AiShaderEvalParamFlt(p_specular1Ior);
   eta = 1.0f / ior;

   // slightly wasteful doing this here, btu it keeps the code simpler
   ior2 = AiShaderEvalParamFlt(p_specular2Ior);
   eta2 = 1.0f / ior2;

   float roughness = AiShaderEvalParamFlt(p_specular1Roughness);
   roughness *= roughness;
   float transmissionRoughness;
   float transmissionIor;
   bool transmissionLinkToSpecular1 =
       AiShaderEvalParamBool(p_transmissionLinkToSpecular1);
   if (transmissionLinkToSpecular1)
   {
      transmissionRoughness = roughness;
      transmissionIor = ior;
   }
   else
   {
      transmissionRoughness = AiShaderEvalParamFlt(p_transmissionRoughness);
      transmissionRoughness *= transmissionRoughness;
      transmissionIor = AiShaderEvalParamFlt(p_transmissionIor);
   }
   AtRGB transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) *
                             AiShaderEvalParamFlt(p_transmissionStrength);
   AtRGB diffuseColor = AiShaderEvalParamRGB(p_diffuseColor) *
                        AiShaderEvalParamFlt(p_diffuseStrength);
   AtRGB backlightColor = AiShaderEvalParamRGB(p_backlightColor) *
                          AiShaderEvalParamFlt(p_backlightStrength);

   // balance diffuse, transmission and backlight
   float t_sum =
       maxh(transmissionColor) + maxh(diffuseColor) + maxh(backlightColor);
   if (t_sum > 1.0f)
   {
      float t_sum_inv = 1.0f / t_sum;
      transmissionColor *= t_sum_inv;
      diffuseColor *= t_sum_inv;
      backlightColor *= t_sum_inv;
   }

   AtRGB ssScattering = AiShaderEvalParamRGB(p_ssScattering);
   AtRGB ssAbsorption = AiShaderEvalParamRGB(p_ssAbsorption);
   float ssDensityScale = AiShaderEvalParamFlt(p_ssDensityScale);
   float ssDirection = AiShaderEvalParamFlt(p_ssDirection);
   float ssInScatteringStrength =
       AiShaderEvalParamFlt(p_ssInScatteringStrength);
   AtRGB ssAttenuationColor = AiShaderEvalParamRGB(p_ssAttenuationColor);
   // if the attenuation color is anything but white, make sure we don't have
   // some channels left at pure white
   if (minh(ssAttenuationColor) < 1.0f)
   {
      ssAttenuationColor = clamp(ssAttenuationColor, AI_RGB_BLACK, rgb(0.999f));
   }
   bool ssSpecifyCoefficients = AiShaderEvalParamBool(p_ssSpecifyCoefficients);

   AtRGB opacity = AiShaderEvalParamRGB(p_opacity);

   // precalculate scattering coefficients as we'll need them for shadows etc.
   AtRGB sigma_t = AI_RGB_BLACK;
   AtRGB sigma_s = AI_RGB_BLACK;
   AtRGB sigma_a = AI_RGB_BLACK;
   AtRGB sigma_s_prime = AI_RGB_BLACK;
   AtRGB sigma_t_prime = AI_RGB_BLACK;
   bool do_attenuation = false;
   bool do_scattering = false;
   if (minh(ssAttenuationColor) < 1 || ssInScatteringStrength > 0 ||
       ssSpecifyCoefficients)
   {
      if (ssSpecifyCoefficients)
      {
         sigma_s = ssScattering * ssDensityScale;
         sigma_a = ssAbsorption * ssDensityScale;
         sigma_t = sigma_s + sigma_a;
         if (maxh(sigma_t) > 0) do_attenuation = true;
         if (maxh(sigma_s) > 0) do_scattering = true;
      }
      else
      {
         sigma_a = (AI_RGB_WHITE - ssAttenuationColor) * ssDensityScale;
         // if we're doing a scattering medium, then make sure absorption is not
         // 0
         if (ssInScatteringStrength > 0) sigma_a = max(sigma_a, rgb(0.01f));
         sigma_s = ssAttenuationColor * ssInScatteringStrength * ssDensityScale;
         sigma_t = sigma_s + sigma_a;
         if (maxh(sigma_a) > 0) do_attenuation = true;
         if (maxh(sigma_s) > 0) do_scattering = true;
      }
      sigma_s_prime = sigma_s * (1.0f - ssDirection);
      sigma_t_prime = (sigma_s_prime + sigma_a);
   }

   // check custom ray type
   int als_raytype = ALS_RAY_UNDEFINED;
   AiStateGetMsgInt("als_raytype", &als_raytype);

   // build a local frame for sampling
   AtVector U, V;
   if (!AiV3isZero(sg->dPdu) && AiV3Exists(sg->dPdu))
   {
      // we have valid a valid dPdu derivative, construct V
      AtVector Utmp = AiV3Normalize(sg->dPdu);
      V = AiV3Normalize(AiV3Cross(sg->Nf, Utmp));
      U = AiV3Cross(V, sg->Nf);
   }
   else
   {
      AiBuildLocalFramePolar(&U, &V, &sg->Nf);
   }

   // if this is a camera ray, prepare the temporary storage for deep groups
   AtRGB* deepGroupPtr = NULL;
   AtRGB result_directGroup[NUM_LIGHT_GROUPS];
   for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
      result_directGroup[i] = AI_RGB_BLACK;
   bool doDeepGroups = !data->standardAovs;
   bool transmitAovs =
       data->transmitAovs;

   if (doDeepGroups && (sg->Rt & AI_RAY_CAMERA))
   {
      // if this is a camera ray allocate the group storage
      deepGroupPtr = (AtRGB*)AiShaderGlobalsQuickAlloc(
          sg, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
      memset(deepGroupPtr, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
      AiStateSetMsgPtr("als_deepGroupPtr", deepGroupPtr);
   }
   else if (doDeepGroups)
   {
      // secondary ray hit - get the pointer from the state
      // if the pointer hasn't been set we're being called from a BSDF that
      // doesn't have deep group support
      // so don't try and do it or we'll be in (crashy) trouble!
      if (!AiStateGetMsgPtr("als_deepGroupPtr", (void**)&deepGroupPtr))
         doDeepGroups = false;
   }

   // If we're in a camera ray, try and allocate the array for transmitting the AOVs, if transmitAovs is enabled
   AtRGB* transmittedAovPtr = NULL;
   if ((sg->Rt & AI_RAY_CAMERA))
   {
     if (transmitAovs)
     {
         transmittedAovPtr =
             (AtRGB*)AiShaderGlobalsQuickAlloc(sg, sizeof(AtRGB) * NUM_AOVs);
         memset(transmittedAovPtr, 0, sizeof(AtRGB) * NUM_AOVs);
         AiStateSetMsgPtr("als_transmittedAovPtr", transmittedAovPtr);
     }
   }
   else // otherwise, we're in a secondary ray so turn on transmitAovs if it has been turned on upstream
   {
      transmitAovs = AiStateGetMsgPtr("als_transmittedAovPtr",
                                      (void**)&transmittedAovPtr);
   }

   // get path throughput so far
   AtRGB path_throughput = AI_RGB_WHITE;
   if (data->do_rr && sg->Rr > 0)
      AiStateGetMsgRGB("als_throughput", &path_throughput);

   DirectionalMessageData* diffusion_msgdata = NULL;
   AiStateGetMsgPtr("als_dmd", (void**)&diffusion_msgdata);

   float sssMix = AiShaderEvalParamFlt(p_sssMix);
   if (als_raytype == ALS_RAY_SSS)
   {
      // compute the diffusion sample
      assert(diffusion_msgdata);
      alsIrradiateSample(
          sg, diffusion_msgdata, data->diffuse_sampler, U, V, data->lightGroups,
          path_throughput, data->trace_set_sss.c_str(),
          data->trace_set_sss_enabled, data->trace_set_sss_inclusive, sssMix);
      sg->out_opacity = AI_RGB_WHITE;
      return;
   }
   // if it's a shadow ray, handle shadow colouring through absorption
   // algorithm based heavily on the example in Kettle
   else if (sg->Rt & AI_RAY_SHADOW)
   {
      // if the object is transmissive and
      AtRGB outOpacity = AI_RGB_WHITE;
      if (maxh(transmissionColor))
      {
         // check transmission through the surface
         float costheta = AiV3Dot(sg->Nf, -sg->Rd);
         float kt = 0.0f;
         if (!data->specular1FresnelMode)
         {
            kt = 1.0f - fresnel(costheta, eta);
         }

         if (kt >= IMPORTANCE_EPS)  // else surface is fully reflective
         {
            if (maxh(sigma_t) > 0.0f)
            {
               AtPoint alsPreviousIntersection;
               AtRGB als_sigma_t = sigma_t_prime;
               if (AiStateGetMsgPnt("alsPreviousIntersection",
                                    &alsPreviousIntersection))
               {
                  AiStateGetMsgRGB("alsPrevious_sigma_t", &als_sigma_t);
                  bool doExtinction = false;
                  if (AiV3Dot(sg->N, sg->Rd) < 0.0f)
                  {
                     // ray is entering a closed volume
                     bool alsInside;
                     AiStateGetMsgBool("alsInside", &alsInside);
                     if (alsInside)
                     {
                        // ray is entering an embedded volume
                        doExtinction = true;
                     }
                     else
                     {
                        // shouldn't get here
                     }
                  }
                  else
                  {
                     // ray is exiting a closed volume
                     doExtinction = true;
                  }

                  if (doExtinction)
                  {
                     float z = AiV3Dist(sg->P, alsPreviousIntersection);
                     outOpacity.r = fast_exp(-z * als_sigma_t.r);
                     outOpacity.g = fast_exp(-z * als_sigma_t.g);
                     outOpacity.b = fast_exp(-z * als_sigma_t.b);
                     outOpacity = -log(outOpacity);
                  }

                  // 'add' the extra opacity due to transmission colour and
                  // fresnel
                  outOpacity =
                      1.0f - ((1.0f - outOpacity) * kt * transmissionColor);
               }
               else
               {
                  // first intersection
                  // tell the next shader invocation that we're now inside the
                  // surface and what our extinction
                  // coefficient is
                  AiStateSetMsgRGB("alsPrevious_sigma_t", sigma_t_prime);
                  AiStateSetMsgBool("alsInside", true);
                  outOpacity = 1.0f - kt * transmissionColor;
               }
            }
            else  // no extinction, shadows are fresnel only.
            {
               AiStateSetMsgRGB("alsPrevious_sigma_t", AI_RGB_BLACK);
               outOpacity = 1.0f - kt * transmissionColor;
            }
         }
      }

      // store intersection position
      AiStateSetMsgPnt("alsPreviousIntersection", sg->P);
      sg->out_opacity = outOpacity * opacity;
      return;
   }
   else if (!diffusion_msgdata)
   {
      // allocate diffusion sample storage
      diffusion_msgdata = (DirectionalMessageData*)AiShaderGlobalsQuickAlloc(
          sg, sizeof(DirectionalMessageData));
      // printf("dmd: %p\n", diffusion_msgdata);
      memset(diffusion_msgdata, 0, sizeof(DirectionalMessageData));
      AiStateSetMsgPtr("als_dmd", diffusion_msgdata);
   }

   if (maxh(opacity) < IMPORTANCE_EPS)
   {
      sg->out.RGB = AI_RGB_BLACK;
      sg->out_opacity = AI_RGB_BLACK;
   }
   sg->out_opacity = opacity;

   // do cryptomattes before early out for fully transparent but after opacity is set
   data->cryptomatte->do_cryptomattes(sg, node, p_crypto_asset_override, 
                                      p_crypto_object_override, 
                                      p_crypto_material_override);

   // early out if we're fully transparent or the object is matte
   if (AiColorIsZero(opacity) ||
       (AiShaderGlobalsIsObjectMatte(sg) && (sg->Rt & AI_RAY_CAMERA)))
      return;

   // Evaluate bump;
   // AtRGB bump = AiShaderEvalParamRGB(p_bump);

   // reset ray type just to be safe
   AiStateSetMsgInt("als_raytype", ALS_RAY_UNDEFINED);

   // Initialize parameter temporaries
   // TODO: reorganize this so we're not evaluating upstream when we don't need
   // the parameters, e.g. in shadow rays
   float diffuseRoughness = AiShaderEvalParamFlt(p_diffuseRoughness);
   AtRGB emissionColor = AiShaderEvalParamRGB(p_emissionColor) *
                         AiShaderEvalParamFlt(p_emissionStrength);

   AtRGB sssRadiusColor = AiShaderEvalParamRGB(p_sssRadiusColor);
   float sssRadius = AiShaderEvalParamFlt(p_sssRadius);
   float sssWeight1 = AiShaderEvalParamFlt(p_sssWeight1);
   AtRGB sssRadiusColor2 = AiShaderEvalParamRGB(p_sssRadiusColor2);
   float sssRadius2 = AiShaderEvalParamFlt(p_sssRadius2);
   float sssWeight2 = AiShaderEvalParamFlt(p_sssWeight2);
   AtRGB sssRadiusColor3 = AiShaderEvalParamRGB(p_sssRadiusColor3);
   float sssRadius3 = AiShaderEvalParamFlt(p_sssRadius3);
   float sssWeight3 = AiShaderEvalParamFlt(p_sssWeight3);
   float sssDensityScale = AiShaderEvalParamFlt(p_sssDensityScale);

   // normalize weights
   float normweight = 1.0f / (sssWeight1 + sssWeight2 + sssWeight3);
   sssWeight1 *= normweight;
   sssWeight2 *= normweight;
   sssWeight3 *= normweight;

   AtRGB specular1Color = AiShaderEvalParamRGB(p_specular1Color) *
                          AiShaderEvalParamFlt(p_specular1Strength);
   AtRGB specular2Color = AiShaderEvalParamRGB(p_specular2Color) *
                          AiShaderEvalParamFlt(p_specular2Strength);
   AtVector specular1Normal = sg->Nf;
   if (data->specular1NormalConnected)
   {
      specular1Normal = AiV3Normalize(AiShaderEvalParamVec(p_specular1Normal));
   }

   AtVector specular2Normal = sg->Nf;
   if (data->specular2NormalConnected)
   {
      specular2Normal = AiV3Normalize(AiShaderEvalParamVec(p_specular2Normal));
   }

   AtVector diffuseNormal = sg->Nf;
   if (data->diffuseNormalConnected)
   {
      diffuseNormal = AiV3Normalize(AiShaderEvalParamVec(p_diffuseNormal));
   }

   AtVector transmissionNormal = sg->Nf;
   if (data->transmissionNormalConnected)
   {
      transmissionNormal =
          AiV3Normalize(AiShaderEvalParamVec(p_transmissionNormal));
   }

   float roughness2 = AiShaderEvalParamFlt(p_specular2Roughness);
   roughness2 *= roughness2;

   float specular1RoughnessDepthScale =
       AiShaderEvalParamFlt(p_specular1RoughnessDepthScale);
   float specular2RoughnessDepthScale =
       AiShaderEvalParamFlt(p_specular2RoughnessDepthScale);
   float transmissionRoughnessDepthScale =
       AiShaderEvalParamFlt(p_transmissionRoughnessDepthScale);

   bool transmissionEnableCaustics =
       AiShaderEvalParamBool(p_transmissionEnableCaustics);

   // adapt the roughness values for anisotropy
   float specular1Anisotropy =
       CLAMP(AiShaderEvalParamFlt(p_specular1Anisotropy), 0.0f, 1.0f);
   float roughness_x = roughness;
   float roughness_y = roughness;
   if (specular1Anisotropy != 0.5f)
   {
      float aniso_t = SQR(fabsf(2.0f * specular1Anisotropy - 1.0f));
      roughness_x = (specular1Anisotropy >= 0.5f)
                        ? roughness
                        : LERP(aniso_t, MAX(roughness, 0.0025f), 1.0f);
      roughness_y = (specular1Anisotropy <= 0.5f)
                        ? roughness
                        : LERP(aniso_t, MAX(roughness, 0.0025f), 1.0f);
   }
   float specular2Anisotropy =
       CLAMP(AiShaderEvalParamFlt(p_specular2Anisotropy), 0.0f, 1.0f);
   float roughness2_x = roughness2;
   float roughness2_y = roughness2;
   if (specular2Anisotropy != 0.5f)
   {
      float aniso_t = SQR(fabsf(2.0f * specular2Anisotropy - 1.0f));
      roughness2_x = (specular2Anisotropy >= 0.5f)
                         ? roughness2
                         : LERP(aniso_t, MAX(roughness2, 0.0025f), 1.0f);
      roughness2_y = (specular2Anisotropy <= 0.5f)
                         ? roughness2
                         : LERP(aniso_t, MAX(roughness2, 0.0025f), 1.0f);
   }
   // Grab the roughness from the previous surface and make sure we're slightly
   // rougher than it to avoid glossy-glossy fireflies
   float alsPreviousRoughness = 0.0f;
   if (sg->Rr > 0)
   {
      AiStateGetMsgFlt("alsPreviousRoughness", &alsPreviousRoughness);
      if (data->do_rr)
      {
         roughness_x = std::max(
             roughness_x, alsPreviousRoughness * specular1RoughnessDepthScale);
         roughness_y = std::max(
             roughness_y, alsPreviousRoughness * specular1RoughnessDepthScale);
         roughness2_x = std::max(
             roughness2_x, alsPreviousRoughness * specular2RoughnessDepthScale);
         roughness2_y = std::max(
             roughness2_y, alsPreviousRoughness * specular2RoughnessDepthScale);
         transmissionRoughness =
             std::max(transmissionRoughness,
                      alsPreviousRoughness * transmissionRoughnessDepthScale);
      }
   }

   // clamp roughnesses
   // roughness = std::max(0.000001f, roughness);
   roughness2 = std::max(0.000001f, roughness2);
   // transmissionRoughness = std::max(0.000001f, transmissionRoughness);

   // rotated frames for anisotropy
   AtVector U1 = U, U2 = U;
   AtVector V1 = V, V2 = V;

   if (specular1Normal != sg->Nf)
   {
      if (!AiV3isZero(sg->dPdu) && AiV3Exists(sg->dPdu))
      {
         // we have valid a valid dPdu derivative, construct V
         AtVector Utmp = AiV3Normalize(sg->dPdu);
         V1 = AiV3Normalize(AiV3Cross(specular1Normal, Utmp));
         U1 = AiV3Cross(V1, specular1Normal);
      }
      else
      {
         AiBuildLocalFramePolar(&U1, &V1, &specular1Normal);
      }
   }

   if (specular2Normal != sg->Nf)
   {
      if (!AiV3isZero(sg->dPdu) && AiV3Exists(sg->dPdu))
      {
         // we have valid a valid dPdu derivative, construct V
         AtVector Utmp = AiV3Normalize(sg->dPdu);
         V2 = AiV3Normalize(AiV3Cross(specular2Normal, Utmp));
         U2 = AiV3Cross(V2, specular2Normal);
      }
      else
      {
         AiBuildLocalFramePolar(&U2, &V2, &specular2Normal);
      }
   }

   if (specular1Anisotropy != 0.5f)
   {
      float specular1Rotation = AiShaderEvalParamFlt(p_specular1Rotation);
      const float cos_phi = cosf(specular1Rotation * AI_PI);
      const float sin_phi = sinf(specular1Rotation * AI_PI);
      U1 = cos_phi * U1 - sin_phi * V1;
      V1 = cos_phi * V1 + sin_phi * U1;
   }
   if (specular2Anisotropy != 0.5f)
   {
      float specular2Rotation = AiShaderEvalParamFlt(p_specular2Rotation);
      const float cos_phi = cosf(specular2Rotation * AI_PI);
      const float sin_phi = sinf(specular2Rotation * AI_PI);
      U2 = cos_phi * U2 - sin_phi * V2;
      V2 = cos_phi * V2 + sin_phi * U2;
   }

   // View direction, omega_o
   AtVector wo = -sg->Rd;

   float diffuseIndirectStrength =
       AiShaderEvalParamFlt(p_diffuseIndirectStrength);
   float backlightIndirectStrength =
       AiShaderEvalParamFlt(p_backlightIndirectStrength);
   float specular1IndirectStrength =
       AiShaderEvalParamFlt(p_specular1IndirectStrength);
   float specular2IndirectStrength =
       AiShaderEvalParamFlt(p_specular2IndirectStrength);

   // Initialize result temporaries
   AtRGB result_diffuseDirect = AI_RGB_BLACK;
   AtRGB result_diffuseDirectRaw = AI_RGB_BLACK;
   AtRGB result_backlightDirect = AI_RGB_BLACK;
   AtRGB result_glossyDirect = AI_RGB_BLACK;
   AtRGB result_glossy2Direct = AI_RGB_BLACK;
   AtRGB result_diffuseIndirect = AI_RGB_BLACK;
   AtRGB result_diffuseIndirectRaw = AI_RGB_BLACK;
   AtRGB result_backlightIndirect = AI_RGB_BLACK;
   AtRGB result_glossyIndirect = AI_RGB_BLACK;
   AtRGB result_glossy2Indirect = AI_RGB_BLACK;
   AtRGB result_sss = AI_RGB_BLACK;
   AtRGB result_ss = AI_RGB_BLACK;
   AtColor result_transmission = AI_RGB_BLACK;
   AtColor result_transmissionDirect = AI_RGB_BLACK;
   AtColor result_emission = AI_RGB_BLACK;

   // Set up flags to early out of calculations based on where we are in the ray
   // tree
   bool do_diffuse = true;
   bool do_backlight = true;
   bool do_glossy = true;
   bool do_glossy2 = true;
   bool do_ss = true;
   bool do_sss = true;
   bool do_transmission = true;
   int glossy_samples = data->GI_glossy_samples;
   int diffuse_samples = data->GI_diffuse_samples;

   bool specular1InternalDirect =
       AiShaderEvalParamBool(p_specular1InternalDirect);
   bool specular2InternalDirect =
       AiShaderEvalParamBool(p_specular2InternalDirect);

   float dummy;
   if (sg->Rr_diff > 0 || sg->Rr_gloss > 1 || sssMix < 0.01f ||
       AiStateGetMsgFlt("als_hairNumIntersections", &dummy) ||
       als_raytype == ALS_RAY_HAIR)
   {
      do_sss = false;
      sssMix = 0.0f;
   }

   if (maxh(diffuseColor) <
           IMPORTANCE_EPS  // disable diffuse if contribution is small
       ||
       sssMix == 1.0f)  // disable diffuse if sss mix is at 100%
   {
      do_diffuse = false;
   }

   if (maxh(backlightColor) <
           IMPORTANCE_EPS  // disable backlight if contribution is small
       ||
       sssMix == 1.0f)  // disable backlight if sss mix is at 100%
   {
      do_backlight = false;
   }

   if ((sg->Rr_diff > 0 &&
        !data->specular1CausticPaths)  // disable glossy->diffuse caustics
       ||
       maxh(specular1Color) <
           IMPORTANCE_EPS  // disable glossy if contribution is small
       ||
       (sg->Rr_refr > 0 &&
        !transmissionEnableCaustics)  // disable glossy->transmitted caustics
       ||
       roughness > 1.0f || als_raytype == ALS_RAY_HAIR ||
       (alsPreviousRoughness > 0.0f &&
        roughness ==
            0.0f))  // kill glossy if roughness has been scaled up too far
   {
      do_glossy = false;
   }
   // if the user has chosen not to do internal direct reflections, turn off the
   // direct glossy calculation if we're backfacing
   bool do_glossy_direct =
       do_glossy && ((AiV3Dot(sg->Nf, sg->Ng) > 0) || specular1InternalDirect);

   if ((sg->Rr_diff > 0 &&
        !data->specular2CausticPaths)  // disable glossy->diffuse caustics
       ||
       maxh(specular2Color) <
           IMPORTANCE_EPS  // disable glossy2 if contribution is small
       ||
       (sg->Rr_refr > 0 &&
        !transmissionEnableCaustics)  // disable glossy->transmitted caustics
       ||
       roughness2 > 1.0f ||
       als_raytype ==
           ALS_RAY_HAIR)  // kill glossy if roughness has been scaled up too far
   {
      do_glossy2 = false;
   }
   // if the user has chosen not to do internal direct reflections, turn off the
   // direct glossy calculation if we're backfacing
   bool do_glossy2_direct =
       do_glossy2 && (AiV3Dot(sg->Nf, sg->Ng) > 0 || specular2InternalDirect);

   if ((sg->Rr_diff > 0 &&
        !data->transmissionCausticPaths)  // disable transmitted caustics
       ||
       maxh(transmissionColor) < IMPORTANCE_EPS ||
       als_raytype ==
           ALS_RAY_HAIR)  // disable transmission if contribution is small
   {
      do_transmission = false;
   }

   // prepare temporaries for light group calculation
   AtRGB lightGroupsDirect[NUM_LIGHT_GROUPS];
   memset(lightGroupsDirect, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);

   // Decide whether to calculate shadow groups or not.
   bool doShadowGroups = (sg->Rt & AI_RAY_CAMERA);

   // Accumulator for transmission integrated according to the specular1 brdf.
   // Will be used to attenuate diffuse,
   // glossy2, sss and transmission
   float kti = 1.0f;
   float kti2 = 1.0f;

   // storage for all deepgroup contributions
   AtRGB deepGroupsGlossy[NUM_LIGHT_GROUPS];
   AtRGB deepGroupsGlossy2[NUM_LIGHT_GROUPS];
   AtRGB deepGroupsDiffuse[NUM_LIGHT_GROUPS];
   AtRGB deepGroupsSss[NUM_LIGHT_GROUPS];
   AtRGB deepGroupsTransmission[NUM_LIGHT_GROUPS];
   AtRGB deepGroupsBacklight[NUM_LIGHT_GROUPS];
   memset(deepGroupsGlossy, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
   memset(deepGroupsGlossy2, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
   memset(deepGroupsDiffuse, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
   memset(deepGroupsSss, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
   memset(deepGroupsTransmission, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
   memset(deepGroupsBacklight, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
   int count = 0;

   // Begin illumination calculation
   int specular1Distribution = AiShaderEvalParamInt(p_specular1Distribution);
   int specular2Distribution = AiShaderEvalParamInt(p_specular2Distribution);

   // Create the BRDF data structures for MIS
   // {
   AtVector Nold = sg->N;
   AtVector Nfold = sg->Nf;
   sg->N = sg->Nf = specular1Normal;
   void* mis;
   mis = GlossyMISCreateData(sg, specular1Distribution, &U1, 0,
                             std::max(roughness_x, 1e-3f),
                             std::max(roughness_y, 1e-3f));
   BrdfData_wrap brdfw;
   brdfw.brdf_data = mis;
   brdfw.sg = sg;

   if (data->specular1FresnelMode == FM_METALLIC)
   {
      AtRGB r = AiShaderEvalParamRGB(p_specular1Reflectivity);
      AtRGB g = AiShaderEvalParamRGB(p_specular1EdgeTint);
      FresnelConductor* fr = (FresnelConductor*)AiShaderGlobalsQuickAlloc(
          sg, sizeof(FresnelConductor));
      new (fr) FresnelConductor(kCustom, r, g);
      brdfw.fr = fr;
   }
   else
   {
      FresnelDielectric* fr = (FresnelDielectric*)AiShaderGlobalsQuickAlloc(
          sg, sizeof(FresnelDielectric));
      new (fr) FresnelDielectric(eta);
      brdfw.fr = fr;
   }
   brdfw.eta = eta;
   brdfw.V = wo;
   brdfw.N = specular1Normal;
   brdfw.kr = 0.0f;

   sg->N = sg->Nf = specular2Normal;
   void* mis2;
   mis2 = GlossyMISCreateData(sg, specular2Distribution, &U2, 0,
                              std::max(roughness2_x, 1e-3f),
                              std::max(roughness2_y, 1e-3f));
   BrdfData_wrap brdfw2;
   brdfw2.brdf_data = mis2;
   brdfw2.sg = sg;

   if (data->specular2FresnelMode == FM_METALLIC)
   {
      AtRGB r = AiShaderEvalParamRGB(p_specular2Reflectivity);
      AtRGB g = AiShaderEvalParamRGB(p_specular2EdgeTint);
      FresnelConductor* fr = (FresnelConductor*)AiShaderGlobalsQuickAlloc(
          sg, sizeof(FresnelConductor));
      new (fr) FresnelConductor(kCustom, r, g);
      brdfw2.fr = fr;
   }
   else
   {
      FresnelDielectric* fr = (FresnelDielectric*)AiShaderGlobalsQuickAlloc(
          sg, sizeof(FresnelDielectric));
      new (fr) FresnelDielectric(eta2);
      brdfw2.fr = fr;
   }

   brdfw2.eta = eta2;
   brdfw2.V = wo;
   brdfw2.N = specular2Normal;
   brdfw2.kr = 0.0f;

   sg->N = Nold;
   sg->Nf = Nfold;

   //sg->Nf = sg->Nf = diffuseNormal;
   sg->Nf = diffuseNormal; // or should that be: sg->N = sg->Nf = diffuseNormal?
   void* dmis = AiOrenNayarMISCreateData(sg, diffuseRoughness);

   if (do_backlight) sg->fhemi = false;
   flipNormals(sg);
   void* bmis = AiOrenNayarMISCreateData(sg, diffuseRoughness);
   flipNormals(sg);
   sg->fhemi = true;

   sg->Nf = Nfold;
   sg->N = Nold;

   AtRGBA shadowGroups[NUM_LIGHT_GROUPS];
   memset(shadowGroups, 0, sizeof(AtRGBA) * NUM_LIGHT_GROUPS);

   // set the global trace set if it's defined
   // this will potentially be overriden by each component as we go through the
   // shader
   if (data->trace_set_all_enabled)
   {
      AiShaderGlobalsSetTraceSet(sg, data->trace_set_all.c_str(),
                                 data->trace_set_all_inclusive);
   }

   // set the shadows trace set if it's defined
   if (data->trace_set_shadows_enabled)
   {
      AiShaderGlobalsSetTraceSet(sg, data->trace_set_shadows.c_str(),
                                 data->trace_set_shadows_inclusive);
   }

   AtVector Nforig = sg->Nf;
   AtVector Norig = sg->N;

   if (do_glossy_direct)
   {
      AtRGB LspecularDirect = AI_RGB_BLACK;

      if (data->specular1NormalConnected)
      {
         sg->N = sg->Nf = specular1Normal;
      }

      AiLightsPrepare(sg);
      while (AiLightsGetSample(sg))
      {
         if (AiLightGetAffectSpecular(sg->Lp))
         {
            // get the group assigned to this light from the hash table using
            // the light's pointer
            int lightGroup = data->lightGroups[sg->Lp];

            // per-light specular and diffuse strength multipliers
            float specular_strength = AiLightGetSpecular(sg->Lp);

            brdfw.ibs = false;
            LspecularDirect =
                AiEvaluateLightSample(sg, &brdfw, GlossyMISSample_wrap,
                                      GlossyMISBRDF_wrap, GlossyMISPDF_wrap) *
                specular_strength * specular1Color;
            // if the light is assigned a valid group number, add this
            // sample's
            // contribution to that light group
            if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
            {
               lightGroupsDirect[lightGroup] += LspecularDirect;
            }
            // accumulate the result
            result_glossyDirect += LspecularDirect;
         }
      }

      sg->Nf = Nforig;
      sg->N = Norig;

      if (data->specular1NormalConnected) AiLightsResetCache(sg);

      if (brdfw.ns > 0.0f)
      {
         kti = 1.0f - (maxh(brdfw.kr_int) / brdfw.ns * maxh(specular1Color));
         assert(AiIsFinite(kti));
      }
   }

   if (do_glossy2_direct)
   {
      if (data->specular2NormalConnected)
      {
         sg->N = sg->Nf = specular2Normal;
      }
      AiLightsPrepare(sg);
      AtRGB Lspecular2Direct = AI_RGB_BLACK;
      while (AiLightsGetSample(sg))
      {
         if (AiLightGetAffectSpecular(sg->Lp))
         {
            // get the group assigned to this light from the hash table using
            // the light's pointer
            int lightGroup = data->lightGroups[sg->Lp];

            // per-light specular and diffuse strength multipliers
            float specular_strength = AiLightGetSpecular(sg->Lp);

            sg->Nf = specular2Normal;
            brdfw2.ibs = false;
            Lspecular2Direct =
                AiEvaluateLightSample(sg, &brdfw2, GlossyMISSample_wrap,
                                      GlossyMISBRDF_wrap, GlossyMISPDF_wrap) *
                kti * specular_strength * specular2Color;
            if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
            {
               lightGroupsDirect[lightGroup] += Lspecular2Direct;
            }
            result_glossy2Direct += Lspecular2Direct;
         }
      }
      sg->Nf = Nforig;
      sg->N = Norig;
      if (data->specular2NormalConnected) AiLightsResetCache(sg);
      if (brdfw2.ns > 0.0f)
         kti *= 1.0f - (maxh(brdfw2.kr_int) / brdfw2.ns * maxh(specular2Color));
   }

   if (do_diffuse)
   {
      if (data->diffuseNormalConnected)
      {
         sg->N = sg->Nf = diffuseNormal;
      }
      AiLightsPrepare(sg);
      AtRGB LdiffuseDirect = AI_RGB_BLACK;
      while (AiLightsGetSample(sg))
      {
         if (AiLightGetAffectDiffuse(sg->Lp))
         {
            // get the group assigned to this light from the hash table using
            // the light's pointer
            int lightGroup = data->lightGroups[sg->Lp];
            // get diffuse strength multiplier
            float diffuse_strength = AiLightGetDiffuse(sg->Lp);

            LdiffuseDirect =
                AiEvaluateLightSample(sg, dmis, AiOrenNayarMISSample,
                                      AiOrenNayarMISBRDF, AiOrenNayarMISPDF) *
                diffuse_strength * kti * (1 - sssMix);
            if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
            {
               lightGroupsDirect[lightGroup] += LdiffuseDirect * diffuseColor;
               shadowGroups[lightGroup].rgb() +=
                   sg->Liu * sg->we * std::max(0.0f, AiV3Dot(sg->Nf, sg->Ld)) *
                   sg->Lo * kti * diffuseColor * AI_ONEOVERPI;
               shadowGroups[lightGroup].a += maxh(sg->Lo) * sg->we;
            }
            result_diffuseDirect += LdiffuseDirect;
            assert(AiIsFinite(result_diffuseDirect));
         }
      }
      sg->Nf = Nforig;
      sg->N = Norig;
   }

   if (do_backlight)
   {
      sg->fhemi = false;
      flipNormals(sg);
      AiLightsPrepare(sg);
      AtRGB LbacklightDirect;
      while (AiLightsGetSample(sg))
      {
         if (AiLightGetAffectDiffuse(sg->Lp))
         {
            // get the group assigned to this light from the hash table using
            // the light's pointer
            int lightGroup = data->lightGroups[sg->Lp];
            float diffuse_strength = AiLightGetDiffuse(sg->Lp);

            LbacklightDirect =
                AiEvaluateLightSample(sg, bmis, AiOrenNayarMISSample,
                                      AiOrenNayarMISBRDF, AiOrenNayarMISPDF) *
                diffuse_strength * kti * backlightColor * (1 - sssMix);
            if (doDeepGroups || sg->Rt & AI_RAY_CAMERA)
            {
               if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
               {
                  lightGroupsDirect[lightGroup] += LbacklightDirect;
               }
            }
            result_backlightDirect += LbacklightDirect;
         }
      }
      flipNormals(sg);
      AiLightsResetCache(sg);
      sg->fhemi = true;
   }

   // only do direct transmission if it's enabled and the ray is leaving a
   // volume
   if (do_transmission && data->transmissionDoDirect &&
       AiV3Dot(sg->Rd, sg->N) > 0.0f)
   {
      sg->fhemi = false;
      AiLightsPrepare(sg);
      AtRGB LtransmissionDirect;
      float t_eta = transmissionIor;
      if (AiV3Dot(sg->N, sg->Rd) > 0.0f) t_eta = 1.0f / t_eta;
      MicrofacetTransmission* mft = MicrofacetTransmission::create(
          sg, transmissionRoughness, transmissionRoughness, t_eta, sg->Nf, U,
          V);
      while (AiLightsGetSample(sg))
      {
         if (AiLightGetAffectSpecular(sg->Lp))
         {
            // get the group assigned to this light from the hash table using
            // the light's pointer
            int lightGroup = data->lightGroups[sg->Lp];
            float specular_strength = AiLightGetSpecular(sg->Lp);

            LtransmissionDirect =
                AiEvaluateLightSample(sg, mft, MicrofacetTransmission::Sample,
                                      MicrofacetTransmission::BTDF,
                                      MicrofacetTransmission::PDF) *
                transmissionColor * specular_strength * kti;

            if (doDeepGroups || sg->Rt & AI_RAY_CAMERA)
            {
               if (lightGroup >= 0 && lightGroup < NUM_LIGHT_GROUPS)
               {
                  lightGroupsDirect[lightGroup] += LtransmissionDirect;
               }
            }
            result_transmissionDirect += LtransmissionDirect;
         }
      }
      AiLightsResetCache(sg);
      sg->fhemi = true;
   }

   // unset the shadows trace set
   if (data->trace_set_shadows_enabled)
   {
      // if we defined a global trace set, re-set this, otherwise, unset
      if (data->trace_set_all_enabled)
      {
         AiShaderGlobalsSetTraceSet(sg, data->trace_set_all.c_str(),
                                    data->trace_set_all_inclusive);
      }
      else
      {
         AiShaderGlobalsUnsetTraceSet(sg);
      }
   }

   // Multiply by the colors
   result_diffuseDirectRaw = result_diffuseDirect;
   result_diffuseDirect *= diffuseColor;

   // Sample BRDFS
   float samples[2];
   AtRay wi_ray;
   AtVector wi;
   AtScrSample scrs;
   AtVector H;
   float kr = 1, kt = 1;
   kti = 1.0f;
   // figure out whether to choose glossy or transmission for russian roulette
   // TODO: unify all the IOR calculations
   bool inside = false;
   if (AiV3Dot(sg->N, sg->Rd) > 0.0f) inside = true;

   float n1 = 1.0f;
   float n2 = ior;

   if (inside)
   {
      n1 = ior;
      n2 = 1.0f;
   }
   AiMakeRay(&wi_ray, AI_RAY_REFRACTED, &sg->P, NULL, AI_BIG, sg);
   bool tir = (!AiRefractRay(&wi_ray, &sg->Nf, n1, n2, sg)) && inside;
   bool rr_transmission =
       (do_glossy && do_transmission && sg->Rr >= data->rrTransmissionDepth &&
        !tir && roughness == 0.0f);
   bool rr_glossy = false;
   if (rr_transmission)
   {
      kr = fresnel(AiV3Dot(-sg->Rd, sg->Nf), eta);

      // get a permuted, stratified random number
      float u = (float(data->perm_table[sg->Rr * data->AA_samples + sg->si]) +
                 sampleTEAFloat(sg->Rr * data->AA_samples + sg->si,
                                TEA_STREAM_ALSURFACE_RR_JITTER)) *
                data->AA_samples_inv;
      // offset based on pixel
      float offset = sampleTEAFloat(sg->y * data->xres + sg->x,
                                    TEA_STREAM_ALSURFACE_RR_OFFSET);
      u = fmodf(u + offset, 1.0f);

      if (u < kr)
      {
         do_glossy &= true;
         do_transmission &= false;
      }
      else
      {
         do_glossy &= false;
         do_transmission &= true;
      }
   }

   bool rr_backlight = do_diffuse && do_backlight && sg->Rr_diff >= 2;
   float bd_prob = 1.0f;
   if (rr_backlight)
   {
      // get a permuted, stratified random number
      float u = (float(data->perm_table[sg->Rr * data->AA_samples + sg->si]) +
                 sampleTEAFloat(sg->Rr * data->AA_samples + sg->si,
                                TEA_STREAM_ALSURFACE_RR_JITTER)) *
                data->AA_samples_inv;
      // offset based on pixel
      float offset = sampleTEAFloat(sg->y * data->xres + sg->x,
                                    TEA_STREAM_ALSURFACE_RR2_OFFSET);
      u = fmodf(u + offset, 1.0f);

      float bd_sum = maxh(diffuseColor) + maxh(backlightColor);
      bd_prob = maxh(backlightColor) / bd_sum;
      if (u < bd_prob)
      {
         do_backlight = true;
         do_diffuse = false;
      }
      else
      {
         do_diffuse = true;
         do_backlight = false;
         bd_prob = 1.0f - bd_prob;
      }
      bd_prob = 1.0f / bd_prob;
   }

   // indirect_specular
   // -----------------
   if (do_glossy && specular1IndirectStrength > 0.0f)
   {
      // set the specular1 trace set if it's defined
      if (data->trace_set_specular1_enabled)
      {
         AiShaderGlobalsSetTraceSet(sg, data->trace_set_specular1.c_str(),
                                    data->trace_set_specular1_inclusive);
      }

      AtSamplerIterator* sampit = AiSamplerIterator(data->glossy_sampler, sg);
      // if we have perfect specular reflection, fall back to a single sample
      // along the reflection direction
      if (roughness == 0.0f)
      {
         if (!tir)
         {
            AiStateSetMsgFlt("alsPreviousRoughness",
                             std::max(roughness_x, roughness_y));
            sg->Nf = specular1Normal;
            AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
            AiReflectRay(&wi_ray, &sg->Nf, sg);
            AtRGB kr;
            if (!rr_transmission)
            {
               kr = brdfw.fr->kr(std::max(0.0f, AiV3Dot(wi_ray.dir, sg->Nf)));
               kti = maxh(kr);
            }
            else
            {
               kr = 1.0f;
               kti = 1.0f;
            }
            // Previously we pulled the sampler here as an optimization. This
            // nets us about a 10-30%
            // speedup in the case of pure dielectrics, but severely fucks up
            // sss, both on the surface
            // being cast, and in reflected surfaces.
            // Remove this for now until we can figure out exactly what's going
            // on
            // AiSamplerGetSample(sampit, samples);

            AtRGB f = kr * specular1Color * specular1IndirectStrength;
            bool cont = true;
            AtRGB throughput = path_throughput * f;
            float rr_p = 1.0f;
#ifdef RR_BOUNCES
            rr_p =
                std::min(1.0f, sqrtf(maxh(throughput) / maxh(path_throughput)));
            if (data->do_rr && sg->Rr > 0)
            {
               cont = false;
               // get a permuted, stratified random number
               int pb = data->AA_samples * data->glossy_samples2;
               int pi = data->glossy_samples2 * sg->si;
               float u =
                   (float(data->perm_table_spec1[sg->Rr * pb + pi]) +
                    sampleTEAFloat(sg->Rr * data->AA_samples + sg->si,
                                   TEA_STREAM_ALSURFACE_RR_SPEC1_JITTER)) /
                   float(pb);
               // offset based on pixel
               float offset =
                   sampleTEAFloat(sg->y * data->xres + sg->x,
                                  TEA_STREAM_ALSURFACE_RR_SPEC1_OFFSET);
               u = fmodf(u + offset, 1.0f);

               if (u < rr_p)
               {
                  cont = true;
                  rr_p = 1.0f / rr_p;
                  throughput *= rr_p;
                  f *= rr_p;
               }
            }
#endif
            if (cont)
            {
               AiStateSetMsgRGB("als_throughput", throughput);
               if (maxh(kr) > IMPORTANCE_EPS)
               {
                  bool hit = AiTrace(&wi_ray, &scrs);
                  result_glossyIndirect =
                      min(scrs.color * f, rgb(data->specular1IndirectClamp));
                  if (doDeepGroups && hit)
                  {

                     for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
                     {
                        deepGroupsGlossy[i] +=
                            min(deepGroupPtr[i] * f,
                                rgb(data->specular1IndirectClamp));
                     }
                  }
               }
            }
            sg->Nf = Nfold;
            kti = 1.0f - kti * maxh(specular1Color) * specular1IndirectStrength;
         }
      }
      else
      {
         AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
         kti = 0.0f;
         AiStateSetMsgFlt("alsPreviousRoughness",
                          std::max(roughness_x, roughness_y));
         sg->Nf = specular1Normal;
         AtRGB kr;
         int ssi = 0;
         while (AiSamplerGetSample(sampit, samples))
         {
            wi = GlossyMISSample(mis, float(samples[0]), float(samples[1]));
            if (AiV3Dot(wi, specular1Normal) > 0.0f)
            {
               // get half-angle vector for fresnel
               wi_ray.dir = wi;
               AiV3Normalize(H, wi + brdfw.V);
               kr = brdfw.fr->kr(std::max(0.0f, AiV3Dot(H, wi)));
               kti += maxh(kr);
               if (maxh(kr) >
                   IMPORTANCE_EPS)  // only trace a ray if it's going to matter
               {
                  AtRGB brdf = GlossyMISBRDF(mis, &wi);
                  float pdf = GlossyMISPDF(mis, &wi);
                  if (pdf <= 0.0f) continue;
                  AtRGB f = brdf / pdf * kr;
                  AtRGB throughput = path_throughput * f * specular1Color *
                                     specular1IndirectStrength;
                  bool cont = true;
                  float rr_p = 1.0f;
#ifdef RR_BOUNCES
                  rr_p = std::min(
                      1.0f, sqrtf(maxh(throughput) / maxh(path_throughput)));
                  if (data->do_rr && sg->Rr > 0)
                  {
                     cont = false;
                     // get a permuted, stratified random number
                     int pb = data->AA_samples * data->glossy_samples2;
                     int pi = data->glossy_samples2 * sg->si + ssi;
                     int idx = sg->Rr * pb + pi;
                     float u =
                         (float(data->perm_table_spec1[idx]) +
                          sampleTEAFloat(
                              idx, TEA_STREAM_ALSURFACE_RR_SPEC1_JITTER)) /
                         float(pb);
                     // offset based on pixel
                     float offset =
                         sampleTEAFloat(sg->y * data->xres + sg->x,
                                        TEA_STREAM_ALSURFACE_RR_SPEC1_OFFSET);
                     u = fmodf(u + offset, 1.0f);

                     if (u < rr_p)
                     {
                        cont = true;
                        rr_p = 1.0f / rr_p;
                        throughput *= rr_p;
                        f *= rr_p;
                     }
                  }
#endif
                  AiStateSetMsgRGB("als_throughput", throughput);
                  // if we're in a camera ray, pass the sample index down to the
                  // child SG
                  if (cont)
                  {
                     bool hit = AiTrace(&wi_ray, &scrs);
                     f *= specular1Color * specular1IndirectStrength;
                     result_glossyIndirect +=
                         min(scrs.color * f, rgb(data->specular1IndirectClamp));

                     // accumulate the lightgroup contributions calculated by
                     // the child shader
                     if (doDeepGroups && hit)
                     {
                        for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
                        {
                           deepGroupsGlossy[i] +=
                               min(deepGroupPtr[i] * f,
                                   rgb(data->specular1IndirectClamp));
                        }
                     }
                  }
               }
            }
            ssi++;
         }  // END while(samples)
         sg->Nf = Nfold;
         result_glossyIndirect *= AiSamplerGetSampleInvCount(sampit);
         kti *= AiSamplerGetSampleInvCount(sampit);
         kti = 1.0f - kti * maxh(specular1Color);

         if (doDeepGroups)
         {
            for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
            {

               deepGroupsGlossy[i] *= AiSamplerGetSampleInvCount(sampit);
            }
         }
      }

      // unset the specular1 trace set
      if (data->trace_set_specular1_enabled)
      {
         // if we defined a global trace set, re-set this, otherwise, unset
         if (data->trace_set_all_enabled)
         {
            AiShaderGlobalsSetTraceSet(sg, data->trace_set_all.c_str(),
                                       data->trace_set_all_inclusive);
         }
         else
         {
            AiShaderGlobalsUnsetTraceSet(sg);
         }
      }
   }  // if (do_glossy)

   // indirect_specular2
   // ------------------
   if (do_glossy2)
   {
      // set the specular2 trace set if it's defined
      if (data->trace_set_specular2_enabled)
      {
         AiShaderGlobalsSetTraceSet(sg, data->trace_set_specular2.c_str(),
                                    data->trace_set_specular2_inclusive);
      }

      AtSamplerIterator* sampit = AiSamplerIterator(data->glossy2_sampler, sg);
      AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
      kti2 = 0.0f;
      AtRGB kr;
      AiStateSetMsgFlt("alsPreviousRoughness",
                       std::max(roughness2_x, roughness2_y));
      sg->Nf = specular2Normal;
      int ssi = 0;
      while (AiSamplerGetSample(sampit, samples))
      {
         wi = GlossyMISSample(mis2, float(samples[0]), float(samples[1]));
         if (AiV3Dot(wi, specular2Normal) > 0.0f)
         {
            wi_ray.dir = wi;
            AiV3Normalize(H, wi + brdfw2.V);
            // add the fresnel for this layer
            kr = brdfw2.fr->kr(std::max(0.0f, AiV3Dot(H, wi)));
            kti2 += maxh(kr);
            AtRGB brdf = GlossyMISBRDF(mis2, &wi);
            float pdf = GlossyMISPDF(mis2, &wi);
            if (pdf <= 0.0f) continue;
            AtRGB f = brdf / pdf * kr * kti;
            AtRGB throughput = path_throughput * f * specular2Color *
                               specular2IndirectStrength;
            AiStateSetMsgRGB("als_throughput", throughput);
            bool cont = true;
            float rr_p = 1.0f;
#ifdef RR_BOUNCES
            rr_p =
                std::min(1.0f, sqrtf(maxh(throughput) / maxh(path_throughput)));
            if (data->do_rr && sg->Rr > 0)
            {
               cont = false;
               // get a permuted, stratified random number
               int pb = data->AA_samples * data->glossy2_samples2;
               int pi = data->glossy2_samples2 * sg->si + ssi;
               int idx = sg->Rr * pb + pi;
               float u =
                   (float(data->perm_table_spec2[idx]) +
                    sampleTEAFloat(idx, TEA_STREAM_ALSURFACE_RR_SPEC2_JITTER)) /
                   float(pb);
               // offset based on pixel
               float offset =
                   sampleTEAFloat(sg->y * data->xres + sg->x,
                                  TEA_STREAM_ALSURFACE_RR_SPEC2_OFFSET);
               u = fmodf(u + offset, 1.0f);

               if (u < rr_p)
               {
                  cont = true;
                  rr_p = 1.0f / rr_p;
                  throughput *= rr_p;
                  f *= rr_p;
               }
            }
#endif
            if (cont &&
                maxh(kr) >
                    IMPORTANCE_EPS)  // only trace a ray if it's going to matter
            {
               bool hit = AiTrace(&wi_ray, &scrs);

               f *= specular2Color * specular2IndirectStrength;
               result_glossy2Indirect +=
                   min(scrs.color * f, rgb(data->specular2IndirectClamp));

               // accumulate the lightgroup contributions calculated by the
               // child shader
               if (doDeepGroups && hit)
               {
                  for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
                  {
                     deepGroupsGlossy2[i] +=
                         min(deepGroupPtr[i] * f,
                             rgb(data->specular1IndirectClamp));
                  }
               }
            }
         }

         ssi++;
      }
      sg->Nf = Nfold;
      result_glossy2Indirect *= AiSamplerGetSampleInvCount(sampit);
      kti2 *= AiSamplerGetSampleInvCount(sampit);
      kti2 = 1.0f - kti2 * maxh(specular2Color);

      if (doDeepGroups)
      {
         for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
         {
            deepGroupsGlossy2[i] *= AiSamplerGetSampleInvCount(sampit);
         }
      }

      // unset the specular2 trace set
      if (data->trace_set_specular2_enabled)
      {
         // if we defined a global trace set, re-set this, otherwise, unset
         if (data->trace_set_all_enabled)
         {
            AiShaderGlobalsSetTraceSet(sg, data->trace_set_all.c_str(),
                                       data->trace_set_all_inclusive);
         }
         else
         {
            AiShaderGlobalsUnsetTraceSet(sg);
         }
      }

   }  // if (do_glossy2)

   // indirect_diffuse
   // ----------------
   if (do_diffuse &&
       kti * kti2 * maxh(diffuseColor) * diffuseIndirectStrength >
           IMPORTANCE_EPS)
   {
      // set the diffuse trace set if it's defined
      if (data->trace_set_diffuse_enabled)
      {
         AiShaderGlobalsSetTraceSet(sg, data->trace_set_diffuse.c_str(),
                                    data->trace_set_diffuse_inclusive);
      }

      float kr = kti * kti2;
      AtSamplerIterator* sampit = AiSamplerIterator(data->diffuse_sampler, sg);
      AiMakeRay(&wi_ray, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg);
      int ssi = 0;
      while (AiSamplerGetSample(sampit, samples))
      {
         // cosine hemisphere sampling as O-N sampling does not work outside of
         // a light loop
         float stheta = sqrtf(float(samples[0]));
         float phi = float(AI_PITIMES2 * samples[1]);
         wi.x = stheta * cosf(phi);
         wi.y = stheta * sinf(phi);
         wi.z = sqrtf(1.0f - float(samples[0]));
         AiV3RotateToFrame(wi, U, V, sg->Nf);

         float cos_theta = AiV3Dot(wi, sg->Nf);
         if (cos_theta <= 0.0f) continue;

         float p = cos_theta * float(AI_ONEOVERPI);

         // trace the ray
         wi_ray.dir = wi;
         AtRGB f = kr * AiOrenNayarMISBRDF(dmis, &wi) * bd_prob / p;
         AtRGB throughput =
             path_throughput * f * diffuseColor * diffuseIndirectStrength;
         float rr_p = 1.0f;
         bool cont = true;
#ifdef RR_BOUNCES
         if (data->do_rr && sg->Rr > 0)
         {
            cont = false;
            // get a permuted, stratified random number
            int pb = data->AA_samples * data->diffuse_samples2;
            int pi = data->diffuse_samples2 * sg->si + ssi;
            int idx = sg->Rr * pb + pi;
            float u =
                (float(data->perm_table_diffuse[idx]) +
                 sampleTEAFloat(idx, TEA_STREAM_ALSURFACE_RR_DIFF_JITTER)) /
                float(pb);
            // offset based on pixel
            float offset = sampleTEAFloat(sg->y * data->xres + sg->x,
                                          TEA_STREAM_ALSURFACE_RR_DIFF_OFFSET);
            u = fmodf(u + offset, 1.0f);

            rr_p =
                std::min(1.0f, sqrtf(maxh(throughput) / maxh(path_throughput)));
            if (u < rr_p)
            {
               cont = true;
               rr_p = 1.0f / rr_p;
               throughput *= rr_p;
               f *= rr_p;
            }
         }
#endif
         if (cont)
         {
            AiStateSetMsgRGB("als_throughput", throughput);
            bool hit = AiTrace(&wi_ray, &scrs);

            result_diffuseIndirectRaw +=
                min(scrs.color * f, rgb(data->diffuseIndirectClamp));

            // accumulate the lightgroup contributions calculated by the child
            // shader
            if (doDeepGroups && hit)
            {
               for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
               {
                  deepGroupsDiffuse[i] +=
                      min(deepGroupPtr[i] * f, rgb(data->diffuseIndirectClamp));
               }
            }
         }

         ssi++;
      }
      float invns = AiSamplerGetSampleInvCount(sampit);
      result_diffuseIndirectRaw *=
          invns * diffuseIndirectStrength * (1 - sssMix);
      result_diffuseIndirect = result_diffuseIndirectRaw * diffuseColor;

      if (doDeepGroups)
      {
         for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
         {
            deepGroupsDiffuse[i] *=
                invns * diffuseColor * diffuseIndirectStrength * (1 - sssMix);
         }
      }

      // unset the diffuse trace set
      if (data->trace_set_diffuse_enabled)
      {
         // if we defined a global trace set, re-set this, otherwise, unset
         if (data->trace_set_all_enabled)
         {
            AiShaderGlobalsSetTraceSet(sg, data->trace_set_all.c_str(),
                                       data->trace_set_all_inclusive);
         }
         else
         {
            AiShaderGlobalsUnsetTraceSet(sg);
         }
      }

   }  // if (do_diffuse)

   // refraction
   // ----------
   AtRGB childAovs[NUM_AOVs];
   memset(childAovs, 0, sizeof(AtRGB) * NUM_AOVs);
   if (do_transmission)
   {
      float samples[2];
      float kt;
      AtRay wi_ray;
      if (data->transmissionNormalConnected)
      {
          if (AiV3Dot(sg->N, transmissionNormal) < 0.0f)
            sg->N = -transmissionNormal;
          else
            sg->N = transmissionNormal;

          if (AiV3Dot(sg->Nf, transmissionNormal) < 0.0f)
            sg->Nf = -transmissionNormal;
          else
            sg->Nf = transmissionNormal;
      }
      AiMakeRay(&wi_ray, AI_RAY_REFRACTED, &sg->P, NULL, AI_BIG, sg);
      AtVector wi, R;
      AtScrSample sample;

      AtRGB mfp = AI_RGB_WHITE / sigma_t_prime;

      // set the transmission trace set if it's defined
      if (data->trace_set_transmission_enabled)
      {
         AiShaderGlobalsSetTraceSet(sg, data->trace_set_transmission.c_str(),
                                    data->trace_set_transmission_inclusive);
      }

      float inv_ns = 1.0f;
      AtSamplerIterator* sampit =
          AiSamplerIterator(data->refraction_sampler, sg);
      if (transmissionRoughness == 0.0f)
      {
         if (!rr_transmission)
         {
            kt = 1.0f - fresnel(transmissionIor, sg->N, wo, R, wi, inside);
         }
         else
         {
            kt = 1.0f;
         }
         float n1, n2;
         if (inside)
         {
            n1 = transmissionIor;
            n2 = 1.0f;
         }
         else
         {
            n1 = 1.0f;
            n2 = transmissionIor;
         }
         bool refraction = AiRefractRay(&wi_ray, &sg->Nf, n1, n2, sg);
         if (refraction)
         {
            // commented this sampler pull out to match the spec. This will make
            // stuff slower, but less noisy. Need more testing before we
            // release.
            // AiSamplerGetSample(sampit, samples);
            AtRGB throughput = path_throughput * kti;
            AiStateSetMsgRGB("als_throughput", throughput);

            if (sg->Rr_refr < data->GI_refraction_depth)
            {
               if (kti * kti2 > IMPORTANCE_EPS)
               {
                  AiStateSetMsgFlt("alsPreviousRoughness", 0.0f);
                  AiTrace(&wi_ray, &sample);
                  AtRGB transmittance = AI_RGB_WHITE;
                  bool hit = !AiColorIsZero(sample.opacity);

                  if (maxh(sigma_t) > 0.0f && !inside)
                  {
                     transmittance.r = fast_exp(float(-sample.z) * sigma_t.r);
                     transmittance.g = fast_exp(float(-sample.z) * sigma_t.g);
                     transmittance.b = fast_exp(float(-sample.z) * sigma_t.b);
                  }

                  AtRGB f = transmittance;
                  result_transmission +=
                      min(sample.color * f, rgb(data->transmissionClamp));
                  assert(AiIsFinite(result_transmission));
                  // accumulate the lightgroup contributions calculated by the
                  // child shader
                  if (doDeepGroups && hit)
                  {
                     for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
                     {
                        deepGroupsTransmission[i] += min(
                            deepGroupPtr[i] * f, rgb(data->transmissionClamp));
                     }
                  }

                  if (transmitAovs)
                  {
                     for (int i = 0; i < NUM_AOVs; ++i)
                     {
                        childAovs[i] += transmittedAovPtr[i] * f;
                     }
                  }

                  // single scattering
                  if (do_attenuation && maxh(sigma_s_prime) > 0.0f && !inside &&
                      do_scattering)
                  {
                     result_ss += AiSSSTraceSingleScatter(
                         sg, bssrdfbrdf(sigma_s_prime / sigma_t_prime), mfp,
                         ssDirection, transmissionIor);
                  }

                  if (minh(sample.opacity) < 1.0f)
                  {
                     AiTraceBackground(&wi_ray, &sample);
                     result_transmission +=
                         min(sample.color, rgb(data->transmissionClamp)) *
                         (AI_RGB_WHITE - sample.opacity);
                     assert(AiIsFinite(result_transmission));
                  }
               }
            }
            else  // trace the background if we've hit the refraction depth
                  // limit
            {
               AiTraceBackground(&wi_ray, &sample);
               result_transmission +=
                   min(sample.color, rgb(data->transmissionClamp));
            }
         }
         else  // total internal reflection
         {
            if (sg->Rr_refr < data->GI_refraction_depth)
            {
                AtRGB throughput = path_throughput * kti;
                AiStateSetMsgRGB("als_throughput", throughput);
                bool hit = AiTrace(&wi_ray, &sample);

                AtRGB transmittance = AI_RGB_WHITE;
                if (maxh(sigma_t) > 0.0f && !inside)
                {
                  transmittance.r = fast_exp(float(-sample.z) * sigma_t.r);
                  transmittance.g = fast_exp(float(-sample.z) * sigma_t.g);
                  transmittance.b = fast_exp(float(-sample.z) * sigma_t.b);
                }
                result_transmission +=
                    min(sample.color * transmittance, rgb(data->transmissionClamp));
                assert(AiIsFinite(result_transmission));
                // accumulate the lightgroup contributions calculated by the child
                // shader
                if (doDeepGroups && hit)
                {
                  for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
                  {
                      deepGroupsTransmission[i] +=
                          min(deepGroupPtr[i] * transmittance,
                              rgb(data->transmissionClamp));
                  }
                }

                if (transmitAovs)
                {
                  for (int i = 0; i < NUM_AOVs; ++i)
                  {
                      childAovs[i] += transmittedAovPtr[i] * transmittance;
                  }
                }
            }
            else  // trace the background if we've hit the refraction depth
                  // limit
              {
                AiTraceBackground(&wi_ray, &sample);
                result_transmission +=
                  min(sample.color, rgb(data->transmissionClamp));
              }
         }
      }
      else
      {

         float t_eta = transmissionIor;
         tir = false;
         if (AiV3Dot(sg->N, sg->Rd) > 0.0f) t_eta = 1.0f / t_eta;
         AtVector U, V;
         AiBuildLocalFramePolar(&U, &V, &sg->N);
         MicrofacetTransmission* mft = MicrofacetTransmission::create(
             sg, transmissionRoughness, transmissionRoughness, t_eta, sg->Nf, U,
             V);

         while (AiSamplerGetSample(sampit, samples))
         {
            AtVector m = mft->sampleMicrofacetNormal(samples[0], samples[1]);
            bool refracted = refraction(sg->Rd, m, t_eta, wi);
            wi_ray.dir = wi;
            AtRGB brdf;
            float pdf;
            if (refracted)
            {
               brdf = mft->btdf(wi_ray.dir);
               pdf = mft->pdf(wi_ray.dir);
            }
            else
            {
               brdf = AI_RGB_WHITE;
               pdf = 1.0f;
            }

            if (AiIsFinite(wi_ray.dir) && pdf > 1e-7f)
            {

               AtRGB f = brdf / pdf;

               AtRGB throughput = path_throughput * kti * f;
               AiStateSetMsgRGB("als_throughput", throughput);
               if (sg->Rr_refr < data->GI_refraction_depth)
               {
                  AiStateSetMsgFlt("alsPreviousRoughness",
                                   transmissionRoughness);
                  AiTrace(&wi_ray, &sample);
                  bool hit = !AiColorIsZero(sample.opacity);
                  AtRGB transmittance = AI_RGB_WHITE;
                  if (maxh(sigma_t) > 0.0f && !inside)
                  {
                     transmittance.r = fast_exp(float(-sample.z) * sigma_t.r);
                     transmittance.g = fast_exp(float(-sample.z) * sigma_t.g);
                     transmittance.b = fast_exp(float(-sample.z) * sigma_t.b);
                     f *= transmittance;
                  }

                  result_transmission +=
                      min(sample.color * f, rgb(data->transmissionClamp));
                  assert(AiIsFinite(result_transmission));
                  // accumulate the lightgroup contributions calculated by the
                  // child shader
                  if (doDeepGroups && hit)
                  {
                     for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
                     {
                        deepGroupsTransmission[i] += min(
                            deepGroupPtr[i] * f, rgb(data->transmissionClamp));
                     }
                  }

                  if (transmitAovs)
                  {
                     for (int i = 0; i < NUM_AOVs; ++i)
                     {
                        childAovs[i] += transmittedAovPtr[i] * f;
                     }
                  }

                  // single scattering
                  if (do_attenuation && maxh(sigma_s_prime) > 0.0f && !inside &&
                      do_scattering)
                  {
                     AtVector N = sg->N;
                     sg->N = m;
                     result_ss += AiSSSTraceSingleScatter(
                         sg, bssrdfbrdf(sigma_s_prime / sigma_t_prime), mfp,
                         ssDirection, transmissionIor);
                     sg->N = N;
                  }

                  if (minh(sample.opacity) < 1.0f)
                  {
                     AiTraceBackground(&wi_ray, &sample);
                     result_transmission +=
                         min(sample.color * f, rgb(data->transmissionClamp)) *
                         (AI_RGB_WHITE - sample.opacity);
                     assert(AiIsFinite(result_transmission));
                  }
               }
               else  // trace the background if we've hit nothing
               {
                  AiTraceBackground(&wi_ray, &sample);
                  result_transmission +=
                      min(sample.color * f, rgb(data->transmissionClamp));
                  assert(AiIsFinite(result_transmission));
               }
            }
         }

         inv_ns = AiSamplerGetSampleInvCount(sampit);
      }

      if (rr_transmission || tir)
      {
         kti = kti2 = 1.0f;
      }

      result_transmission *= inv_ns * transmissionColor * kti * kti2;
      assert(AiIsFinite(result_transmission));
      result_ss *= inv_ns * transmissionColor * kti * kti2;

      if (doDeepGroups)
      {
         for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
         {
            deepGroupsTransmission[i] *=
                inv_ns * transmissionColor * kti * kti2;
         }
      }

      if (transmitAovs)
      {
         for (int i = 0; i < NUM_AOVs; ++i)
         {
            childAovs[i] *= inv_ns * transmissionColor * kti * kti2;
         }
      }

      // unset the transmission trace set
      if (data->trace_set_transmission_enabled)
      {
         // if we defined a global trace set, re-set this, otherwise, unset
         if (data->trace_set_all_enabled)
         {
            AiShaderGlobalsSetTraceSet(sg, data->trace_set_all.c_str(),
                                       data->trace_set_all_inclusive);
         }
         else
         {
            AiShaderGlobalsUnsetTraceSet(sg);
         }
      }

      // reset the normal
      sg->N = Norig;
      sg->Nf = Nforig;

   }  // if (do_transmission)
   result_transmission += result_transmissionDirect;
   assert(AiIsFinite(result_transmission));

   // backlight
   // ---------
   if (do_backlight &&
       kti * kti2 * maxh(backlightColor) * backlightIndirectStrength >
           IMPORTANCE_EPS)
   {
      // set the backlight trace set if it's defined
      if (data->trace_set_backlight_enabled)
      {
         AiShaderGlobalsSetTraceSet(sg, data->trace_set_backlight.c_str(),
                                    data->trace_set_backlight_inclusive);
      }

      flipNormals(sg);
      float kr = kti * kti2;
      AtSamplerIterator* sampit =
          AiSamplerIterator(data->backlight_sampler, sg);
      AiMakeRay(&wi_ray, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg);
      int ssi = 0;
      while (AiSamplerGetSample(sampit, samples))
      {
         // cosine hemisphere sampling as O-N sampling does not work outside of
         // a light loop
         float stheta = sqrtf(float(samples[0]));
         float phi = float(AI_PITIMES2 * samples[1]);
         wi.x = stheta * cosf(phi);
         wi.y = stheta * sinf(phi);
         wi.z = sqrtf(1.0f - float(samples[0]));
         AiV3RotateToFrame(wi, U, V, sg->Nf);

         float cos_theta = AiV3Dot(wi, sg->Nf);
         if (cos_theta <= 0.0f) continue;

         float p = cos_theta * float(AI_ONEOVERPI);

         // trace the ray
         wi_ray.dir = wi;
         AtRGB f = kr * AiOrenNayarMISBRDF(bmis, &wi) * bd_prob / p;
         AtRGB throughput = path_throughput * f;
         AiStateSetMsgRGB("als_throughput", throughput);
         bool cont = true;
         float rr_p = 1.0f;
#ifdef RR_BOUNCES
         rr_p = std::min(1.0f, sqrtf(maxh(throughput) / maxh(path_throughput)));
         if (sg->Rt > 0)
         {
            cont = false;
            // get a permuted, stratified random number
            int pb = data->AA_samples * data->diffuse_samples2;
            int pi = data->diffuse_samples2 * sg->si + ssi;
            int idx = sg->Rr * pb + pi;
            float u = (float(data->perm_table_backlight[idx]) +
                       sampleTEAFloat(
                           idx, TEA_STREAM_ALSURFACE_RR_BACKLIGHT_JITTER)) /
                      float(pb);
            // offset based on pixel
            float offset =
                sampleTEAFloat(sg->y * data->xres + sg->x,
                               TEA_STREAM_ALSURFACE_RR_BACKLIGHT_OFFSET);
            u = fmodf(u + offset, 1.0f);

            if (u < rr_p)
            {
               cont = true;
               rr_p = 1.0f / rr_p;
               throughput *= rr_p;
               f *= rr_p;
            }
         }
#endif
         if (cont)
         {
            bool hit = AiTrace(&wi_ray, &scrs);
            result_backlightIndirect += scrs.color * f;

            // accumulate the lightgroup contributions calculated by the child
            // shader
            if (doDeepGroups && hit)
            {
               for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
               {
                  deepGroupsBacklight[i] += deepGroupPtr[i] * f;
               }
            }
         }

         ssi++;
      }
      result_backlightIndirect *= AiSamplerGetSampleInvCount(sampit) *
                                  backlightColor * backlightIndirectStrength *
                                  (1 - sssMix);

      if (doDeepGroups)
      {
         for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
         {
            deepGroupsBacklight[i] *= AiSamplerGetSampleInvCount(sampit) *
                                      backlightColor *
                                      backlightIndirectStrength;
         }
      }
      flipNormals(sg);

      // unset the backlight trace set
      if (data->trace_set_backlight_enabled)
      {
         // if we defined a global trace set, re-set this, otherwise, unset
         if (data->trace_set_all_enabled)
         {
            AiShaderGlobalsSetTraceSet(sg, data->trace_set_all.c_str(),
                                       data->trace_set_all_inclusive);
         }
         else
         {
            AiShaderGlobalsUnsetTraceSet(sg);
         }
      }

   }  // if (do_backlight)

   // Emission
   result_emission = emissionColor * kti * kti2;

   // Diffusion multiple scattering
   if (do_sss)
   {
      if (data->sssMode == SSSMODE_CUBIC)
      {
         AtRGB radius =
             max(rgb(0.0001), sssRadius * sssRadiusColor / sssDensityScale);
#if AI_VERSION_MAJOR_NUM > 0
         // if the user has only specified one layer (default) then just use
         // that
         if (sssWeight2 == 0.0f && sssWeight3 == 0.0f)
         {
            AtRGB weights[3] = {AI_RGB_RED, AI_RGB_GREEN, AI_RGB_BLUE};
            float r[3] = {radius.r, radius.g, radius.b};
            result_sss = AiBSSRDFCubic(sg, r, weights, 3);
         }
         else
         {
            // AtRGB r1 = sssRadius*sssRadiusColor/sssDensityScale;
            AtRGB r2 = max(rgb(0.0001),
                           sssRadius2 * sssRadiusColor2 / sssDensityScale);
            AtRGB r3 = max(rgb(0.0001),
                           sssRadius3 * sssRadiusColor3 / sssDensityScale);
            AtRGB weights[9] = {
                AI_RGB_RED * sssWeight1,   AI_RGB_GREEN * sssWeight1,
                AI_RGB_BLUE * sssWeight1,  AI_RGB_RED * sssWeight2,
                AI_RGB_GREEN * sssWeight2, AI_RGB_BLUE * sssWeight2,
                AI_RGB_RED * sssWeight3,   AI_RGB_GREEN * sssWeight3,
                AI_RGB_BLUE * sssWeight3};
            float r[9] = {radius.r, radius.g, radius.b, r2.r, r2.g,
                          r2.b,     r3.r,     r3.g,     r3.b};
            result_sss = AiBSSRDFCubic(sg, r, weights, 9);
         }
#else
         result_sss = AiSSSPointCloudLookupCubic(sg, radius);
#endif
      }
      else if (data->sssMode == SSSMODE_EMPIRICAL)
      {
         AtRGB radius = max(rgb(0.0001),
                            sssRadius * sssRadiusColor / (sssDensityScale * 7));

         AtRGB result_sss_direct = AI_RGB_BLACK;
         AtRGB result_sss_indirect = AI_RGB_BLACK;
         if (sssWeight2 == 0.0f && sssWeight3 == 0.0f)
         {
            AtRGB weights[3] = {AI_RGB_RED, AI_RGB_GREEN, AI_RGB_BLUE};
            AtRGB c = AI_RGB_WHITE;
            AiBSSRDFEmpirical(sg, result_sss_direct, result_sss_indirect,
                              &radius.r, &c.r, weights, 3);
            result_sss = result_sss_direct + result_sss_indirect;
         }
         else 
         {
            AtRGB r2 = max(rgb(0.0001),
                           sssRadius2 * sssRadiusColor2 / (sssDensityScale * 7));
            AtRGB r3 = max(rgb(0.0001),
                           sssRadius3 * sssRadiusColor3 / (sssDensityScale * 7));
            AtRGB weights[9] = {
                AI_RGB_RED * sssWeight1,   AI_RGB_GREEN * sssWeight1,
                AI_RGB_BLUE * sssWeight1,  AI_RGB_RED * sssWeight2,
                AI_RGB_GREEN * sssWeight2, AI_RGB_BLUE * sssWeight2,
                AI_RGB_RED * sssWeight3,   AI_RGB_GREEN * sssWeight3,
                AI_RGB_BLUE * sssWeight3};
            float r[9] = {radius.r, radius.g, radius.b, r2.r, r2.g,
                          r2.b,     r3.r,     r3.g,     r3.b};
            AtRGB c[3] = {AI_RGB_WHITE, AI_RGB_WHITE, AI_RGB_WHITE};
            AiBSSRDFEmpirical(sg, result_sss_direct, result_sss_indirect,
                              r, &c[0].r, weights, 9);
            result_sss = result_sss_direct + result_sss_indirect;
         }
      }
      else
      {
         // memset(diffusion_msgdata, 0, sizeof(DirectionalMessageData));
         int nc = 3;
         if (sssWeight2 > 0.0f) nc = 6;
         if (sssWeight3 > 0.0f) nc = 9;

         float Rd[9] = {
             sssRadiusColor.r,  sssRadiusColor.g,  sssRadiusColor.b,
             sssRadiusColor2.r, sssRadiusColor2.g, sssRadiusColor2.b,
             sssRadiusColor3.r, sssRadiusColor3.g, sssRadiusColor3.b};
         float radii[9] = {sssRadius,  sssRadius,  sssRadius,
                           sssRadius2, sssRadius2, sssRadius2,
                           sssRadius3, sssRadius3, sssRadius3};
         AtRGB weights[9] = {
             AI_RGB_RED * sssWeight1,   AI_RGB_GREEN * sssWeight1,
             AI_RGB_BLUE * sssWeight1,  AI_RGB_RED * sssWeight2,
             AI_RGB_GREEN * sssWeight2, AI_RGB_BLUE * sssWeight2,
             AI_RGB_RED * sssWeight3,   AI_RGB_GREEN * sssWeight3,
             AI_RGB_BLUE * sssWeight3};
         memcpy(diffusion_msgdata->weights, weights, sizeof(AtRGB) * nc);
         for (int i = 0; i < nc; ++i)
         {
            diffusion_msgdata->sp[i] = ScatteringProfileDirectional(
                MAX(Rd[i], 0.001f), sssDensityScale / radii[i]);
         }

         /*
         float g = 0.0f;
         // skin2
         const float sigma_s_prime[3] = {1.09, 1.59, 1.79};
         const float sigma_a[3] = {0.013, 0.07, 0.145};
         // marble
         const float sigma_s_prime[3] = {2.19, 2.62, 3.00};
         const float sigma_a[3] = {0.0021, 0.0041, 0.0071};
         // wholemilk
         const float sigma_s_prime[3] = {2.55, 3.21, 3.77};
         const float sigma_a[3] = {0.0011, 0.0024, 0.014};
         // ketchup
         const float sigma_s_prime[3] = {0.18, 0.07, 0.03};
         const float sigma_a[3] = {0.061, 0.97, 1.45};
         sp[0] = ScatteringProfileDirectional(sigma_s_prime[0] / (1-g) *
         sssDensityScale*10, sigma_a[0] * sssDensityScale*10, g);
         sp[1] = ScatteringProfileDirectional(sigma_s_prime[1] / (1-g) *
         sssDensityScale*10, sigma_a[1] * sssDensityScale*10, g);
         sp[2] = ScatteringProfileDirectional(sigma_s_prime[2] / (1-g) *
         sssDensityScale*10, sigma_a[2] * sssDensityScale*10, g);
         */
         AtRGB result_sss_direct;
         AtRGB result_sss_indirect;
         AtRGB lightGroupsSss[NUM_LIGHT_GROUPS];
         memset(lightGroupsSss, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
         result_sss = alsDiffusion(
             sg, diffusion_msgdata, data->sss_sampler,
             data->sssMode == SSSMODE_DIRECTIONAL, nc, result_sss_direct,
             result_sss_indirect, lightGroupsSss, deepGroupsSss, deepGroupPtr,
             data->trace_set_sss.c_str(), data->trace_set_sss_enabled,
             data->trace_set_sss_inclusive);
         for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
         {
            lightGroupsDirect[i] +=
                lightGroupsSss[i] * diffuseColor * sssMix * kti * kti2;
            ;
            deepGroupsSss[i] *= diffuseColor * sssMix * kti * kti2;
            ;
         }
      }
      result_sss *= diffuseColor * sssMix * kti * kti2;
   }

   // blend sss and diffuse
   result_diffuseIndirect *= (1 - sssMix);
   result_backlightIndirect *= (1 - sssMix);

   // Now accumulate the deep group brdf results onto the relevant samples
   if (sg->Rt & AI_RAY_CAMERA)
   {
      // stick some AOVs in the state for an alCel to pick up
      AiStateSetMsgRGB("als_diffuse_color", diffuseColor);
      AiStateSetMsgRGB("als_direct_diffuse_raw", result_diffuseDirectRaw);
      AiStateSetMsgRGB("als_direct_specular", result_glossyDirect);
      AiStateSetMsgRGB("als_indirect_diffuse", result_diffuseIndirect);

      if (data->standardAovs)
      {
         AtRGB tmp;
         tmp = result_diffuseDirect + result_backlightDirect;
         if (tmp != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_direct_diffuse].c_str(), tmp);
         tmp = result_diffuseIndirect + result_backlightIndirect;
         if (tmp != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_indirect_diffuse].c_str(), tmp);
         tmp = result_glossyDirect + result_glossy2Direct;
         if (tmp != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_direct_specular].c_str(), tmp);
         tmp = result_glossyIndirect + result_glossy2Indirect;
         if (tmp != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_indirect_specular].c_str(), tmp);
         tmp = result_transmission + result_ss;
         if (tmp != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_refraction].c_str(), tmp);

         if (result_sss != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_sss].c_str(), result_sss);
         if (result_emission != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_emission].c_str(), result_emission);
      }
      else if (transmitAovs && do_transmission)
      {
         for (int i = 0; i < NUM_AOVs; ++i)
         {
            if (i == k_refraction)
            {
               if (result_transmission != AI_RGB_BLACK)
               {
                  AiAOVSetRGB(sg, data->aovs[k_refraction].c_str(),
                              result_transmission);
               }
               continue;
            }
            else if (i == k_direct_specular)
            {
               AtRGB tmp =
                   result_glossyDirect + transmittedAovPtr[k_direct_specular];
               if (tmp != AI_RGB_BLACK)
               {
                  AiAOVSetRGB(sg, data->aovs[k_direct_specular].c_str(), tmp);
               }
               continue;
            }
            else if (i == k_indirect_specular)
            {
               AtRGB tmp = result_glossyIndirect +
                           transmittedAovPtr[k_indirect_specular];
               if (tmp != AI_RGB_BLACK)
               {
                  AiAOVSetRGB(sg, data->aovs[k_indirect_specular].c_str(), tmp);
               }
               continue;
            }

            if (transmittedAovPtr[i] != AI_RGB_BLACK)
               AiAOVSetRGB(sg, data->aovs[i].c_str(), transmittedAovPtr[i]);
         }
         if (doDeepGroups)
           {
             AtRGB deepGroups[NUM_LIGHT_GROUPS];
             memset(deepGroups, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
             for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
               {
                 deepGroups[i] = deepGroupsDiffuse[i] + deepGroupsGlossy[i] +
                   deepGroupsGlossy2[i] +
                   deepGroupsTransmission[i] +
                   deepGroupsBacklight[i] + deepGroupsSss[i] +
                   lightGroupsDirect[i];

                 if (deepGroups[i] != AI_RGB_BLACK)
                   {
                     deepGroups[i] =
                       min(deepGroups[i], rgb(data->aov_light_group_clamp[i]));
                     AiAOVSetRGB(sg, data->aovs[k_light_group_1 + i].c_str(),
                                 deepGroups[i]);
                   }
               }
           }
      }
      else
      {
         if (doDeepGroups)
         {
            AtRGB deepGroups[NUM_LIGHT_GROUPS];
            memset(deepGroups, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
            for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
            {
               deepGroups[i] = deepGroupsDiffuse[i] + deepGroupsGlossy[i] +
                               deepGroupsGlossy2[i] +
                               deepGroupsTransmission[i] +
                               deepGroupsBacklight[i] + deepGroupsSss[i] +
                               lightGroupsDirect[i];

               if (deepGroups[i] != AI_RGB_BLACK)
               {
                  deepGroups[i] =
                      min(deepGroups[i], rgb(data->aov_light_group_clamp[i]));
                  AiAOVSetRGB(sg, data->aovs[k_light_group_1 + i].c_str(),
                              deepGroups[i]);
               }
            }
         }
         else
         {
            for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
            {
               if (lightGroupsDirect[i] != AI_RGB_BLACK)
               {
                  lightGroupsDirect[i] =
                      min(lightGroupsDirect[i],
                          rgb(data->aov_light_group_clamp[i]));
                  AiAOVSetRGB(sg, data->aovs[k_light_group_1 + i].c_str(),
                              lightGroupsDirect[i]);
               }
            }
         }

         for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
         {
            if (shadowGroups[i] != AI_RGBA_BLACK)
               AiAOVSetRGBA(sg, data->aovs_rgba[k_shadow_group_1 + i].c_str(),
                            shadowGroups[i]);
         }

         diffuseColor = min(diffuseColor, rgb(data->aov_diffuse_color_clamp));
         result_diffuseDirect =
             min(result_diffuseDirect, rgb(data->aov_direct_diffuse_clamp));
         result_diffuseDirectRaw = min(result_diffuseDirectRaw,
                                       rgb(data->aov_direct_diffuse_raw_clamp));
         result_backlightDirect =
             min(result_backlightDirect, rgb(data->aov_direct_backlight_clamp));
         result_sss = min(result_sss, rgb(data->aov_sss_clamp));
         result_glossyDirect =
             min(result_glossyDirect, rgb(data->aov_direct_specular_clamp));
         result_glossy2Direct =
             min(result_glossy2Direct, rgb(data->aov_direct_specular_2_clamp));
         result_diffuseIndirect =
             min(result_diffuseIndirect, rgb(data->aov_indirect_diffuse_clamp));
         result_diffuseIndirectRaw =
             min(result_diffuseIndirectRaw,
                 rgb(data->aov_indirect_diffuse_raw_clamp));
         result_backlightIndirect = min(
             result_backlightIndirect, rgb(data->aov_indirect_backlight_clamp));
         result_glossyIndirect =
             min(result_glossyIndirect, rgb(data->aov_indirect_specular_clamp));
         result_glossy2Indirect = min(result_glossy2Indirect,
                                      rgb(data->aov_indirect_specular_2_clamp));
         result_ss = min(result_ss, rgb(data->aov_single_scatter_clamp));
         result_transmission =
             min(result_transmission, rgb(data->aov_refraction_clamp));
         result_emission = min(result_emission, rgb(data->aov_emission_clamp));

         if (diffuseColor != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_diffuse_color].c_str(), diffuseColor);
         if (result_diffuseDirect != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_direct_diffuse].c_str(),
                        result_diffuseDirect);
         if (result_diffuseDirectRaw != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_direct_diffuse_raw].c_str(),
                        result_diffuseDirectRaw);
         if (result_backlightDirect != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_direct_backlight].c_str(),
                        result_backlightDirect);
         if (result_sss != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_sss].c_str(), result_sss);
         if (result_glossyDirect != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_direct_specular].c_str(),
                        result_glossyDirect);
         if (result_glossy2Direct != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_direct_specular_2].c_str(),
                        result_glossy2Direct);
         if (result_diffuseIndirect != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_indirect_diffuse].c_str(),
                        result_diffuseIndirect);
         if (result_diffuseIndirectRaw != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_indirect_diffuse_raw].c_str(),
                        result_diffuseIndirectRaw);
         if (result_backlightIndirect != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_indirect_backlight].c_str(),
                        result_backlightIndirect);
         if (result_glossyIndirect != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_indirect_specular].c_str(),
                        result_glossyIndirect);
         if (result_glossy2Indirect != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_indirect_specular_2].c_str(),
                        result_glossy2Indirect);
         if (result_ss != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_single_scatter].c_str(), result_ss);
         if (result_transmission != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_refraction].c_str(),
                        result_transmission);
         if (result_emission != AI_RGB_BLACK)
            AiAOVSetRGB(sg, data->aovs[k_emission].c_str(), result_emission);

         // write IDs
         for (int i = 0; i < NUM_ID_AOVS; ++i)
         {
            AtRGB tmp;
            // check if output is enabled first in case we have an expensive
            // network upstream
            if (AiAOVEnabled(data->aovs[k_id_1 + i].c_str(), AI_TYPE_RGB))
            {
               tmp = AiShaderEvalParamRGB(p_id1 + i);

               // check if we're overriding it with a per-object id
               if (AiNodeLookUpUserParameter(sg->Op, id_names[i]))
               {
                  tmp = AiNodeGetRGB(sg->Op, id_names[i]);
               }

               if (tmp != AI_RGB_BLACK)
                  AiAOVSetRGB(sg, data->aovs[k_id_1 + i].c_str(), tmp);
            }
         }

         // write data AOVs
         AtRGB uv = AtRGB(sg->u, sg->v, 0.0f);
         AiAOVSetRGB(sg, data->aovs[k_uv].c_str(), uv);
         AtRGB depth =
                 AtRGB(float(sg->Rl), AiV3Dot(sg->Nf, wo), sg->P.y);
         AiAOVSetRGB(sg, data->aovs[k_depth].c_str(), depth);
      }
   }
   else  // we're in a secondary ray //
   {
      if (doDeepGroups)
      {
         for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
         {
            deepGroupPtr[i] = deepGroupsDiffuse[i] + deepGroupsGlossy[i] +
                              deepGroupsGlossy2[i] + deepGroupsTransmission[i] +
                              deepGroupsBacklight[i] + deepGroupsSss[i] +
                              lightGroupsDirect[i];
         }
      }
      if (transmitAovs)
      {
         if (do_transmission)
         {
            for (int i = 0; i < NUM_AOVs; ++i)
            {
               transmittedAovPtr[i] = childAovs[i];
            }
         }
         else
         {
            transmittedAovPtr[k_diffuse_color] = diffuseColor;
            transmittedAovPtr[k_direct_diffuse] = result_diffuseDirect;
            transmittedAovPtr[k_direct_diffuse_raw] = result_diffuseDirectRaw;
            transmittedAovPtr[k_indirect_diffuse] = result_diffuseIndirect;
            transmittedAovPtr[k_indirect_diffuse_raw] =
                result_diffuseIndirectRaw;
            transmittedAovPtr[k_direct_backlight] = result_backlightDirect;
            transmittedAovPtr[k_indirect_backlight] = result_backlightIndirect;
            transmittedAovPtr[k_direct_specular] = result_glossyDirect;
            transmittedAovPtr[k_indirect_specular] = result_glossyIndirect;
            transmittedAovPtr[k_direct_specular_2] = result_glossy2Direct;
            transmittedAovPtr[k_indirect_specular_2] = result_glossy2Indirect;
            transmittedAovPtr[k_single_scatter] = result_ss;
            transmittedAovPtr[k_sss] = result_sss;
            transmittedAovPtr[k_refraction] = result_transmission;
            transmittedAovPtr[k_emission] = result_emission;
            transmittedAovPtr[k_uv] = rgb(sg->u, sg->v, 0.0f);
            transmittedAovPtr[k_depth] =
                rgb(float(sg->Rl), AiV3Dot(sg->Nf, wo), 0.0f);
            for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
            {
               transmittedAovPtr[k_light_group_1 + i] = lightGroupsDirect[i];
            }
            for (int i = 0; i < NUM_ID_AOVS; ++i)
            {
               AtRGB tmp = AiShaderEvalParamRGB(p_id1 + i);

               // check if we're overriding it with a per-object id
               if (AiNodeLookUpUserParameter(sg->Op, id_names[i]))
               {
                  tmp = AiNodeGetRGB(sg->Op, id_names[i]);
               }

               if (tmp != AI_RGB_BLACK) transmittedAovPtr[k_id_1 + i] = tmp;
            }
         }
      }
   }

   assert(AiIsFinite(result_diffuseDirect));
   assert(AiIsFinite(result_backlightDirect));
   assert(AiIsFinite(result_sss));
   assert(AiIsFinite(result_glossyDirect));
   assert(AiIsFinite(result_glossy2Direct));
   assert(AiIsFinite(result_diffuseIndirect));
   assert(AiIsFinite(result_glossyIndirect));
   assert(AiIsFinite(result_glossy2Indirect));
   assert(AiIsFinite(result_ss));
   assert(AiIsFinite(result_transmission));
   assert(AiIsFinite(result_emission));

   // Sum final result from temporaries
   //
   sg->out.RGB = result_diffuseDirect + result_backlightDirect + result_sss +
                 result_glossyDirect + result_glossy2Direct +
                 result_diffuseIndirect + result_backlightIndirect +
                 result_glossyIndirect + result_glossy2Indirect + result_ss +
                 result_transmission + result_emission;

   assert(AiIsFinite(sg->out.RGB));
}
