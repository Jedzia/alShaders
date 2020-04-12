#include "alUtil.h"
#include <ai.h>
#include <vector>
#include <string>
#include <cassert>
#include "cryptomatte/cryptomatte.h"

#include "aovs.h"

AI_SHADER_NODE_EXPORT_METHODS(alLayer)

struct ShaderData
{
   // AOV names
   std::vector<std::string> aovs;
   std::vector<std::string> aovs_rgba;
   bool standardAovs;
   CryptomatteData* cryptomatte;
};

enum alLayerParams
{
   p_layer1,
   p_layer2,
   p_mix,
   p_debug,

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

   p_aiEnableMatte,
   p_aiMatteColor,
   p_aiMatteColorA
};

enum DebugModes
{
   kOff = 0,
   kLayer1,
   kLayer2,
   kMixer
};

static const char* debugModeNames[] = {"off", "layer1", "layer2", "mixer",
                                       NULL};

node_parameters
{
   AiParameterRGB("layer1", 0.0f, 0.0f, 0.0f);
   AiParameterRGB("layer2", 0.0f, 0.0f, 0.0f);
   AiParameterFlt("mix", 0.0f);
   AiParameterEnum("debug", kOff, debugModeNames);

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

   AiParameterBOOL("aiEnableMatte", false);
   AiParameterRGB("aiMatteColor", 0.0f, 0.0f, 0.0f);
   AiParameterFlt("aiMatteColorA", 0.0f);
}

node_loader
{
   if (i > 0) return 0;
   node->methods = alLayer;
   node->output_type = AI_TYPE_RGB;
   node->name = "alLayer";
   node->node_type = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return true;
}

node_initialize
{
   ShaderData* data = new ShaderData;
   data->cryptomatte = new CryptomatteData();
   AiNodeSetLocalData(node, data);
}

node_finish
{
   if (AiNodeGetLocalData(node))
   {
      ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
      if (data->cryptomatte)
         delete data->cryptomatte;
      AiNodeSetLocalData(node, NULL);
      delete data;
   }
}

node_update
{
   ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);

   // set up AOVs
   REGISTER_AOVS
   data->standardAovs = params[p_standardAovs].BOOL;
   data->cryptomatte->setup_all(AiNodeGetStr(node, "aov_crypto_asset"), 
      AiNodeGetStr(node, "aov_crypto_object"), AiNodeGetStr(node, "aov_crypto_material"));
}

shader_evaluate
{
   ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);

   AtRGB result = AI_RGB_BLACK;
   AtRGB result_opacity = AI_RGB_WHITE;

   float mix = AiShaderEvalParamFlt(p_mix);
   int debug = AiShaderEvalParamEnum(p_debug);

   if (debug == kMixer)
   {
      result = AtRGB(mix, mix, mix);
   }
   else
   {
      if (debug == kLayer1)
         mix = 0.0f;
      else if (debug == kLayer2)
         mix = 1.0f;

      int als_raytype = ALS_RAY_UNDEFINED;
      AiStateGetMsgInt("als_raytype", &als_raytype);

      if (als_raytype != ALS_RAY_SSS && mix >= (1.0f - IMPORTANCE_EPS))
      {
         result = AiShaderEvalParamRGB(p_layer2);
         result_opacity = sg->out_opacity;
      }
      else if (als_raytype != ALS_RAY_SSS && mix <= IMPORTANCE_EPS)
      {
         result = AiShaderEvalParamRGB(p_layer1);
         result_opacity = sg->out_opacity;
      }
      else
      {
         if (sg->Rt & AI_RAY_CAMERA)  // handle aovs
         {
            AiStateSetMsgInt("als_context", ALS_CONTEXT_LAYER);
            // RGB AOVs
            AtRGB tmp[NUM_AOVs];
            memset(tmp, 0, sizeof(AtRGB) * NUM_AOVs);
            AtRGB layer1 = AiShaderEvalParamRGB(p_layer1);
            AtRGB layer1_opacity = sg->out_opacity;
            for (size_t i = 0; i < data->aovs.size(); ++i)
            {
               if (!AiAOVGetRGB(sg, data->aovs[i].c_str(), tmp[i]))
               {
                  tmp[i] = AI_RGB_BLACK;
               }
               AiAOVSetRGB(sg, data->aovs[i].c_str(), AI_RGB_BLACK);
            }

            // RGBA AOVs
            AtRGBA tmp_rgba[NUM_AOVs_RGBA];
            memset(tmp_rgba, 0, sizeof(AtRGBA) * NUM_AOVs_RGBA);
            for (size_t i = 0; i < data->aovs_rgba.size(); ++i)
            {
               if (!AiAOVGetRGBA(sg, data->aovs_rgba[i].c_str(), tmp_rgba[i]))
               {
                  tmp_rgba[i] = AI_RGBA_BLACK;
               }
               AiAOVSetRGBA(sg, data->aovs_rgba[i].c_str(), AI_RGBA_BLACK);
            }

            AtRGB layer2 = AiShaderEvalParamRGB(p_layer2);
            result = lerp(layer1, layer2, mix);
            result_opacity = lerp(layer1_opacity, sg->out_opacity, mix);
            for (size_t i = 0; i < data->aovs.size(); ++i)
            {
               AtRGB tmp2;
               if (!AiAOVGetRGB(sg, data->aovs[i].c_str(), tmp2))
               {
                  tmp2 = AI_RGB_BLACK;
               }
               AiAOVSetRGB(sg, data->aovs[i].c_str(), lerp(tmp[i], tmp2, mix));
            }

            for (size_t i = 0; i < data->aovs_rgba.size(); ++i)
            {
               AtRGBA tmp_rgba2;
               if (!AiAOVGetRGBA(sg, data->aovs_rgba[i].c_str(), tmp_rgba2))
               {
                  tmp_rgba2 = AI_RGBA_BLACK;
               }
               AiAOVSetRGBA(sg, data->aovs_rgba[i].c_str(),
                            lerp(tmp_rgba[i], tmp_rgba2, mix));
            }
            AiStateSetMsgInt("als_context", ALS_CONTEXT_NONE);
         }
         else  // just layer the results
         {
            AiStateSetMsgInt("als_context", ALS_CONTEXT_LAYER);
            AtRGB deepGroupTmp1[NUM_LIGHT_GROUPS];
            AtRGB deepGroupTmp2[NUM_LIGHT_GROUPS];
            memset(deepGroupTmp1, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
            memset(deepGroupTmp2, 0, sizeof(AtRGB) * NUM_LIGHT_GROUPS);
            AtRGB* deepGroupPtr = NULL;
            AtRGB layer1 = AiShaderEvalParamRGB(p_layer1);
            AtRGB layer1_opacity = sg->out_opacity;
            if (AiStateGetMsgPtr("als_deepGroupPtr", (void**)&deepGroupPtr))
            {
               memcpy(deepGroupTmp1, deepGroupPtr,
                      sizeof(AtRGB) * NUM_LIGHT_GROUPS);
            }
            AtRGB layer2 = AiShaderEvalParamRGB(p_layer2);
            result = lerp(layer1, layer2, mix);
            result_opacity = lerp(layer1_opacity, sg->out_opacity, mix);
            if (AiStateGetMsgPtr("als_deepGroupPtr", (void**)&deepGroupPtr))
            {
               memcpy(deepGroupTmp2, deepGroupPtr,
                      sizeof(AtRGB) * NUM_LIGHT_GROUPS);
               for (int i = 0; i < NUM_LIGHT_GROUPS; ++i)
               {
                  deepGroupPtr[i] =
                      lerp(deepGroupTmp1[i], deepGroupTmp2[i], mix);
               }
            }
            AiStateSetMsgInt("als_context", ALS_CONTEXT_NONE);
         }
      }
   }

   sg->out.RGB = result;
   sg->out_opacity = result_opacity;

   data->cryptomatte->do_cryptomattes(sg, node, -1, -1, -1);
}
