#pragma once
// aovs.h

#define NUM_LIGHT_GROUPS 16
#define NUM_AOVs 41
#define NUM_AOVs_RGBA 16

#define REGISTER_AOVS \
data->aovs.clear(); \
data->aovs.push_back(AiNodeGetStr(node, "aov_diffuse_color")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_direct_diffuse")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_direct_diffuse_raw")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_indirect_diffuse")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_indirect_diffuse_raw")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_direct_backlight")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_indirect_backlight")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_direct_specular")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_indirect_specular")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_direct_specular_2")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_indirect_specular_2")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_single_scatter")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_sss")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_refraction")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_emission")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_uv")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_depth")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_1")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_2")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_3")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_4")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_5")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_6")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_7")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_8")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_9")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_10")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_11")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_12")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_13")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_14")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_15")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_light_group_16")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_id_1")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_id_2")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_id_3")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_id_4")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_id_5")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_id_6")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_id_7")); \
data->aovs.push_back(AiNodeGetStr(node, "aov_id_8")); \
assert(NUM_AOVs == data->aovs.size() && "NUM_AOVs does not match size of aovs array!"); \
for (size_t i=0; i < data->aovs.size(); ++i) \
    	AiAOVRegister(data->aovs[i].c_str(), AI_TYPE_RGB, AI_AOV_BLEND_OPACITY); \
data->aovs_rgba.clear(); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_1")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_2")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_3")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_4")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_5")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_6")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_7")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_8")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_9")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_10")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_11")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_12")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_13")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_14")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_15")); \
data->aovs_rgba.push_back(AiNodeGetStr(node, "aov_shadow_group_16")); \
assert(NUM_AOVs_RGBA == data->aovs_rgba.size() && "NUM_AOVs_RGBA does not match size of aovs_rgba array!"); \
for (size_t i=0; i < data->aovs_rgba.size(); ++i) \
        AiAOVRegister(data->aovs_rgba[i].c_str(), AI_TYPE_RGBA, AI_AOV_BLEND_OPACITY); \

enum AovIndices
{
	k_diffuse_color = 0,
    k_direct_diffuse,
    k_direct_diffuse_raw,
    k_indirect_diffuse,
    k_indirect_diffuse_raw,
    k_direct_backlight,
    k_indirect_backlight,
    k_direct_specular,
    k_indirect_specular,
    k_direct_specular_2,
    k_indirect_specular_2,
    k_single_scatter,
    k_sss,
    k_refraction,
    k_emission,
    k_uv,
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

enum AovRGBIndices
{
    k_shadow_group_1=0,
    k_shadow_group_2,
    k_shadow_group_3,
    k_shadow_group_4,
    k_shadow_group_5,
    k_shadow_group_6,
    k_shadow_group_7,
    k_shadow_group_8,
    k_shadow_group_9,
    k_shadow_group_10,
    k_shadow_group_11,
    k_shadow_group_12,
    k_shadow_group_13,
    k_shadow_group_14,
    k_shadow_group_15,
    k_shadow_group_16
};
