#pragma once

#include <ai.h>
#include <map>
#include <string>
#include <vector>

#include "fresnel.h"
#include "stats.h"

struct ShaderData
{
   ShaderData()
   : sss_samples_taken("sss_samples")
   {

   }
   AtSampler* diffuse_sampler;
   AtSampler* sss_sampler;
   AtSampler* glossy_sampler;
   AtSampler* glossy2_sampler;
   AtSampler* refraction_sampler;
   AtSampler* backlight_sampler;
   int GI_diffuse_depth;
   int GI_reflection_depth;
   int GI_refraction_depth;
   int GI_glossy_depth;
   int GI_diffuse_samples;
   int diffuse_samples2;
   int GI_glossy_samples;
   int glossy_samples2;
   int GI_glossy2_samples;
   int glossy2_samples2;
   int GI_refraction_samples;
   int refraction_samples2;
   int diffuse_sample_offset;
   int glossy_sample_offset;
   int glossy2_sample_offset;
   int refraction_sample_offset;
   int backlight_sample_offset;
   int total_samples;
   int sss_bssrdf_samples;
   int sss_bssrdf_samples2;
   AtCritSec cs;
   std::map<AtNode*, int> lightGroups;
   std::map<AtNode*, float> shadowDensities;
   bool specular1NormalConnected;
   bool specular2NormalConnected;
   bool standardAovs;
   bool transmitAovs;
   int numLights;
   int rrTransmissionDepth;
   bool transmissionDoDirect;

   // data for doing RR
   bool do_rr;
   int AA_samples;
   float AA_samples_inv;
   int total_depth;
   int* perm_table;
   int* perm_table_diffuse;
   int* perm_table_spec1;
   int* perm_table_spec2;
   int* perm_table_backlight;
   int* perm_table_sss;
   int xres;

   float specular1IndirectClamp;
   float specular2IndirectClamp;
   float transmissionClamp;
   float diffuseIndirectClamp;

   // AOV names
   std::vector<std::string> aovs;
   std::vector<std::string> aovs_rgba;

   // Fresnel
   int specular1FresnelMode;
   int specular2FresnelMode;

   std::string trace_set_all;
   bool trace_set_all_enabled;
   bool trace_set_all_inclusive;

   std::string trace_set_shadows;
   bool trace_set_shadows_enabled;
   bool trace_set_shadows_inclusive;

   std::string trace_set_diffuse;
   bool trace_set_diffuse_enabled;
   bool trace_set_diffuse_inclusive;

   std::string trace_set_backlight;
   bool trace_set_backlight_enabled;
   bool trace_set_backlight_inclusive;

   std::string trace_set_specular1;
   bool trace_set_specular1_enabled;
   bool trace_set_specular1_inclusive;

   std::string trace_set_specular2;
   bool trace_set_specular2_enabled;
   bool trace_set_specular2_inclusive;

   std::string trace_set_transmission;
   bool trace_set_transmission_enabled;
   bool trace_set_transmission_inclusive;

   std::string trace_set_sss;
   bool trace_set_sss_enabled;
   bool trace_set_sss_inclusive;

   bool cel_connected;

   int sssMode;
   Range sss_samples_taken;
};

#define RAND_STREAM_ALSURFACE_RR_PERMUTE 0
#define RAND_STREAM_ALSURFACE_RR_DIFF_PERMUTE 10000
#define RAND_STREAM_ALSURFACE_RR_SPEC1_PERMUTE 20000
#define RAND_STREAM_ALSURFACE_RR_SPEC2_PERMUTE 30000
#define RAND_STREAM_ALSURFACE_RR_BACKLIGHT_PERMUTE 40000
#define RAND_STREAM_ALSURFACE_RR_SSS_PERMUTE 50000


#define TEA_STREAM_ALSURFACE_RR_OFFSET 0
#define TEA_STREAM_ALSURFACE_RR_JITTER 1

#define TEA_STREAM_ALSURFACE_RR_DIFF_OFFSET 2
#define TEA_STREAM_ALSURFACE_RR_DIFF_JITTER 3

#define TEA_STREAM_ALSURFACE_RR_SPEC1_OFFSET 4
#define TEA_STREAM_ALSURFACE_RR_SPEC1_JITTER 5

#define TEA_STREAM_ALSURFACE_RR_SPEC2_OFFSET 6
#define TEA_STREAM_ALSURFACE_RR_SPEC2_JITTER 7

#define TEA_STREAM_ALSURFACE_RR_BACKLIGHT_OFFSET 8
#define TEA_STREAM_ALSURFACE_RR_BACKLIGHT_JITTER 9

#define TEA_STREAM_ALSURFACE_RR_SSS_OFFSET 10
#define TEA_STREAM_ALSURFACE_RR_SSS_JITTER 11

#define TEA_STREAM_ALSURFACE_RR2_OFFSET 12
