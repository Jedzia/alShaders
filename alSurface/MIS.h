#pragma once
#include <ai.h>
#include "alUtil.h"
#include "fresnel.h"
#include <cassert>

struct BrdfData_wrap
{
   BrdfData_wrap() : brdf_data(NULL), fr(NULL), kr_int(AI_RGB_BLACK), ns(0.0f)
   {
   }

   void* brdf_data;
   AtShaderGlobals* sg;
   Fresnel* fr;
   float eta;
   AtVector V;
   AtVector N;
   mutable AtRGB kr;
   mutable AtRGB kr_int;
   mutable AtRGB brdf;
   mutable float pdf;
   mutable float ns;
   mutable bool ibs;
};

/* ToDoJed: Fix for Porting >v5, AiEvaluateLightSample removed
AtRGB AiMicrofacetMISBRDF_wrap(const void* brdf_data, const AtVector* indir)
{
   const BrdfData_wrap* brdfw =
       reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   AtRGB result = AI_RGB_BLACK;
   if (brdfw->ibs)
   {
      result = brdfw->kr * brdfw->brdf;
   }
   else
   {
      AtVector H;
      AiV3Normalize(H, (*indir) + brdfw->V);


// old AI_API AI_PURE  AtColor  AiMicrofacetMISBRDF(const void* brdf_data, const AtVector* indir);
// new AI_API AI_DEVICE AtBSDF* AiMicrofacetBSDF(const AtShaderGlobals* sg, const AtRGB& weight,
//                                          int distribution, const AtVector& N, const AtVector *U,
//                                          float eta, float rx, float ry, uint8_t exit_type = 0,
//                                          const AtString label = AtString());
      result = brdfw->fr->kr(std::max(0.0f, AiV3Dot(H, *indir))) *
               AiMicrofacetMISBRDF(brdfw->brdf_data, indir);
   }

   return result;
}

float AiMicrofacetMISPDF_wrap(const void* brdf_data, const AtVector* indir)
{
   const BrdfData_wrap* brdfw =
       reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   if (brdfw->ibs)
   {
      return brdfw->pdf;
   }
   else
   {
      return AiMicrofacetMISPDF(brdfw->brdf_data, indir);
   }
}

AtVector AiMicrofacetMISSample_wrap(const void* brdf_data, float randx,
                                    float randy)
{
   const BrdfData_wrap* brdfw =
       reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   brdfw->ibs = true;
   AtVector indir = AiMicrofacetMISSample(brdfw->brdf_data, randx, randy);
   if (!AiV3IsZero(indir))
   {
      AtVector H;
      AiV3Normalize(H, (indir) + brdfw->V);
      brdfw->kr = brdfw->fr->kr(std::max(0.0f, AiV3Dot(H, indir)));
      brdfw->brdf = AiMicrofacetMISBRDF(brdfw->brdf_data, &indir);
      brdfw->pdf = AiMicrofacetMISPDF(brdfw->brdf_data, &indir);
      if (brdfw->pdf > 0.0f)
      {
         AtRGB w = brdfw->brdf / brdfw->pdf;
         brdfw->kr_int += brdfw->kr * w;
         assert(AiIsFinite(brdfw->kr_int));
         brdfw->ns++;
      }
   }
   return indir;
}

AtRGB AiCookTorranceMISBRDF_wrap(const void* brdf_data, const AtVector* indir)
{
   const BrdfData_wrap* brdfw =
       reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   AtRGB result = AI_RGB_BLACK;
   if (brdfw->ibs)
   {
      result = brdfw->kr * brdfw->brdf;
   }
   else
   {
      AtVector H;
      AiV3Normalize(H, (*indir) + brdfw->V);
      result = brdfw->fr->kr(std::max(0.0f, AiV3Dot(H, *indir))) *
               AiCookTorranceMISBRDF(brdfw->brdf_data, indir);
   }

   return result;
}

float AiCookTorranceMISPDF_wrap(const void* brdf_data, const AtVector* indir)
{
   const BrdfData_wrap* brdfw =
       reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   if (brdfw->ibs)
   {
      return brdfw->pdf;
   }
   else
   {
      return AiCookTorranceMISPDF(brdfw->brdf_data, indir);
   }
}

AtVector AiCookTorranceMISSample_wrap(const void* brdf_data, float randx,
                                      float randy)
{
   const BrdfData_wrap* brdfw =
       reinterpret_cast<const BrdfData_wrap*>(brdf_data);
   brdfw->ibs = true;
   AtVector indir = AiCookTorranceMISSample(brdfw->brdf_data, randx, randy);
   if (!AiV3IsZero(indir))
   {
      AtVector H;
      AiV3Normalize(H, (indir) + brdfw->V);
      brdfw->kr = brdfw->fr->kr(std::max(0.0f, AiV3Dot(H, indir)));
      brdfw->brdf = AiCookTorranceMISBRDF(brdfw->brdf_data, &indir);
      brdfw->pdf = AiCookTorranceMISPDF(brdfw->brdf_data, &indir);
      if (brdfw->pdf > 0.0f)
      {
         AtRGB w = brdfw->brdf / brdfw->pdf;
         brdfw->kr_int += brdfw->kr * w;
         assert(AiIsFinite(brdfw->kr_int));
         brdfw->ns++;
      }
   }
   return indir;
}
*/
