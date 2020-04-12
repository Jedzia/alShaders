#include <ai.h>
#include <cstring>
#include "tea.h"
#include "alUtil.h"

AI_SHADER_NODE_EXPORT_METHODS(alFlake)

enum alFlakeParams
{
   p_amount = 0,
   p_size,
   p_divergence,
   p_P,
   p_space
};

enum FlakeSpace
{
   S_TANGENT = 0,
   S_WORLD
};

const char* space_names[] = 
{
   "tangent",
   "world",
   NULL
};

node_parameters
{
   AiParameterFlt("amount", 0.7f);
   AiParameterFlt("size", 0.01f);
   AiParameterFlt("divergence", 0.5f);
   AiParameterVec("P", 0.0f, 0.0, 0.0);
   AiParameterEnum("space", S_TANGENT, space_names);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alFlake;
   node->output_type = AI_TYPE_VECTOR;
   node->name        = "alFlake";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return true;
}

node_initialize
{

}

node_finish
{

}

node_update
{

}

shader_evaluate
{
   AtVector result(0.0f, 0.0f, 1.0f);

   float amount = clamp(AiShaderEvalParamFlt(p_amount), 0.0f, 1.0f);
   float size = AiShaderEvalParamFlt(p_size);
   float divergence = clamp(AiShaderEvalParamFlt(p_divergence), 0.0f, 1.0f);
   int space = AiShaderEvalParamInt(p_space);

   AtPoint P = sg->P;
   if (AiNodeIsLinked(node, "P"))
   {
      P = AiShaderEvalParamPnt(p_P);
   }

   P /= size;

   // get a cellular id
    uint32_t id;
   AiCellular(P, 1, 1.0f, 2.0f, 1.0f, NULL, NULL, &id);

   // get two random numbers based on the id
   float u1 = sampleTEAFloat(id, 0);
   if (u1 < amount)
   {
      u1 /= amount;
      float u2 = sampleTEAFloat(id, 1);

      // get a random vector in the hemisphere
      AtVector d = uniformSampleHemisphere(u1, u2);
      // swap y and z to make a normal map
      std::swap(d.y, d.z);

      // blend it in 
      result = lerp(result, d, divergence);
   }

   if (space == S_WORLD)
   {
      // build a local tangent frame to transform the normals
      AtVector U, V;
      if (!AiV3IsSmall(sg->dPdu) && AiV3Exists(sg->dPdu))
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

      sg->out.VEC() = U * result.x + V * result.y + sg->Nf * result.z;
   }
   else
   {
      sg->out.VEC() = result * .5f + AtVector(0.5f, 0.5f, 0.5f);
   }
}


