#include <ai.h>
#include <cstring>
#include "Remap.h"

AI_SHADER_NODE_EXPORT_METHODS(alInputScalar)

enum alInputVectorParams
{
   p_input,
   p_userName,
   REMAP_FLOAT_PARAM_ENUM
};

enum Inputs
{
   IN_facingratio = 0,
   IN_area,
   IN_fi,
   IN_Rl,
   IN_Rr,
   IN_USER
};

static const char* InputNames[] = {
   "facing-ratio",
   "area",
   "face-index",
   "ray-length",
   "ray-depth",
   "User",
   NULL
};

node_parameters
{
   AiParameterEnum("input", IN_facingratio, InputNames);
   AiParameterStr("userName", "");
   REMAP_FLOAT_PARAM_DECLARE;
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alInputScalar;
   node->output_type = AI_TYPE_FLOAT;
   node->name        = "alInputScalar";
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
   int input = AiShaderEvalParamInt(p_input);
   const char* userName = AiShaderEvalParamStr(p_userName);
   float result = 0.0f;

   // first select the input vector to use
   switch(input)
   {
   case IN_facingratio:
      result = AiV3Dot(sg->Nf, -sg->Rd);
      break;
   case IN_area:
      result = AiShaderGlobalsArea(sg);
      break;
   case IN_fi:
      result = float(sg->fi);
      break;
   case IN_Rl:
      result = sg->Rl;
      break;
   case IN_Rr:
       result = float(sg->bounces);
      break;
   case IN_USER:
      AiUDataGetFlt(userName, &result);
      break;
   default:
      break;
   }

   RemapFloat r = REMAP_FLOAT_CREATE;
   result = r.remap(result);

   sg->out.FLT() = result;
}


