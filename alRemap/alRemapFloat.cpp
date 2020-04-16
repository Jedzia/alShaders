#include "alUtil.h"
#include "Remap.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alRemapFloatMtd)

enum alRemapParams
{
	p_input,
	REMAP_FLOAT_PARAM_ENUM,
   p_mask
};

node_parameters
{
	AiParameterFlt("input", 0.0f);
	REMAP_FLOAT_PARAM_DECLARE;
    AiParameterFlt("mask", 1.0f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alRemapFloatMtd;
   node->output_type = AI_TYPE_FLOAT;
   node->name        = "alRemapFloat";
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
	float input = AiShaderEvalParamFlt(p_input);
   float mask = AiShaderEvalParamFlt(p_mask);

   float result = input;
   if (mask > 0.0f)
   {
	  RemapFloat r = REMAP_FLOAT_CREATE;
     result = lerp(input, r.remap(input), mask);
   }

	sg->out.FLT() = result;
}


