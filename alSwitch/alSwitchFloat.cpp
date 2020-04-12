#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alSwitchFloatMtd)

enum alSwitchParams
{
    p_inputA,
    p_inputB,
    p_inputC,
    p_inputD,
    p_inputE,
    p_inputF,
    p_inputG,
    p_inputH,
    p_mix,
    p_threshold
};

node_parameters
{
    AiParameterFlt("inputA", 0.0f);
    AiParameterFlt("inputB", 1.0f);
    AiParameterFlt("inputC", .15f);
    AiParameterFlt("inputD", .30f);
    AiParameterFlt("inputE", .45f);
    AiParameterFlt("inputF", .60f);
    AiParameterFlt("inputG", .75f);
    AiParameterFlt("inputH", .90f);
    AiParameterFlt("mix", 1.0f);
    AiParameterFlt("threshold", 0.5f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alSwitchFloatMtd;
   node->output_type = AI_TYPE_FLOAT;
   node->name        = "alSwitchFloat";
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
    float mix = AiShaderEvalParamFlt(p_mix);
    float threshold = AiShaderEvalParamFlt(p_threshold);

    int input = floorf(mix);
    if (mix-input >= threshold) input++;
    input = clamp(input, 0, 7);

    float result = 0.0f;

    switch(input)
    {
    case 0:
        result = AiShaderEvalParamFlt(p_inputA);
        break;
    case 1:
        result = AiShaderEvalParamFlt(p_inputB);
        break;
    case 2:
        result = AiShaderEvalParamFlt(p_inputC);
        break;
    case 3:
        result = AiShaderEvalParamFlt(p_inputD);
        break;
    case 4:
        result = AiShaderEvalParamFlt(p_inputE);
        break;
    case 5:
        result = AiShaderEvalParamFlt(p_inputF);
        break;
    case 6:
        result = AiShaderEvalParamFlt(p_inputG);
        break;
    case 7:
        result = AiShaderEvalParamFlt(p_inputH);
        break;
    default:
        // should never get here
        result = 0.0f;
        break;
    }

    sg->out.FLT = result;
}


