#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alSwitchColorMtd)

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
    AiParameterRGB("inputA", 0.f, 0.f, 0.f);
    AiParameterRGB("inputB", 1.f, 1.f, 1.f);
    AiParameterRGB("inputC", 1.f, 0.f, 0.f);
    AiParameterRGB("inputD", 0.f, 1.f, 0.f);
    AiParameterRGB("inputE", 0.f, 0.f, 1.f);
    AiParameterRGB("inputF", 1.f, 1.f, 0.f);
    AiParameterRGB("inputG", 1.f, 0.f, 1.f);
    AiParameterRGB("inputH", 0.f, 1.f, 1.f);
    AiParameterFLT("mix", 1.0f);
    AiParameterFLT("threshold", 0.5f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alSwitchColorMtd;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alSwitchColor";
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

    AtRGB result = AI_RGB_BLACK;

    switch(input)
    {
    case 0:
        result = AiShaderEvalParamRGB(p_inputA);
        break;
    case 1:
        result = AiShaderEvalParamRGB(p_inputB);
        break;
    case 2:
        result = AiShaderEvalParamRGB(p_inputC);
        break;
    case 3:
        result = AiShaderEvalParamRGB(p_inputD);
        break;
    case 4:
        result = AiShaderEvalParamRGB(p_inputE);
        break;
    case 5:
        result = AiShaderEvalParamRGB(p_inputF);
        break;
    case 6:
        result = AiShaderEvalParamRGB(p_inputG);
        break;
    case 7:
        result = AiShaderEvalParamRGB(p_inputH);
        break;
    default:
        // should never get here
        result = AI_RGB_BLACK;
        break;
    }

    sg->out.RGB() = result;
}


