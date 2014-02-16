#include <stdlib.h>

#include "Remap.h"

#include <ai.h>
AI_SHADER_NODE_EXPORT_METHODS(alJitterColorMtd)


enum alTriplanarParams
{
    p_input,
    p_minSaturation,
    p_maxSaturation,
    p_minHueOffset,
    p_maxHueOffset,
    p_minGain,
    p_maxGain,
    p_clamp,
    p_signal
};

node_parameters
{
    AiParameterRGB("input", 1.0, 1.0, 1.0);
    AiParameterFLT("minSaturation", 0.0);
    AiParameterFLT("maxSaturation", 1.0);
    AiParameterFLT("minHueOffset", -.1);
    AiParameterFLT("maxHueOffset", 0.1);
    AiParameterFLT("minGain", 0.5);
    AiParameterFLT("maxGain", 1.5);
    AiParameterBOOL("clamp", true);
    AiParameterFLT("signal", 0.f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alJitterColorMtd;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alJitterColor";
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
    AtRGB input = AiShaderEvalParamRGB(p_input);

    float signal = AiShaderEvalParamFlt(p_signal);
    float minSat = AiShaderEvalParamFlt(p_minSaturation);
    float maxSat = AiShaderEvalParamFlt(p_maxSaturation);
    float minHue = AiShaderEvalParamFlt(p_minHueOffset);
    float maxHue = AiShaderEvalParamFlt(p_maxHueOffset);
    float minGain = AiShaderEvalParamFlt(p_minGain);
    float maxGain = AiShaderEvalParamFlt(p_maxGain);

    bool clamp = AiShaderEvalParamBool(p_clamp);

    AtRGB result = input;

    float saturation = fmodf((signal + 0.9311) * 51731.13215, 1.f);
    saturation = lerp(minSat, maxSat, saturation);
    if (saturation != 1.0f)
    {
        float l = luminance(result);
        result = lerp(rgb(l), result, saturation);
    }

    // hue
    float hueOffset = fmodf((signal + 1.3311) * 173.1231, 1.f);
    hueOffset = lerp(minHue, maxHue, hueOffset);
    if (hueOffset != 0.0f)
    {
        AtRGB hsv = rgb2hsv(result);
        hsv.r += hueOffset;
        result = hsv2rgb(hsv);
    }

    float gain = fmodf((signal + 0.65416) * 413.7254, 1.f);
    gain = lerp(minGain, maxGain, gain);
    result *= gain;

    if(clamp){
        result = AiColorClamp(result, 0.f, 1.f);
    }

    sg->out.RGB = result;
}
