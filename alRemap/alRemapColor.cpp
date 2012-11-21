#include "alUtil.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alRemapColorMtd)

enum alRemapParams
{
	p_input,
	p_gamma,
	p_saturation,
	p_hueOffset,
	p_contrast,
	p_contrastPivot,
	p_contrastSoftClip,
	p_scale,
	p_exposure,
	p_mask,
};

node_parameters
{
	AiParameterRGB("input", 0.18f, 0.18f, 0.18f);
	AiParameterFLT("gamma", 1.0f);
	AiParameterFLT("saturation", 1.0f);
	AiParameterFLT("hueOffset", 0.0f);
	AiParameterFLT("contrast", 1.0f);
	AiParameterFLT("contrastPivot", 0.18f);
	AiParameterFLT("contrastSoftClip", 0.0f);
	AiParameterFLT("scale", 1.0f);
	AiParameterFLT("exposure", 0.f);
	AiParameterFLT("mask", 1.0f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alRemapColorMtd;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alRemapColor";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return TRUE;
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
	AtFloat gamma = AiShaderEvalParamFlt(p_gamma);
	AtFloat saturation = AiShaderEvalParamFlt(p_saturation);
	AtFloat hueOffset = AiShaderEvalParamFlt(p_hueOffset);
	AtFloat contrastVal = AiShaderEvalParamFlt(p_contrast);
	AtFloat contrastPivot = AiShaderEvalParamFlt(p_contrastPivot);
	AtFloat contrastSoftClip = AiShaderEvalParamFlt(p_contrastSoftClip);
	AtFloat scale = AiShaderEvalParamFlt(p_scale);
	AtFloat exposure = AiShaderEvalParamFlt(p_exposure);
	AtFloat mask = AiShaderEvalParamFlt(p_mask);

	AtRGB result = input;
	if (mask > 0.0f)
	{
		// gamma
		result = pow(input, 1.0f/gamma);

		// saturation
		if (saturation != 1.0f)
		{
			float l = luminance(result);
			result = lerp(rgb(l), result, saturation);
		}

		// hue
		if (hueOffset != 0.0f)
		{
			AtRGB hsv = rgb2hsv(result);
			hsv.r += hueOffset;
			hsv.r = fmod(hueOffset, 1.0f);
			result = hsv2rgb(hsv);
		}

		// contrast
		if (contrastVal != 1.0f)
		{
			result = contrast(result, contrastVal, contrastPivot, contrastSoftClip);
		}

		// gain and exposure
		result = result * powf(2.0f, exposure) * scale;

		// mask
		if (mask < 1.0f)
		{
			result = lerp(input, result, mask);
		}
	}
	sg->out.RGB = result;
}


