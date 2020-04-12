#include "Remap.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alFractal);

enum NoiseSpaceEnum
{
	NS_WORLD = 0,
	NS_OBJECT,
	NS_PREF,
	NS_UV
};

static const char* noiseSpaceNames[] =
{
		"world",
		"object",
		"Pref",
		"UV",
		NULL
};

enum NoiseMode
{
	NM_SCALAR=0,
	NM_VECTOR
};

static const char* noiseModeNames[] = 
{
	"scalar",
	"vector",
	NULL
};

enum alNoiseParams
{
	p_mode,
	p_space,
	p_frequency,
	p_scale,
	p_time,
	p_octaves,
	p_distortion,
	p_lacunarity,
	p_gain,
	p_turbulent,
	REMAP_FLOAT_PARAM_ENUM,
	p_color1,
	p_color2,
	p_P
};

node_parameters
{
	AiParameterENUM("mode", 0, noiseModeNames);
	AiParameterENUM("space", 0, noiseSpaceNames);
	AiParameterFLT("frequency", 1.0f);
	AiParameterPNT("scale", 1.0f, 1.0f, 1.0f);
	AiParameterFLT("time", 0.0f);
	AiParameterINT("octaves", 8);
	AiParameterFLT("distortion", 0.0f);
	AiParameterFLT("lacunarity", 2.121f);
	AiParameterFLT("gain", 0.5f);
	AiParameterBOOL("turbulent", false);
	REMAP_FLOAT_PARAM_DECLARE;
	AiParameterRGB("color1", 0.0f, 0.0f, 0.0f);
	AiParameterRGB("color2", 1.0f, 1.0f, 1.0f);
	AiParameterPnt("P", 0.0f, 0.0f, 0.0f);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alFractal;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alFractal";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return true;
}

struct ShaderData
{
	int mode;
	int space;
	int octaves;
	bool turbulent;
	bool ridged;
};

node_initialize
{
	ShaderData* data = new ShaderData;
	AiNodeSetLocalData(node, data);
}

node_finish
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	delete data;
}

node_update
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	data->mode = params[p_mode].INT;
	data->space = params[p_space].INT;
	data->octaves = params[p_octaves].INT;
	data->turbulent = params[p_turbulent].BOOL;
}

shader_evaluate
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	float frequency = AiShaderEvalParamFlt(p_frequency);
	AtPoint scale = AiShaderEvalParamVec(p_scale);
	float gain = AiShaderEvalParamFlt(p_gain);
	float lacunarity = AiShaderEvalParamFlt(p_lacunarity);
	float distortion = AiShaderEvalParamFlt(p_distortion);
	float time = AiShaderEvalParamFlt(p_time);
	AtRGB color1 = AiShaderEvalParamRGB(p_color1);
	AtRGB color2 = AiShaderEvalParamRGB(p_color2);
	AtPoint Pin = AiShaderEvalParamPnt(p_P);

	// choose what space we want to calculate in
	AtPoint P;
	if (AiNodeIsLinked(node, "P"))
	{
		P = Pin;
	}
	else
	{
		switch (data->space)
		{
		case NS_OBJECT:
			P = sg->Po;
			break;
		case NS_UV:
			P.x = sg->u;
			P.y = sg->v;
			P.z = 0.0f;
			break;
		case NS_PREF:
			if (!AiUDataGetPnt("Pref", &P))
				P = sg->Po;
			break;
		default:
			P = sg->P;
			break;
		}
	}

	P *= frequency;

	P *= scale;

	if (data->mode == NM_SCALAR)
	{
		float n = 0.0f;
		float amp = 1.0f;
		float weight = 1;
		float v;
		float nrm = 0.0f;
		for (int i=0; i < data->octaves; ++i)
		{
			AtPoint PP = P;
			if (distortion != 0.0f)
				PP += distortion * AiVNoise3(P, 1, 0, 0);
			v = AiPerlin4(PP, time);
			if (data->turbulent) v = fabs(v);
			n += v * amp;
			amp *= gain;

			P *= lacunarity;
		}
		
		RemapFloat r = REMAP_FLOAT_CREATE;
		n = r.remap(n);

		sg->out.RGB() = AiColorLerp(n, color1, color2);
	}
	else
	{
		AtRGB n= rgb(0, 0, 0);
		float amp = 1.0f;
		AtRGB weight = rgb(0, 0, 0);
		AtRGB v;
		for (int i=0; i < data->octaves; ++i)
		{
			AtPoint PP = P;
			if (distortion != 0.0f)
				PP += distortion * AiVNoise3(P, 1, 0, 0);
			v = rgb(AiVNoise4(PP, time, 1, 0, 0));
			if (data->turbulent) v = fabs(v);
			n += v * amp;
			amp *= gain;

			P *= lacunarity;
		}
        RemapFloat r = REMAP_FLOAT_CREATE;
        n.r = r.remap(n.r);
        n.g = r.remap(n.g);
        n.b = r.remap(n.b);

		sg->out.RGB() = n;
	}
}


