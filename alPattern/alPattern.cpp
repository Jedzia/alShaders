#include "Remap.h"
#include <ai.h>


AI_SHADER_NODE_EXPORT_METHODS(alPattern)

enum PatternSpaceEnum
{
	NS_WORLD = 0,
	NS_OBJECT,
	NS_PREF,
	NS_UV
};

static const char* patternSpaceNames[] =
{
	"world",
	"object",
	"Pref",
	"UV",
	NULL
};

enum PatternEnum
{
	PN_SINE = 0,
	PN_SQUARE,
	PN_SAW
};

static const char* patternNames[] =
{
	"sine",
	"square",
	"saw",
	NULL
};

static const char* axisNames[] =
{
	"X",
	"Y",
	"Z",
	NULL
};

enum alPatternParams
{
	p_space,
	p_axis,
	p_shape,
	p_frequency,
	p_offset,

	p_color1,
	p_color2,
	p_P,
	REMAP_FLOAT_PARAM_ENUM
};

node_parameters
{
	AiParameterEnum("space", 0, patternSpaceNames);
	AiParameterEnum("axis", 0, axisNames);
	AiParameterEnum("shape", 0, patternNames);
	AiParameterFlt("frequency", 5.0f);
	AiParameterFlt("offset", 0.0f);
	AiParameterRGB("color1", 0.0f, 0.0f, 0.0f);
	AiParameterRGB("color2", 1.0f, 1.0f, 1.0f);
	AiParameterVec("P", 0.0f, 0.0f, 0.0f);
	REMAP_FLOAT_PARAM_DECLARE;
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alPattern;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alPattern";
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
	int space = AiShaderEvalParamInt(p_space);
	int axis = AiShaderEvalParamInt(p_axis);
	int shape = AiShaderEvalParamInt(p_shape);
	float frequency = AiShaderEvalParamFlt(p_frequency);
	float offset = AiShaderEvalParamFlt(p_offset);
	AtRGB color1 = AiShaderEvalParamRGB(p_color1);
	AtRGB color2 = AiShaderEvalParamRGB(p_color2);

	// choose what space we want to calculate in
	AtPoint P;
	if (AiNodeIsLinked(node, "P"))
	{
		P = AiShaderEvalParamVec(p_P);
	}
	else
	{
		switch (space)
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
			if (!AiUDataGetVec(AtString("Pref"), P))
				P = sg->Po;
			break;
		default:
			P = sg->P;
			break;
		}
	}
	// scale the space
	P *= frequency;

	float result = 0.0f;
	float x = P.x;
	if (axis == 1) x = P.y;
	else if (axis == 2) x = P.z;
	x += offset;
	switch (shape)
	{
	case PN_SINE:
		result = sinf(x);
		break;
	case PN_SQUARE:
		result = sinf(x) > 0.0f ? 1.0f : 0.0f;
		break;
	case PN_SAW:
		result = modulo(x*0.25f, 1.0f);
		break;
	default:
		result = 0.0f;
		break;
	}

	RemapFloat r = REMAP_FLOAT_CREATE;
	result = r.remap(result);

	sg->out.RGB() = AiLerp(result, color1, color2);
}


