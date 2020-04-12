#include <ai.h>
#include <cstring>

#include "../common/Color.h"

AI_SHADER_NODE_EXPORT_METHODS(alColorSpace)

enum alColorSpaceParams
{
	p_input,
	p_sourceSpace
};

enum SourceSpace
{
	SS_SRGB = 0,
	SS_CINEON,
	SS_LOGC
};

static const char* SourceSpaceNames[] = {
	"sRGB",
	"Cineon",
	"LogC",
	NULL
};

node_parameters
{
	AiParameterRGB("input", 0.0f, 0.0f, 0.0f);
	AiParameterEnum("sourceSpace", SS_SRGB, SourceSpaceNames);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alColorSpace;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alColorSpace";
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
	int sourceSpace = AiShaderEvalParamInt(p_sourceSpace);
	
	AtRGB result = input;

	switch(sourceSpace)
	{
	case SS_SRGB:
		result = sRgbToLin(input);
		break;
	case SS_CINEON:
		result = cineonToLin(input);
		break;
	case SS_LOGC:
		result = logCToLin(input);
	default:
		break;
	}

	sg->out.RGB() = result;
}


