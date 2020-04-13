#include <ai.h>
#include <alUtil.h>

#define NUM_LAYERS 8

AI_SHADER_NODE_EXPORT_METHODS(alLayerFloat)

const char* param_names[] = 
{
	"layer1",
	"layer1a",
	"layer1name",
	"layer1enabled",
	"layer2",
	"layer2a",
	"layer2name",
	"layer2enabled",
	"layer3",
	"layer3a",
	"layer3name",
	"layer3enabled",
	"layer4",
	"layer4a",
	"layer4name",
	"layer4enabled",
	"layer5",
	"layer5a",
	"layer5name",
	"layer5enabled",
	"layer6",
	"layer6a",
	"layer6name",
	"layer6enabled",
	"layer7",
	"layer7a",
	"layer7name",
	"layer7enabled",
	"layer8",
	"layer8a",
	"layer8name",
	"layer8enabled"
};

node_parameters
{
	AiParameterFlt("layer1", 0.0f);
	AiParameterFlt("layer1a", 0.0f);
	AiParameterStr("layer1name", "");
	AiParameterBool("layer1enabled", true);
	AiParameterFlt("layer2", 0.0f);
	AiParameterFlt("layer2a", 0.0f);
	AiParameterStr("layer2name", "");
	AiParameterBool("layer2enabled", true);
	AiParameterFlt("layer3", 0.0f);
	AiParameterFlt("layer3a", 0.0f);
	AiParameterStr("layer3name", "");
	AiParameterBool("layer3enabled", true);
	AiParameterFlt("layer4", 0.0f);
	AiParameterFlt("layer4a", 0.0f);
	AiParameterStr("layer4name", "");
	AiParameterBool("layer4enabled", true);
	AiParameterFlt("layer5", 0.0f);
	AiParameterFlt("layer5a", 0.0f);
	AiParameterStr("layer5name", "");
	AiParameterBool("layer5enabled", true);
	AiParameterFlt("layer6", 0.0f);
	AiParameterFlt("layer6a", 0.0f);
	AiParameterStr("layer6name", "");
	AiParameterBool("layer6enabled", true);
	AiParameterFlt("layer7", 0.0f);
	AiParameterFlt("layer7a", 0.0f);
	AiParameterStr("layer7name", "");
	AiParameterBool("layer7enabled", true);
	AiParameterFlt("layer8", 0.0f);
	AiParameterFlt("layer8a", 0.0f);
	AiParameterStr("layer8name", "");
	AiParameterBool("layer8enabled", true);
}

node_loader
{
	if (i>0) return false;
	node->methods     = alLayerFloat;
	node->output_type = AI_TYPE_FLOAT;
	node->name        = "alLayerFloat";
	node->node_type   = AI_NODE_SHADER;
	strcpy(node->version, AI_VERSION);
	return true;
}

struct ShaderData
{
	int max_layer;
};

node_initialize
{
	ShaderData* data = (ShaderData*) AiMalloc(sizeof(ShaderData));
	AiNodeSetLocalData(node, data);
}

node_finish
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	AiFree(data);
	AiNodeSetLocalData(node, NULL);
}

#define PARAMS_PER_LAYER 4

node_update
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);

	data->max_layer = -1;
	// check to see what the top layer should be
	for (int i=0; i < NUM_LAYERS; ++i)
	{
		// if the alpha is either linked or non-zero, layer is active
		if (AiNodeIsLinked(node, param_names[i*PARAMS_PER_LAYER+1]) ||
                AiNodeGetFlt(node, param_names[i*PARAMS_PER_LAYER+1]) != 0.0f)
		{
			data->max_layer = i;
		}
	}

}

shader_evaluate
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);

	float result = 0.0f;

	for (int i=0; i <= data->max_layer; ++i)
	{
		if (!AiShaderEvalParamBool(i*PARAMS_PER_LAYER+3)) continue;
		float layerVal = AiShaderEvalParamFlt(i*PARAMS_PER_LAYER);
		float layerAlpha = AiShaderEvalParamFlt(i*PARAMS_PER_LAYER+1);
		result = lerp(result, layerVal, layerAlpha);
	}

	sg->out.FLT() = result;
}
