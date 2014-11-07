// All PS blend modes (plus some more) nicked from here:
// http://mouaif.wordpress.com/2009/01/05/photoshop-math-with-glsl-shaders/

#include <ai.h>
#include <alUtil.h>

#define NUM_LAYERS 8

AI_SHADER_NODE_EXPORT_METHODS(alLayerColor)

const char* param_names[] = 
{
	"layer1",
	"layer1a",
   "layer1blend",
	"layer2",
	"layer2a",
   "layer2blend",
	"layer3",
	"layer3a",
   "layer3blend",
	"layer4",
	"layer4a",
   "layer4blend",
	"layer5",
	"layer5a",
   "layer5blend",
	"layer6",
	"layer6a",
   "layer6blend",
	"layer7",
	"layer7a",
   "layer7blend",
	"layer8",
	"layer8a",
   "layer8blend"
};

enum BlendModes
{
   BM_NORMAL=0,
   BM_LIGHTEN,
   BM_DARKEN,
   BM_MULTIPLY,
   BM_AVERAGE,
   BM_ADD,
   BM_SUBTRACT,
   BM_DIFFERENCE,
   BM_NEGATION,
   BM_EXCLUSION,
   BM_SCREEN,
   BM_OVERLAY,
   BM_SOFTLIGHT,
   BM_HARDLIGHT,
   BM_COLORDODGE,
   BM_COLORBURN,
   BM_LINEARDODGE,
   BM_LINEARBURN,
   BM_LINEARLIGHT,
   BM_VIVIDLIGHT,
   BM_PINLIGHT,
   BM_HARDMIX,
   BM_REFLECT,
   BM_GLOW,
   BM_PHOENIX

};

const char* BlendModeNames[] = 
{
   "Normal",
   "Lighten",
   "Darken",
   "Multiply",
   "Average",
   "Add",
   "Subtract",
   "Difference",
   "Negation",
   "Exclusion",
   "Screen",
   "Overlay",
   "Soft Light",
   "Hard Light",
   "Color Dodge",
   "Color Burn",
   "Linear Dodge",
   "Linear Burn",
   "Linear Light",
   "Vivid Light",
   "Pin Light",
   "Hard Mix",
   "Reflect",
   "Glow",
   "Phoenix",
   NULL
};

node_parameters
{
	AiParameterRGB("layer1", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer1a", 0.0f);
   AiParameterEnum("layer1blend", 0, BlendModeNames);
	AiParameterRGB("layer2", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer2a", 0.0f);
   AiParameterEnum("layer2blend", 0, BlendModeNames);
	AiParameterRGB("layer3", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer3a", 0.0f);
   AiParameterEnum("layer3blend", 0, BlendModeNames);
	AiParameterRGB("layer4", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer4a", 0.0f);
   AiParameterEnum("layer4blend", 0, BlendModeNames);
	AiParameterRGB("layer5", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer5a", 0.0f);
   AiParameterEnum("layer5blend", 0, BlendModeNames);
	AiParameterRGB("layer6", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer6a", 0.0f);
   AiParameterEnum("layer6blend", 0, BlendModeNames);
	AiParameterRGB("layer7", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer7a", 0.0f);
   AiParameterEnum("layer7blend", 0, BlendModeNames);
	AiParameterRGB("layer8", 0.0f, 0.0f, 0.0f);
	AiParameterFlt("layer8a", 0.0f);
   AiParameterEnum("layer8blend", 0, BlendModeNames);
}

node_loader
{
	if (i>0) return false;
	node->methods     = alLayerColor;
	node->output_type = AI_TYPE_RGB;
	node->name        = "alLayerColor";
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

node_update
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);
	
	data->max_layer = -1;
	// check to see what the top layer should be
	for (int i=0; i < NUM_LAYERS; ++i)
	{
		// if the alpha is either linked or non-zero, layer is active
		if (AiNodeIsLinked(node, param_names[i*3+1]) ||
			params[i*3+1].FLT != 0.0f)
		{
			data->max_layer = i;
		}
	}

}

inline AtRGB overlay(const AtRGB& l1, const AtRGB& l2)
{
   AtRGB result;
   for (int i=0; i < 3; ++i)
   {
      l1[i] < 0.5f ? result[i] = 2.0f * l1[i] * l2[i] : result[i] = (1.0f - 2.0f * (1.0f-l1[i]) * (1.0f-l2[i]));
   }

   return result;
}

inline AtRGB softlight(const AtRGB& l1, const AtRGB& l2)
{
   AtRGB result;
   for (int i=0; i < 3; ++i)
   {
      if (l2[i] < 0.5f)
         result[i] = 2.0f * l1[i] * l2[i] + SQR(l1[i]) * (1.0f - 2.0f * l2[i]);
      else 
         result[i] = sqrtf(l1[i]) * (2.0f * l2[i] - 1.0f) + 2.0f * l1[i] * (1.0f - l2[i]);
   }

   return result;
}

inline float colordodgef(float l1, float l2)
{
   if (l2 == 1.0f)
      return l2;
   else 
      return std::min(l1 / (1.0f - l2), 1.0f);
}

inline AtRGB colordodge(const AtRGB& l1, const AtRGB& l2)
{
   AtRGB result;
   for (int i=0; i < 3; ++i)
   {
      result[i] = colordodgef(l1[i], l2[i]);
   }

   return result;
}

inline float colorburnf(float l1, float l2)
{
   if (l2 == 0.0f)
      return l2;
   else 
      return std::max(1.0f - (1.0f - l1) / l2, 0.0f);
}

inline AtRGB colorburn(const AtRGB& l1, const AtRGB& l2)
{
   AtRGB result;
   for (int i=0; i < 3; ++i)
   {
      result[i] = colorburnf(l1[i], l2[i]);
   }

   return result;
}

inline AtRGB linearlight(const AtRGB& l1, const AtRGB& l2)
{
   AtRGB result;
   for (int i=0; i < 3; ++i)
   {
      if (l2[i] < 0.5f)
         result[i] = std::max(l1[i] + 2.0f * l2[i] - 1.0f, 0.0f);
      else
         result[i] = std::min(l1[i] + 2.0f * (l2[i] - 0.5f), 1.0f);
   }

   return result;
}

inline AtRGB vividlight(const AtRGB& l1, const AtRGB& l2)
{
   AtRGB result;
   for (int i=0; i < 3; ++i)
   {
      if (l2[i] < 0.5f)
         result[i] = colorburnf(l1[i], 2.0f * l2[i]);
      else
         result[i] = colordodgef(l1[i], 2.0f * (l2[i] - 0.5f));
   }

   return result;
}

inline AtRGB pinlight(const AtRGB& l1, const AtRGB& l2)
{
   AtRGB result;
   for (int i=0; i < 3; ++i)
   {
      if (l2[i] < 0.5f)
         result[i] = std::min(l1[i], 2.0f * l2[i]);
      else
         result[i] = std::max(l1[i], 2.0f * (l2[i] - 0.5f));
   }

   return result;
}

inline AtRGB hardmix(const AtRGB& l1, const AtRGB& l2)
{
   AtRGB result = vividlight(l1, l2);
   for (int i=0; i < 3; ++i)
   {
      result[i] < 0.5f ? result[i] = 0.0f : result[i] = 1.0f;
   }

   return result;
}

inline AtRGB reflect(const AtRGB& l1, const AtRGB& l2)
{
   AtRGB result;
   for (int i=0; i < 3; ++i)
   {
      if (l2[i] == 1.0f)
         result[i] = l2[i];
      else
         result[i] = std::min(SQR(l1[i]) / (1.0f - l2[i]), 1.0f);
   }

   return result;
}

AtRGB blend(const AtRGB& l1, const AtRGB& l2, float a, int mode)
{
   AtRGB result;
   switch(mode)
   {
   case BM_LIGHTEN:
      result = max(l1, l2);
      break;
   case BM_DARKEN:
      result = min(l1, l2);
      break;
   case BM_MULTIPLY:
      result = l1 * l2;
      break;
   case BM_AVERAGE:
      result = (l1 + l2) * 0.5f;
      break;
   case BM_ADD:
   case BM_LINEARDODGE:
      result = min(l1 + l2, AI_RGB_WHITE);
      break;
   case BM_SUBTRACT:
   case BM_LINEARBURN:
      result = max(l1 + l2 - AI_RGB_WHITE, AI_RGB_BLACK);
      break;
   case BM_DIFFERENCE:
      result = fabs(l1 - l2);
      break;
   case BM_NEGATION:
      result = AI_RGB_WHITE - fabs(AI_RGB_WHITE - l1 - l2);
      break;
   case BM_EXCLUSION:
      result = l1 + l2 - (2.0f * l1 * l2);
      break;
   case BM_SCREEN:
      result = AI_RGB_WHITE - ((AI_RGB_WHITE-l1) * (AI_RGB_WHITE-l2));
      break;
   case BM_OVERLAY:
      result = overlay(l1, l2);
      break;
   case BM_SOFTLIGHT:
      result = softlight(l1, l2);
      break;
   case BM_HARDLIGHT:
      result = overlay(l2, l1);
      break;
   case BM_COLORDODGE:
      result = colordodge(l1, l2);
      break;
   case BM_COLORBURN:
      result = colorburn(l1, l2);
      break;
   case BM_LINEARLIGHT:
      result = linearlight(l1, l2);
      break;
   case BM_VIVIDLIGHT:
      result = vividlight(l1, l2);
      break;
   case BM_PINLIGHT:
      result = pinlight(l1, l2);
      break;
   case BM_HARDMIX:
      result = hardmix(l1, l2);
      break;
   case BM_REFLECT:
      result = reflect(l1, l2);
      break;
   case BM_GLOW:
      result = reflect(l2, l1);
      break;
   case BM_PHOENIX:
      result = min(l1, l2) - max(l1, l2) + AI_RGB_WHITE;
      break;
   default:
      result = l2;
      break;
   }

   return lerp(l1, result, a);
}

shader_evaluate
{
	ShaderData* data = (ShaderData*)AiNodeGetLocalData(node);

	AtRGB result = AI_RGB_BLACK;

	for (int i=0; i <= data->max_layer; ++i)
	{
		AtRGB layerVal = AiShaderEvalParamRGB(i*3);
		float layerAlpha = AiShaderEvalParamFlt(i*3+1);
      int layerBlend = AiShaderEvalParamInt(i*3+2);
		result = blend(result, layerVal, layerAlpha, layerBlend);
	}

	sg->out.RGB = result;
}