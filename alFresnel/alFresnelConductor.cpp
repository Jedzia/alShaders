#include <ai.h>
#include <cstring>
#include "fresnel.h"

AI_SHADER_NODE_EXPORT_METHODS(alFresnelConductor)

node_loader
{
    if (i>0) return 0;
    node->methods     = alFresnelConductor;
    node->output_type = AI_TYPE_RGB;
    node->name        = "alFresnelConductor";
    node->node_type   = AI_NODE_SHADER;
    strcpy(node->version, AI_VERSION);
    return true;
}

enum alFresnelConductorParams
{
	p_material=0,
	p_normalize,
	p_reflectivity,
	p_edgetint
};

static const char* alFresnelConductorMaterialNames[] = 
{
    "aluminium",
    "chrome",
    "copper",
    "gold",
    "silver",
    "platinum",
    "titanium",
    "tungsten",
    "custom",
    NULL
};

node_parameters
{
	AiParameterEnum("material", 0, alFresnelConductorMaterialNames);
	AiParameterBool("normalize", false);
	AiParameterRGB("reflectivity", .94, .78, .37);
    AiMetaDataSetBool(nentry, "reflectivity", "always_linear", true);  // no inverse-gamma correction
	AiParameterRGB("edgetint", 1.0, 0.98, 0.73);
    AiMetaDataSetBool(nentry, "edgetint", "always_linear", true);  // no inverse-gamma correction
}

node_initialize
{

}

node_update
{

}

node_finish
{

}

shader_evaluate
{
	int material = AiShaderEvalParamInt(p_material);
	AtRGB r = AiShaderEvalParamRGB(p_reflectivity);
	AtRGB g = AiShaderEvalParamRGB(p_edgetint);
	FresnelConductor fr(material, r, g);
   AiStateSetMsgRGB("als_fr_r", fr._r);
   AiStateSetMsgRGB("als_fr_g", fr._g);
	sg->out.RGB() = fr.kr(AiV3Dot(-sg->Rd, sg->Nf));
}