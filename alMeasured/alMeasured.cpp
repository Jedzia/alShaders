#include <ai.h>
#include <cstring>
#include <string>

#include "../common/alUtil.h"
#include "brdf_read.h"

/*
Shader to read the MERL BRDF Database from
"A Data-Driven Reflectance Model",
Wojciech Matusik, Hanspeter Pfister, Matt Brand and Leonard McMillan,
ACM Transactions on Graphics 22, 3(2003), 759-769.
BibTeX:
@article {Matusik:2003, 
	author = "Wojciech Matusik and Hanspeter Pfister and Matt Brand and Leonard McMillan",
	title = "A Data-Driven Reflectance Model",
	journal = "ACM Transactions on Graphics",
	year = "2003",
	month = jul,
	volume = "22",
	number = "3",
	pages = "759-769"
}
*/

AI_SHADER_NODE_EXPORT_METHODS(alMeasured)

struct ShaderData
{
	ShaderData() : brdf(NULL), ok(false) {}
	std::string filename;
	double* brdf;
	bool ok;
};

enum alMeasuredParams
{
	p_filename
};

node_parameters
{
	AiParameterStr("filename", "");
}


node_initialize
{
	ShaderData* data = new ShaderData;
    AiNodeSetLocalData(node,data);
}

node_finish
{
	if (AiNodeGetLocalData(node))
    {
        ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);
        AiNodeSetLocalData(node, NULL);
        delete data;
    }
}

node_update
{
	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
	std::string newfn = std::string(params[p_filename].STR);
	// if the filename has changed, reload the data
	if (newfn != data->filename)
	{
		data->filename = newfn;
		if (data->brdf) free(data->brdf);
		data->ok = read_brdf(newfn.c_str(), data->brdf);
		if (!data->ok) 
		{
			data->brdf = NULL;
			std::cerr << "ERROR loading brdf file \"" << newfn << "\"" << std::endl;
		}
		else
		{
			std::cout << "OK loading brdf file \"" << newfn << "\"" << std::endl;	
		}
	}
}

shader_evaluate
{
	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

	AtRGB result = rgb(0.0f, 1.0f, 1.0f);

	if (data->ok)
	{
		AtVector U, V;
        AiV3BuildLocalFramePolar(U, V, sg->Nf);
		AtVector wo = -sg->Rd;
		float theta_out = acosf(AiV3Dot(wo, sg->Nf));
		float phi_out = sphericalPhi(wo, U, V);

		result = AI_RGB_BLACK;
		AiLightsPrepare(sg);
        AtLightSample ls;
        while (AiLightsGetSample(sg, ls))
		{
			float cos_theta_in = AiV3Dot(sg->Ld, sg->Nf);
			float theta_in = acosf(cos_theta_in);
			float phi_in = sphericalPhi(sg->Ld, U, V);
			double r, g, b;
			lookup_brdf_val(data->brdf, theta_in, phi_in, theta_out, phi_out, r, g, b);
			result += sg->Li * JedPortGetShaderWeight(sg) * rgb(float(r), float(g), float(b)) * cos_theta_in;
			//result += sg->Li * sg->we * AiV3Dot(sg->Ld, sg->Nf) * AI_ONEOVERPI * 0.18f;
		}
	}

	sg->out.RGB() = result;
}

node_loader
{
	if (i>0) return 0;
   	node->methods     = alMeasured;
   	node->output_type = AI_TYPE_RGB;
   	node->name        = "alMeasured";
   	node->node_type   = AI_NODE_SHADER;
   	strcpy(node->version, AI_VERSION);
   	return true;
}