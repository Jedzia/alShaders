#include <ai.h>
#include <cstring>

AI_SHADER_NODE_EXPORT_METHODS(alInputVector)

enum alInputVectorParams
{
	p_input,
	p_type,
	p_userName,
	p_vector,
	p_matrix,
	p_coordinates
};

enum Inputs
{
	IN_P = 0,
	IN_PO,
	IN_N,
	IN_Nf,
	IN_Ng,
	IN_Ngf,
	IN_Ns,
	IN_dPdu,
	IN_dPdv,
	IN_Ld,
	IN_Rd,
	IN_UV,
	IN_USER,
	IN_CUSTOM
};

static const char* InputNames[] = {
	"P",
	"Po",
	"N",
	"Nf",
	"Ng",
	"Ngf",
	"Ns",
	"dPdu",
	"dPdv",
	"Ld",
	"Rd",
	"uv",
	"User",
	"Custom",
	NULL
};

enum Types
{
	T_POINT = 0,
	T_VECTOR = 1
};

static const char* TypeNames[] =
{
	"Point",
	"Vector",
	NULL
};

enum Coordinates
{
	CS_CARTESIAN = 0,
	CS_SPHERICAL,
	CS_NORM_SPHERICAL
};

static const char* coordinatesNames[] = 
{
	"cartesian",
	"spherical",
	"normalized spherical",
	NULL
};

node_parameters
{
	AiParameterEnum("input", IN_P, InputNames);
	AiParameterEnum("type", T_POINT, TypeNames);
	AiParameterStr("userName", "");
	AiParameterVec("vector", 0.0f, 0.0f, 0.0f);
	auto mtx = AiM4Identity();
	AiParameterMtx("matrix", mtx);
	AiParameterEnum("coordinates", CS_CARTESIAN, coordinatesNames);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alInputVector;
   node->output_type = AI_TYPE_VECTOR;
   node->name        = "alInputVector";
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
	int input = AiShaderEvalParamInt(p_input);
	int type = AiShaderEvalParamInt(p_type);
	AtString userName = AiShaderEvalParamStr(p_userName);
	AtVector vector = AiShaderEvalParamVec(p_vector);
	AtMatrix& mtx = *(AiShaderEvalParamMtx(p_matrix));
	int coordinates = AiShaderEvalParamInt(p_coordinates);
	AtVector result;

	// first select the input vector to use
	switch(input)
	{
	case IN_P:
		vector = sg->P;
		break;
	case IN_PO:
		vector = sg->Po;
		break;
	case IN_N:
		vector = sg->N;
		break;
	case IN_Nf:
		vector = sg->Nf;
		break;
	case IN_Ng:
		vector = sg->Ng;
		break;
	case IN_Ngf:
		vector = sg->Ngf;
		break;
	case IN_Ns:
		vector = sg->Ns;
		break;
	case IN_dPdu:
		vector = sg->dPdu;
		break;
	case IN_dPdv:
		vector = sg->dPdv;
		break;
	case IN_Ld:
		{
			AtLightSample ls;
			AiLightsGetSample(sg, ls);
			vector = ls.Ld;
		}
		break;
	case IN_Rd:
		vector = sg->Rd;
		break;
	case IN_UV:
		vector = AtVector(sg->u, sg->v, 0);
		break;
	case IN_USER:
		AiUDataGetVec(userName, vector);
		break;
	default:
		break;
	}

	// now perform the appropriate xform based on the type
	if (type == T_POINT)
	{
		result = AiM4PointByMatrixMult(mtx, vector);
	}
	else
	{
		result = AiM4VectorByMatrixMult(mtx, vector);
	}

	// convert to spherical coordinates
	if (coordinates > CS_CARTESIAN)
	{
		float theta = acosf(result.z);
		float phi = atan2f(result.z, result.x);

		if (coordinates == CS_NORM_SPHERICAL)
		{
			theta /= AI_PI;
			phi = phi / AI_PITIMES2 + 0.5f;
		}

		result.x = phi;
		result.y = theta;
		result.z = 0.0f;
	}

	sg->out.VEC() = result;
}


