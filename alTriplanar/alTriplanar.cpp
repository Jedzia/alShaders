#include "Remap.h"
#include <ai.h>

AI_SHADER_NODE_EXPORT_METHODS(alTriplanar)

enum TriplanarSpaceEnum
{
	NS_WORLD = 0,
	NS_OBJECT,
	NS_PREF,
};

static const char* triplanarSpaceNames[] =
{
	"world",
	"object",
	"Pref",
	NULL
};

enum TriplanarNormalEnum
{
    N_GEOMETRIC = 0,
    N_SMOOTH,
    N_SMOOTHNOBUMP,
};

static const char* triplanarNormalNames[] =
{
    "Geometric",
    "Smooth",
    "SmoothNoBump",
    NULL
};

enum TriplanarTiling
{
    TM_REGULAR = 0,
    TM_CELLNOISE
};

static const char* triplanarTilingNames[] =
{
    "regular",
    "cellnoise",
    NULL
};

enum alTriplanarParams
{
    p_input,
    p_space,
    p_normal,
    p_tiling,
    p_frequency,
    p_texture,
    p_blendSoftness,
    p_cellSoftness,
    p_scalex,
    p_scaley,
    p_scalez,
    p_offsetx,
    p_offsety,
    p_offsetz,
    p_rotx,
    p_roty,
    p_rotz,
    p_rotjitterx,
    p_rotjittery,
    p_rotjitterz,
    p_mipMapBias
};

node_parameters
{
	AiParameterRGB("input", 0.0f, 0.0f, 0.0f);
    AiParameterEnum("space", 0, triplanarSpaceNames);
    AiParameterEnum("normal", 0, triplanarNormalNames);
    AiParameterEnum("tiling", 0, triplanarTilingNames);
    AiParameterFlt("frequency", 1.0f);
    AiParameterStr("texture", "");
    AiParameterFlt("blendSoftness", 0.1);
    AiParameterFlt("cellSoftness", 0.1);
    AiParameterFlt("scalex", 1.0f);
    AiParameterFlt("scaley", 1.0f);
    AiParameterFlt("scalez", 1.0f);
    AiParameterFlt("offsetx", 0.0f);
    AiParameterFlt("offsety", 0.0f);
    AiParameterFlt("offsetz", 0.0f);
    AiParameterFlt("rotx", 0.0f);
    AiParameterFlt("roty", 0.0f);
    AiParameterFlt("rotz", 0.0f);
    AiParameterFlt("rotjitterx", 1.0f);
    AiParameterFlt("rotjittery", 1.0f);
    AiParameterFlt("rotjitterz", 1.0f);
    AiParameterInt("mipMapBias", 0);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alTriplanar;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alTriplanar";
   node->node_type   = AI_NODE_SHADER;
   strcpy(node->version, AI_VERSION);
   return true;
}

struct ShaderData
{
    AtTextureHandle* texturehandle;
    AtTextureParams *textureparams;
};

node_initialize
{
    ShaderData *data = new ShaderData;
    const char *texname = params[p_texture].STR;
    data->texturehandle = AiTextureHandleCreate(texname);
    data->textureparams = new AtTextureParams;
    AiTextureParamsSetDefaults(data->textureparams);
    AiNodeSetLocalData(node, data);
}

node_finish
{
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
    AiTextureHandleDestroy(data->texturehandle);
    delete data->textureparams;
    delete data;
}

node_update
{
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
    data->textureparams->mipmap_bias = params[p_mipMapBias].INT;

}

struct SGCache{
	
	void initCache(const AtShaderGlobals *sg){
		u = sg->u;
		v = sg->v;
		dudx = sg->dudx;
		dudy = sg->dudy;
		dvdx = sg->dvdx;
		dvdy = sg->dvdy;
	}

	void restoreSG(AtShaderGlobals *sg){
		sg->u = u;
		sg->v = v;
		sg->dudx = dudx;
		sg->dudy = dudy;
		sg->dvdx = dvdx;
		sg->dvdy = dvdy;
	}

	float u;
	float v;
	float dudx;
	float dudy;
	float dvdx;
	float dvdy;
};

void getProjectionGeometry(const AtNode* node, const AtShaderGlobals *sg, int space, int normal, AtPoint *P, AtVector *N, AtVector *dPdx, AtVector *dPdy){
    AtVector baseN;
    switch (normal)
    {
    case N_GEOMETRIC:
        baseN = sg->Ng;
        break;
    case N_SMOOTH:
        baseN = sg->N;
        break;
    case N_SMOOTHNOBUMP:
        baseN = sg->Ns;
        break;
    default:
        baseN = sg->N;
        break;
    }

    switch (space)
    {
	case NS_WORLD:
		*P = sg->P;
		*N = baseN;
		*dPdx = sg->dPdx;
		*dPdy = sg->dPdy;
		break;
    case NS_OBJECT:
        *P = sg->Po;
        *N = AiShaderGlobalsTransformNormal(sg, baseN, AI_WORLD_TO_OBJECT);
		*dPdx = AiShaderGlobalsTransformVector(sg, sg->dPdx, AI_WORLD_TO_OBJECT);
		*dPdy = AiShaderGlobalsTransformVector(sg, sg->dPdy, AI_WORLD_TO_OBJECT);
        break;
    case NS_PREF:
        if (!AiUDataGetPnt("Pref", P)){
            AiMsgWarning("[alTriplanar] could not get Pref");
            // TODO: Output warning about not finding the correct data.
            *P = sg->Po;
            *N = AiShaderGlobalsTransformNormal(sg, baseN, AI_WORLD_TO_OBJECT);
			*dPdx = AiShaderGlobalsTransformVector(sg, sg->dPdx, AI_WORLD_TO_OBJECT);
			*dPdy = AiShaderGlobalsTransformVector(sg, sg->dPdy, AI_WORLD_TO_OBJECT);	
        } else {
            AiMsgWarning("[alTriplanar got Pref]");
			AiUDataGetDxyDerivativesPnt("Pref", dPdx, dPdy);
            *N = AiV3Normalize(AiV3Cross(*dPdx, *dPdy));
		}
        break;
    default:
        *P = sg->P;
        *N = baseN;
   		*dPdx = sg->dPdx;
		*dPdy = sg->dPdy;
        break;
    }
}

void computeBlendWeights(const AtVector N, int space, float blendSoftness, float *weights){
    weights[0] = fabsf(N.x);
    weights[1] = fabsf(N.y);
    weights[2] = fabsf(N.z);
    float weightsum = 0.f;
    for(int i=0; i<3; ++i){
        weights[i] = weights[i] - (1.f-blendSoftness)/2.f;
        weights[i] = MAX(weights[i], 0.00f);
        weightsum += weights[i];
    }
    if(weightsum){
        for(int i=0; i<3; ++i){
            weights[i] /= weightsum;
        }
    }
}

void rotateUVs(AtPoint &P, float degrees){
    AtVector orientVectorX;
    const double d2r = 1. / 360. * AI_PI * 2;
    double phi = d2r * degrees;
    orientVectorX.x = cosf(phi);
    orientVectorX.y = sinf(phi);
    orientVectorX.z = 0.f;

    AtVector orientVectorZ;
    orientVectorZ.x = 0.f;
    orientVectorZ.y = 0.f;
    orientVectorZ.z = 1.f;

    AtVector orientVectorY = AiV3Cross(orientVectorX, orientVectorZ);

    AiV3RotateToFrame(P, orientVectorX, orientVectorY, orientVectorZ);
}

AtRGBA tileRegular(const AtPoint &P, const AtVector &dPdx, const AtVector dPdy,
						 const AtPoint &scale, const AtPoint &offset,
                         float *weights, const AtPoint &rot, AtShaderGlobals *sg,
                         AtTextureHandle *handle, AtTextureParams *params){
    AtRGBA textureResult[3];
    bool textureAccessX = false;
    bool textureAccessY = false;
    bool textureAccessZ = false;

    // lookup X
    AtPoint ProjP;
    ProjP.x = (P.z + 123.94 + offset.x) * scale.x ;
    ProjP.y = (P.y + 87.22 + offset.x) * scale.x;
    ProjP.z = 0.;
    rotateUVs(ProjP, rot.x);

    sg->u = ProjP.x;
    sg->v = ProjP.y;
	sg->dudx = dPdx.z * scale.x;
	sg->dudy = dPdy.z * scale.x;
	sg->dvdx = dPdx.y * scale.x;
	sg->dvdy = dPdy.y * scale.x;

    if(weights[0] > 0.){
        textureResult[0] = AiTextureHandleAccess(sg, handle, params, &textureAccessX);
    } else {
        textureResult[0] = AI_RGBA_BLACK;
        textureAccessX = true;
    }

    // lookup Y
    ProjP.x = (P.x + 74.1 + offset.y) * scale.y;
    ProjP.y = (P.z + 9.2 + offset.y) * scale.y;
    ProjP.z = 0.;
    rotateUVs(ProjP, rot.y);

    sg->u = ProjP.x;
    sg->v = ProjP.y;
	sg->dudx = dPdx.x * scale.y;
	sg->dudy = dPdy.x * scale.y;
	sg->dvdx = dPdx.z * scale.y;
	sg->dvdy = dPdy.z * scale.y;

    if(weights[1] > 0.){
        textureResult[1] = AiTextureHandleAccess(sg, handle, params, &textureAccessY);
    } else {
        textureResult[1] = AI_RGBA_BLACK;
        textureAccessY = true;
    }

    // lookup Z
    ProjP.x = (P.x + 123.94 + offset.z) * scale.z;
    ProjP.y = (P.y + 87.22 + offset.z) * scale.z;
    ProjP.z = 0.;
    rotateUVs(ProjP, rot.z);

    sg->u = ProjP.x;
    sg->v = ProjP.y;
 	sg->dudx = dPdx.x * scale.z;
	sg->dudy = dPdy.x * scale.z;
	sg->dvdx = dPdx.y * scale.z;
	sg->dvdy = dPdy.y * scale.z;

    if(weights[2] > 0.){
        textureResult[2] = AiTextureHandleAccess(sg, handle, params, &textureAccessZ);
    } else {
        textureResult[2] = AI_RGBA_BLACK;
        textureAccessZ = true;
    }

    if(textureAccessX && textureAccessY && textureAccessZ){
        AtRGBA result = AI_RGBA_BLACK;
        result += textureResult[0] * weights[0];
        result += textureResult[1] * weights[1];
        result += textureResult[2] * weights[2];
        return result;
    }
    else {
        // Something went wrong during lookup.
        // TODO: Log the error
        return AI_RGBA_RED;
    }
}

bool lookupCellNoise(float u, float v, float dudx, float dudy, float dvdx, float dvdy,
							const float cellSoftness, float rot, float rotjitter,
                            AtShaderGlobals *sg, AtTextureHandle *handle,
                            AtTextureParams *params, AtRGBA *textureResult){
    AtPoint P;
    P.x = u;
    P.y = v;
    P.z = 0.f;

    int samples = (cellSoftness == 0) ? 1 : 3;
    // run cellnoise
    float weights[3];
    float f[3];
    AtVector delta[3];
    uint32_t id[3];
    AiCellular(P, samples, 1, 1.92, 1, f, delta, id);

    if(samples == 1){
        weights[0] = 1.f;
    } else {
        // find closest cell
        float closestDistance = 100000.f;
        float distances[3];
        for(int i=0; i<samples; ++i){
            distances[i] = AiV3Length(delta[i]);
            closestDistance = MIN(distances[i], closestDistance);
        }

        float weightsum = 0.f;
        for(int i=0; i<samples; ++i){
            float diff = distances[i] - closestDistance;
            weights[i] = cellSoftness - diff;
            weights[i] = MAX(0.f, weights[i]);
            weightsum += weights[i];
        }
        if(weightsum){
            for(int i=0; i<samples; ++i){
                weights[i] /= weightsum;
            }
        }
    }

    bool success = false;
    *textureResult = AI_RGBA_BLACK;
	sg->dudx = dudx;
	sg->dudy = dudy;
	sg->dvdx = dvdx;
	sg->dvdy = dvdy;
    for(int i=0; i<samples; ++i){
        if(weights[i] > 0.f){
            // pick direction for orientation
            AtVector orientVectorX;
            double jitter = (random(id[i])-0.5) * rotjitter;
            double phi = modulo(rot/360. + jitter, 1.f) * AI_PI * 2.;
            orientVectorX.x = cosf(phi);
            orientVectorX.y = sinf(phi);
            orientVectorX.z = 0.f;

            AtVector orientVectorZ;
            orientVectorZ.x = 0.f;
            orientVectorZ.y = 0.f;
            orientVectorZ.z = 1.f;

            AtVector orientVectorY = AiV3Cross(orientVectorX, orientVectorZ);

            AiV3RotateToFrame(delta[i], orientVectorX, orientVectorY, orientVectorZ);

            // find new uv coordinates, set the center of the cell to be 0.5/0.5;
            sg->u = delta[i].x * 0.75 - 0.5;
            sg->v = delta[i].y * 0.75 - 0.5;

            // texture lookup
            bool currentSuccess = false;
            *textureResult += AiTextureHandleAccess(sg, handle, params, &success) * weights[i];
            success |= currentSuccess;
        }
    }

    return success;
}

AtRGBA tileCellnoise(const AtPoint &P, const AtVector &dPdx, const AtVector &dPdy,
						   const AtPoint &scale, const AtPoint &offset,
                           float *weights, float cellSoftness,
                           const AtPoint &rot, const AtPoint &rotjitter, AtShaderGlobals *sg,
                           AtTextureHandle *handle, AtTextureParams *params){
    AtRGBA textureResult[3];

    bool textureAccessX = true;
    if(weights[0] > 0.){
        textureAccessX = lookupCellNoise((P.y + offset.x) * scale.x,
                                          (P.z + offset.x) * scale.x,
										  dPdx.y * scale.x,
										  dPdy.y * scale.x,
										  dPdx.z * scale.x,
										  dPdy.z * scale.x,
                                          cellSoftness,
                                          rot.x,
                                          rotjitter.x,
                                          sg,
                                          handle,
                                          params,
                                          &textureResult[0]
                                          );
    } else {
        textureResult[0] = AI_RGBA_BLACK;
    }


    bool textureAccessY = true;
    if(weights[1] > 0.){
        textureAccessY = lookupCellNoise((P.x + offset.y) * scale.y,
                                          (P.z + offset.y) * scale.y,
    									  dPdx.x * scale.y,
										  dPdy.x * scale.y,
										  dPdx.z * scale.y,
										  dPdy.z * scale.y,
                                          cellSoftness,
                                          rot.y,
                                          rotjitter.y,
                                          sg,
                                          handle,
                                          params,
                                          &textureResult[1]
                                          );
    } else {
        textureResult[1] = AI_RGBA_BLACK;
    }

    bool textureAccessZ = true;
    if(weights[2] > 0.){
        textureAccessZ = lookupCellNoise((P.y + offset.z) * scale.z,
                                          (P.x + offset.z) * scale.z,
        								  dPdx.y * scale.z,
										  dPdy.y * scale.z,
										  dPdx.x * scale.z,
										  dPdy.x * scale.z,
                                      	  cellSoftness,
                                          rot.z,
                                          rotjitter.z,
                                          sg,
                                          handle,
                                          params,
                                          &textureResult[2]
                                          );
    } else {
        textureResult[2] = AI_RGBA_BLACK;
    }

    if(textureAccessX && textureAccessY && textureAccessZ){
        AtRGBA result = AI_RGBA_BLACK;
        result += textureResult[0] * weights[0];
        result += textureResult[1] * weights[1];
        result += textureResult[2] * weights[2];
        return result;
    }
    else {
        // Something went wrong during lookup.
        // TODO: Log the error
        return AI_RGBA_RED;
    }
}

shader_evaluate
{
    // get shader parameters

    AtRGB input = AiShaderEvalParamRGB(p_input);

	int space = AiShaderEvalParamInt(p_space);

    int normal = AiShaderEvalParamInt(p_normal);

    int tiling = AiShaderEvalParamInt(p_tiling);

    float frequency = AiShaderEvalParamFlt(p_frequency);

    float blendSoftness = AiShaderEvalParamFlt(p_blendSoftness);
    blendSoftness = CLAMP(blendSoftness, 0.f, 1.f);

    float cellSoftness = AiShaderEvalParamFlt(p_cellSoftness);
    cellSoftness = CLAMP(cellSoftness, 0.f, 1.f);

    AtPoint scale;
    scale.x = 1.f/AiShaderEvalParamFlt(p_scalex);
    scale.y = 1.f/AiShaderEvalParamFlt(p_scaley);
    scale.z = 1.f/AiShaderEvalParamFlt(p_scalez);

    AtPoint offset;
    offset.x = AiShaderEvalParamFlt(p_offsetx);
    offset.y = AiShaderEvalParamFlt(p_offsety);
    offset.z = AiShaderEvalParamFlt(p_offsetz);

    AtPoint rot;
    rot.x = AiShaderEvalParamFlt(p_rotx);
    rot.y = AiShaderEvalParamFlt(p_roty);
    rot.z = AiShaderEvalParamFlt(p_rotz);

    AtPoint rotjitter;
    rotjitter.x = AiShaderEvalParamFlt(p_rotjitterx);
    rotjitter.y = AiShaderEvalParamFlt(p_rotjittery);
    rotjitter.z = AiShaderEvalParamFlt(p_rotjitterz);

        // get local data
    ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

        // set up P and blend weights
    AtPoint P;
    AtPoint N;
	AtVector dPdx;
	AtVector dPdy;
	SGCache SGC;
	SGC.initCache(sg);
    getProjectionGeometry(node, sg, space, normal, &P, &N, &dPdx, &dPdy);
    float weights[3];
    computeBlendWeights(N, space, blendSoftness, weights);

    P *= frequency;
        // compute texture values
    AtRGBA result = AI_RGBA_RED;
    switch(tiling){
        case TM_CELLNOISE:
            result = tileCellnoise(P, dPdx, dPdy, scale, offset, weights, cellSoftness, rot, rotjitter, sg, data->texturehandle, data->textureparams);
            break;
        case TM_REGULAR:
            result = tileRegular(P, dPdx, dPdy, scale, offset, weights, rot, sg, data->texturehandle, data->textureparams);
            break;
        default:
            // TODO: We should never end up here. Log the error to inform the shader writer.
            result = AI_RGBA_BLUE;
            break;
    }
    sg->out.RGB() = lerp(input, result.rgb(), result.a);

		// clean up after ourselves
	SGC.restoreSG(sg);
}


