#include "alUtil.h"
#include <ai.h>
#include <string.h>
#include <sstream>

/*
    This shader is pretty stupid, but I really needed a way to grab per-frame textures in mtoa, so here goes..
*/

AI_SHADER_NODE_EXPORT_METHODS(alImageMtd)

enum alImageParams
{
    p_filename,
    p_frame
};

node_parameters
{
    AiParameterStr("filename", "");
    AiParameterInt("frame", 0);
}

node_loader
{
   if (i>0) return 0;
   node->methods     = alImageMtd;
   node->output_type = AI_TYPE_RGB;
   node->name        = "alImage";
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

    int framenum = params[p_frame].INT;
    std::stringstream out;
    out << framenum;
    std::string frameString;
    if(framenum > 999)
        frameString = out.str();
    else if(framenum > 99)
        frameString = out.str();
    else if(framenum > 9)
        frameString = out.str();
    else
        frameString = out.str();

    std::string texname = std::string(params[p_filename].STR);
    size_t found = texname.find("####");
    if (found!=std::string::npos){
        texname.replace(found,4,frameString);
    }

    data->texturehandle = AiTextureHandleCreate(texname.c_str());
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

}

shader_evaluate
{
   ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

   bool success;
   sg->out.RGB() = AiTextureHandleAccess(sg, data->texturehandle, *data->textureparams, &success).rgb();
}


