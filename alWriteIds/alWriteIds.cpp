#include <cassert>
#include <string.h>

#include <ai.h>

#include <ids.h>

using namespace WriteIds;

AI_SHADER_NODE_EXPORT_METHODS(alWriteIds);

struct ShaderData;
typedef void (*WriteIdsForAtStringsToFloatAOVsFunc)(AtShaderGlobals*, const ShaderData*);
struct ShaderData
{
    AtString aov_names[WI_NUM_AOVS];
    WriteIdsForAtStringsToFloatAOVsFunc writeIdsFunc;
};

void writeAtUInt64IdForAtStringToFloatAOVs(AtUInt64 id64,
                                           const AtString str,
                                           AtShaderGlobals* sg,
                                           const ShaderData* data)
{
    void* string = const_cast<void*>(static_cast<const void*>(str.c_str()));
    AiAOVSetPtr(sg, data->aov_names[WI_AOV_STR_INDEX], string);
    uint32_t* id32 = reinterpret_cast<uint32_t*>(&id64);
    for (int f = WI_AOV_FLT1_INDEX; f <= WI_AOV_FLT2_INDEX; f++)
    {
        AiAOVSetFlt(sg, data->aov_names[f], *reinterpret_cast<float*>(id32));
        id32++;
    }
}

bool hasParamEntry(const AtNode *node, const AtString param, AtByte type)
{
    const AtParamEntry *pentry = AiNodeEntryLookUpParameterAtString(AiNodeGetNodeEntry(node), param );
    return (pentry && (AiParamGetType(pentry) == type));
}

void writeIdsFromShaderNodeNameAtStringHashes(AtShaderGlobals* sg, const ShaderData* data)
{
    assert(hasParamEntry(sg->Op, shaderParamName, AI_TYPE_ARRAY) && "there should always be a shader array param" );
    const AtArray* shaderArray = AiNodeGetArray(sg->Op, shaderParamName);
    assert(shaderArray->type == AI_TYPE_NODE && "the shader array param should hold node types");
    assert(shaderArray->nelements >= 1 && "the shader array param should have one or more elements");
    for(uint32_t i = 0; i < shaderArray->nelements; i++)
    {
        const AtNode *shaderNode = static_cast<AtNode *>(AiArrayGetPtr(shaderArray, i));
        if (AiNodeEntryGetType(AiNodeGetNodeEntry(shaderNode)) == AI_NODE_SHADER)
        {
            const AtString shaderName = AiNodeGetStr(shaderNode, nameParamName);
            assert(sizeof(AtUInt64) == sizeof(size_t));
            AtUInt64 id = (AtUInt64) shaderName.hash();
            writeAtUInt64IdForAtStringToFloatAOVs(id, shaderName, sg, data);
        }
    }
}

void writeIdFromShapeNodeNameAtStringHash(AtShaderGlobals* sg,const ShaderData* data)
{
    const AtString shapeName = AiNodeGetStr(sg->Op, nameParamName);
    assert(sizeof(AtUInt64) == sizeof(size_t));
    AtUInt64 id = (AtUInt64) shapeName.hash();
    writeAtUInt64IdForAtStringToFloatAOVs(id, shapeName, sg, data);
}

void writeIdFromProceduralNodeNameAtStringHash(AtShaderGlobals* sg, const ShaderData* data)
{
    if (sg->proc == NULL) { return; } // There may or may not be a procedural node
    const AtString procName = AiNodeGetStr(sg->proc, nameParamName);
    assert(sizeof(AtUInt64) == sizeof(size_t));
    AtUInt64 id = (AtUInt64) procName.hash();
    writeAtUInt64IdForAtStringToFloatAOVs(id, procName, sg, data);
}

void writeIdsToFloatAOVs(AtShaderGlobals* sg, const ShaderData* data)
{
    for (int i = 0; i < WI_NUM_AOVS; i++)
    {
        if (data->aov_names[i].empty()) { return; }
        if (!AiAOVEnabledAtString(data->aov_names[i], aov_types[i])) { return; }
    }
    data->writeIdsFunc(sg, data);
}

enum WriteIdsParams
{
    p_passthrough,
    p_stringOption,
    p_aov_string,
    p_aov_idFloat1,
    p_aov_idFloat2,
};

node_parameters // static void Parameters(AtList* params, AtMetaDataStore* mds)
{
    AiParameterRGB("passthrough", 0.0f, 0.0f, 0.0f);
    AiParameterENUM("string", WI_SHAPE_NAME, stringOptionsNames);
    for (int i = 0; i < WI_NUM_AOVS; i++)
    {
        AiParameterStr(aov_paramNames[i], aov_namesDefaultValues[i]);
        AiMetaDataSetBool(mds, aov_paramNames[i], "linkable", false);
    }
}

node_initialize // static void Initialize(AtNode* node, AtParamValue* params)
{
    ShaderData* data = new ShaderData();
    AiNodeSetLocalData(node, data);
}

node_update // static void Update(AtNode* node, AtParamValue* params)
{
    ShaderData *data = reinterpret_cast<ShaderData*>(AiNodeGetLocalData(node));
    switch (params[p_stringOption].INT)
    {
        case WI_SHAPE_NAME:
            data->writeIdsFunc = writeIdFromShapeNodeNameAtStringHash;
            break;
        case WI_PROCEDURAL_NAME:
            data->writeIdsFunc = writeIdFromProceduralNodeNameAtStringHash;
            break;
        case WI_SHADER_NAME:
            data->writeIdsFunc = writeIdsFromShaderNodeNameAtStringHashes;
            break;
        default:
            AiMsgWarning("[alWriteIds] Unknown string option, defaulting to shape node name");
            data->writeIdsFunc = writeIdFromShapeNodeNameAtStringHash;
            break;
    }
    for (int i = 0; i < WI_NUM_AOVS; i++)
    {
        data->aov_names[i] = AiNodeGetStr(node, aov_paramNames[i]);
        AiAOVRegister(data->aov_names[i].c_str(), aov_types[i], AI_AOV_BLEND_NONE);
    }
}

node_finish // static void Finish(AtNode* node)
{

    ShaderData* data = reinterpret_cast<ShaderData*>(AiNodeGetLocalData(node));
    delete data;
}

shader_evaluate // static void Evaluate(AtNode* node, AtShaderGlobals* sg)
{
    sg->out.RGB() = AiShaderEvalParamRGB(p_passthrough);
    if (!(sg->Rt & AI_RAY_CAMERA)) { return; }
    ShaderData* data = reinterpret_cast<ShaderData*>(AiNodeGetLocalData(node));
    writeIdsToFloatAOVs(sg, data);
}

node_loader
{
    if (i > 0)
        return false;

    node->methods     = alWriteIds;
    node->output_type = AI_TYPE_RGB;
    node->name        = "alWriteIds";
    node->node_type   = AI_NODE_SHADER;
    strcpy(node->version, AI_VERSION);
    return true;
}
