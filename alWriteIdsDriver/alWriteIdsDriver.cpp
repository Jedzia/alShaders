#include <cassert>
#include <string.h>
#include <fstream>
#include <map>
#include <vector>

#include <ai.h>

#include <ids.h>

using namespace WriteIds;

AI_DRIVER_NODE_EXPORT_METHODS(alWriteIdsDriver);


struct FloatAOVs
{
    float values [WI_NUM_FLT_AOVS];
};

typedef std::map<const char*, FloatAOVs> CstringsToFloatAOVsMap;
typedef std::vector< std::pair<const char*, int> > NameAndTypePairVector;

typedef struct {
    bool wrongOutputs;
    NameAndTypePairVector outputs;
    CstringsToFloatAOVsMap cStringToFloatAOVsMap;
} DriverData;

node_loader
{
    if (i>0)
      return false;
    node->methods = alWriteIdsDriver;
    node->output_type = AI_TYPE_NONE;
    node->name = "driver_ids";
    node->node_type = AI_NODE_DRIVER;
    strcpy(node->version, AI_VERSION);
    return true;
}

node_parameters
{
    AiParameterStr("filename", "shapeName.ids");
}

node_initialize
{
    DriverData *data = new DriverData();
    static const char *required_aovs[] = { "FLOAT Z", "FLOAT A", NULL };
    AiRawDriverInitialize(node, required_aovs, true,  data);
}

node_update
{
}

driver_supports_pixel_type
{
    // This function is not needed for a raw driver
    return true;
}

driver_open
{
    DriverData * data = reinterpret_cast<DriverData *>(AiDriverGetLocalData(node));

    // Get the AOV outputs we expect
    const char *name = NULL;
    int type = AI_TYPE_UNDEFINED;
    int nptr = 0;
    int nflt = 0;
    while (AiOutputIteratorGetNext(iterator, &name, &type, 0))
    {
        switch(type)
        {
            case AI_TYPE_POINTER :
            {
                if (nptr < WI_NUM_STR_AOVS)
                {
                    std::pair<const char*, int> output(name, type);
                    data->outputs.push_back(output);
                    nptr++;
                }
                break;
            }
            case AI_TYPE_FLOAT :
            {
                if (nflt < WI_NUM_FLT_AOVS)
                {
                    std::pair<const char*, int> output(name, type);
                    data->outputs.push_back(output);
                    nflt++;
                }
                break;
            }
        }
    }

    // Sanity check the AOV outputs
    if ((data->outputs.size() == WI_NUM_AOVS))
    {
        data->wrongOutputs = false;

        std::string infoStr = "";
        for (NameAndTypePairVector::iterator it = data->outputs.begin(); it != data->outputs.end(); ++it)
        {
            infoStr += "\"";
            infoStr += AiParamGetTypeName(it->second);
            infoStr += "\" ";
            infoStr += it->first;
            infoStr += " ";
        }
        AiMsgInfo("[alWriteIdsDriver] correct aov outputs found in this order -- %s", infoStr.c_str());
    }
    else
    {
        AiMsgWarning("[alWriteIdsDriver] exactly one pointer and two float aov outputs were not specified - ignoring all outputs");
        data->wrongOutputs = true;
    }
}

driver_extension
{
   static const char *extensions[] = { "ids", NULL };
   return extensions;
}

driver_needs_bucket
{
   return true;
}

driver_prepare_bucket
{
}

driver_process_bucket
{
}

bool AOVSampleIteratorHasAllAOVValues(const AtAOVSampleIterator* iter, std::vector< std::pair<const char*, int> >& outputs)
{
    bool result = false;
    for (NameAndTypePairVector::iterator it = outputs.begin(); it != outputs.end(); ++it)
    {
        if (!AiAOVSampleIteratorHasAOVValue(iter, it->first, it->second))
        {
            result = false;
            break;
        }
        result = true;
    }
    return result;
}

driver_write_bucket
{
    DriverData * data = reinterpret_cast<DriverData *>(AiDriverGetLocalData(node));
    if (data->wrongOutputs) { return; }

    for (int y = bucket_yo; y < bucket_yo + bucket_size_y; y++)
    {
        for (int x = bucket_xo; x < bucket_xo + bucket_size_x; x++)
        {
            // Iterator for samples on this pixel
            AiAOVSampleIteratorInitPixel(sample_iterator, x, y);
            while (AiAOVSampleIteratorGetNext(sample_iterator))
            {
                while (AiAOVSampleIteratorGetNextDepth(sample_iterator))
                {
                    if (AOVSampleIteratorHasAllAOVValues(sample_iterator, data->outputs))
                    {
                        int nflt = 0;
                        std::pair<const char*, FloatAOVs> mapping;
                        for (NameAndTypePairVector::iterator it = data->outputs.begin(); it != data->outputs.end(); ++it)
                        {
                            switch(it->second)
                            {
                                case AI_TYPE_POINTER :
                                {
                                    mapping.first =  static_cast<const char*>(AiAOVSampleIteratorGetAOVPtr(sample_iterator, it->first));
                                    break;
                                }
                                case AI_TYPE_FLOAT :
                                {
                                    // Retain the order in which the floats were specified in the outputs
                                    assert(nflt < WI_NUM_FLT_AOVS);
                                    mapping.second.values[nflt] = AiAOVSampleIteratorGetAOVFlt(sample_iterator, it->first);
                                    nflt++;
                                    break;
                                }
                            }
                        }
                        // Inserted only if string key is not equivalent to the key of any other element already in the container
                        data->cStringToFloatAOVsMap.insert(mapping);
                    }
                }
            }
        }
    }
}

driver_close
{
    DriverData * data = reinterpret_cast<DriverData *>(AiDriverGetLocalData(node));
    if (data->wrongOutputs) { return; }

    std::ofstream idsOutFileStream;
    AtString filename = AiNodeGetStr(node, filenameParamName);
    try
    {
        // Arnold should have already checked that the dir is writeable by this point
        idsOutFileStream.open(filename.c_str(), std::ios::out | std::ios::trunc);
    }
    catch(std::ofstream::failure &writeErr)
    {
        AiMsgWarning("[alWriteIdsDriver] failed to open %s (%s)", filename.c_str(), writeErr.what());
    }
    assert(idsOutFileStream.is_open() && "File is not opened");

    if (idsOutFileStream.is_open())
    {
        for ( CstringsToFloatAOVsMap::iterator it = data->cStringToFloatAOVsMap.begin(); it != data->cStringToFloatAOVsMap.end(); ++it )
        {
            idsOutFileStream  << it->first << " ";
            for(int i = 0; i < WI_NUM_FLT_AOVS; i++)
            {
                idsOutFileStream << *reinterpret_cast<const uint32_t*>(&it->second.values[i]) << " ";
            }
            idsOutFileStream <<  std::endl;
        }
        idsOutFileStream.close();
        AiMsgInfo("[alWriteIdsDriver] wrote id mappings to file %s", filename.c_str());
    }
}

node_finish
{
    DriverData * data = reinterpret_cast<DriverData *>(AiDriverGetLocalData(node));
    // Free local data
    delete data;
    AiDriverDestroy(node);
}
