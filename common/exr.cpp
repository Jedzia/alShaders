#include <OpenEXR/ImfHeader.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>

void writeRGBEXR(const char* filename, const float* pixels, int w, int h)
{
    Imf::Header header(w, h);
    header.channels().insert("r", Imf::Channel(Imf::FLOAT));
    header.channels().insert("g", Imf::Channel(Imf::FLOAT));
    header.channels().insert("b", Imf::Channel(Imf::FLOAT));

    Imf::OutputFile file(filename, header);

    Imf::FrameBuffer fb;

    fb.insert("r", Imf::Slice(Imf::FLOAT, (char*)pixels, sizeof(float)*3, sizeof(float)*3*w));
    fb.insert("g", Imf::Slice(Imf::FLOAT, (char*)(pixels+1), sizeof(float)*3, sizeof(float)*3*w));
    fb.insert("b", Imf::Slice(Imf::FLOAT, (char*)(pixels+2), sizeof(float)*3, sizeof(float)*3*w));

    file.setFrameBuffer(fb);
    file.writePixels(h);
}

#define THICKNESS 64
void writeThickRGBEXR(const char* filename, const float* pixels, int w, int h)
{
    float* tmp = new float[3 * w * THICKNESS];
    float* tmp2 = tmp;
    for (int i=0; i < THICKNESS; ++i)
    {   
        memcpy(tmp2, pixels, sizeof(float)*3*w);
        tmp2 += w*3;
    }
    writeRGBEXR(filename, tmp, w, THICKNESS);
    delete[] tmp;
}