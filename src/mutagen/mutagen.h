class WaveFile {

struct WavInfo
{
    FILE* pHandle;
    unsigned short wNumChannels;
    unsigned int lSampleRate;
    unsigned short wBytesPerSample;
    unsigned int lNumFrames;
    unsigned int lNumDataBytes;
};

int wavOpenRead(char* pszPath, WavInfo* pWavInfo);
FILE* wavOpenWrite(char* pszPath);
int wavRead(WavInfo* pWavInfo, float* pwData);
int wavWrite(WavInfo* pWavInfo, float* pwData);

}